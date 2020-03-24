r"""Support routines for running the optimization in parallel across the
objectives.

The time-propagation that is the main numerical effort in an optimization with
Krotov's method can naturally be performed in parallel for the different
objectives. There are three time-propagations that happen inside
:func:`.optimize_pulses`:

1. A forward propagation of the :attr:`~.Objective.initial_state` of each
   objective under the initial guess pulse.

2. A backward propagation of the states $\ket{\chi_k}$ constructed by the
   `chi_constructor` routine that is passed to :func:`.optimize_pulses`, where
   the number of states is the same as the number of objectives.

3. A forward propagation of the :attr:`~.Objective.initial_state` of each
   objective under the optimized pulse in each iteration. This can only
   be parallelized *per time step*, as the propagated states from each time
   step collectively determine the pulse update for the next time step, which
   is then used for the next propagation step. (In this sense Krotov's method
   is "sequential")


The :func:`.optimize_pulses` routine has a parameter `parallel_map` that can
receive a tuple of three "map" functions to enable parallelization,
corresponding to the three propagation listed above. If not given,
:func:`qutip.parallel.serial_map` is used for all three propagations, running
in serial. Any alternative "map" must have the same interface as
:func:`qutip.parallel.serial_map`.

It would be natural to assume that :func:`qutip.parallel.parallel_map`
respectively the slightly improved :func:`parallel_map` provided in this module
would be a good choice for parallel execution, using multiple CPUs on the same
machine.  However, this function is only a good choice for the propagation (1)
and (2): these run in parallel over the entire time grid without any
communication, and thus minimal overhead.  However, this is not true for the
propagation (3), which must synchronize after each time step. In that case, the
"naive" use of :func:`qutip.parallel.parallel_map` results in a communication
overhead that completely dominates the propagation, and actually makes the
optimization slower (potentially by more than an order of magnitude).

The function :func:`parallel_map_fw_prop_step` provided in this module is an
appropriate alternative implementation that uses long-running processes,
internal caching, and minimal inter-process communication to eliminate the
communication overhead as much as possible. However, the internal caching is
valid only under the assumption that the `propagate` function does not have
side effects.

In general,

.. code-block:: python

    parallel_map=(
        krotov.parallelization.parallel_map,
        krotov.parallelization.parallel_map,
        krotov.parallelization.parallel_map_fw_prop_step,
    )

is a decent choice for enabling parallelization for a typical multi-objective
optimization (but don't expect wonders: general pure-python parallelization is
an unsolved problem.)

You may implement your own "map" functions to exploit parallelization paradigms
other than Python's built-in :mod:`multiprocessing`, provided here. This
includes distributed propagation, e.g. through ipyparallel_ clusters. To write
your own `parallel_map` functions, review the source code of
:func:`.optimize_pulses` in detail.

In most cases, it will be difficult to obtain a linear speedup from
parallelization: even with carefully tuned manual interprocess communication,
the communication overhead can be substantial. For best results, it would be
necessary to use `parallel_map` functions implemented in Cython, where the GIL
can be released and the entire propagation (and storage of propagated states)
can be done in shared-memory with no overhead.

Also note that the overhead of multi-process parallelization is
platform-dependent. On Linux, subprocesses are "forked" which causes them to
inherit the current state of the parent process without any explicit (and
expensive) inter-process communication (IPC). On other platforms, most notably
Windows and the combination of macOS with Python 3.8, subprocesses are
"spawned" instead of "forked": The subprocesses start from a clean slate, and
all objects must be transfered from the parent process via IPC. This is very
slow, and you should not expect to be able to achieve any speedup from
parallelization on such platforms.

Another caveat on platforms using "spawn" is that certain objects by default
cannot be transferred via IPC, due to limitations of the :mod:`pickle`
protocol. This affects :ref:`lambda` and functions defined in
Jupyter notebooks, in particular. The third-party :mod:`loky` library provides
an alternative implementation for multi-processes parallelization that does not
have these restrictions, but causes even more overhead.

You may attempt to use the various options to :func:`set_parallelization` in
order to find a combination of settings that minimizes the runtime in your
particular environment.

.. _ipyparallel: https://ipyparallel.readthedocs.io/en/latest/
"""
import contextlib
import multiprocessing
from concurrent.futures import ProcessPoolExecutor
from functools import partial

from qutip.parallel import serial_map
from qutip.ui.progressbar import BaseProgressBar, TextProgressBar
from threadpoolctl import threadpool_limits

from .conversions import plug_in_pulse_values


try:

    import loky
    from loky import get_reusable_executor as LokyReusableExecutor
    from loky.process_executor import (
        ProcessPoolExecutor as LokyProcessPoolExecutor,
    )

    _HAS_LOKY = True

except ImportError:

    _HAS_LOKY = False

USE_LOKY = False
"""Whether to use :mod:`loky` instead of :mod:`multiprocessing`.

Set by :func:`set_parallelization`.
"""

USE_THREADPOOL_LIMITS = True
"""Whether to limit the number of low-level BLAS/OpenMP threads.

When using multi-process parallelization, *nested parallelization* must be
avoided. That is, low-level numerical routines e.g. in :mod:`numpy` should not
be allowed to use multiple threads. This would lead to over-subscribing CPUs
and can slow down the entire program by orders of magnitude.

If True, threadpoolctl_ will be used internally to attempt to eliminate any
nested threads.

Set by :func:`set_parallelization`.

.. note::

    Alternatively (or in addition), you may want to consider setting the
    following environment variables in your shell:

    .. code-block:: shell

        export MKL_NUM_THREADS=1
        export NUMEXPR_NUM_THREADS=1
        export OMP_NUM_THREADS=1

.. _threadpoolctl: https://github.com/joblib/threadpoolctl
"""


__all__ = [
    'set_parallelization',
    'parallel_map',
    'serial_map',
    'Consumer',
    'FwPropStepTask',
    'parallel_map_fw_prop_step',
]


@contextlib.contextmanager
def _no_threadpool_limits(*args, **kwargs):  # pragma: nocover
    """No-op replacement for :func:`threadpool_limits`."""
    # this is just an empty context manager that does nothing.
    yield None


def set_parallelization(
    use_loky=False,
    start_method=None,
    loky_pickler=None,
    use_threadpool_limits=True,
):  # pragma: nocover
    """Configure multi-process parallelization.

    Args:
        use_loky (bool): Value for :obj:`USE_LOKY`.
        start_method (None or str): One of 'fork', 'spawn', and 'forkserver',
            see :func:`multiprocessing.set_start_method`. If ``use_loky=True``,
            also 'loky' and 'loky_int_main', see :mod:`loky`. If None, a
            platform and version-dependent default will be chosen automatically
            (e.g., 'fork' on Linux, 'spawn' on Windows, 'loky' if
            ``use_loky=True``)
        loky_pickler (None or str): Serialization module to use for
            :mod:`loky`. One of 'cloudpickle', 'pickle'. This forces the
            serialization for *all* objects. The default value None chooses the
            serialization automatically depending of the type of object. Using
            'cloudpickle' is signficiantly slower than 'pickle' (but 'pickle'
            cannot serialize all objects, such as lambda functions or functions
            defined in a Jupyter notebook).
        use_threadpool_limits (bool): Value for :obj:`USE_THREADPOOL_LIMITS`.

    Raises:
        ImportError: if ``use_loky=True`` but :mod:`loky` is not installed.

    Note:
        When working in Jupyter notebooks on systems that use the 'spawn'
        `start_method` (Windows, or macOS with Python >= 3.8), you may have to
        use :mod:`loky` (``use_loky=True``). This will incur a signficiant
        increase in multi-processing overhead. Use Linux if you can.

    Warning:
        This function should only be called once per script/notebook, at its
        very beginning. The :obj:`USE_LOKY` and :obj:`USE_THREADPOOL_LIMITS`
        variables may be set at any time.
    """
    global USE_LOKY
    start_methods = ['fork', 'spawn', 'forkserver']
    if use_loky:
        start_methods.extend(['loky', 'loky_int_main'])
    if start_method is not None:
        if start_method not in start_methods:
            raise ValueError("start_method not in %s" % str(start_methods))
    if use_loky:
        if not _HAS_LOKY:
            raise ImportError("The loky library is not installed.")
        USE_LOKY = True
        loky.backend.context.set_start_method(start_method)
        if loky_pickler is not None:
            loky.set_loky_pickler(loky_pickler)
    else:
        multiprocessing.set_start_method(start_method)


def parallel_map(
    task,
    values,
    task_args=None,
    task_kwargs=None,
    num_cpus=None,
    progress_bar=None,
):
    """Map function `task` onto `values`, in parallel.

    This function's interface is identical to
    :func:`qutip.parallel.parallel_map` as of QuTiP 4.5.0, but has the option
    of using :mod:`loky` as a backend (see :func:`set_parallelization`). It
    also eliminates internal threads, according to
    :obj:`USE_THREADPOOL_LIMITS`.
    """
    # TODO: if QuTiP's parallel_map catches up, we can remove this function,
    # and put QuTiP's parallel_map into __all__ to maintain krotov's interface.
    if task_args is None:
        task_args = ()
    if task_kwargs is None:
        task_kwargs = {}

    if num_cpus is None:
        num_cpus = multiprocessing.cpu_count()

    if progress_bar is None:
        progress_bar = BaseProgressBar()
    if progress_bar is True:
        progress_bar = TextProgressBar()

    progress_bar.start(len(values))
    nfinished = [0]

    def _update_progress_bar(x):
        nfinished[0] += 1
        progress_bar.update(nfinished[0])

    if USE_LOKY:
        Executor = LokyReusableExecutor
        if USE_THREADPOOL_LIMITS:
            Executor = partial(
                LokyReusableExecutor,
                initializer=_process_threadpool_limits_initializier,
            )
    else:
        Executor = ProcessPoolExecutor

    _threadpool_limits = _no_threadpool_limits
    if USE_THREADPOOL_LIMITS:
        _threadpool_limits = threadpool_limits

    with _threadpool_limits(limits=1):
        with Executor(max_workers=num_cpus) as executor:
            jobs = []
            try:
                for value in values:
                    args = (value,) + tuple(task_args)
                    job = executor.submit(task, *args, **task_kwargs)
                    job.add_done_callback(_update_progress_bar)
                    jobs.append(job)
                res = [job.result() for job in jobs]
            except KeyboardInterrupt as e:
                raise e

    progress_bar.finished()
    return res


def _process_threadpool_limits_initializier():
    """Initializer for settings threadpool limits.

    This is an initializer for :mod:`loky` Executors that deactivates threads
    in the spawned sub-processes.
    """
    import numpy  # required for loky's autodetection
    from threadpoolctl import threadpool_limits

    threadpool_limits(limits=1)


class Consumer(multiprocessing.Process):
    """A process-based task consumer.

    This is for internal use in :func:`parallel_map_fw_prop_step`.

    Args:
        task_queue (multiprocessing.JoinableQueue): A queue from which to read
            tasks.
        result_queue (multiprocessing.Queue): A queue where to put the results
            of a task
        data: cached (in-process) data that will be passed to each task
    """

    def __init__(self, task_queue, result_queue, data):
        super().__init__()
        self.task_queue = task_queue
        self.result_queue = result_queue
        self.data = data

    def run(self):
        """Execute all tasks on the `task_queue`.

        Each task must be a callable that takes `data` as its only argument.
        The return value of the task will be put on the `result_queue`. A None
        value on the `task_queue` acts as a "poison pill", causing the
        :class:`Consumer` process to shut down.
        """
        _threadpool_limits = _no_threadpool_limits
        if USE_THREADPOOL_LIMITS:
            _threadpool_limits = threadpool_limits

        while True:
            next_task = self.task_queue.get()
            if next_task is None:
                # Poison pill means shutdown
                self.task_queue.task_done()
                break
            with _threadpool_limits(limits=1):
                answer = next_task(self.data)
            self.task_queue.task_done()
            self.result_queue.put(answer)


class FwPropStepTask:
    """A task that performs a single forward-propagation step.

    The task object is a callable, receiving the single tuple of the same
    form as `task_args` in :func:`parallel_map_fw_prop_step` as input. This
    `data` is internally cached by the :class:`Consumer` that will execute
    the task.

    This is for internal use in :func:`parallel_map_fw_prop_step`.

    Args:
        i_state (int): The index of the state to propagation. That is, the
            index of the objective from whose :attr:`~.Objective.initial_state`
            the propagation started
        pulse_vals (list[float]): the values of the pulses at `time_index` to
            use.
        time_index (int): the index of the interval on the time grid covered by
            the propagation step

    The passed arguments update the internal state (`data`) of the
    :class:`Consumer` executing the task; they are the minimal information that
    must be passed via inter-process communication to enable the forward
    propagation (assuming `propagate` in :func:`.optimize_pulses` has no
    side-effects)
    """

    def __init__(self, i_state, pulse_vals, time_index):
        # The task object itself gets send to the Consumer via IPC (pickling).
        # Since this is only a few scalar values, we largely avoid
        # communication overhead
        # warnings.warn(
        #    "FwPropStepTask is deprecated and will be removed in version 2.0",
        #    warnings.DeprecationWarning,
        # )
        self.i_state = i_state
        self.pulse_vals = pulse_vals
        self.time_index = time_index

    def __call__(self, data):
        (
            states,
            objectives,
            pulses,
            pulses_mapping,
            tlist,
            _,  # time_index
            propagators,
        ) = data
        # the data is passed by the Consumer, and is cached locally inside of
        # each process. Thus, it does not contribute to the IPC communication
        # overhead
        time_index = self.time_index
        pulse_vals = self.pulse_vals
        i_state = self.i_state
        for (pulse, pulse_val) in zip(pulses, pulse_vals):
            pulse[time_index] = pulse_val

        for (pulse, pulse_val) in zip(pulses, pulse_vals):
            pulse[time_index] = pulse_val
        state = states[i_state]
        mapping = pulses_mapping[i_state]
        obj = objectives[i_state]
        H = plug_in_pulse_values(obj.H, pulses, mapping[0], time_index)
        c_ops = [
            plug_in_pulse_values(c_op, pulses, mapping[ic + 1], time_index)
            for (ic, c_op) in enumerate(obj.c_ops)
        ]
        dt = tlist[time_index + 1] - tlist[time_index]
        states[i_state] = propagators[i_state](H, state, dt, c_ops)
        # While there is no significant IPC-communication overhead associated
        # with the *input* of the task, the resulting state returned here still
        # must go through the `result_queue` of the Consumer. This is the main
        # bottleneck of this implementation.
        return states[i_state]


def parallel_map_fw_prop_step(shared, values, task_args):
    """`parallel_map` function for the forward-propagation by one time step.

    Args:
        shared: A global object to which we can attach attributes for sharing
            data between different calls to :func:`parallel_map_fw_prop_step`,
            allowing us to have long-running :class:`Consumer` processes,
            avoiding process-management overhead.
            This happens to be a callable (the original internal routine for
            performing a forward-propagation), but here, it is (ab-)used as a
            storage object only.
        values (list): a list 0..(N-1) where N is the number of objectives
        task_args (tuple): A tuple of 7 components:

            1. A list of states to propagate, one for each objective.
            2. The list of objectives
            3. The list of optimized pulses (updated up to `time_index`)
            4. The "pulses mapping", cf :func:`.extract_controls_mapping`
            5. The list of time grid points
            6. The index of the interval on the time grid over which to
               propagate
            7. A list of `propagate` callables, as passed to
               :func:`.optimize_pulses`.  The propagators must not have
               side-effects in order for :func:`parallel_map_fw_prop_step` to
               work correctly.
    """
    # `shared` is the original task function
    # (krotov.optimize._forward_propagation_step), but here we abuse it
    # as a shared namespace, between calls to `my_map`, by setting custom
    # data attributes on it
    if USE_LOKY:
        return _parallel_map_fw_prop_step_loky(shared, values, task_args)
    tlist = task_args[4]
    pulses = task_args[2]
    time_index = task_args[5]
    n = len(values)
    if time_index == 0:
        # set up a Consumer process that stays active for the entire
        # propagation (multiple calls of `parallel_map_fw_prop_step`)
        shared.tasks = [multiprocessing.JoinableQueue() for i in range(n)]
        shared.results = [multiprocessing.Queue() for i in range(n)]
        shared.consumers = [
            Consumer(shared.tasks[i], shared.results[i], data=task_args)
            for i in range(n)
        ]
        for consumer in shared.consumers:
            consumer.start()
    # Assign tasks
    pulse_vals = [pulses[i][time_index] for i in range(len(pulses))]
    for i_state in values:
        shared.tasks[i_state].put(
            FwPropStepTask(i_state, pulse_vals, time_index)
        )
    for i_state in values:
        shared.tasks[i_state].join()  # wait to finish
    # collect results
    res = [shared.results[i].get() for i in range(n)]

    if time_index == len(tlist) - 2:  # end of time grid
        for i in range(n):
            shared.tasks[i].put(None)  # add poison pill
            shared.tasks[i].join()  # wait to finish
    return res


def _parallel_map_fw_prop_step_loky(shared, values, task_args):
    """Loky-based implementation of :func:`parallel_map_fw_prop_step`."""
    tlist = task_args[4]
    pulses = task_args[2]
    time_index = task_args[5]
    n = len(values)
    if time_index == 0:
        # we only send the full task_args through IPC once, for the first time
        # step. Subsequent time steps will reuse the data
        shared.executors = [
            LokyProcessPoolExecutor(
                max_workers=1,
                initializer=partial(
                    _pmfw_initializer, limit_thread_pool=USE_THREADPOOL_LIMITS
                ),
                initargs=(
                    state_index,
                    task_args[0][state_index],  # initial_state
                    task_args[1][state_index],  # objective
                    task_args[2],  # pulses
                    task_args[3],  # pulses_mapping
                    task_args[4],  # tlist
                    task_args[6][state_index],  # propagator
                ),
            )
            for state_index in range(n)
        ]

    pulse_vals = [pulses[i][time_index] for i in range(len(pulses))]

    res = []
    jobs = []
    for i_state in values:
        jobs.append(
            shared.executors[i_state].submit(
                _pmfw_forward_prop_step, pulse_vals, time_index
            )
        )
    res = [job.result() for job in jobs]

    if time_index == len(tlist) - 2:  # end of time grid
        del shared.executors

    return res


def _pmfw_initializer(
    state_index,
    initial_state,
    objective,
    pulses,
    pulses_mapping,
    tlist,
    propagator,
    limit_thread_pool=True,
):
    """Copy `task_args` into a process-local global variable."""
    # for internal use in _parallel_map_fw_prop_step_loky
    global _pmfw_data

    if limit_thread_pool:
        import numpy  # required for loky's autodetection
        from threadpoolctl import threadpool_limits

        threadpool_limits(limits=1)

    _pmfw_data = dict(
        state_index=state_index,
        state=initial_state,
        objective=objective,
        pulses=pulses,
        pulses_mapping=pulses_mapping,
        tlist=tlist,
        propagator=propagator,
    )


def _pmfw_forward_prop_step(pulse_vals, time_index):
    """Propagate a single step, with minimal non-local data.

    Most of the data taken from the process-local global storage set by the
    initializer.
    """
    # for internal use in _parallel_map_fw_prop_step_loky
    global _pmfw_data
    _threadpool_limits = _no_threadpool_limits
    if USE_THREADPOOL_LIMITS:
        _threadpool_limits = threadpool_limits
    i_state = _pmfw_data['state_index']
    pulses = _pmfw_data['pulses']
    for (pulse, pulse_val) in zip(pulses, pulse_vals):
        pulse[time_index] = pulse_val
    mapping = _pmfw_data['pulses_mapping'][i_state]
    obj = _pmfw_data['objective']
    H = plug_in_pulse_values(obj.H, pulses, mapping[0], time_index)
    c_ops = [
        plug_in_pulse_values(c_op, pulses, mapping[ic + 1], time_index)
        for (ic, c_op) in enumerate(obj.c_ops)
    ]
    tlist = _pmfw_data['tlist']
    dt = tlist[time_index + 1] - tlist[time_index]
    with _threadpool_limits(limits=1):
        next_state = _pmfw_data['propagator'](
            H, _pmfw_data['state'], dt, c_ops
        )
    _pmfw_data['state'] = next_state
    return next_state
