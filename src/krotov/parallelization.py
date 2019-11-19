r"""Support routines for running the optimization in parallel across the
objectives

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
:func:`qutip.parallel.serial_map` is used for all three propations, running in
serial. Any alternative "map" must have the same interface as
:func:`qutip.parallel.serial_map`.

It would be natural to assume that :func:`qutip.parallel.parallel_map` would be
a good choice for parallel execution, using multiple CPUs on the same machine.
However, this function is only a good choice for the propagation (1) and (2):
these run in parallel over the entire time grid without any communication, and
thus minimal overhead.  However, this is not true for the propagation (3),
which must synchronize after each time step. In that case, the "naive" use of
:func:`qutip.parallel.parallel_map` results in a communication overhead that
completely dominates the propagation, and actually makes the optimization
slower (potentially by more than an order of magnitude).

The function :func:`parallel_map_fw_prop_step` provided in this module is an
appropriate alternative implementation that uses long-running processes,
internal caching, and minimal inter-process communication to eliminate the
communication overhead as much as possible. However, the internal caching is
valid only under the assumption that the `propagate` function does not have
side effects.

In general,

.. code-block:: python

    parallel_map=(
        qutip.parallel_map,
        qutip.parallel_map,
        krotov.parallelization.parallel_map_fw_prop_step,
    )

is a decent choice for enabling parallelization for a typical multi-objective
optimization.

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

.. _ipyparallel: https://ipyparallel.readthedocs.io/en/latest/
"""
import multiprocessing

from .conversions import plug_in_pulse_values


__all__ = ['Consumer', 'FwPropStepTask', 'parallel_map_fw_prop_step']


class Consumer(multiprocessing.Process):
    """A process-based task consumer

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
        while True:
            next_task = self.task_queue.get()
            if next_task is None:
                # Poison pill means shutdown
                self.task_queue.task_done()
                break
            answer = next_task(self.data)
            self.task_queue.task_done()
            self.result_queue.put(answer)


class FwPropStepTask:
    """A task that performs a single forward-propagation step

    The task object is a callable, receiving the single tuple of the same
    form as `task_args` in :func:`parallel_map_fw_prop_step` as input. This
    `data` is internally cached by the :class:`Consumer` that will execute
    the task.

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
    """`parallel_map` function for the forward-propagation by one time step

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
    # `shared` is the original task function, but here we abuse it
    # as a shared namespace, between calls to `my_map`, by setting custom
    # data attributes on it
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
