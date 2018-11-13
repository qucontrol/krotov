"""Routines to aid in running in parallel"""

__all__ = 'serial_map'


def serial_map(task, values, task_args=tuple(), task_kwargs=None, **kwargs):
    """Apply function `task` to all `values`

    Equivalent to::

        [task(value, *task_args, **task_kwargs) for value in values]
    """
    if task_kwargs is None:
        task_kwargs = {}
    return [task(value, *task_args, **task_kwargs) for value in values]
