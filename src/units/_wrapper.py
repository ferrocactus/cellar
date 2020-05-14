from . import *


def wrap(step, method):
    """
    Wrapper function that takes a pipeline step and method
    and returns the corresponding object.

    Args:
        step (string): The step in the pipeline.
        method (string): Method to use in the given step.
    Returns:
        object (Unit): Object of the right type.
    """
    if step not in translation_dict:
        raise NotImplementedError("{0} step not implemented.".format(step))
    if method not in translation_dict[step]:
        raise NotImplementedError("{0} method not implemented.".format(method))
    return translation_dict[step][method]
