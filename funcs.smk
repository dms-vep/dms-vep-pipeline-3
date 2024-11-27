"""Utility functions."""


import io


def yaml_str(o):
    """Return str with YAML representation of Python object `o`."""
    with io.StringIO() as stream:
        y = yaml.YAML().dump(o, stream)
        return stream.getvalue()
