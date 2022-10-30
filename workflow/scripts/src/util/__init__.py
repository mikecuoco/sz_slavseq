#!/usr/bin/env python
__author__ = "Apuã Paquola"


def static_var(varname, value):
    def decorate(func):
        setattr(func, varname, value)
        return func
    return decorate


@static_var("complements", str.maketrans('acgtACGT', 'tgcaTGCA'))
def rc(seq):
    """ Reverse complement
    """
    return seq.translate(rc.complements)[::-1]

