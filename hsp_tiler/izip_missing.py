#Author: Eli Draizen
#Date: 5-9-14
#File: izip_missing.py

"""Module to fix itertools.izip: izip_missing 
"""

def izip_missing(iterA, iterB, **kwds):
    """Iterate through two iterables, while making sure they are in the same
    order. If there are missing values, you can skip the value entirely or 
    return only the iterator with the value and a special fill value.

    Parameters:
    ___________
    iterA : the first iterator
    iterB : the second iterator
    key : function that returns items to compare. Must return strings, ints, or
        an object with the __lt__, __gt__, and __eq__ methods. Optional.
    fillvalue : The value to return if the item is missing. Optional.

    Returns:
    ________
    A : item from first iterator, or fillValue
    B : item from second iterator, or fillValue
    """
    #Get the comparison functions
    key = kwds.get("key", lambda x: x)
    keyA = kwds.get("keyA", key)
    keyB = kwds.get("keyB", key)

    useMissing = "fillvalue" in kwds
    fillvalue = kwds.get("fillvalue")

    verbose = kwds.get("verbose", False)

    #Start both iterators
    A = iterA.next()
    B = iterB.next()
    try:
        while True:
            if keyA(A) == keyB(B):
                if verbose:
                    print keyA(A), "==", keyB(B)
                yield A, B
                A = iterA.next()
                B = iterB.next()
            elif keyA(A) < keyB(B):
                if verbose:
                    print keyA(A), "<", keyB(B)
                if useMissing:
                    yield A, fillvalue
                A = iterA.next()
            elif keyA(A) > keyB(B):
                if verbose:
                    print keyA(A), ">", keyB(B)
                if useMissing:
                    yield fillvalue, B
                B = iterB.next()
            else:
                raise RuntimeError("Invalid compartor")
    except StopIteration:
        pass
