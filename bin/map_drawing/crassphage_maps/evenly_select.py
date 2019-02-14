"""
Select some events
"""
import numpy as np


def evenly_select(n, m):
    """
    Evenly select M elements from a list of length N.
    This returns a list of True/False (actually 0,1)
    See https://stackoverflow.com/questions/46494029/nearly-evenly-select-items-from-a-list

    Then use itertools.compress to create the new list

    :param n: The length of the list
    :param m: The number of elements to return
    :return: A list of 0/1 where 1 should be selected
    """
    if n == m:
        return np.ones(n, dtype=int)
    assert n > m
    if m > n/2:
        cut = np.ones(n, dtype=int)
        q, r = divmod(n, n - m)
        indices = [q * i + min(i, r) for i in range(n - m)]
        cut[indices] = False
    else:
        cut = np.zeros(n, dtype=int)
        q, r = divmod(n, m)
        indices = [q * i + min(i, r) for i in range(m)]
        cut[indices] = True

    return cut
