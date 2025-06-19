# FILE: deriv2D.py

# Import packages
import numpy as np

def deriv2D(ima, 
            x=False, 
            y=False):
    """
    Compute 1D derivative of a 2D image array along x or y direction.

    Parameters
    ----------
    ima : 2D np.ndarray
        Input image or 2D array.
    x : bool
        If True, compute derivative along x-axis (columns).
    y : bool
        If True, compute derivative along y-axis (rows).

    Returns
    -------
    result : 2D np.ndarray
        Array of same shape as ima with derivatives computed.
    """
    result = np.zeros_like(ima)

    if y:
        for i in range(ima.shape[0]):
            result[i, :] = deriv1d(ima[i, :])
    else:
        for i in range(ima.shape[1]):
            result[:, i] = deriv1d(ima[:, i])

    return result

def deriv1d(arr):
    """
    Approximate 1D derivative using central differences.
    """
    d = np.zeros_like(arr)
    d[1:-1] = (arr[2:] - arr[:-2]) / 2.0
    d[0] = arr[1] - arr[0]       # forward difference at start
    d[-1] = arr[-1] - arr[-2]    # backward difference at end
    return d