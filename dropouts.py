"""
Dropout selection
It applies a dropout redshift selection for galaxies from the BoRG
survey at z~10.
Mangitudes for this case:
    # mag f160w = mag[0][]
    # mag f098m = mag[1][]
    # mag f125w = mag[2][]
    # mag f606w = mag[3][]

This file can be modified by the user.
"""

import numpy as np

def main(mag_iso, mag_auto, sn, st):
    """
    Dropout selection for galaxies from BoRG at z~10
    Args:
        mag_iso (float array) = ISO magnitudes for each band and object.
        mag_auto (float array) = AUTO magnitudes for each band and object.
        sn (float array) = Signal to noise for each band and object.
        st  (int) = Indicates the identification status of the object
                (non detected, blended, detected).
    Returns:
        drops (int) = Number of artificial sources that pass the redshift
                      selection criteria.
        drops_array (int array) = array with the dropout status, 1 for objects
                                  that pass the criteria and 0 for objects that
                                  don't.
    """

    drops = 0  # Initialise dropouts count.
    drops_array = np.zeros(sn.shape[0])  # Initialise dropouts status array.

    for i in xrange(sn.shape[0]):
        if (st[i] >= 0 and (mag_iso[i, 2] - mag_iso[i, 0]) > 1.5 and
           sn[i, 3] < 1.5 and sn[i, 1] < 1.5):
            drops = drops + 1
            drops_array[i] = 1
    return drops, drops_array

if __name__ == "__main__":
    main()
