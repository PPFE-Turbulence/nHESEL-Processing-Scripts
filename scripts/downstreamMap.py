# Imports
import numpy as np

### Packaging of downstream map
def downstreamPosition(x):
    """
    Function used to map an upstream variable (x) to a downstream variable (s) on the diverter.

    This is done by doing a linear interpolation on (x,s) data obtained through mapping flux lines from the midplane (upstream) to their corresponding positions on the diverter (downstream).

    The valid range for x is [0, 09327518], and therefore not all the way out to the wall. This is due to a zero-correction.
    """
    
    # First define the data obtained through suffering
    upstreamPositionsFiltered = np.array([-0.01606274, -0.00912414, -0.00246905,  0.0037867 ,  0.00964887,
        0.01512325,  0.0202156 ,  0.02493166,  0.02927706,  0.03325733,
        0.03687789,  0.04014398,  0.04306067,  0.04563285,  0.0478652 ,
        0.04976216,  0.05118268,  0.05303391,  0.05532461,  0.05807933,
        0.06090076,  0.06343153,  0.06569142,  0.06791227,  0.07020983,
        0.07241228,  0.07430134,  0.07587983,  0.07726298,  0.07852373,
        0.07968339,  0.08074674,  0.08171469,  0.08258883,  0.08337003,
        0.08406285,  0.08468587,  0.08529435,  0.08593923,  0.08651237,
        0.08677069,  0.0867216 ,  0.08669119,  0.08700075,  0.08776468,
        0.08877431,  0.08976889,  0.09098582,  0.09289251,  0.09523871,
        0.09722371,  0.09863326,  0.09975609,  0.1004806 ,  0.10068002])
    downstreamPositions = np.array([-0.05266045, -0.04076938, -0.02038469,  0.        ,  0.02038469,
        0.04076938,  0.06115407,  0.08153875,  0.10192344,  0.12230813,
        0.14269282,  0.16307751,  0.1834622 ,  0.20384689,  0.22423157,
        0.24461626,  0.26500095,  0.28538564,  0.30577033,  0.32936976,
        0.34126083,  0.36164551,  0.37353658,  0.39392127,  0.4175207 ,
        0.42941177,  0.44979646,  0.46168752,  0.48207221,  0.49396328,
        0.51434797,  0.52623904,  0.54662373,  0.5585148 ,  0.57889948,
        0.59079055,  0.61117524,  0.62306631,  0.643451  ,  0.65534207,
        0.67894149,  0.69932618,  0.71121725,  0.73160194,  0.74349301,
        0.76709244,  0.7789835 ,  0.79936819,  0.81125926,  0.82315033,
        0.8350414 ,  0.85542609,  0.86731716,  0.87920822,  0.89959291])

    # Then check the input
    if x > upstreamPositionsFiltered[-1]:
        raise Exception(f"x to map is currently {x}, but the mapping only goes up to {upstreamPositionsFiltered[-1]}")

    # Then figure out which indices it lies between on the upstream values
    for i in range(len(upstreamPositionsFiltered)):
        pos = upstreamPositionsFiltered[i]
        if pos >= x:
            # Handling the first position
            if i == 0:
                return downstreamPositions[0]
            # Handling the last position
            elif i == (len(upstreamPositionsFiltered)-1):
                return downstreamPositions[-1]
            # And the interpolation of everything in-between
            elif i != 0:
                lowerIndex = i-1
                upperIndex = i
                lowerPos = upstreamPositionsFiltered[lowerIndex]
                upperPos = upstreamPositionsFiltered[upperIndex]
                lowerPosS = downstreamPositions[lowerIndex]
                upperPosS = downstreamPositions[upperIndex]
            break
        else:
            # Handling the last position
            if i == (len(upstreamPositionsFiltered)-1):
                return downstreamPositions[-1]
            else:
                pass

    # Fuck it wikipedia: https://en.wikipedia.org/wiki/Linear_interpolation
    downstreamPos = lowerPosS + (x - lowerPos)*((upperPosS - lowerPosS)/(upperPos - lowerPos))

    return downstreamPos

def downstreamMap(xPositions):
    """
    Function used to map an iterable (list, array, etc.) of upstream variables (x) to downstream variables (s) on the diverter.

    Returns a numpy array of mapped downstream variables.
    """
    
    # Do a simply list comprehension, calling downstreamPosition for each value
    mappedPositions = np.array([downstreamPosition(x) for x in xPositions])

    return mappedPositions

