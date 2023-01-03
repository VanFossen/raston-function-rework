# import matplotlib.pyplot as plt
# import numpy as np


def __param_check(data):
    """
    Checks user given parameters and makes sure they are valid.

        Parameters:
            data (dict): the user given parameters

        Returns:
            True if params are good. Else, returns False
    """

    # check if number of parameters is correct
    if len(data) != 11:
        print("  not enough params. total params: %s" % (len(data)))
        return False

    # check if parameter names are correct
    valid_params = [
        "minWave",
        "maxWave",
        "molecule",
        "pressure",
        "resolution",
        "numScan",
        "zeroFill",
        "source",
        "beamsplitter",
        "cellWindow",
        "detector",
    ]

    for key, value in data.items():
        if (key in valid_params) and (data[key] is not None):
            print(f"  {key}: {value}")
        else:
            print(f"  error with key: {key}. Value is: {value}")
            return False

    return True


def __calc_wstep(resolution, zero_fill):
    """
    Calculates the apropriate wstep for a spectrum based on the given resolution and zero fill.

        Parameters:
            resolution (int): the given resolution
            xero_fill (int): the given zero fill

        Returns:
            The calculated wstep
    """

    match resolution:
        case 1:
            if zero_fill == 0:
                wstep = 0.481927711
            elif zero_fill == 1:
                wstep = 0.240963855
            elif zero_fill == 2:
                wstep = 0.120481928

        case 0.5:
            if zero_fill == 0:
                wstep = 0.240963855
            elif zero_fill == 1:
                wstep = 0.120481928
            elif zero_fill == 2:
                wstep = 0.060240964

        case 0.25:
            if zero_fill == 0:
                wstep = 0.120481928
            elif zero_fill == 1:
                wstep = 0.060240964
            elif zero_fill == 2:
                wstep = 0.030120482

        case 0.125:
            if zero_fill == 0:
                wstep = 0.060240964
            elif zero_fill == 1:
                wstep = 0.030120482
            elif zero_fill == 2:
                wstep = 0.015060241

        case 0.0625:
            if zero_fill == 0:
                wstep = 0.030120482
            elif zero_fill == 1:
                wstep = 0.015060241
            elif zero_fill == 2:
                wstep = 0.00753012

    return wstep


# def __loadData(s):
#     """
#     Takes a spectrum and converts its data to a dictionary of x and y values.

#             Parameters:
#                 s (Spectrum): The spectrum to convert

#             Returns:
#                 A dictionary of x (key) and y (value) values
#     """

#     data = {}

#     for key, val in zip(s[0], s[1]):
#         data[float(key)] = float(val)

#     return data


# def __graph(spectrum):
#     xs = []
#     ys = []
#     for key in spectrum:
#         xs.append(float(key))
#         ys.append(float(spectrum[key]))
#     plt.plot(np.array(xs), np.array(ys), "blue")
#     plt.show()


# def __zeroY(data):
#     """
#     Returns a function of y = 1 for generating background samples.

#             Parameters:
#                 data (int): the range of x-values for the function

#             Returns:
#                 The approximated function
#     """
#     return (data * 0) + 1
