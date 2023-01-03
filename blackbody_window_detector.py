from radis import calc_spectrum, Spectrum, SerialSlabs
from radis.spectrum.operations import add_array, multiply

from specutils import SpectralRegion
from specutils.manipulation import noise_region_uncertainty
from specutils.fitting import find_lines_threshold, find_lines_derivative

import astropy.units as u
import numpy as np

# -------------------------------------
# ------------- blackbody -------------
# -------------------------------------
def __sPlanck(spectrum, temp):
    """
    Returns a function that generates a Blackbody Spectrum.

            Parameters:
                data (int): the range of x-values for the spectrum

            Returns:
                The Blackbody Spectrum
    """

    H = 6.62606957e-34
    C = 2.99792458e8
    K_B = 1.3806488e-23

    return ((0.2 * H * (C**2)) / (((spectrum * (10**-9)) ** 4) * spectrum)) * (
        1 / (np.exp((H * C) / ((spectrum * (10**-9)) * K_B * temp)) - 1)
    )


# --------------------------------------
# --------------- window ---------------
# --------------------------------------
def __CaF2(data):
    """
    Returns a function that approximates a CaF2 Cell Window.

            Parameters:
                data (int): the range of x-values for the function

            Returns:
                The approximated function
    """

    return (0.93091) / (1 + (11.12929 / (data / 1000)) ** -12.43933) ** 4.32574


def __ZnSe(data):
    """
    Returns a function that approximates a ZnSe Cell Window.

            Parameters:
                data (int): the range of x-values for the function

            Returns:
                The approximated function
    """

    x_um = data / 1000
    return (0.71015) / ((1 + (20.99353 / x_um) ** -19.31355) ** 1.44348) + -0.13265 / (
        2.25051 * np.sqrt(np.pi / (4 * np.log(2)))
    ) * np.exp(-4 * np.log(2) * ((x_um - 16.75) ** 2) / (2.25051**2))


def __sapphire(data):
    """
    Returns a function that approximates a Sapphire Window.

            Parameters:
                data (int): the range of x-values for the function

            Returns:
                The approximated function
    """

    return 0.78928 / (1 + (11.9544 / (data / 1000)) ** -12.07226) ** 6903.57039


def __AR_ZnSe(data):
    """
    Returns a function that approximates a AR_ZnSe beamsplitter.

            Parameters:
                data (int): the range of x-values for the function

            Returns:
                The approximated function
    """

    x_um = data / 1000
    return (
        (0.82609) / ((1 + ((34.63971 / x_um) ** -8.56269)) ** 186.34792)
        + -0.47
        / (0.55 * np.sqrt(np.pi / (4 * np.log(2))))
        * np.exp(-4 * np.log(2) * ((x_um - 1.47) ** 2) / (0.55**2))
        + -0.03456
        / (0.4 * np.sqrt(np.pi / (4 * np.log(2))))
        * np.exp(-4 * np.log(2) * ((x_um - 2.88) ** 2) / (0.4**2))
        + -0.009
        / (0.3 * np.sqrt(np.pi / (4 * np.log(2))))
        * np.exp(-4 * np.log(2) * ((x_um - 6.16) ** 2) / (0.3**2))
        + -0.09
        / (1 * np.sqrt(np.pi / (4 * np.log(2))))
        * np.exp(-4 * np.log(2) * ((x_um - 16.2) ** 2) / (1**2))
        + -0.08
        / (1 * np.sqrt(np.pi / (4 * np.log(2))))
        * np.exp(-4 * np.log(2) * ((x_um - 17.4) ** 2) / (1**2))
        + 1.12
        / (8 * np.sqrt(np.pi / (4 * np.log(2))))
        * np.exp(-4 * np.log(2) * ((x_um - 9.5) ** 2) / (8**2))
        + 0.11546
        / (2 * np.sqrt(np.pi / (4 * np.log(2))))
        * np.exp(-4 * np.log(2) * ((x_um - 4.9) ** 2) / (2**2))
        + 0.21751
        / (2 * np.sqrt(np.pi / (4 * np.log(2))))
        * np.exp(-4 * np.log(2) * ((x_um - 2.6) ** 2) / (2**2))
        + -0.05
        / (0.07 * np.sqrt(np.pi / (4 * np.log(2))))
        * np.exp(-4 * np.log(2) * ((x_um - 0.8) ** 2) / (0.07**2))
    )


def __AR_CaF2(data):
    """
    Returns a function that approximates a AR_CaF2 beamsplitter.

            Parameters:
                data (int): the range of x-values for the function

            Returns:
                The approximated function
    """

    x_um = data / 1000
    return (
        (0.9795) / ((1 + ((18.77617 / x_um) ** -6.94246)) ** 91.98745)
        + -0.06
        / (0.08 * np.sqrt(np.pi / (4 * np.log(2))))
        * np.exp(-4 * np.log(2) * ((x_um - 0.76) ** 2) / (0.08**2))
        + -0.06
        / (0.2 * np.sqrt(np.pi / (4 * np.log(2))))
        * np.exp(-4 * np.log(2) * (x_um - 1.06) ** 2 / 0.20**2)
        + -0.6
        / (3.0 * np.sqrt(np.pi / (4 * np.log(2))))
        * np.exp(-4 * np.log(2) * ((x_um - 4.85) ** 2) / (3.0**2))
        + -0.35
        / (1.0 * np.sqrt(np.pi / (4 * np.log(2))))
        * np.exp(-4 * np.log(2) * ((x_um - 9.40) ** 2) / (1.00**2))
        + 0.05
        / (0.8 * np.sqrt(np.pi / (4 * np.log(2))))
        * np.exp(-4 * np.log(2) * ((x_um - 2.60) ** 2) / (0.8**2))
        + 0.04
        / (0.5 * np.sqrt(np.pi / (4 * np.log(2))))
        * np.exp(-4 * np.log(2) * ((x_um - 7.75) ** 2) / (0.50**2))
        + -0.01
        / (0.6 * np.sqrt(np.pi / (4 * np.log(2))))
        * np.exp(-4 * np.log(2) * ((x_um - 6.55) ** 2) / (0.6**2))
        + 0.01
        / (0.5 * np.sqrt(np.pi / (4 * np.log(2))))
        * np.exp(-4 * np.log(2) * ((x_um - 1.82) ** 2) / (0.5**2))
    )


# --------------------------------------
# -------------- detector --------------
# --------------------------------------
def __InSb(data):
    """
    Returns a function that approximates a InSb dectector.

            Parameters:
                data (int): the range of x-values for the function

            Returns:
                The approximated function
    """

    x_um = data / 1000
    return 1.97163e11 * (1 / (1 + np.exp(-(x_um - 5.3939) / 1.6624))) * (
        1 - 1 / (1 + np.exp(-(x_um - 5.3939) / 0.11925))
    ) + (3.3e10) / (2.44977 * np.sqrt(np.pi / (4 * np.log(2)))) * np.exp(
        -4 * np.log(2) * ((x_um - 5) ** 2) / (2.44977**2)
    )


def __MCT(data):
    """
    Returns a function that approximates a MCT dectector.

            Parameters:
                data (int): the range of x-values for the function

            Returns:
                The approximated function
    """

    x_um = data / 1000
    return (
        (1.98748 * (10**9))
        + (2.10252 * (10**10))
        * (1 / (1 + np.exp(-(x_um - 20.15819) / 5.73688)))
        * (1 - 1 / (1 + np.exp(-(x_um - 20.15819) / 1.11659)))
        + (1.3 * (10**9))
        / (2 * np.sqrt(np.pi / (4 * np.log(2))))
        * np.exp(-4 * np.log(2) * ((x_um - 18.6) ** 2) / (2**2))
    )
