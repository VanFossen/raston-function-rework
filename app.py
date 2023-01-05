#%%
import json

import matplotlib.pyplot as plt
import numpy as np
from radis import SerialSlabs, Spectrum, calc_spectrum
from radis.spectrum.operations import add_array

from blackbody_window_detector import (
    __MCT,
    __AR_CaF2,
    __AR_ZnSe,
    __CaF2,
    __InSb,
    __sapphire,
    __sPlanck,
    __ZnSe,
)
from helper import __calc_wstep, __param_check

DEBUG = False


def __log(string, spectrum=False, title=False):
    if DEBUG:
        print(string)
        if spectrum is not False:
            spectrum.plot("transmittance_noslit")
            plt.title(title)


def main():
    # read local data file into a dictionary
    __log("reading 'data.json'...")
    with open("./data.json", "r") as data_json:
        data = json.load(data_json)

    # verify user input is valid
    __log("verifying user parameters...")
    if not __param_check(data):
        print("parameter check failed")
        return

    # resolution of wavenumber grid (cm^-1)
    #   https://radis.readthedocs.io/en/latest/source/radis.lbl.calc.html#radis.lbl.calc.calc_spectrum:~:text=wstep%20(float%20(,%27auto%27)
    __log("calculating wstep...")
    wstep = __calc_wstep(data["resolution"], data["zeroFill"])

    try:
        # ----- a.) transmission spectrum of gas sample -----
        #   https://radis.readthedocs.io/en/latest/source/radis.lbl.calc.html#radis.lbl.calc.calc_spectrum
        __log("performing 'calc_spectrum()'...")
        s = calc_spectrum(
            data["minWave"],
            data["maxWave"],
            molecule=data["molecule"],
            isotope="1,2,3",
            pressure=data["pressure"],
            Tgas=294.15,
            path_length=10,
            wstep=wstep,
            databank="hitran",
            verbose=False,
            warnings={"AccuracyError": "ignore"},
        )
        __log("plotting result of 'calc_spectrum'...", s, "raw spectrum")
    except Exception as e:
        print(f"ERROR: {e}")
        return False

    # ----- Pre-processing -----
    # Generate the necessary spectra for each component of the following processing steps. The spectra are generated based on the function provided in the call to the Spectrum constructor

    # return wavelength in defined medium.
    #   https://radis.readthedocs.io/en/latest/source/radis.spectrum.spectrum.html#radis.spectrum.spectrum.Spectrum.get_wavelength
    __log("performing 'get_wavenumber()'...")
    w = s.get_wavenumber()

    # NOTE I believe the following can be simplified into a new function or an existing in `blackbody_window_detector.py`
    # -------------------------------------------------------------------------------------------------------------------
    __log("performing 'Spectrum()' for sPlanck...")
    spec_sPlanck = Spectrum(
        {"wavenumber": w, "transmittance_noslit": __sPlanck(w, data["source"])},
        wunit="cm",
        units={"transmittance_noslit": ""},
        name="sPlanck",
    )
    __log("plotting result of 'spec_sPlanck'...", spec_sPlanck, "raw sPlanck")

    __log("performing 'Spectrum()' for AR_ZnSe...")
    spec_AR_ZnSe = Spectrum(
        {"wavenumber": w, "transmittance_noslit": __AR_ZnSe(w)},
        wunit="cm",
        units={"transmittance_noslit": ""},
        name="AR_ZnSe",
    )
    __log("plotting result of 'spec_AR_ZnSe'...", spec_AR_ZnSe, "raw AR_ZnSe")

    __log("performing 'Spectrum()' for AR_CaF2...")
    spec_AR_CaF2 = Spectrum(
        {"wavenumber": w, "transmittance_noslit": __AR_CaF2(w)},
        wunit="cm",
        units={"transmittance_noslit": ""},
        name="AR_CaF2",
    )
    __log("plotting result of 'spec_AR_CaF2'...", spec_AR_CaF2, "raw AR_CaF2")

    __log("performing 'Spectrum()' for CaF2...")
    spec_CaF2 = Spectrum(
        {"wavenumber": w, "transmittance_noslit": __CaF2(w)},
        wunit="cm",
        units={"transmittance_noslit": ""},
        name="CaF2",
    )
    __log("plotting result of 'spec_CaF2'...", spec_CaF2, "raw CaF2")

    __log("performing 'Spectrum()' for ZnSe...")
    spec_ZnSe = Spectrum(
        {"wavenumber": w, "transmittance_noslit": __ZnSe(w)},
        wunit="cm",
        units={"transmittance_noslit": ""},
        name="ZnSe",
    )
    __log("plotting result of 'spec_ZnSe'...", spec_ZnSe, "raw ZnSe")

    __log("performing 'Spectrum()' for sapphire...")
    spec_sapphire = Spectrum(
        {"wavenumber": w, "transmittance_noslit": __sapphire(w)},
        wunit="cm",
        units={"transmittance_noslit": ""},
        name="sapphire",
    )
    __log("plotting result of 'spec_sapphire'...", spec_sapphire, "raw sapphire")

    __log("performing 'Spectrum()' for MCT...")
    spec_MCT = Spectrum(
        {"wavenumber": w, "transmittance_noslit": __MCT(w)},
        wunit="cm",
        units={"transmittance_noslit": ""},
        name="MCT",
    )
    __log("plotting result of 'spec_MCT'...", spec_MCT, "raw MCT")

    __log("performing 'Spectrum()' for InSb...")
    spec_InSb = Spectrum(
        {"wavenumber": w, "transmittance_noslit": __InSb(w)},
        wunit="cm",
        units={"transmittance_noslit": ""},
        name="InSb",
    )
    __log("plotting result of 'spec_InSb'...", spec_InSb, "raw InSb")

    # ----- b.) blackbody spectrum of source -----
    # SerialSlabs() multiplies the transmittance values (y-values) of the provided spectrum, s, with the corresponding transmittance values of the spec_sPlanck spectrum
    #   https://radis.readthedocs.io/en/latest/source/radis.los.slabs.html#radis.los.slabs.SerialSlabs
    __log("combining spectrum with blackbody ('sPlanck()')...")
    spectrum = SerialSlabs(s, spec_sPlanck)
    __log("plotting new combined spectrum...", spectrum, "raw + blackbody")

    # This loop simulates scans and runs as many times as the user indicated in "Number of Scans"
    for x in range(data["numScan"]):
        # ----- c.) transmission spectrum of windows/beamsplitter -----
        # ----- c.1) Beamsplitter -----
        __log("combining spectrum with beamsplitter ('AR_ZnSe()' or 'AR_CaF2()')...")
        match data["beamsplitter"]:
            case "AR_ZnSe":
                spectrum = SerialSlabs(spectrum, spec_AR_ZnSe)
            case "AR_CaF2":
                spectrum = SerialSlabs(spectrum, spec_AR_CaF2)
        __log(
            "plotting beamsplitter result...",
            spectrum,
            "raw + blackbody + beamsplitter",
        )

        # ----- c.2) cell windows -----
        __log("combining spectrum with cell windows ('CaF2()' or 'ZnSe()')...")
        match data["cellWindow"]:
            case "CaF2":
                spectrum = SerialSlabs(spectrum, spec_CaF2)
                spectrum = SerialSlabs(spectrum, spec_CaF2)
            case "ZnSe":
                spectrum = SerialSlabs(spectrum, spec_ZnSe)
                spectrum = SerialSlabs(spectrum, spec_ZnSe)
        __log(
            "plotting cell windows result...",
            spectrum,
            "raw + blackbody + beamsplitter + cell window",
        )

        # ----- d.) detector response spectrum -----
        __log(
            "combining spectrum with detector ('ZnSe()' & 'MCT()' or 'sapphire()' & 'InSb()')..."
        )
        match data["detector"]:
            case "MCT":
                spectrum = SerialSlabs(spectrum, spec_ZnSe)
                spectrum = SerialSlabs(spectrum, spec_MCT)
            case "InSb":
                spectrum = SerialSlabs(spectrum, spec_sapphire)
                spectrum = SerialSlabs(spectrum, spec_InSb)
        __log(
            "plotting detector result...",
            spectrum,
            "raw + blackbody + beamsplitter + cell window + detector",
        )

        # add random noise to spectrum
        __log("adding noise to spectrum...")
        spectrum = add_array(
            spectrum,
            np.random.normal(0, 200000000, len(spectrum.get_wavenumber())),
            var="transmittance_noslit",
        )
        __log(
            "plotting noise result...",
            spectrum,
            "raw + blackbody + beamsplitter + cell window + detector + noise",
        )

    # NOTE additional code that I am unsure if needed or needs to be rewritten
    # ------------------------------------------------------------------------
    # ------------------------------------------------------------------------
    # # Normalize data
    # # numbers = __loadData(
    # #     spectrum.get("transmittance_noslit", wunit="nm", Iunit="default")
    # # )
    # # factor = 1 / sum(numbers.values())

    # # spectrum = multiply(spectrum, factor, var="transmittance_noslit")

    # # Post-processing - Find Peaks
    # # Not done on background samples
    # # https://radis.readthedocs.io/en/latest/auto_examples/plot_specutils_processing.html#sphx-glr-auto-examples-plot-specutils-processing-py
    # if find_peaks:
    #     find_peaks = spectrum.to_specutils()
    #     noise_region = SpectralRegion(
    #         (1 / data["minWave"]) / u.cm, (1 / data["maxWave"]) / u.cm
    #     )
    #     find_peaks = noise_region_uncertainty(find_peaks, noise_region)
    #     lines = find_lines_threshold(find_peaks, noise_factor=6)
    #     print()
    #     print(lines)


if __name__ == "__main__":
    main()

# %%
