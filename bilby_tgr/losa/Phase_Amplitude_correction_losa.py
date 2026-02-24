import lal
import numpy as np

LAL_PI = lal.PI
MTSUN_SI = lal.MTSUN_SI
MPC_SI = lal.PC_SI * 10**6
c_SI = lal.C_SI
EC = 0.577216  # Euler constant


# Point particle contributions to the phase

def pp_contr_2PN_to_LO(eta):
    return 743 / 126 + 22 / 3 * eta


def pp_contr_3PN_to_LO(eta):
    return - 64 * LAL_PI / 5


def pp_contr_4PN_to_LO(eta):
    return (1755623 / 84672 + 32633 / 756 * eta + 367 / 12 * eta * eta)


def pp_contr_5PN_to_LO(eta):
    return (- 20807 / 210 - 574 / 15 * eta) * LAL_PI


def pp_contr_6PN_to_LO(eta):
    return (-28907482848623 / 35206617600 + 9472 / 75 * LAL_PI * LAL_PI + 13696 / 105 * EC
            + (3311653861 / 1524096 - 451 / 6 * LAL_PI * LAL_PI) *
            eta + 2030687 / 18144 * eta * eta
            + 66287 / 648 * eta * eta * eta)


def pp_contr_6PNLOG_to_LO(_):
    return 13696 / 105


def pp_contr_7PN_to_LO(eta):
    return (- 158992529 / 317520 - 1015907 / 1890 *
            eta + 419 / 945 * eta * eta) * LAL_PI


# Tidal correction contributions to the phase

def tidal_contr_10PN_to_LO(mbyM):
    return 576 - 528 * mbyM


def tidal_contr_12PN_to_LO(mbyM):
    return (46029 / 14 + 4021 / 42 * mbyM - 5191 * mbyM * mbyM + 2066 * mbyM * mbyM * mbyM)


def tidal_contr_14PN_to_LO(mbyM):
    return (226452487 / 15876 + 48108089 / 7056 * mbyM - 576641 / 28 * mbyM * mbyM - 3358153 / 252 * mbyM * mbyM * mbyM
            + 56924 / 3 * mbyM * mbyM * mbyM * mbyM - 15220 / 3 * mbyM * mbyM * mbyM * mbyM * mbyM)


def tidal_contr_15PN_to_LO(mbyM):
    return (- 172581 / 5 + 838063 / 35 * mbyM + 55458 / 7 *
            mbyM * mbyM - 4856 / 35 * mbyM * mbyM * mbyM) * LAL_PI


# Aligned Spin Contributions to the phase

def as_contr_3PN_to_LO(eta, delta, chi_s, chi_a):
    return 452 / 15 * delta * chi_a + (452 / 15 - 304 / 15 * eta) * chi_s


def as_contr_4PN_to_LO(eta, delta, chi_s, chi_a):
    return ((80 * chi_a * chi_a + chi_s * chi_s) * eta - 81 / 4 *
            (chi_a * chi_a + chi_s * chi_s + 2 * delta * chi_a * chi_s))


def as_contr_5PN_to_LO(eta, delta, chi_s, chi_a):
    return ((300301 / 630 + 5812 / 45 * eta) * delta * chi_a +
            (300301 / 630 - 97232 / 315 * eta - 5384 / 45 * eta * eta) * chi_s)


def as_contr_6PN_to_LO(eta, delta, chi_s, chi_a):
    return (- 12388 / 25 * LAL_PI * delta * chi_a - (12388 / 25 - 25328 / 75 * eta) * LAL_PI * chi_s
            + (84452 / 225 - 515783 / 450 * eta) * delta * chi_a * chi_s +
            (42226 / 225 - 2106779 / 2100 * eta + 936 * eta * eta) * chi_a * chi_a
            + (42226 / 225 - 5629937 / 6300 * eta + 60679 / 225 * eta * eta) * chi_s * chi_s)


def as_contr_7PN_to_LO(eta, delta, chi_s, chi_a):
    return ((30198206713 / 3810240 - 39553501 / 22680 * eta + 194267 / 540 * eta * eta) * delta * chi_a
            + (2928 / 5 - 2304 * eta) * LAL_PI * chi_a * chi_a -
            (4193 / 5 - 9872 / 3 * eta) * delta * chi_a * chi_a * chi_a
            + (30198206713 / 3810240 - 10028546561 / 952560 * eta - 12880727 / 11340 * eta * eta
                - 61921 / 135 * eta * eta * eta) * chi_s + 5856 / 5 * LAL_PI * delta * chi_a * chi_s
            - (12579 / 5 - 31970 / 3 * eta + 7744 / 3 * eta * eta) * chi_a * chi_a * chi_s
            + (2928 - 192 * eta) / 5 * LAL_PI * chi_s * chi_s -
            (12579 / 5 - 20672 / 15 * eta) * delta * chi_a * chi_s * chi_s
            - (4193 / 5 - 10814 / 15 * eta + 712 / 15 * eta * eta) * chi_s * chi_s * chi_s)


pp_cor_dict = {
    '2PN': pp_contr_2PN_to_LO, '3PN': pp_contr_3PN_to_LO, '4PN': pp_contr_4PN_to_LO, '5PN': pp_contr_5PN_to_LO,
    '6PN': pp_contr_6PN_to_LO, '6PNLOG': pp_contr_6PNLOG_to_LO, '7PN': pp_contr_7PN_to_LO
}


as_cor_dict = {
    '3PN': as_contr_3PN_to_LO, '4PN': as_contr_4PN_to_LO, '5PN': as_contr_5PN_to_LO,
    '6PN': as_contr_6PN_to_LO, '7PN': as_contr_7PN_to_LO
}


tidal_cor_dict = {
    '10PN': tidal_contr_10PN_to_LO, '12PN': tidal_contr_12PN_to_LO, '14PN': tidal_contr_14PN_to_LO,
    '15PN': tidal_contr_15PN_to_LO
}

# Phase Correction


def phase_correction_losa(frequency_array,
                          mass_1,
                          mass_2,
                          a_by_c,
                          a_1=0,
                          a_2=0,
                          tilt_1=0,
                          tilt_2=0,
                          lambda_1=0,
                          lambda_2=0,
                          pn_losa_order=-1
                          ):
    '''
    Phase correction due to the line-of-sight acceleration of the center of mass of the binary. Adds -4 PN
    correction term of arXiv:2302.09651 to the waveform. It also includes the aligned spin and tidal
    correction terms (mass-type quadrupolar tidal corrections).

    Parameters
    ----------
    frequency_array: array_like
        The frequencies at which we want to calculate the strain

    mass_1: float
        The mass of the heavier object in solar masses
    mass_2: float
        The mass of the lighter object in solar masses
    a_by_c: float
        The line-of-sight acceleration of the center of the binary in s^{-1} (divided by the speed of light)
    a_1: float
        Dimensionless primary spin magnitude
        Default is 0
    a_2: float
        Dimensionless secondary spin magnitude
        Default is 0
    tilt_1: float
        Primary tilt angle
        Default is 0
    tilt_2: float
        Secondary tilt angle
        Default is 0
    lambda_1: float
        Dimensionless tidal deformability of mass_1
        Default is 0
    lambda_2: float
        Dimensionless tidal deformability of mass_2
        Default is 0
    pn_losa_order: int, optional
        The Next to leading order post-newtonian correction to phase due to LOSA.
        Default is -1 (all terms up to 3.5PN order in point-particle, aligned-spin,
        and 7.5PN in tidal order are considered).

        Note: pn_losa_order should be an integer. Set it to a value 2 times the required PN order.
              For example: pn_losa_order = 3 will represent 1.5PN order.

    Returns
    -------
    array: an array of phase corrections of the size same as that of the frequency_array.
    '''

    total_mass = mass_1 + mass_2
    X1 = mass_1 / total_mass
    X2 = mass_2 / total_mass

    eta = X1 * X2
    delta = X1 - X2

    v = (LAL_PI * total_mass * MTSUN_SI * frequency_array) ** (1 / 3)

    v2 = v * v
    v5 = v2 * v2 * v
    v6 = v5 * v
    v7 = v6 * v
    v13 = v6 * v7

    v_dict = {
        '2PN': v2, '3PN': v2 * v, '4PN': v2 * v2, '5PN': v5, '6PN': v6, '6PNLOG': np.log(4 * v) * v6, '7PN': v7,
        '10PN': v5 * v5, '12PN': v6 * v6, '14PN': v7 * v7, '15PN': v13 * v2
    }

    LO_corr = (25 * a_by_c / 65536 / eta / eta) * \
        (total_mass * MTSUN_SI) / v13         # LO: leading order

    SHO = 1                                  # SHO: Sum of higher order corrections

    # Adding point particle higher order corrections

    if (pn_losa_order != -1):
        SHO += sum(pp_cor_dict[key](eta) * v_dict[key]
                   for key in pp_cor_dict.keys() if key <= f'{pn_losa_order}' + 'PN')
        SHO += pp_cor_dict['6PNLOG'](eta) * \
            v_dict['6PNLOG'] if pn_losa_order == 6 else 0

    else:
        SHO += sum(pp_cor_dict[key](eta) * v_dict[key]
                   for key in pp_cor_dict.keys())

    # Adding aligned spin contributions

    if (a_1 != 0 or a_2 != 0):
        spin_1z = a_1 * np.cos(tilt_1)
        spin_2z = a_2 * np.cos(tilt_2)
        chi_s = (spin_1z + spin_2z) / 2
        chi_a = (spin_1z - spin_2z) / 2

        if (pn_losa_order != -1):
            SHO += sum(as_cor_dict[key](eta, delta, chi_s, chi_a) * v_dict[key]
                       for key in as_cor_dict.keys() if key <= f'{pn_losa_order}' + 'PN')

        else:
            SHO += sum(as_cor_dict[key](eta, delta, chi_s, chi_a)
                       * v_dict[key] for key in as_cor_dict.keys())

    # Adding tidal corrections

    if (lambda_1 != 0 or lambda_2 != 0):
        X1_4 = X1 * X1 * X1 * X1
        X2_4 = X2 * X2 * X2 * X2

        if (pn_losa_order != -1):
            SHO += sum(
                (lambda_1 *
                 X1_4 *
                 tidal_cor_dict[key](X1) +
                 lambda_2 *
                 X2_4 *
                 tidal_cor_dict[key](X2)) *
                v_dict[key] for key in tidal_cor_dict.keys() if key <= f'{10 + pn_losa_order}' + 'PN')

        else:
            SHO += sum((lambda_1 *
                        X1_4 *
                        tidal_cor_dict[key](X1) +
                        lambda_2 *
                        X2_4 *
                        tidal_cor_dict[key](X2)) *
                       v_dict[key] for key in tidal_cor_dict.keys())

    return LO_corr * SHO


# Amplitude correction

def amplitude_correction_losa(frequency_array, mass_1, mass_2, a_by_c):
    '''
    Amplitude correction due to the line-of-sight acceleration of the center of mass of the binary. Adds -4 PN
    correction term of arXiv:2302.09651 to the waveform. This function does not take PN order of the amplitude
    because it assumes the amplitude to be newtonian i.e. 0 PN.

    Parameters
    ----------
    frequency_array: array_like
        The frequencies at which we want to calculate the strain
    mass_1: float
        The mass of the heavier object in solar masses
    mass_2: float
        The mass of the lighter object in solar masses
    a_by_c: float
        The line-of-sight acceleration of the center of the binary in s^{-1} (divided by the speed of light)
    luminosity_distance: float
        The luminosity distance in Mpc.

    Returns
    -------
    array: an array of relative amplitude corrections (A/A_lead - 1) of the size same as that of the frequency_array.
    '''

    total_mass = mass_1 + mass_2
    eta = mass_1 * mass_2 / total_mass**2

    v = (LAL_PI * total_mass * MTSUN_SI * frequency_array) ** (1 / 3)

    v2 = v * v
    v4 = v2 * v2

    cor = 65 / 2048 * a_by_c * (total_mass * MTSUN_SI) / eta / (v4 * v4)

    return cor
