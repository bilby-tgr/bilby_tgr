import bilby.gw.prior as prior
import os
from bilby.core.utils import logger

DEFAULT_PRIOR_DIR = os.path.join(os.path.dirname(__file__), 'prior_files')


class BBHPriorDict(prior.BBHPriorDict):
    def __init__(self, dictionary=None, filename=None, aligned_spin=False,
                 conversion_function=None, v4_version=False):
        """ Initialises a Prior set for Binary Black holes with pSEOB parameter.
        Adapted from bilby.

        Parameters
        ==========
        dictionary: dict, optional
            See superclass
        filename: str, optional
            See superclass
        conversion_function: func
            Function to convert between sampled parameters and constraints.
            By default this generates many additional parameters, see
            BBHPriorDict.default_conversion_function
        """
        if dictionary is None and filename is None:
            try:
                prior_files = {
                    (True, False): ('aligned_spins_bbh.prior', 'Using aligned-spin prior for pSEOBNRv5HM'),
                    (True, True): ('aligned_spins_bbh_v4.prior', 'Using aligned-spin prior for pSEOBNRv4HM_PA'),
                    (False, False): ('precessing_spins_bbh.prior', 'Using precessing-spin prior for pSEOBNRv5PHM'),
                }
                fname, log_message = prior_files[(aligned_spin, v4_version)]
                logger.info(log_message)

            except KeyError:
                raise NotImplementedError("A precessing-spin prior is not supported for pSEOBNRv4HM_PA")

            filename = os.path.join(DEFAULT_PRIOR_DIR, fname)
            logger.info(f'No prior given, using default BBH priors in {filename}')
        elif filename is not None:
            if not os.path.isfile(filename):
                filename = os.path.join(DEFAULT_PRIOR_DIR, filename)
        super(prior.BBHPriorDict, self).__init__(dictionary=dictionary, filename=filename,
                                                 conversion_function=conversion_function)


# (Adapted from bilby)
prior.Prior._default_latex_labels = {
    'mass_1': '$m_1$',
    'mass_2': '$m_2$',
    'total_mass': '$M$',
    'chirp_mass': r'$\mathcal{M}$',
    'mass_ratio': '$q$',
    'symmetric_mass_ratio': r'$\eta$',
    'a_1': '$a_1$',
    'a_2': '$a_2$',
    'tilt_1': r'$\theta_1$',
    'tilt_2': r'$\theta_2$',
    'cos_tilt_1': r'$\cos\theta_1$',
    'cos_tilt_2': r'$\cos\theta_2$',
    'phi_12': r'$\Delta\phi$',
    'phi_jl': r'$\phi_{JL}$',
    'luminosity_distance': '$d_L$',
    'dec': r'$\mathrm{DEC}$',
    'ra': r'$\mathrm{RA}$',
    'iota': r'$\iota$',
    'cos_iota': r'$\cos\iota$',
    'theta_jn': r'$\theta_{JN}$',
    'cos_theta_jn': r'$\cos\theta_{JN}$',
    'psi': r'$\psi$',
    'phase': r'$\phi$',
    'geocent_time': '$t_c$',
    'time_jitter': '$t_j$',
    'lambda_1': r'$\Lambda_1$',
    'lambda_2': r'$\Lambda_2$',
    'lambda_tilde': r'$\tilde{\Lambda}$',
    'delta_lambda_tilde': r'$\delta\tilde{\Lambda}$',
    'chi_1': r'$\chi_1$',
    'chi_2': r'$\chi_2$',
    'chi_1_in_plane': r'$\chi_{1, \perp}$',
    'chi_2_in_plane': r'$\chi_{2, \perp}$',
    'dchi_Minus2': r'$\delta\varphi_{-2}$',
    'domega220': r'$\delta f_{220}$',
    'domega210': r'$\delta f_{210}$',
    'domega330': r'$\delta f_{330}$',
    'domega320': r'$\delta f_{320}$',
    'domega440': r'$\delta f_{440}$',
    'domega430': r'$\delta f_{430}$',
    'domega550': r'$\delta f_{550}$',
    'dtau220': r'$\delta \tau_{220}$',
    'dtau210': r'$\delta \tau_{210}$',
    'dtau330': r'$\delta \tau_{330}$',
    'dtau320': r'$\delta \tau_{320}$',
    'dtau440': r'$\delta \tau_{440}$',
    'dtau430': r'$\delta \tau_{430}$',
    'dtau550': r'$\delta \tau_{550}$',
    'dA22': r'$\delta A_{22}$',
    'dA21': r'$\delta A_{21}$',
    'dA33': r'$\delta A_{33}$',
    'dA32': r'$\delta A_{32}$',
    'dA44': r'$\delta A_{44}$',
    'dA43': r'$\delta A_{43}$',
    'dA55': r'$\delta A_{55}$',
    'dw22': r'$\delta \omega_{22}$',
    'dw21': r'$\delta \omega_{21}$',
    'dw33': r'$\delta \omega_{33}$',
    'dw32': r'$\delta \omega_{32}$',
    'dw44': r'$\delta \omega_{44}$',
    'dw43': r'$\delta \omega_{43}$',
    'dw55': r'$\delta \omega_{55}$',
    'da6': r'$\delta a_{6}$',
    'ddSO': r'$\delta d_{\mathrm{SO}}$',
    'dTpeak': r'$\delta \Delta_t_{\mathrm{peak}}$',
}
