import bilby.gw.prior as prior
import os
from bilby.core.utils import logger


DEFAULT_PRIOR_DIR = os.path.join(os.path.dirname(__file__), 'prior_files')


class BBHPriorDict(prior.CBCPriorDict):
    def __init__(self, dictionary=None, filename=None, aligned_spin=True,
                 conversion_function=None):
        """
        (Adapted from bilby)
        Initialises a Prior set for Binary Black holes with FTI parameters.

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
            if aligned_spin:
                fname = 'aligned_spins_bbh.prior'
                logger.info('Using aligned spin prior')
            else:
                fname = 'precessing_spins_bbh.prior'
            filename = os.path.join(DEFAULT_PRIOR_DIR, fname)
            logger.info('No prior given, using default BBH priors in {}.'.format(filename))
        elif filename is not None:
            if not os.path.isfile(filename):
                filename = os.path.join(DEFAULT_PRIOR_DIR, filename)
        super(BBHPriorDict, self).__init__(dictionary=dictionary, filename=filename,
                                           conversion_function=conversion_function)

    def test_redundancy(self, key, disable_logging=False):
        """
        (Adapted from bilby)
        Test whether adding the key would be redundant.
        Already existing keys return True.

        Parameters
        ==========
        key: str
            The key to test.
        disable_logging: bool, optional
            Disable logging in this function call. Default is False.

        Returns
        =======
        redundant: bool
            Whether the key is redundant or not
        """
        if key in self:
            logger.debug('{} already in prior'.format(key))
            return True

        sampling_parameters = {key for key in self if not isinstance(
            self[key], prior.Constraint)}

        mass_parameters = {'mass_1', 'mass_2', 'chirp_mass', 'total_mass', 'mass_ratio', 'symmetric_mass_ratio'}
        spin_tilt_1_parameters = {'tilt_1', 'cos_tilt_1'}
        spin_tilt_2_parameters = {'tilt_2', 'cos_tilt_2'}
        spin_azimuth_parameters = {'phi_1', 'phi_2', 'phi_12', 'phi_jl'}
        inclination_parameters = {'theta_jn', 'cos_theta_jn'}
        distance_parameters = {'luminosity_distance', 'comoving_distance', 'redshift'}
        siqm_parameters = {'dchi_kappaS', 'dchi_kappaA', 'dchi_kappa1', 'dchi_kappa2'}

        if key in {'dchi_kappaS', 'dchi_kappaA'}:
            if 'dchi_kappa1' in sampling_parameters or 'dchi_kappa2' in sampling_parameters:
                logger.disabled = False
                logger.warning("Incompatible dkappas provided. This may lead to unexpected behaviour. \
                                Provide (dchi_kappaS and dchi_kappaA) xor (dchi_kappa1 and dchi_kappa2) instead.")
                return True
        if key in {'dchi_kappa1', 'dchi_kappa2'}:
            if 'dchi_kappaS' in sampling_parameters or 'dchi_kappaA' in sampling_parameters:
                logger.disabled = False
                logger.warning("Incompatible dkappas provided. This may lead to unexpected behaviour. \
                                Provide (dchi_kappaS and dchi_kappaA) xor (dchi_kappa1 and dchi_kappa2) instead.")
                return True

        for independent_parameters, parameter_set in \
                zip([2, 2, 1, 1, 1, 1, 2],
                    [mass_parameters, spin_azimuth_parameters,
                     spin_tilt_1_parameters, spin_tilt_2_parameters,
                     inclination_parameters, distance_parameters,
                     siqm_parameters]):
            if key in parameter_set:
                if len(parameter_set.intersection(
                        sampling_parameters)) >= independent_parameters:
                    logger.disabled = disable_logging
                    logger.warning('{} already in prior. '
                                   'This may lead to unexpected behaviour.'
                                   .format(parameter_set.intersection(self)))
                    logger.disabled = False
                    return True

        return False


class BNSPriorDict(prior.CBCPriorDict):
    def __init__(self, dictionary=None, filename=None, aligned_spin=True,
                 conversion_function=None):
        """
        (Adapted from bilby)
        Initialises a Prior set for Binary Neutron Stars with FTI parameters.

        Parameters
        ==========
        dictionary: dict, optional
            See superclass
        filename: str, optional
            See superclass
        conversion_function: func
            Function to convert between sampled parameters and constraints.
            By default this generates many additional parameters, see
            BNSPriorDict.default_conversion_function
        """
        if aligned_spin:
            default_file = 'aligned_spins_bns_tides_on.prior'
        else:
            default_file = 'precessing_spins_bns_tides_on.prior'
        if dictionary is None and filename is None:
            filename = os.path.join(DEFAULT_PRIOR_DIR, default_file)
            logger.info('No prior given, using default BNS priors in {}.'.format(filename))
        elif filename is not None:
            if not os.path.isfile(filename):
                filename = os.path.join(DEFAULT_PRIOR_DIR, filename)
        super(BNSPriorDict, self).__init__(dictionary=dictionary, filename=filename,
                                           conversion_function=conversion_function)

    def test_redundancy(self, key, disable_logging=False):
        logger.disabled = disable_logging
        logger.info("Performing redundancy check using BBHPriorDict(self).test_redundancy")
        bbh_redundancy = BBHPriorDict(self).test_redundancy(key)
        logger.disabled = False

        if bbh_redundancy:
            return True
        redundant = False

        sampling_parameters = {key for key in self if not isinstance(
            self[key], (prior.Constraint))}

        tidal_parameters = \
            {'lambda_1', 'lambda_2', 'lambda_tilde', 'delta_lambda_tilde'}

        if key in tidal_parameters:
            if len(tidal_parameters.intersection(sampling_parameters)) > 2:
                redundant = True
                logger.disabled = disable_logging
                logger.warning('{} already in prior. '
                               'This may lead to unexpected behaviour.'
                               .format(tidal_parameters.intersection(self)))
                logger.disabled = False
            elif len(tidal_parameters.intersection(sampling_parameters)) == 2:
                redundant = True
        return redundant

    @property
    def tidal(self):
        """ Return true if priors include phase parameters """
        return self.is_nonempty_intersection("tidal")


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
    'dchi_0': r'$\delta\varphi_{0}$',
    'dchi_1': r'$\delta\varphi_{1}$',
    'dchi_2': r'$\delta\varphi_{2}$',
    'dchi_3': r'$\delta\varphi_{3}$',
    'dchi_3S': r'$\delta\varphi_{3S}$',
    'dchi_3NS': r'$\delta\varphi_{3NS}$',
    'dchi_4': r'$\delta\varphi_{4}$',
    'dchi_4S': r'$\delta\varphi_{4S}$',
    'dchi_4NS': r'$\delta\varphi_{4NS}$',
    'dchi_5l': r'$\delta\varphi_{5l}$',
    'dchi_5lS': r'$\delta\varphi_{5lS}$',
    'dchi_5lNS': r'$\delta\varphi_{5lNS}$',
    'dchi_6': r'$\delta\varphi_{6}$',
    'dchi_6S': r'$\delta\varphi_{6S}$',
    'dchi_6NS': r'$\delta\varphi_{6NS}$',
    'dchi_6l': r'$\delta\varphi_{6l}$',
    'dchi_7': r'$\delta\varphi_{7}$',
    'dchi_7S': r'$\delta\varphi_{7S}$',
    'dchi_7NS': r'$\delta\varphi_{7NS}$',
    'dchi_kappaS': r'$\delta\kappa_s$',
    'dchi_kappaA': r'$\delta\kappa_a$',
    'dchi_kappa1': r'$\delta\kappa_1$',
    'dchi_kappa2': r'$\delta\kappa_2$',
    'f_window_div_f_peak': r'$f^\mathrm{tape} / f_{22}^\mathrm{peak}$',
    'NumCycles': r'$\mathcal{N}_\mathrm{GW}$',
}
