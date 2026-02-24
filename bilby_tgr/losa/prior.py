import bilby.gw.prior as prior
import os
from bilby.core.utils import logger


DEFAULT_PRIOR_DIR = os.path.join(os.path.dirname(__file__), 'prior_files')


class BBHPriorDict(prior.CBCPriorDict):
    def __init__(self, dictionary=None, filename=None, aligned_spin=True,
                 conversion_function=None):
        """
        Initialises a Prior set for Binary Black holes with LOSA parameter.

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
                logger.info('Using precessing spin prior')
            filename = os.path.join(DEFAULT_PRIOR_DIR, fname)
            logger.info('No prior given, using default BBH priors in {}.'.format(filename))
        elif filename is not None:
            if not os.path.isfile(filename):
                filename = os.path.join(DEFAULT_PRIOR_DIR, filename)
        super(BBHPriorDict, self).__init__(dictionary=dictionary, filename=filename,
                                           conversion_function=conversion_function)

    def test_redundancy(self, key, disable_logging=False):
        """
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
        losa_parameter = {'a_by_c'}

        for independent_parameters, parameter_set in \
                zip([2, 2, 1, 1, 1, 1, 1],
                    [mass_parameters, spin_azimuth_parameters,
                     spin_tilt_1_parameters, spin_tilt_2_parameters,
                     inclination_parameters, distance_parameters,
                     losa_parameter]):
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
        if dictionary is None and filename is None:
            if aligned_spin:
                fname = 'aligned_spins_bns.prior'
                logger.info('Using aligned spin prior')
            else:
                fname = 'precessing_spins_bns.prior'
                logger.info('Using precessing spin prior')
            filename = os.path.join(DEFAULT_PRIOR_DIR, fname)
            logger.info('No prior given, using default BNS priors in {}.'.format(filename))
        elif filename is not None:
            if not os.path.isfile(filename):
                filename = os.path.join(DEFAULT_PRIOR_DIR, filename)
        super(BBHPriorDict, self).__init__(dictionary=dictionary, filename=filename,
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
    'a_by_c': r'$a/c$'
}
