import bilby.gw.prior as prior
import os
from bilby.core.utils import logger

DEFAULT_PRIOR_DIR = os.path.join(os.path.dirname(__file__), 'prior_files')


class BBHPriorDict(prior.BBHPriorDict):
    def __init__(self, dictionary=None, filename=None, aligned_spin=False,
                 conversion_function=None):
        """ Initialises a Prior set for Binary Black holes with LIV parameter.
        Copied from bilby.

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
        super(prior.BBHPriorDict, self).__init__(dictionary=dictionary, filename=filename,
                                                 conversion_function=conversion_function)


class BNSPriorDict(prior.BNSPriorDict):
    def __init__(self, dictionary=None, filename=None, aligned_spin=True,
                 conversion_function=None):
        """ Initialises a Prior set for Binary Neutron Stars.
        Copied from bilby.

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
        super(prior.BNSPriorDict, self).__init__(dictionary=dictionary, filename=filename,
                                                 conversion_function=conversion_function)
