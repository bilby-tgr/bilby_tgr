import cbcflow
from asimov import logger


class Applicator:
    """
    Wrapper to the general TGR applicator to expand the copied PE priors
    """

    def __init__(self, ledger):
        hook_data = ledger.data["hooks"]["applicator"]["mdrflow"]
        self.ledger = ledger
        self.library = cbcflow.core.database.LocalLibraryDatabase(
            hook_data["library location"]
        )
        self.library.git_pull_from_remote(automated=True)
        # in case of disparity between data product in cbcflow and config,
        # choose which one to prefer (config by default)
        if 'data preference' in hook_data:
            if hook_data['data preference'] == 'cbcflow':
                self.prefer_config = False
            elif hook_data['data preference'] == 'config':
                self.prefer_config = True
            else:
                logger.warning("Unrecognized data preference applicator option. Setting to config")
                self.prefer_config = True
        else:
            self.prefer_config = True

    def run(self, sid=None):
        """
        Apply information from CBCFlow to an asimov event based on PE illustrative result.
        Changes some settings on event level to be compatible with FTI.
        Applies FTI selection criteria and adds FTI analyses.

        Parameters
        ==========
        sid: str
            S-name of event
        """

        from asimov.cli.application import apply_via_plugin

        # Get event settings from PE illustrative result using tgrflow
        apply_via_plugin(sid, 'tgrflow')
        event = self.ledger.get_event(sid)[0]

        priors = event.meta['priors']

        priors['chirp mass']['minimum'] = priors['chirp mass']['minimum'] * 0.8
        priors['chirp mass']['maximum'] = priors['chirp mass']['maximum'] * 1.2

        self.ledger.update_event(event)
