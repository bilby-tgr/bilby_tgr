def sort_mdr_by_subtype(analysis):
    """
    Sort MDR analysis into the proper subcategory based on alpha used

    input:
    analysis - asimov production for given event

    output:
    uid - string determining the analysis subcategory for given gr test
    """
    try:
        uid = "alpha_" + str(analysis.meta["priors"]["alpha"]["peak"]).replace('.', 'p').replace('-', 'm')
    except KeyError:
        uid = 'uncategorized'
    return uid


def fill_in_mdr_specific_metadata(analysis, corresponding_analysis):
    """
    For mdr analysis, fill in fields in metadata other than result list.

    input:
    analysis - asimov production for given event
    corresponding_analysis - equivalent subanalysis as stored in cbcflow

    output:
    analysis_output - dictionary used to update cbcflow with new information
    """
    analysis_output = {}
    try:
        analysis_output["Alpha"] = analysis.meta["priors"]['alpha']['peak']
    except KeyError:
        pass
    analysis_output["AnalysisSoftware"] = str(analysis.pipeline)
    if corresponding_analysis is None or "Description" not in corresponding_analysis:
        try:
            text = analysis_output["Alpha"]
            analysis_output["Description"] = f"MDR analyses for deviation exponent of {text}"
        except KeyError:
            pass
    if corresponding_analysis is None or "MaximumAEff" not in corresponding_analysis:
        try:
            alpha = analysis.meta["priors"]['alpha']['peak']
        except KeyError:
            alpha = None
        alphas = [0, 0.5, 1.0, 1.5, 2.5, 3.0, 3.5, 4.0, -1.0, -2.0, -3.0]
        bounds = [3e-19, 2e-18, 1e-18, 6e-18, 2e-18, 9e-19, 7e-19, 1e-18, 2.5e-20, 2e-21, 1.7e-22]
        for alpha_listed, bound in zip(alphas, bounds):
            if alpha == alpha_listed:
                analysis_output["MaximumAEff"] = bound
                analysis_output["MinimumAEff"] = -bound
    return analysis_output
