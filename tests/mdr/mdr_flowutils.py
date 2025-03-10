import unittest
from bilby_tgr.mdr.asimov import flowutils


class example_event:
    def __init__(self, alpha):
        self.meta = {}
        self.meta['priors'] = {'alpha': {'peak': alpha}}
        self.pipeline = 'bilbymdr'


class TestMDRFlowUtils(unittest.TestCase):
    "Tests if the MDR specific tgr_flow logic behaves correctly."
    def test_sort_mdr(self):
        "Test if different MDR tests are correctly sorted."
        alphas = [-3.0, -2.0, -1.0, 0.0, 0.5, 1.5, 2.5, 3.0, 3.5, 4.0]
        keys = [
            'alpha_m3p0', 'alpha_m2p0', 'alpha_m1p0', 'alpha_0p0', 'alpha_0p5',
            'alpha_1p5', 'alpha_2p5', 'alpha_3p0', 'alpha_3p5', 'alpha_4p0']
        for alpha, key in zip(alphas, keys):
            predicted_key = flowutils.sort_mdr_by_subtype(example_event(alpha))
            self.assertEqual(key, predicted_key)

    def test_fill_in_mdr_metadata(self):
        "Test if mdr specific metadata is correctly produced."
        predicted_metadata = {
            'Alpha': 0.0,
            'AnalysisSoftware': 'bilbymdr',
            'Description': 'MDR analyses for deviation exponent of 0.0',
            'MaximumAEff': 3e-19,
            'MinimumAEff': -3e-19}
        alpha = 0.0
        metadata = flowutils.fill_in_mdr_specific_metadata(example_event(alpha), None)
        self.assertEqual(metadata, predicted_metadata)


if __name__ == "__main__":
    unittest.main()
