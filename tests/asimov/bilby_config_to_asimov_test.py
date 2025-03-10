import unittest
import bilby_tgr
import os
import json


class TestConfigParsing(unittest.TestCase):
    def setUp(self):
        self.location = os.path.dirname(__file__)

    def tearDown(self):
        del self.location

    def compare_input_output(self, input, output):
        parsed_config = bilby_tgr.asimov.utils.bilby_config_to_asimov(input)
        with open(output, 'r') as file:
            predicted_ledger = json.load(file)
        self.assertEqual(parsed_config, predicted_ledger)

    def test_prior_dict_config(self):
        """Tests config parsing with prior dictionary specified"""
        input = os.path.join(self.location, 'config_prior_dict.ini')
        output = os.path.join(self.location, 'ledger_prior_dict.json')
        self.compare_input_output(input, output)


if __name__ == "__main__":
    unittest.main()
