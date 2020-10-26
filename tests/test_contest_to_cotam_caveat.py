import unittest
from mock import mock_open, patch

from contest_to_contam_caveat import extract_qc_call


class TestMethods(unittest.TestCase):
    def test_extract_contam_qc_ok1(self):
        test_data = \
"""
name\tpopulation\tpopulation_fit\tcontamination\tconfidence_interval_95_width\tconfidence_interval_95_low\tconfidence_interval_95_high\tsites
META\tCEU\tn/a\t0.1\t0.1\t0.1\t0.2\t875
"""

        with patch('extract_coverage_caveat.extract_qc_call.open',
                   mock_open(read_data=test_data), create=True) as test_input_file:
            self.assertEquals(extract_qc_call(test_input_file.return_value, 1), "OK")

    def test_extract_contam_qc_fail1(self):
        test_data = \
"""
name\tpopulation\tpopulation_fit\tcontamination\tconfidence_interval_95_width\tconfidence_interval_95_low\tconfidence_interval_95_high\tsites
META\tCEU\tn/a\t2\t0.1\t0.1\t0.2\t875
"""

        with patch('extract_coverage_caveat.extract_qc_call.open',
                   mock_open(read_data=test_data), create=True) as test_input_file:
            self.assertEquals(extract_qc_call(test_input_file.return_value, 1), "FAIL")

    def test_extract_contam_qc_ok2(self):
        test_data = \
"""
name\tpopulation\tpopulation_fit\tcontamination\tconfidence_interval_95_width\tconfidence_interval_95_low\tconfidence_interval_95_high\tsites
"""

        with patch('extract_coverage_caveat.extract_qc_call.open',
                   mock_open(read_data=test_data), create=True) as test_input_file:
            self.assertEquals(extract_qc_call(test_input_file.return_value, 1), "OK")

    def test_extract_contam_qc_ok3(self):
        test_data = \
"""
Warning: We're throwing out lane META since it has fewer than 500 read bases at genotyped positions
name\tpopulation\tpopulation_fit\tcontamination\tconfidence_interval_95_width\tconfidence_interval_95_low\tconfidence_interval_95_high\tsites
"""

        with patch('extract_coverage_caveat.extract_qc_call.open',
                   mock_open(read_data=test_data), create=True) as test_input_file:
            self.assertEquals(extract_qc_call(test_input_file.return_value, 1), "OK")