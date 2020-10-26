import unittest
from mock import mock_open, patch

from extract_coverage_caveat import extract_qc_call


class TestMethods(unittest.TestCase):
    def test_extract_coverage_caveat_ok1(self):
        test_data = \
"""
# target_coverage_histogram, bam: path.bam
# coverage\tbases_at_coverage\ttotal_bases\tfraction_bases_at_coverage
110\t3\t548\t0.04
111\t7\t548\t0.21
112\t7\t548\t0.25
113\t10\t548\t0.25
120\t7\t548\t0.25
"""

        with patch('extract_coverage_caveat.extract_qc_call.open',
                   mock_open(read_data=test_data), create=True) as test_input_file:
            self.assertEquals(extract_qc_call(test_input_file.return_value, 0.95, 100, 0.95, 50), "OK")

    def test_extract_coverage_caveat_ok2(self):
        test_data = \
"""
# target_coverage_histogram, bam: path.bam
# coverage\tbases_at_coverage\ttotal_bases\tfraction_bases_at_coverage
40\t3\t548\t0.04
101\t7\t548\t0.21
112\t7\t548\t0.25
113\t10\t548\t0.25
120\t7\t548\t0.25
"""

        with patch('extract_coverage_caveat.extract_qc_call.open',
                   mock_open(read_data=test_data), create=True) as test_input_file:
            self.assertEquals(extract_qc_call(test_input_file.return_value, 0.95, 100, 0.95, 50), "OK")

    def test_extract_coverage_caveat_warn1(self):
        test_data = \
"""
# target_coverage_histogram, bam: path.bam
# coverage\tbases_at_coverage\ttotal_bases\tfraction_bases_at_coverage
40\t3\t548\t0.04
80\t7\t548\t0.21
112\t7\t548\t0.25
113\t10\t548\t0.25
120\t7\t548\t0.25
"""

        with patch('extract_coverage_caveat.extract_qc_call.open',
                   mock_open(read_data=test_data), create=True) as test_input_file:
            self.assertEquals(extract_qc_call(test_input_file.return_value, 0.95, 100, 0.95, 50), "WARN")

    def test_extract_coverage_caveat_fail1(self):
        test_data = \
"""
# target_coverage_histogram, bam: path.bam
# coverage\tbases_at_coverage\ttotal_bases\tfraction_bases_at_coverage
40\t3\t548\t0.04
45\t7\t548\t0.21
112\t7\t548\t0.25
113\t10\t548\t0.25
120\t7\t548\t0.25
"""

        with patch('extract_coverage_caveat.extract_qc_call.open',
                   mock_open(read_data=test_data), create=True) as test_input_file:
            self.assertEquals(extract_qc_call(test_input_file.return_value, 0.95, 100, 0.95, 50), "FAIL")

    def test_extract_coverage_empty_file(self):
        test_data = \
"""
# target_coverage_histogram, bam: path.bam
# coverage\tbases_at_coverage\ttotal_bases\tfraction_bases_at_coverage
"""

        with patch('extract_coverage_caveat.extract_qc_call.open',
                   mock_open(read_data=test_data), create=True) as test_input_file:
            self.assertEquals(extract_qc_call(test_input_file.return_value, 0.95, 100, 0.95, 50), "FAIL")
