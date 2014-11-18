import unittest
from treetoreads import TreeToReads

class PrimesTestCase(unittest.TestCase):
    """Tests for `primes.py`."""

    def test_is_five_prime(self):
        """Is five successfully determined to be prime?"""
        self.assertTrue(TreeToReads(configfi='tests/input/test1.cfg',run=0))

if __name__ == '__main__':
    unittest.main()