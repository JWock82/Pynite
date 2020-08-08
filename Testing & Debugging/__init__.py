# -*- coding: utf-8 -*-
"""
MIT License

Copyright (c) 2020 D. Craig Brinck, SE; tamalone1
"""

import unittest

if __name__ == '__main__':

    # Run the tests in this module
    # Use warnings flag to suppress the PendingDeprecationWarning 
    # from numpy.matrix
    # unittest.main(warnings='ignore')
    loader = unittest.TestLoader()
    testSuite = loader.discover("Testing & Debugging", pattern='test_*.py')
    testRunner = unittest.TextTestRunner(warnings='ignore')
    testRunner.run(testSuite)