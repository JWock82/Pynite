# -*- coding: utf-8 -*-
"""
MIT License

Copyright (c) 2020 D. Craig Brinck, SE; tamalone1
"""

import unittest

# Run the tests in this module
# Use warnings flag to suppress the PendingDeprecationWarning 
# from numpy.matrix
# unittest.main(warnings='ignore')
test_suite = unittest.TestLoader().discover("Testing", pattern='test_*.py')

# `TextTestRunner` does not exit the module. CI will get confused unless we save the result
# and send the proper exit code.
result = unittest.TextTestRunner().run(test_suite)

# Send the proper exit code for Travis-CI to read
if result.wasSuccessful():
    exit(0)
else:
    exit(1)
