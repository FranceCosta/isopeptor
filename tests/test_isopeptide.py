#! /usr/env/python
# -*- coding: utf-8 -*-

import unittest
import doctest
import isopeptor.isopeptide

# Load doctests
def load_tests(loader, tests, ignore):
    tests.addTests(
        doctest.DocTestSuite(
            isopeptor.isopeptide,
            optionflags=doctest.NORMALIZE_WHITESPACE
        )
    )
    return tests

if __name__ == "__main__":
    unittest.main()