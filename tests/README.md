# Tests for module Data_Reduction

Files which start with "test_" are unittest-based.  Other files are user
directed tests.

## Contents

### `check_ndfile.py`

This checks the base classes `Observation` and `Map` using a `t`-file but any CSV will do.

### `test_malargue.py`

This checks the subclasses `Observation` and `Map` in submodule `Malargue`.  It also checks the `Channel` subclass.

### `test-tfiles.py`

This will test the subclasses `Observation` and `Map` in submodule `GAVRT`.

## Notes

* I wonder if `GAVRT.Mysql` should become `GAVRT` and `GAVRT` be renamed `OldGAVRT`.