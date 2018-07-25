#!/bin/bash
echo Did you remember to increment the version number?
python setup.py sdist bdist_wheel
twine upload dist/*
