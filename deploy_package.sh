# build dist
python3 setup.py sdist

# test on test pypi
python3 -m twine upload --repository-url https://test.pypi.org/legacy/ dist/dracoon-0.0.1.tar.gz

# upload to pypi
python3 -m twine upload dist/dracoon-0.0.1.tar.gz
