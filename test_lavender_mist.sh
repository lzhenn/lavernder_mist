#pip uninstall lavender_mist 
#rm -rf build dist *egg-info
python setup.py sdist
pip install .
cd test_case
python test.py
cd ..
