#!/usr/bin/env bash
#For executing this script you need root privileges

# run building
python3 setup.py build

# run installation
python3 setup.py install

# remove test.py in Python packages
cd /usr/local/lib/python3.4/dist-packages/MirrorTree/

rm test.py

echo "Installation DONE! :)"
