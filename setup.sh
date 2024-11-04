#! /bin/sh

# Unpack HOOMD-blue and reconfigure with `TwoStepVicsek`.
tar -xf hoomd-v2.8.2.tar.gz
cd hoomd-v2.8.2
cp ../TwoStepColoredVicsek/* hoomd/md/

# Compile HOOMD-blue (build using 4 CPUs).
mkdir build && cd build
cmake ..
make -j4

# Create virtual environment and install required libraries.
pip3 install --upgrade pip
pip3 install virtualenv --upgrade
python3 -m venv v-env --without-pip --system-site-packages
v-env/bin/python3 -m pip install -r requirements.txt
