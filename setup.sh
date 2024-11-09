#! /bin/sh

# Unpack HOOMD-blue and reconfigure with `TwoStepVicsek`.
tar -xf hoomd-v2.8.2.tar.gz
cd hoomd-v2.8.2
cp ../TwoStepColoredVicsek/* hoomd/md/

# Compile HOOMD-blue (build using 4 CPUs).
mkdir build && cd build
cmake ..
make -j4
cd ./../../

# Create virtual environment and install required libraries.
pip3.6 install --upgrade pip
pip3.6 install virtualenv --upgrade
python3.6 -m venv v-env --without-pip --system-site-packages
v-env/bin/python3.6 -m pip install -r requirements.txt

