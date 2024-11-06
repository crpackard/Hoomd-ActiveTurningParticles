#! /bin/sh

cp ./run_simulation.py ./hoomd-v2.8.2/build/run_simulation.py
cd ./hoomd-v2.8.2/build/

python3.9 ./run_simulation.py --L 100 --rho 2.0 --beta 0.0 --tau 100 --sym 2
#python3.9 ./run_simulation.py
