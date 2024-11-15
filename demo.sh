#! /bin/sh

# Activate the virtual environment created by setup.sh.
source ./v-env/bin/activate

# Create copy of python script in the /build directory.
cp ./run_simulation.py ./hoomd-v2.8.2/build/run_simulation.py
cd ./hoomd-v2.8.2/build/

# Simulation: Polar vortex lattice (Chiral).
L=300     # Linear system size.
rho=2.0   # Global particle density.
v=1.0     # Self-propulsion speed.
sym=1     # Orientational alignment interaction symmetry.
alpha=0.1 # Orientational alignment coupling strength.
beta=1.0  # Rotational frequency coupling strength.
eta=0.01  # Rotational frequency noise strength.
w0=0.0    # Rotational frequency bias.
tau=100.0 # Persistence/memory time.
python3.6 ./run_simulation.py $L $rho $v $sym $alpha $beta $eta $w0 $tau

# Simulation: Polar vortex lattice (High-Density).
L=200     # Linear system size.
rho=4.0   # Global particle density.
v=1.0     # Self-propulsion speed.
sym=1     # Orientational alignment interaction symmetry.
alpha=0.1 # Orientational alignment coupling strength.
beta=0.0  # Rotational frequency coupling strength.
eta=0.01  # Rotational frequency noise strength.
w0=0.0    # Rotational frequency bias.
tau=100.0 # Persistence/memory time.
python3.6 ./run_simulation.py $L $rho $v $sym $alpha $beta $eta $w0 $tau

# Simulation: Polar vortex lattice (Medium-Density).
L=200     # Linear system size.
rho=2.0   # Global particle density.
v=1.0     # Self-propulsion speed.
sym=1     # Orientational alignment interaction symmetry.
alpha=0.1 # Orientational alignment coupling strength.
beta=0.0  # Rotational frequency coupling strength.
eta=0.01  # Rotational frequency noise strength.
w0=0.0    # Rotational frequency bias.
tau=100.0 # Persistence/memory time.
python3.6 ./run_simulation.py $L $rho $v $sym $alpha $beta $eta $w0 $tau

# Simulation: Nematic vortex lattice.
L=200     # Linear system size.
rho=2.0   # Global particle density.
v=1.0     # Self-propulsion speed.
sym=2     # Orientational alignment interaction symmetry.
alpha=0.1 # Orientational alignment coupling strength.
beta=0.0  # Rotational frequency coupling strength.
eta=0.01  # Rotational frequency noise strength.
w0=0.0    # Rotational frequency bias.
tau=100.0 # Persistence/memory time.
python3.6 ./run_simulation.py $L $rho $v $sym $alpha $beta $eta $w0 $tau
