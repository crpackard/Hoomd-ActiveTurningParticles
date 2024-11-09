
# Import these modules for file/directory management.
import sys, os, random

# .
from typing import List

# These imports are needed for creating initial conditions of simulation.
import numpy as np
from math import pi as π

# These imports are needed for running simulations (path structures assumes code is run from /hoomd-v2.8.2/build directory).
import hoomd, gsd.hoomd
from hoomd.md import nlist as nl

# These imports are used to visualize simulation dynamics.
import imageio
from tqdm import tqdm
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
rc = {'font.size': 12,
  'font.family' : 'serif',
  'mathtext.fontset' : 'stix',
  'savefig.dpi': 300}
plt.rcParams.update(rc)
plt.rcParams['font.serif'] = ['Times New Roman'] + plt.rcParams['font.serif']

# Location where simulation snapshots will be saved to.
output_dir: str = os.path.join('..', '..', 'data')
assert os.path.exists(output_dir), 'This script is intended to be run from the directory: ./hoomd-v2.8.2/build/'

# Script takes parameter input at command-line.
def main(
  L: int = None,
  ρ: float = None,
  v: float = None,
  α: float = None,
  β: float = None,
  η: float = None,
  τ: float = None,
  ω0: float = None,
  sym: int = 1
) -> None:

  assert (sym==1) or (sym==2), f'particle interaction symmetry must be either `1` or `2`'

  # Create sub-directory based on interaction symmetry and rotational frequency coupling.
  subdir_name: str = f'm={sym}_b={β:.3f}_w0={ω0:.3f}'
  subdir_path: str = os.path.join(output_dir, subdir_name)
  if not os.path.exists(subdir_path):
    os.mkdir(subdir_path)

  # Create another sub-directory based on length and time scales.
  subsubdir_name: str = f'L={L}_d={ρ:.3f}_v={v:.3f}_τ={τ:.3f}'
  subsubdir_path: str = os.path.join(subdir_path, subsubdir_name)
  if not os.path.exists(subsubdir_path):
    os.mkdir(subsubdir_path)

  # Create a simulation filename based on Vicsek-like coupling and noise strength.
  simfile_name: str = f'a={α:.3f}_e={η:.3f}.gsd'
  simfile_path: str = os.path.join(subsubdir_path, simfile_name)
  if not os.path.exists(simfile_path):
    print(f'\nCreating new simulation datafile: {simfile_path}')
  else:
    print(f'\nSimulation datafile already exists: {simfile_path}')
    return True

  # Initialize HOOMD-Blue simulation.
  sim = hoomd.context.SimulationContext()
  with sim:
    #print(L); import sys; sys.exit()
    hoomd.context.initialize()

    # Number of particles in simulation.
    N = int(ρ * L**2)

    # Define initial conditions (randomized positions, orientations, and rotational frequencies).
    xs = list(np.random.uniform(-L/2, L/2, N))
    ys = list(np.random.uniform(-L/2, L/2, N))
    θs = list(np.random.uniform(-π, π, N))
    ωs = list(np.random.uniform(-π, π, N))
    unitcell = hoomd.lattice.unitcell(
      N=N,
      a1=[L, 0, 0],
      a2=[0, L, 0],
      a3=[0, 0, 1],
      dimensions=2,
      position=[[xs[ii], ys[ii], 0] for ii in range(N)],
      type_name=['A']*N,
      orientation=[[θs[ii], 0, 0, 1] for ii in range(N)],
      moment_inertia=[[np.cos(ω), np.sin(ω), 0] for ω in ωs])
    hoomd.init.create_lattice(unitcell, n=1)
    print(f'\nSimulation initialized with N={N} particles')

    # Setup the HOOMD-Blue integrator.
    nl_cell = nl.cell(r_buff=0)
    all = hoomd.group.all()
    hoomd.md.integrate.mode_standard(dt=1)
    integrator = hoomd.md.integrate.vicsek(
      group=all,
      nlist=nl_cell,
      v0=v,
      alpha=α,
      beta=β,
      delta=η,
      tau=τ,
      bias=ω0,
      sym=sym,
      seed=random.randint(1,10e5)
    )

    # Run the simulation for `dt` time-steps and save a snapshot of the very last one.
    def evolve(dt: int) -> None:
      dump = hoomd.dump.gsd(simfile_path, period=dt, group=all, overwrite=False)
      hoomd.run(dt)
      dump.disable()

    # Save snapshots at logarithmically spaced time-step intervals.
    evolve(dt=int(1e0)-0)
    evolve(dt=int(1e1)-int(1e0))
    evolve(dt=int(1e2)-int(1e1))
    evolve(dt=int(1e3)-int(1e2))
    evolve(dt=int(1e4)-int(1e3))
    if False:
      evolve(dt=int(1e5)-int(1e4))
      evolve(dt=int(1e6)-int(1e5))
    tf = int(1e6)

    # Create a movie of the simulation.
    movie_file_path = simfile_path.replace('.gsd', '.gif')
    tmp_gsd_path: str = simfile_path.replace('.gsd', '_tmp.gsd')
    with imageio.get_writer(movie_file_path, mode='I', fps=30) as writer:
      for tIdx in tqdm(range(100)):
        # Evolve the simulation a single time-step forward, and load the data.
        hoomd.dump.gsd(tmp_gsd_path, period=1, group=all, overwrite=True)
        hoomd.run(1)
        hoomd_snapshots = gsd.hoomd.open(tmp_gsd_path, mode='rb')
        snapshot = hoomd_snapshots[0]
        # Configure a new matplotlib figure.
        fig = plt.figure()
        ax = fig.add_subplot()
        ax.set_title(f'$t={tf+tIdx}$')
        ax.set_xlim(-L/2, L/2)
        ax.set_ylim(-L/2, L/2)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_facecolor('k')
        # Plot the position of each particle, and its instantaneous rotational frequency.
        rx: np.ndarray = np.array([r[0] for r in snapshot.particles.position])
        ry: np.ndarray = np.array([r[1] for r in snapshot.particles.position])
        ωx: np.ndarray = np.array([ω[0] for ω in snapshot.particles.moment_inertia])
        ωy: np.ndarray = np.array([ω[1] for ω in snapshot.particles.moment_inertia])
        ωθ: np.ndarray = np.arctan2(ωy, ωx)
        if (sym==1):
          cmap = 'coolwarm'
          z, zlbl, zmin, zmax = ωθ, r'$\omega$', -0.2, 0.2
        elif (sym==2):
          cmap = 'coolwarm'
          z, zlbl, zmin, zmax = ωθ, r'$\omega$', -0.2, 0.2
        else:
          pass
        ax.scatter(rx, ry, c=z, cmap=cmap, vmin=zmin, vmax=zmax, s=0.10)
        # Add a color bar legend (align it tight to the plot window).
        fake_data = np.asarray([[0]])
        cmap = plt.cm.get_cmap(cmap)
        norm = plt.Normalize(zmin, zmax)
        im = ax.imshow(fake_data, cmap=cmap, vmin=zmin, vmax=zmax)
        divider = make_axes_locatable(ax)
        cax = divider.append_axes('right', size='5%', pad=0.05)
        color_bar = plt.colorbar(im,
          ax=ax,
          cax=cax,
          orientation='vertical',
          ticks=[zmin, 0, zmax],
        )
        color_bar.ax.set_title(zlbl)
        # Save the plot, load the image into the movie, and delete the temporary .gsd and .png files.
        tmp_png_path = os.path.join('.', 'tmp.png')
        plt.savefig(tmp_png_path, bbox_inches='tight', dpi=150)
        plt.close()
        writer.append_data(imageio.imread(tmp_png_path))
        os.remove(tmp_png_path)
      os.remove(tmp_gsd_path)

if __name__ == "__main__":
  main(
    L=int(sys.argv[1]),    # Linear system size.
    ρ=float(sys.argv[2]),  # Global particle density.
    v=float(sys.argv[3]),  # Self-propulsion speed.
    sym=int(sys.argv[4]),  # Orientational alignment interaction symmetry.
    α=float(sys.argv[5]),  # Orientational alignment coupling strength.
    β=float(sys.argv[6]),  # Rotational frequency coupling strength.
    η=float(sys.argv[7]),  # Rotational frequency noise strength.
    ω0=float(sys.argv[8]), # Rotational frequency bias.
    τ=float(sys.argv[9]),  # Persistence/memory time.
  )
