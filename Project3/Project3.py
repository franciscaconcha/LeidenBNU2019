import numpy
from amuse.lab import *
from amuse.ext.galactics_model import new_galactics_model

M_galaxy = 1.0e12 | units.MSun  # Mass of the galaxy
R_galaxy = 10 | units.kpc  # Radius of the galaxy
n_halo = 20000  # Number of particles for halo
n_bulge = 10000 # Number of particles for bulge
n_disk = 10000 # Number of particles for disk

# Converter: used to turn n-body system units to the physical units we are interested in
converter=nbody_system.nbody_to_si(M_galaxy, R_galaxy)
galaxy1 = new_galactics_model(n_halo,
                               converter,
                               bulge_number_of_particles=n_bulge,
                               disk_number_of_particles=n_disk)
    
# rotate function from the galactics_model module
galaxy1.rotate(0., numpy.pi/2, numpy.pi/4)

# Set position and velocity for the galaxy
galaxy1.position += [100.0, 100, 0] | units.kpc
galaxy1.velocity += [-10.0, 0.0, -10.0] | units.km/units.s

# Start SPH code
dynamics = Gadget2(converter, number_of_workers=4)
dynamics.parameters.epsilon_squared = (100 | units.parsec)**2  # Softening length

set1 = dynamics.particles.add_particles(galaxy1)
dynamics.particles.move_to_center()

t_end = 200 | units.Myr

# Evolve galactic model
dynamics.evolve_model(t_end)

