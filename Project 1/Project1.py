import numpy as np
import matplotlib.pyplot as plt

import amuse.lab as al
import amuse.ext as ae


star1 = al.Particle(mass=10.  | al.units.MSun)
star2 = al.Particle(mass=0.08 | al.units.MSun)


print 'Starting EVtwin'

stellar = al.EVtwin()
#stellar = al.MESA()

print 'Started EVtwin'

stellar.particles.add_particle( star1 )
stellar.particles.add_particle( star2 )

print 'Added particles'

stellar.evolve_model( 1. | al.units.Myr )

print 'Finished initial evolve'


N_sph = 1000

star1_sph = ae.star_to_sph.convert_stellar_model_to_SPH(stellar.particles[0], N_sph, seed=874325)
star2_sph = ae.star_to_sph.convert_stellar_model_to_SPH(stellar.particles[1], N_sph, seed=874325)

print 'Converted stellar models'

stellar.stop()


distance = 10. | al.units.RSun

star2_sph.x += distance
star2_sph.vx -= np.sqrt(2.*al.units.G*stellar.particles.mass.sum()/distance)


sph_particles = al.Particles()
sph_particles.add_particles(star1_sph)
sph_particles.add_particles(star2_sph)
sph_particles.move_to_center()


converter = al.nbody_system.nbody_to_si(1. | al.units.hour, 1. | al.units.RSun)

hydro = al.Gadget2(converter)
hydro.gas_particles.add_particles(sph_particles)

hydro.evolve_model( 10. | al.units.yr )
hydro.gas_particles.copy_values_of_attributes_to(['x', 'y', 'z', 'vx', 'vy', 'vz', 'density', 'u', 'pressure'], sph_particles)

hydro.stop()


sph_particles.move_to_center()

star3 = al.convert_SPH_to_stellar_model(sph_particles)

print 'Converted SPH'

stellar = al.EVtwin()
stellar.new_particle_from_model(star3, 0. | al.units.Myr)

for i in range(5):

	print "Starting step " + str(i+1)

	stellar.evolve_model( 1. | al.units.Myr )

	print stellar.particles[0].mass
	print stellar.particles[0].radius
	print stellar.particles[0].temperature
	print stellar.particles[0].stellar_type

stellar.stop()
