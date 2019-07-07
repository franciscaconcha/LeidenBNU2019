import numpy as np
import matplotlib.pyplot as plt

import amuse.lab as al
import amuse.ext.bridge as bridge


def plot_cluster(X, Y, num=1):

	fig = plt.figure(num)
	ax = fig.add_subplot(111, aspect='equal')

	ext = 1.2*np.max(np.abs([X, Y]))

	ax.set_xlim(-ext, ext)
	ax.set_ylim(-ext, ext)

	ax.scatter(X, Y, s=5.)


def R_virial_to_plummer (Rv):
	# Convert the virial radius of a plummer sphere to its plummer radius

	a = 3.*np.pi/16.*Rv

	return a


class GalacticPotential:
	# class defining the gravitational forces exerted by the galaxy
	# crucial is the get_gravity_at_point function, which gives the gravitational acceleration of a particle at some location

	def __init__ (self, R, M, alpha):
		# Enclosed mass within r: M(r) = M*(r/R)^alpha
		# So M is the mass inside R, alpha gives scaling

		self.radius = R
		self.mass = M
		self.alpha = alpha


	def get_gravity_at_point(self, eps, x, y, z):

		r2 = x**2. + y**2. + z**2.
		r = r2**0.5

		m = self.mass*(r/self.radius)**self.alpha

		# Newton: if a mass distribution is spherically symmetric, shells outside your position don't count, 
		# and those inside are effectively a central point mass
		fr = al.constants.G*m/r2

		ax = -fr*x/r
		ay = -fr*y/r
		az = -fr*z/r

		return ax, ay, az


	def circular_velocity(self, r):

		m = self.mass*(r/self.radius)**self.alpha
		vc = (al.constants.G*m/r)**0.5
		return vc


N = 100						# Stars in the cluster
Rv = 1. | al.units.parsec	# virial radius of cluster

RSun = 50 | al.units.parsec	# Position of cluster in the galaxy

M = al.new_kroupa_mass_distribution(N, 100. | al.units.MSun)	# N samples from Kroupa IMF, with each max. 100 MSun

# Contains the plummer radius and total mass of the plummer distribution
converter = al.nbody_system.nbody_to_si( R_virial_to_plummer(Rv), M.sum() )

bodies = al.new_plummer_model(N, converter)
bodies.mass = M




t_end = 2.5 | al.units.Myr	# End time
dt = 0.01 | al.units.Myr	# Bridge timestep


# Galaxy model with 1.6x10^10 MSun within 1 kpc, and power law scaling with index 1.2
galaxy_gravity = GalacticPotential(1. | al.units.kpc, 1.6E10 | al.units.MSun, 1.2)


# Place cluster in galaxy
bodies.x += RSun
bodies.vy += galaxy_gravity.circular_velocity(RSun)

cluster_gravity = al.ph4(converter)
cluster_gravity.particles.add_particles(bodies)


# Bridge couples a number of gravitational codes to another code. In this case, 
# the cluster evolves under its own gravity and under that of the galaxy.
# The galaxy is represented by just a smooth potential, and its forces are updated
# every timestep dt. 
bridge_gravity = bridge.bridge()
bridge_gravity.add_system(cluster_gravity, (galaxy_gravity,))	# The codes in the bracket (1 in this case) influence the first code


# Plot initial conditions
plot_cluster(bridge_gravity.particles.x.value_in(al.units.kpc), bridge_gravity.particles.y.value_in(al.units.kpc), num=1)


bridge_gravity.evolve_model(t_end, timestep=dt)


# Plot outcome
plot_cluster(bridge_gravity.particles.x.value_in(al.units.kpc), bridge_gravity.particles.y.value_in(al.units.kpc), num=2)




cluster_gravity.stop()

plt.show()
