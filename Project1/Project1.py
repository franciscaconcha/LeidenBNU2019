import numpy
from amuse.community.fi.interface import Fi
from amuse.lab import *
from amuse.units import units
from amuse.units import nbody_system
from amuse.ext.protodisk import ProtoPlanetaryDisk

def rotate(disk, a, b, c):
    """ Function to rotate a disk.
        
        :disk: ProtoPlanetaryDisk particles to be rotated
        :param a: angle in the x axis
        :param b: angle in the y axis
        :param c: angle in the z axis"""
          
    a = numpy.radians(a)
    ca, sa = numpy.cos(a), numpy.sin(a)
    Rx = numpy.array([[1, 0, 0], [0, ca, -sa], [0, sa, ca]])

    b = numpy.radians(b)
    cb, sb = numpy.cos(b), numpy.sin(b)
    Ry = numpy.array([[cb, 0, sb], [0, 1, 0], [-sb, 0, cb]])

    c = numpy.radians(c)
    cc, sc = numpy.cos(c), numpy.sin(c)
    Rz = numpy.array([[cc, -sc, 0], [sc, cc, 0], [0, 0, 1]])

    R = numpy.dot(Rz, numpy.dot(Ry, Rx))

    for p in disk:
        r = numpy.array([p.x, p.y, p.z])
        v = numpy.array([p.vx, p.vy, p.vz])

        nr = numpy.dot(R, r)
        nv = numpy.dot(R, v)

        p.x, p.y, p.z = nr[0], nr[1], nr[2]
        p.vx, p.vy, p.vz = nv[0], nv[1], nv[2]

    return disk


# Create star with 1 solar mass
stars = Particles(mass=1. | units.MSun)

# Solar radius
stars[0].radius = 700000. | units.km

# Locate the star at a certain position
stars[0].x, stars[0].y, stars[0].z = 100. | units.au, 0. | units.au, 0. | units.au

# Give it a certain velocity
stars[0].vx, stars[0].vy, stars[0].vz = 20. | (units.au / units.yr), 0. | (units.au / units.yr), 0. | (units.au / units.yr)

# Create a converter. In this way we turn n-body system units to the physical units we are interested in
converter = nbody_system.nbody_to_si(1. | units.MSun, 1. | units.au)

# Create protoplanetary disk
N_particles = 1000  # Number of SPH particles
disk = ProtoPlanetaryDisk(N_particles, convert_nbody=converter, densitypower=1.5, Rmin=0.5, Rmax=100, q_out=1.)  # Rmin and Rmax in AU

disk_particles = disk.result  # Disk particles
disk_particles.h_smooth = 0.06 | units.au  # Smoothing length

# Rotate the disk 90 degrees around x axis
disk = rotate(disk_particles, 90, 0, 0)

# Start SPH code and add particles
sph = Fi(converter)
sph.gas_particles.add_particles(disk_particles)  # Disk particles are added as gas particles
sph.dm_particles.add_particles(stars[0])  # Star is added as dark matter particle

# Evolve code in time steps
t_end = 100 | units.yr
dt = 1 | units.yr
t = 0 | units.yr

while t < t_end:
    t += dt
    print t
    sph.evolve_model(t)
