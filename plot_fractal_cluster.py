import numpy
from amuse.lab import *
from amuse.ext.orbital_elements import orbital_elements_from_binary
from amuse.community.fractalcluster.interface import new_fractal_cluster_model
from matplotlib import pyplot
from distinct_colours import get_distinct

def plot_single_image(particles, lim=1):
    """ Script to plot a single frame of a set of particles.

        :param particles: set of particles to plot
        :param lim: limits for the plot (parsec)
    """
    colors = get_distinct(2)

    fig = pyplot.figure(figsize=(10,10))

    pyplot.xlabel("X [pc]")
    pyplot.ylabel("Y [pc]")

    x = particles.x.value_in(units.parsec)
    y = particles.y.value_in(units.parsec)

    pyplot.scatter(x, y, lw=0, c=colors[0], s=8)  # Plot all stars
    pyplot.scatter(x[x < 0.5], y[x < 0.5], lw=0, c=colors[1], s=8)  # Plot stars with x<0.5 in a different color

    pyplot.xlim((-lim, lim))
    pyplot.ylim((-lim, lim))

    save_file = 'FractalCluster.png'
    pyplot.savefig(save_file)
    print '\nSaved figure in file', save_file,'\n'
    pyplot.show()


N = 2000  # Number of particles (stars) for the cluster
Rvir = 0.5 | units.parsec  # Virial radius of the cluster
Qvir = 0.5  # Virial ratio of the cluster (Q=0.5 == equilibrium)
Fd = 1.6  # Fractal dimension
numpy.random.seed(1234)  # Random seed for star positions

masses = new_kroupa_mass_distribution(N)  # Stellar masses
converter = nbody_system.nbody_to_si(masses.sum(), Rvir)  # Converter, turn n-body system units into relevant physical units

# Create new fractal cluster
bodies = new_fractal_cluster_model(N=N, fractal_dimension=Fd,
                                   random_seed=1234,
                                   convert_nbody=converter)

bodies.mass = masses  # Assign stellar masses to cluster particles
bodies.move_to_center()
bodies.scale_to_standard(converter, virial_ratio=Qvir)

plot_single_image(bodies, lim=1)
