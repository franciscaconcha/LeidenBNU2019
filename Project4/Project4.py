from amuse.lab import *
from amuse.ext.orbital_elements import orbital_elements_from_binary
from amuse.community.fractalcluster.interface import new_fractal_cluster_model


N = 2500  # Number of stars
Rvir = 0.5  # Virial radius
Qvir = 0.5  # Viral ratio
t_end = 1 | units.Myr

masses = new_kroupa_mass_distribution(N, 100|units.MSun)  # Kroupa mass distribution with maximum mass 100 MSun

# Convert n-body system units to our relevant physical units
converter=nbody_system.nbody_to_si(masses.sum(),Rvir)  

# New Plummer sphere distribution
bodies = new_plummer_model(N, converter)
bodies.scale_to_standard(converter, virial_ratio=Qvir)

# Asign the Kroupa distribution as stellar masses
bodies.stellar_mass = masses

# Total mass of the particles in the code is stellar mass + disk mass
bodies.disk_mass = 0.01 * bodies.stellar_mass
bodies.mass = bodies.stellar_mass + bodies.disk_mass

# Assigning disk radius
bodies.disk_radius = 400 | units.AU

# Assigning collisional radius for particle
bodies.radius = 10 * bodies.disk_radius

# Start gravity code
gravity = ph4(converter, number_of_workers=2)
gravity.parameters.epsilon_squared = (100|units.AU)**2

# Add particles to gravity code
gravity.particles.add_particles(bodies)

# Communication channels
channel_from_gravity = gravity.particles.new_channel_to(bodies)
channel_to_gravity = bodies.new_channel_to(gravity.particles)

# Enable the "collision detection" stopping condition
stopping_condition = gravity.stopping_conditions.collision_detection
stopping_condition.enable()

# Save data
write_set_to_file(bodies.savepoint(0|units.Myr), filename, 'hdf5',
                      append_to_file=False)
    
# Energy of the system
Etot_init = gravity.kinetic_energy + gravity.potential_energy
Etot_prev = Etot_init

# Evolve system
dt = t_end/10.
time = 0 | units.yr

while gravity.model_time < t_end:
    time += dt
    gravity.evolve_model(time)

    while stopping_condition.is_set():
        channel_from_gravity.copy()
        Ek_enc = gravity.kinetic_energy 
        Ep_enc = gravity.potential_energy
        for ci in range(len(stopping_condition.particles(0))): 
            bodies_in_enc \
                = Particles(particles=[stopping_condition.particles(0)[ci],
                                       stopping_condition.particles(1)[ci]])
            local_bodies_in_enc \
                = bodies_in_enc.get_intersecting_subset_in(bodies)

            # This is the function that updated truncated disk radii and masses
            resolve_close_encounter(gravity.model_time, local_bodies_in_enc)
            
            print "At time=", gravity.model_time.value_in(units.Myr), \
                "Rdisk=", local_bodies_in_enc.disk_radius.in_(units.AU)
            channel_to_gravity.copy_attributes(["radius"])

    gravity.evolve_model(time)
    channel_to_gravity.copy_attributes(["mass"])
        Etot = gravity.kinetic_energy + gravity.potential_energy
        print "T=", gravity.model_time, 
        print "E= ", Etot, "Q= ", \
              gravity.kinetic_energy/gravity.potential_energy
        print "dE=", (Etot-Etot_init)/Etot, "ddE=", (Etot-Etot_prev)/Etot
        Etot_init -= (Etot_prev-Etot)
        Etot_prev = Etot

    gravity.stop()

