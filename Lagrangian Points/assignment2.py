"""
Homework assingment 2 CAP - an orbit through the
Lagrange point of two bodies

Lennart van Sluijs & Marco TvdB
s1376861 & s1246984
"""

# import packages
import matplotlib.pyplot as plt
from amuse.plot import xlabel, ylabel, effective_iso_potential_plot
import numpy as np

# import relevant amuse modules
from amuse.lab import *
from amuse.units import units, constants, nbody_system
from amuse.community.hermite0.interface import Hermite
from amuse.datamodel import Particles

def Earth_Moon_system():
	"""
	Create the Earth - Moon system.
	"""

	# parameters of the Earth Moon system
	distance = 384400 # distance between two main bodies in km
	Mearth = 5.97219 * 10**24 #Earth mass in kg
	M1 = Mearth # Earth mass
	M2 = 0.012*Mearth # Lunar mass
	vorbit = 1.022 # orbital speed of Moon around the Earth in km/s

	particles = Particles(2)
	particles.mass = [M1,  M2] | units.kg
	particles.position = [[0, 0, 0], [distance, 0, 0]] | units.km
	particles.velocity = [[0, 0, 0], [0, vorbit, 0]] | units.kms

	return particles
    
def create_testparticle(position, velocity):
	"""
	Create number test particles with opposite x velocities starting at
	approximately the same point (x0,0).

	DEFAULT PARAMETERS:
	x0 = 349931 the estimated Earth - Moon Lagrange point
	vorbit = lunar orbital velocity
	"""

	particles = Particles(1)
	particles.mass = 0 | units.kg
	particles.radius = 0 | units.m
	particles.position = position
	particles.velocity = velocity
	
	return particles
    
def integrate_orbit(t_tot, gravity):
	"""
	Integrate the orbits. t_tot is total integration time in days, gravity
	is the integrator.
	"""
	
	#declare some variables for the simulation
	time = t_tot  	| units.day #time to be integrated
	t = 0.	  		| units.day	#actual integrated time
	dt = 0.01 		| units.day	#steps of integration

	#save positions to plot orbits later
	x = [[] for i in range(len(gravity.particles))]	
	y = [[] for i in range(len(gravity.particles))]	

	#start simulation
	E_begin = gravity.kinetic_energy + gravity.potential_energy
	
	while t <= time:
		
		gravity.evolve_model(t)
		t += dt
		print 'integrating t = ', t, '/', time		,'\r',

		for i in range(len(gravity.particles)):
			x[i].append(gravity.particles[i].x.value_in(units.km))
			y[i].append(gravity.particles[i].y.value_in(units.km))
	
	E_end = gravity.kinetic_energy + gravity.potential_energy
	
	return x,y, E_begin, E_end

def make_effective_iso_potential_plot(gravity):
	"""
	This function is copied from the amuse documentation.
	The last part has been altered to also plot our Moon Earth system
	and the test particles orbit.
	"""
	
	omega = (constants.G * gravity.particles.total_mass() /
	(5.0|units.AU**3)).sqrt()
	center_of_mass = gravity.particles.center_of_mass()[:]
	plt.rcParams.update({'font.size': 30})
	figure = plt.figure(figsize = (12, 12))
	ax = plt.gca()
	ax.get_yaxis().get_major_formatter().set_useOffset(False)
	ax.minorticks_on() 

	current_axes = plt.subplot(1, 1, 1)
	current_axes.set_aspect("equal", adjustable = "box")
	lim = 6e5
	potential = effective_iso_potential_plot(gravity,
								 omega,
								 xlim = [-lim, lim] | units.km,
								 ylim = [-lim, lim] | units.km,
								 center_of_rotation = center_of_mass,
								 fraction_screen_filled=0.85)
                                 

	# Simulation towards the Moon
	position = [329931, 0, 0] 	| units.km # position of the lagrange point
	velocity = [0.14, 1.022, 0] | units.kms # these parameters worked well
	gravity.particles.add_particles(create_testparticle(position, velocity))
	x, y, E_begin, E_end = integrate_orbit(1.85, gravity) # towards the moon 1.85 days is sufficent

	#initialize integrator
	converter = nbody_system.nbody_to_si(1.0|units.MSun, 1.0|units.AU)
	gravity = Hermite(converter)
	gravity.particles.add_particles(Earth_Moon_system())
	gravity.particles.add_particles(create_testparticle(position, velocity))
	gravity.particles.velocity = gravity.particles.velocity * -1.
	x2, y2, E_begin2, E_end2 = integrate_orbit(10, gravity) # towards the Earth 10 days is sufficient
	
	# make numpy arrays
	x = np.array(x)
	y = np.array(y)
	x2 = np.array(x2)
	y2 = np.array(y2)

	# Make plot
	plt.scatter([0],[0], color = 'b', s=25, label = 'Earth')
	plt.plot(x[1],y[1],color='r',lw=2, label='Moon')
	plt.plot(x[2],y[2], color = 'k', lw=2,label ='Test particle')
	plt.plot(x2[1],y2[1],color='r',lw=2)
	plt.plot(x2[2],y2[2], color = 'k', lw=2)
	plt.legend(fontsize = 12)

	# Plot langrangian point
	plt.scatter(329931, 0, marker='+', s=30, zorder=2, c='m')
	plt.text(329931, -31306, "$L_1$", size=20, zorder =3, color='m')
	
	# Stop integrator
	gravity.cleanup_code()
	gravity.stop()
	
	# Label and show
	xlabel('x')
	ylabel('y')
	plt.ticklabel_format(style='sci', axis='both', scilimits=(0,0))
	plt.tight_layout()
	plt.savefig("lagrange_points.png")
	plt.show() # show figure
	
# --------------------------- Main -------------------------------------

#initialize integrator
converter = nbody_system.nbody_to_si(1.0|units.MSun, 1.0|units.AU)
gravity = Hermite(converter)

#add Moon and Earth
gravity.particles.add_particles(Earth_Moon_system())
make_effective_iso_potential_plot(gravity)
