"""
CAP - Assignment 4
Lennart van Sluijs, Marco TvdB
"""

import sys
import numpy as np
from matplotlib import pyplot as plt
from amuse.lab import *
from amuse.plot import plot
from amuse.ext.sph_to_star import convert_SPH_to_stellar_model
import matplotlib.animation as animation

from amuse.ext.orbital_elements import orbital_elements_from_binary
from amuse.community.fractalcluster.interface import new_fractal_cluster_model

from distinct_colours import *
from prepare_figure import *

def return_evolved_star_hydro(mass, time, Nsph):
    star =  Particle(mass=mass)
    stellar = EVtwin()
    star = stellar.particles.add_particle(star)
    stellar.evolve_model(time)
    Nsph = Nsph * int(mass.value_in(units.MSun))
    star_in_sph = convert_stellar_model_to_SPH(star, Nsph).gas_particles
    stellar.stop()
    return star_in_sph

def merge_two_stars_sph(Mprim, Msec, Rsec, b, vsec, t_end, dt, Nsph, t_coll):
    
    # primary star - our target
    primary_in_sph = return_evolved_star_hydro(Mprim, t_coll, Nsph)
    primary_in_sph = relax_sph_realization(primary_in_sph)
    
    # neutron star - our bullet
    secondary_in_sph = Particles(1)
    secondary_in_sph.mass = Msec
    secondary_in_sph.radius = Rsec
        
    # combine two a system
    R_primary = primary_in_sph.x.max()
    M = primary_in_sph.mass.sum() + secondary_in_sph.mass.sum()
    secondary_in_sph.x = np.sqrt(4-b**2) * R_primary
    secondary_in_sph.y = b * R_primary
    secondary_in_sph.z = 0. | units.km
    secondary_in_sph.vx = -vsec
    secondary_in_sph.vy = 0. | units.kms
    secondary_in_sph.vz = 0. | units.kms
    secondary_in_sph.u = 0. | units.m**2 / units.s**2
        
    # initalize hydro simulation
    converter=nbody_system.nbody_to_si(Mprim, 1.0|units.AU)
    hydro = Gadget2(converter)
    hydro.gas_particles.add_particles(primary_in_sph)
    hydro.gas_particles.add_particles(secondary_in_sph)
    
    # array to save parameters
    x = [[] for i in range(len(hydro.particles))]
    y = [[] for i in range(len(hydro.particles))]
    vx = [[] for i in range(len(hydro.particles))]
    vy = [[] for i in range(len(hydro.particles))]
    Ekin = [[] for i in range(len(hydro.particles))]
    Epot = [[] for i in range(len(hydro.particles))]
    
    frames = 0
    print ''

	
    # run simulation
    while hydro.model_time<t_end:
        frames += 1
        hydro.evolve_model(hydro.model_time + dt)
        print 'running sph:',hydro.model_time + dt, '\r',
               
        for i in range(len(hydro.particles)):
			x[i].append(hydro.particles[i].x.value_in(units.RSun))
			y[i].append(hydro.particles[i].y.value_in(units.RSun))
			vx[i].append(hydro.particles[i].x.value_in(units.RSun))
			vy[i].append(hydro.particles[i].y.value_in(units.RSun))
			Ekin[i].append(hydro.particles[i].kinetic_energy)
			Epot[i].append(hydro.particles[i].potential_energy)
			
    # check which particles belong to which object
    clumps = find_clumps_with_hop(hydro.particles, converter)
    
    # find the particles that are part of the star
    for clump in clumps:
        clump.scale_to_standard(converter, virial_ratio=0.5)
    plot_single_image(clumps, lim=50)
    
    # caculate relevant parameters
    Nstar = float(len(clumps[0]))
    Nloss = Nsph*Mprim.value_in(units.MSun) - float(len(clumps[0]))
    Mstar = clumps[0].mass.sum()
    Mloss = M  - Mstar
    Rstar = abs(clumps[0].x.max() - clumps[0].x.min())/2.
    xstar = clumps[0].x.mean()
    ystar = clumps[0].y.mean()
    vxstar = clumps[0].vx.mean()
    vystar = clumps[0].vy.mean()
    Ekinstar = clumps[0].kinetic_energy()
    Epotstar = clumps[0].potential_energy()
    Ekin_neutron = -1
    Epot_neutron = -1
    Ekin = hydro.kinetic_energy
    Epot = hydro.potential_energy
    
    # save the parameters

    systempar = np.array([Nstar, Nloss, Mstar.value_in(units.MSun), Mloss.value_in(units.MSun),Rstar.value_in(units.RSun),
    Ekinstar.value_in(units.m**2 * units.kg * units.s**-2), Epotstar.value_in(units.m**2 * units.kg * units.s**-2),
    Ekin.value_in(units.m**2 * units.kg * units.s**-2), Epot.value_in(units.m**2 * units.kg * units.s**-2),
    Ekin_neutron, Epot_neutron,
    x[-1][-1], y[-1][-1], vx[-1][-1], vy[-1][-1], xstar.value_in(units.RSun), ystar.value_in(units.RSun), vxstar.value_in(units.kms), vystar.value_in(units.kms)])
        
    # stop the hydro simulation
    hydro.stop()
    
    return x, y, vx, vy, frames, systempar
    
def plot_single_image(groups_of_particles, lim=5):
    left, width = 0.1, 0.4
    bottom, height = 0.1, 0.4
    bottom_h = left_h = left+width+0.05
    rect_xy = [left, bottom, width, height]
    rect_xz = [left, bottom_h, width, 0.4]
    rect_yz = [left_h, bottom, 0.4, height]

    from distinct_colours import get_distinct
    colors = get_distinct(12)

    fig = plt.figure(figsize=(10,10))

    xy = plt.axes(rect_xy)
    xz = plt.axes(rect_xz)
    yz = plt.axes(rect_yz)
    xy.set_xlabel("X [pc]")
    xy.set_ylabel("Y [pc]")
    xz.set_ylabel("Z [pc]")
    yz.set_xlabel("Z [pc]")

    i = 0
    for group in groups_of_particles:
        x = group.x.value_in(units.RSun)
        y = group.y.value_in(units.RSun)
        z = group.z.value_in(units.RSun)
        xy.scatter(x, y, lw=0, c=colors[min(11, i)], s=8)
        xz.scatter(x, z, lw=0, c=colors[min(11, i)], s=8)
        yz.scatter(z, y, lw=0, c=colors[min(11, i)], s=8)
        i += 1


    save_file = 'FractalClusterHop.png'
    plt.savefig(save_file)
    print '\nSaved figure in file', save_file,'\n'
    #plt.show()
    
def find_clumps_with_hop(particles, unit_converter):

    hop = Hop(unit_converter)
    hop.particles.add_particles(particles)
    hop.calculate_densities()

    mean_density = hop.particles.density.mean() 
    hop.parameters.peak_density_threshold = mean_density
    hop.parameters.saddle_density_threshold = 0.9*mean_density
    hop.parameters.outer_density_threshold = 0.1*mean_density

    hop.do_hop()
    result = [x.get_intersecting_subset_in(particles) for x in hop.groups()]
    hop.stop()

    return result

def relax_sph_realization(sph_star):

    dynamical_timescale = sph_star.dynamical_timescale()
    converter = nbody_system.nbody_to_si(dynamical_timescale, 1|units.RSun)
    hydro = Gadget2(converter, number_of_workers=2)
    hydro.gas_particles.add_particles(sph_star)

    to_hydro = sph_star.new_channel_to(hydro.gas_particles)
    to_framework = hydro.gas_particles.new_channel_to(sph_star)

    ts_factor = 2.5
    t_end = ts_factor * sph_star.dynamical_timescale(mass_fraction=0.9)
    n_steps = ts_factor * 100
    velocity_damp_factor = 1.0 - (ts_factor*2*np.pi)/n_steps
    dt = t_end/float(n_steps)
    time = 0|units.day
    while time < t_end:
        print 'evolving star:', time, '\r',
        time += dt
        hydro.evolve_model(time)
        hydro.gas_particles.velocity = velocity_damp_factor * hydro.gas_particles.velocity
    to_framework.copy()
    hydro.stop()
    return sph_star



def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("-M", unit=units.MSun,
                      dest="Mprim", type="float",default = 10|units.MSun,
                      help="Mass of the primary star [%default] MSun")
    result.add_option("-m", unit=units.MSun,
                      dest="Msec", type="float",default = 1.3|units.MSun, # a neutron star has a mass of about
                      help="Mass of the neutron star [%default] MSun") # 1.3 MSun
                      # see: http://iopscience.iop.org/article/10.1088/0004-637X/757/1/55/pdf
    result.add_option("-r", unit=units.MSun,
                  dest="Rsec", type="float",default = 10 |units.km, # a neutron star has a radius of about
                  help="Radius of the neutron star [%default] km") # 10 km
                  # see: Neutron stars, wikipedia
    result.add_option("-b",
              dest="b", type="float",default = 0, # b is impact parameter, default is head-on collision. In terms of R_primary
              help="Impact paimport matplotlib.animation as animationrameter [%default]")
    result.add_option("-v", unit=units.kms,
              dest="vsec", type="float",default = 1000 |units.kms, # inital velocity of 1000 km/s
              help="Velocity of the secondary star [%default] km/s")
    result.add_option("-T", unit=units.kms,
          dest="t_end", type="float",default = 25 |units.hour, # initally run for 2 hours
          help="Duration of simulation [%default] h")
    result.add_option("-n",
          dest="dt", type="float",default = 0.05 | units.hour,
          help="Amount of time before output is generated.")
    result.add_option("-N", 
                      dest="Nsph", type="int",default = 100,
                      help="Number of sph particles per MSun [%default]")
    result.add_option("-t", unit=units.Myr, 
                      dest="t_coll", type="float", default = 0.01|units.Myr,
                      help="end time of the simulation [%default] Myr")
    return result






def update(t):
	# Optionally. Clear axes and reset limits
	plt.gca().cla() 
	ax.set_xlim(-25,15)
	ax.set_ylim(-25,15)
	
	for i in range(len(x)):
		ax.scatter(x[i][t], y[i][t], marker="o", c='b', s=1)
	ax.scatter(x[-1][t], y[-1][t], marker="o", c='k', s=15)
	
	title = str(t*0.05) + ' hours'
	ax.set_title(title)





if __name__ in ('__main__', '__plot__'):
    
     
    # run a single simulation
    x, y, vx, vy, frames, par = merge_two_stars_sph(Mprim, Msec, Rsec, b, vsec, t_end, dt, Nsph, t_coll) 
  
    print '\nSimulation ended!!'	

    # Set up formatting for the movie files
    Writer = animation.writers['ffmpeg']
    writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)


    fig, ax = plt.subplots()
    ani = animation.FuncAnimation(fig, update, frames=frames, interval=60)
    ani.save('movie_b=0_v=1500kms.mp4', writer=writer)
    plt.show()

    print 'Done!!'


