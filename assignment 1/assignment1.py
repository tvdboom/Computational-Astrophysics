import numpy as np
import matplotlib.pyplot as plt
from amuse.lab import *


#other possible values for v1 and v2
#v1=0.2869236336, v2=0.0791847624
#v1=0.3420307307, v2=0.1809369236
#v1=0.3697718457, v2=0.1910065395

def make_three_bodies(m1=1.0, m2=1.0, m3=0.5, v1=0.2869236336, v2=0.0791847624):
	
	position = [0 for i in range(3)]
	velocity = [0 for i in range(3)]
	position[0] = [-1, 0, 0] | nbody_system.length
	position[1] = [1, 0, 0] | nbody_system.length
	position[2] = [0, 0, 0] |nbody_system.length
	
	velocity[0] = [v1, v2, 0] | nbody_system.length/nbody_system.time
	velocity[1] = [v1, v2, 0] | nbody_system.length/nbody_system.time
	velocity[2] = [-2*v1/m3, -2*v2/m3, 0] | nbody_system.length/nbody_system.time

	
	#create new particles
	bodies = Particles(3)
	bodies.mass = [m1, m2, m3] | nbody_system.mass
	bodies.radius = [0, 0, 0] |nbody_system.length
	bodies.position = position
	bodies.velocity = velocity

	return bodies



#initialize integrator
gravity = Hermite()

#insert particles
bodies = make_three_bodies()
gravity.particles.add_particles(bodies)

	
#declare some variables for the simulation
time = 100 | nbody_system.time  #time to be integrated
t = 0. | nbody_system.time		#actual integrated time
dt = 0.05 | nbody_system.time	#steps of integration

#save positions to plot orbits later
x = [[] for i in range(int(time/dt))] 	
y = [[] for i in range(int(time/dt))] 	

#start simulation
while t <= time:
	
	gravity.evolve_model(t)
	t += dt
	print 'integrating t = -', t, '/', time,			'\r',
	
	for i in range(3):
		x[i].append(gravity.particles[i].x.value_in(nbody_system.length))
		y[i].append(gravity.particles[i].y.value_in(nbody_system.length))
	
	
	
#stop integrator
print '\nDone integrating!\n'
gravity.cleanup_code()
gravity.stop()


# Make multiple plots for a video

color = ['k','b','r']
steps = 400
snap = 5.

for j in range(steps):
	#make plot every snap timesteps
	if j/snap == int(j/snap):		
		fig, ax = plt.subplots(figsize=(40, 40))
		
		print 'Making plots... :', int(j/snap)+1,'/', int(steps/snap),			'\r',
		for i in range(3):
			plt.plot(x[i],y[i],c=color[i], lw=1, alpha=1, zorder=1)
			plt.scatter(x[i][j], y[i][j], c=color[i], s=250, zorder=2)
		
		plt.tight_layout()
		ax.axes.get_xaxis().set_visible(False)
		ax.axes.get_yaxis().set_visible(False)
		plt.savefig('./Plots/p3/image%.4d.png'%int(j/snap))

print '\nDone!'



