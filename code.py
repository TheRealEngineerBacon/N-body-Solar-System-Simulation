import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from mpl_toolkits.mplot3d.art3d import Line3DCollection
import numpy as np
from numba import njit

#SETUP STAGE
N = 5 #number of particles, Sun (N[0]) plus N - 1 planetary bodies
T = 5e9 #total time for simulation
t = 4e4 #time step (s)
G = 6.67430e-11

#Create iterative list for attribute calculations.
#timestep, x, y, z, vx, vy, vz, ax, ay, az, mass
particlelist = np.zeros((int(T/t), 11, N), dtype=np.float64)

#Fill timestep columns for each particle for reference.
for particle in range(0, N):
    for timestep in range(0, int(T/t)):
        particlelist[timestep][0][particle] = t

#Calculate mass values for star and  and N - 1 system bodies.
star_mass = np.random.uniform(2e29, 2e31)
for step in range(0, int(T/t)):
    particlelist[step][10][0] = star_mass
for particles in range(1, N):
    body_mass = np.random.uniform(3.285e20, 1.898e27)
    for step in range (0, int(T/t)):
        particlelist[step][10][particles] = body_mass

#Calculate random x, y, and z positions.
particlelist[0][1][0] = 0
particlelist[0][2][0] = 0
particlelist[0][3][0] = 0
for positions in range(1, N):
    radius = np.random.uniform(5.79e11, 4.4951e12)
    angle = np.random.uniform(0, 2*np.pi)
    z_angle = np.random.uniform((np.pi/2) - 0.122348, (np.pi/2) + 0.122348)
    particlelist[0][1][positions] = np.cos(angle) * radius
    particlelist[0][2][positions] = np.sin(angle) * radius
    particlelist[0][3][positions] = radius * -np.cos(z_angle)

#Calculate random velocity components for x and y axes.
particlelist[0][4][0] = 0
particlelist[0][5][0] = 0
particlelist[0][6][0] = 0
for i in range(1, N):
    distance = np.sqrt(((particlelist[0][2][0] - particlelist[0][2][i])**2 + (particlelist[0][1][0] - particlelist[0][1][i])**2 + (particlelist[0][3][0] - particlelist[0][3][i])**2))
    vel = np.sqrt((G * particlelist[0][10][0]) / distance) * np.random.uniform(1.00, 1.10)
    theta = np.arctan((particlelist[0][2][0] - particlelist[0][2][i])/(particlelist[0][1][0] - particlelist[0][1][i]))
    incl_angle = np.random.uniform(-0.122348, 0.122348)
    if particlelist[0][3][i] > particlelist[0][3][0]:
        if particlelist[0][1][i] < particlelist[0][1][0] and particlelist[0][2][i] > particlelist[0][2][0]:
            particlelist[0][4][i] = vel * np.sin(theta)
            particlelist[0][5][i] = -1 * vel * np.cos(theta)
            particlelist[0][6][i] = np.sqrt(particlelist[0][4][i]**2 + particlelist[0][5][i]**2) * np.tan(incl_angle)
            
        elif particlelist[0][1][i] < particlelist[0][1][0] and particlelist[0][2][i] < particlelist[0][2][0]:
            particlelist[0][4][i] = vel * np.sin(theta)
            particlelist[0][5][i] = -1 * vel * np.cos(theta)
            particlelist[0][6][i] = np.sqrt(particlelist[0][4][i]**2 + particlelist[0][5][i]**2) * np.tan(incl_angle)
        
        elif particlelist[0][1][i] > particlelist[0][1][0] and particlelist[0][2][i] > particlelist[0][2][0]:
            particlelist[0][4][i] = -1 * vel * np.sin(theta)
            particlelist[0][5][i] = vel * np.cos(theta)
            particlelist[0][6][i] = np.sqrt(particlelist[0][4][i]**2 + particlelist[0][5][i]**2) * np.tan(incl_angle)
    
        elif particlelist[0][1][i] > particlelist[0][1][0] and particlelist[0][2][i] < particlelist[0][2][0]:
            particlelist[0][4][i] = -1 * vel * np.sin(theta)
            particlelist[0][5][i] = vel * np.cos(theta)
            particlelist[0][6][i] = np.sqrt(particlelist[0][4][i]**2 + particlelist[0][5][i]**2) * np.tan(incl_angle)
            
    elif particlelist[0][3][i] < particlelist[0][3][0]:
        if particlelist[0][1][i] < particlelist[0][1][0] and particlelist[0][2][i] > particlelist[0][2][0]:
            particlelist[0][4][i] = vel * np.sin(theta)
            particlelist[0][5][i] = -1 * vel * np.cos(theta)
            particlelist[0][6][i] = np.sqrt(particlelist[0][4][i]**2 + particlelist[0][5][i]**2) * np.tan(incl_angle)
            
        elif particlelist[0][1][i] < particlelist[0][1][0] and particlelist[0][2][i] < particlelist[0][2][0]:
            particlelist[0][4][i] = vel * np.sin(theta)
            particlelist[0][5][i] = -1 * vel * np.cos(theta)
            particlelist[0][6][i] = np.sqrt(particlelist[0][4][i]**2 + particlelist[0][5][i]**2) * np.tan(incl_angle)
        
        elif particlelist[0][1][i] > particlelist[0][1][0] and particlelist[0][2][i] > particlelist[0][2][0]:
            particlelist[0][4][i] = -1 * vel * np.sin(theta)
            particlelist[0][5][i] = vel * np.cos(theta)
            particlelist[0][6][i] = np.sqrt(particlelist[0][4][i]**2 + particlelist[0][5][i]**2) * np.tan(incl_angle)
    
        elif particlelist[0][1][i] > particlelist[0][1][0] and particlelist[0][2][i] < particlelist[0][2][0]:
            particlelist[0][4][i] = -1 * vel * np.sin(theta)
            particlelist[0][5][i] = vel * np.cos(theta)
            particlelist[0][6][i] = np.sqrt(particlelist[0][4][i]**2 + particlelist[0][5][i]**2) * np.tan(incl_angle)


#Calculate initial acceleration values for each particle.
for n in range(0, N): #n
    forces = np.ndarray((N, 3))
    for i in range(0, N): #i
        if i == n:
            forces[i, :] = 0
        else:
            distance = np.sqrt((particlelist[0][1][n] - particlelist[0][1][i])**2 + (particlelist[0][2][n] - particlelist[0][2][i])**2 + (particlelist[0][3][n] - particlelist[0][3][i])**2)
            force = (G * particlelist[0][10][i] * particlelist[0][10][n]) / distance**2
            theta = np.arctan((particlelist[0][2][n] - particlelist[0][2][i])/(particlelist[0][1][n] - particlelist[0][1][i]))
            z_comp = abs((particlelist[0][3][i] - particlelist[0][3][n])/ distance)
            step = 0
            if particlelist[step][3][i] > particlelist[step][3][n]:
                if particlelist[step][1][i] < particlelist[step][1][n] and particlelist[step][2][i] > particlelist[step][2][n]:
                    forces[i, :] = [-1 * force * np.cos(theta), -1 * force * np.sin(theta), force * np.cos(z_comp)]
                                    
                elif particlelist[step][1][i] < particlelist[step][1][n] and particlelist[step][2][i] < particlelist[step][2][n]:
                    forces[i, :] = [-1 * force * np.cos(theta), -1 * force * np.sin(theta), force * z_comp]
                    
                elif particlelist[step][1][i] > particlelist[step][1][n] and particlelist[step][2][i] > particlelist[step][2][n]:
                    forces[i, :] = [force * np.cos(theta), force * np.sin(theta), force * z_comp]
                    
                elif particlelist[step][1][i] > particlelist[step][1][n] and particlelist[step][2][i] < particlelist[step][2][n]:
                    forces[i, :] = [force * np.cos(theta), force * np.sin(theta), force * z_comp]
                                        
            elif particlelist[step][3][i] < particlelist[step][3][n]:
                if particlelist[step][1][i] < particlelist[step][1][n] and particlelist[step][2][i] > particlelist[step][2][n]:
                    forces[i, :] = [-1 * force * np.cos(theta), -1 * force * np.sin(theta), -1 * force * z_comp]
                                    
                elif particlelist[step][1][i] < particlelist[step][1][n] and particlelist[step][2][i] < particlelist[step][2][n]:
                    forces[i, :] = [-1 * force * np.cos(theta), -1 * force * np.sin(theta), -1 * force * z_comp]
                                
                elif particlelist[step][1][i] > particlelist[step][1][n] and particlelist[step][2][i] > particlelist[step][2][n]:
                    forces[i, :] = [force * np.cos(theta), force * np.sin(theta), -1 * force * z_comp]
                    
                elif particlelist[step][1][i] > particlelist[step][1][n] and particlelist[step][2][i] < particlelist[step][2][n]:
                    forces[i, :] = [force * np.cos(theta), force * np.sin(theta), -1 * force * z_comp]
                    
    particlelist[0][7][n] = np.sum(forces[:, 0]) / particlelist[0][10][n]
    particlelist[0][8][n] = np.sum(forces[:, 1]) / particlelist[0][10][n]
    particlelist[0][9][n] = np.sum(forces[:, 2]) / particlelist[0][10][n]




#=============================================================================
#Evaluation Stage

@njit
def calc_vel(step, n, i):
    calculated_vel = (particlelist[step-1][i+3][n] * particlelist[step][0][n]) + particlelist[step-1][i][n]
    return calculated_vel

@njit
def calc_poss(step, n, i):
    calculated_poss = particlelist[step-1][i][n] + (particlelist[step-1][i+3][n] * particlelist[step][0][n]) + (0.5 * particlelist[step-1][i+6][n] * (particlelist[step][0][n]**2))
    return calculated_poss

@njit
def sum_forces(step, n):
    forces = np.zeros((N, 3), dtype=np.float64)
    for i in range(0, N): #Calculate force from all points and compute acceleration.
        if i == n:
            forces[i, :] = 0
        else:
            distance = np.sqrt((particlelist[step][1][n] - particlelist[step][1][i])**2 + (particlelist[step][2][n] - particlelist[step][2][i])**2 + (particlelist[step][3][n] - particlelist[step][3][i])**2)
            force = (G * particlelist[step][10][i] * particlelist[step][10][n]) / distance**2
            theta = np.arctan((particlelist[step][2][n] - particlelist[step][2][i])/(particlelist[step][1][n] - particlelist[step][1][i]))
            z_comp = abs((particlelist[step][3][i] - particlelist[step][3][n]) / distance)
                
            if particlelist[step][3][i] > particlelist[step][3][n]:
                if particlelist[step][1][i] < particlelist[step][1][n] and particlelist[step][2][i] > particlelist[step][2][n]:
                    forces[i, :] = [-1 * force * np.cos(theta), -1 * force * np.sin(theta), force * z_comp]
                
                elif particlelist[step][1][i] < particlelist[step][1][n] and particlelist[step][2][i] < particlelist[step][2][n]:
                    forces[i, :] = [-1 * force * np.cos(theta), -1 * force * np.sin(theta), force * z_comp]
                
                elif particlelist[step][1][i] > particlelist[step][1][n] and particlelist[step][2][i] > particlelist[step][2][n]:
                    forces[i, :] = [force * np.cos(theta), force * np.sin(theta), force * z_comp]

                elif particlelist[step][1][i] > particlelist[step][1][n] and particlelist[step][2][i] < particlelist[step][2][n]:
                    forces[i, :] = [force * np.cos(theta), force * np.sin(theta), force * z_comp]
                    
            elif particlelist[step][3][i] < particlelist[step][3][n]:
                if particlelist[step][1][i] < particlelist[step][1][n] and particlelist[step][2][i] > particlelist[step][2][n]:
                    forces[i, :] = [-1 * force * np.cos(theta), -1 * force * np.sin(theta), -1 * force * z_comp]
                
                elif particlelist[step][1][i] < particlelist[step][1][n] and particlelist[step][2][i] < particlelist[step][2][n]:
                    forces[i, :] = [-1 * force * np.cos(theta), -1 * force * np.sin(theta), -1 * force * z_comp]
                
                elif particlelist[step][1][i] > particlelist[step][1][n] and particlelist[step][2][i] > particlelist[step][2][n]:
                    forces[i, :] = [force * np.cos(theta), force * np.sin(theta), -1 * force * z_comp]

                elif particlelist[step][1][i] > particlelist[step][1][n] and particlelist[step][2][i] < particlelist[step][2][n]:
                    forces[i, :] = [force * np.cos(theta), force * np.sin(theta), -1 * force * z_comp]
    return forces

#Over T with interval t, calculate v_f values for x and y.
for step in range(1, int(T/t)):
    for n in range(0, N):
        #Calculate new velocities
        # particlelist[step][4][n] = (particlelist[step-1][7][n] * particlelist[step][0][n]) + particlelist[step-1][4][n]
        # particlelist[step][5][n] = (particlelist[step-1][8][n] * particlelist[step][0][n]) + particlelist[step-1][5][n]
        # particlelist[step][6][n] = (particlelist[step-1][9][n] * particlelist[step][0][n]) + particlelist[step-1][6][n]
        for i in range(4, 7):
            particlelist[step][i][n] = calc_vel(step, n, i)
        
        #Calculate new positions.
        # particlelist[step][1][n] = particlelist[step-1][1][n] + (particlelist[step-1][4][n] * particlelist[step][0][n]) + (0.5 * particlelist[step-1][7][n] * (particlelist[step][0][n]**2))
        # particlelist[step][2][n] = particlelist[step-1][2][n] + (particlelist[step-1][5][n] * particlelist[step][0][n]) + (0.5 * particlelist[step-1][8][n] * (particlelist[step][0][n]**2))
        # particlelist[step][3][n] = particlelist[step-1][3][n] + (particlelist[step-1][6][n] * particlelist[step][0][n]) + (0.5 * particlelist[step-1][9][n] * (particlelist[step][0][n]**2))
        for i in range(1, 4):
            particlelist[step][i][n] = calc_poss(step, n, i)
        
        
    for n in range(0, N):
        # forces = np.ndarray((N, 3), dtype=np.float64)
        # for i in range(0, N): #Calculate force from all points and compute acceleration.
        #     if i == n:
        #         forces[i, :] = 0
        #     else:
        #         distance = np.sqrt((particlelist[step][1][n] - particlelist[step][1][i])**2 + (particlelist[step][2][n] - particlelist[step][2][i])**2 + (particlelist[step][3][n] - particlelist[step][3][i])**2)
        #         force = (G * particlelist[step][10][i] * particlelist[step][10][n]) / distance**2
        #         theta = np.arctan((particlelist[step][2][n] - particlelist[step][2][i])/(particlelist[step][1][n] - particlelist[step][1][i]))
        #         z_comp = abs((particlelist[step][3][i] - particlelist[step][3][n]) / distance)
                    
        #         if particlelist[step][3][i] > particlelist[step][3][n]:
        #             if particlelist[step][1][i] < particlelist[step][1][n] and particlelist[step][2][i] > particlelist[step][2][n]:
        #                 forces[i, :] = [-1 * force * np.cos(theta), -1 * force * np.sin(theta), force * z_comp]
                    
        #             elif particlelist[step][1][i] < particlelist[step][1][n] and particlelist[step][2][i] < particlelist[step][2][n]:
        #                 forces[i, :] = [-1 * force * np.cos(theta), -1 * force * np.sin(theta), force * z_comp]
                    
        #             elif particlelist[step][1][i] > particlelist[step][1][n] and particlelist[step][2][i] > particlelist[step][2][n]:
        #                 forces[i, :] = [force * np.cos(theta), force * np.sin(theta), force * z_comp]
    
        #             elif particlelist[step][1][i] > particlelist[step][1][n] and particlelist[step][2][i] < particlelist[step][2][n]:
        #                 forces[i, :] = [force * np.cos(theta), force * np.sin(theta), force * z_comp]
                        
        #         elif particlelist[step][3][i] < particlelist[step][3][n]:
        #             if particlelist[step][1][i] < particlelist[step][1][n] and particlelist[step][2][i] > particlelist[step][2][n]:
        #                 forces[i, :] = [-1 * force * np.cos(theta), -1 * force * np.sin(theta), -1 * force * z_comp]
                    
        #             elif particlelist[step][1][i] < particlelist[step][1][n] and particlelist[step][2][i] < particlelist[step][2][n]:
        #                 forces[i, :] = [-1 * force * np.cos(theta), -1 * force * np.sin(theta), -1 * force * z_comp]
                    
        #             elif particlelist[step][1][i] > particlelist[step][1][n] and particlelist[step][2][i] > particlelist[step][2][n]:
        #                 forces[i, :] = [force * np.cos(theta), force * np.sin(theta), -1 * force * z_comp]
    
        #             elif particlelist[step][1][i] > particlelist[step][1][n] and particlelist[step][2][i] < particlelist[step][2][n]:
        #                 forces[i, :] = [force * np.cos(theta), force * np.sin(theta), -1 * force * z_comp]
        forces = sum_forces(step, n)
                
        particlelist[step][7][n] = np.sum(forces[:, 0]) / particlelist[0][10][n]
        particlelist[step][8][n] = np.sum(forces[:, 1]) / particlelist[0][10][n]
        particlelist[step][9][n] = np.sum(forces[:, 2]) / particlelist[0][10][n]
    print(int(step/(T/t)*100))

#Uncomment to save particlelist array to npy file.
#np.save(rf'C:\Users\DellPC\Desktop\MAPRENDERS\physicsim\testsave{N, T, t}.npy', particlelist)


#=============================================================================
#Display Stage
#Create a single plot of a single time step.
fig = plt.figure(dpi=100)
ax = fig.add_subplot(111, projection='3d')
ax.scatter(particlelist[0][1][:], particlelist[0][2][:], particlelist[0][3][:], 'o', color='black')

trajectories = list()
percentage = 0.01
for particle in range (0, N):
    addition = np.ndarray((int((T/t) * percentage), 3))
    point = 0
    for step in np.arange(0, int(T/t), int((T/t)/((T/t) * percentage))):
        addition[point][0] = particlelist[step][1][particle]
        addition[point][1] = particlelist[step][2][particle]
        addition[point][2] = particlelist[step][3][particle]
        point += 1
    trajectories.append(addition)

#linecolor = [(234/255, 230/255, 229/255, 1)]
linecolor = [(0/255, 0/255, 0/255, 1)]

line_segments = Line3DCollection(trajectories, linewidths=0.5, colors=linecolor)
ax.add_collection(line_segments)

bounds = 0

for particle in range(1, N):
    for column in range (1, 4):
        for row in range(0, int(T/t)):
            if particlelist[row][column][particle] > bounds:
                bounds = particlelist[row][column][particle]

ax.auto_scale_xyz([-1 * bounds, bounds], [-1 * bounds, bounds], [-1 * bounds, bounds])
#ax.set_axis_off()
#ax.set_facecolor('black')
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')

plt.tight_layout()
plt.show()

# #Create and save plot for each time step.
# plt.switch_backend('Agg')
# print(int(T/t))
# for plots in range(0, int(T/t)):
#     fig = plt.figure(dpi=200)
#     ax = fig.gca()
#     ax.plot(particlelist[plots][1][:], particlelist[plots][2][:], 'o', color='white')
#     ax.set_xlim(-1 * bounds, bounds)
#     ax.set_ylim(-1 * bounds, bounds)
#     ax.set_facecolor('black')
#     #plt.axis('off')
#     plt.tight_layout()
#     #plt.show()
#     plt.savefig(fname=rf'C:\Users\DellPC\Desktop\MAPRENDERS\physicsim\iter{plots}.png')
#     plt.close()
#     #plt.clf()
