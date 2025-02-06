#Python code for Viscek model
#Rupesh Mahore, Gal Ross, Capucine Beraud IISc
#Feb, 2025

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation


N = 300             # Number of particles
L = 10.0            # Domain size (assumed square domain LxL)
v = 0.03            # Particle speed
eta = 0.1            # Noise amplitude (in radians)
r = 0.5             # Interaction radius
dt = 1.0            # Time step
num_steps = 800     # Number of simulation steps


# Defining initial random positions and angles for all the particles
positions = np.random.rand(N, 2) * L
angles = np.random.rand(N) * 2 * np.pi

def apply_periodic_boundary(positions, L):
    """
    Ensures that the particle positions lie between 0 and L.
    """
    return positions % L

def find_neighbors(i, positions, r, L):
    """
    returns the locations of all particles that are less than r away from the current particle.
    
    Params:
        i:          index of current particle
        positions:  list of particle coordinates
        r:          maximum distance
        L:          grid length
        
    Returns:
        neighbors:  list of neighbor locations 
    """
    diff = positions - positions[i]
    diff = (diff + L/2) % L - L/2
    distances = np.sqrt(np.sum(diff**2, axis=1))
    neighbors = np.where(distances < r)[0]
    return neighbors

def update_positions(positions, angles, L, v, eta, r):
    """
    returns the locations of all particles that are less than r awaz from the current particle.
    
    Params:
        i:          index of current particle
        positions:  list of particle coordinates
        r:          maximum distance
        L:          grid length
        
    Returns:
        neighbors:  list of neighbor locations 
    """
    new_positions = np.copy(positions)
    new_angles = alignment(positions, angles, L, v, eta, r)
        
    new_positions[:,0] += v * np.cos(new_angles) * dt
    new_positions[:,1] += v * np.sin(new_angles) * dt
    
    new_positions = apply_periodic_boundary(new_positions, L)
    
    return new_positions, new_angles

def alignment(positions, angles, L, v, eta, r):
    """
    computes the angles after alignment with neighboring particles.
    
    Params:
        positions:  list of particle coordinates
        angles:  list of particle angles
        L:          grid length
        v:          driving speed
        eta:        noise strength
        r:          maximum distance of interaction
        
    Returns:
        neighbors:  list of neighbor locations 
    """
    new_angles = np.zeros_like(angles)
    for i in range(len(positions)):
        neighbors = find_neighbors(i, positions, r, L)
        influence = np.angle(np.sum(np.exp(1.j * angles[neighbors])))
        new_angles[i] = influence
    new_angles += np.random.normal(0,eta, len(positions))
    new_angles %= 2*np.pi
    return new_angles

# Set up the plot for animation using quiver
fig, ax = plt.subplots()
ax.set_xlim(0, L)
ax.set_ylim(0, L)
ax.set_aspect('equal')
ax.set_title('Vicsek Model with Quiver Visualization')

# Initialize quiver plot. U and V components based on initial angles.
U = np.cos(angles)
V = np.sin(angles)
Q = ax.quiver(positions[:,0], positions[:,1], U, V, angles='xy', scale_units='xy', scale=8, color='blue')

# Add a text object to display the frame number
frame_text = ax.text(0.05, 0.95, f'Frame: 0', transform=ax.transAxes, fontsize=12, verticalalignment='top', horizontalalignment='left', color='black')


def animate(frame):
    """
    updates the particle locations / angles then updates the plot at the frame
    
    Params:
        frame:      current timestep
        
    Returns:
        Q:          updated list of direction vectors
        frame_text: text object in animation with current timestep 
    """
    global positions, angles
    prev_positions = positions.copy()
    positions, angles = update_positions(positions, angles, L, v, eta, r)
    # order = 1/(v*dt) * np.linalg.norm(np.mean(positions - prev_positions, axis=0))
    diff = (positions - prev_positions + L/2) % L - L/2
    order = 1/(v*dt) * np.mean(diff, axis=0)
    print(order)
    # Update velocity components for quiver plot
    U = np.cos(angles)
    V = np.sin(angles)
    
    # Update quiver plot data
    Q.set_offsets(positions)
    Q.set_UVC(U, V)

    # Update the frame number text
    frame_text.set_text(f'Frame: {frame+1}')
    
    # print(frame, np.sum(positions[:,0] > 9.5))

    return Q,frame_text

ani = animation.FuncAnimation(fig, animate, frames=num_steps, interval=50, blit=True)
plt.show()

writervideo = animation.FFMpegWriter(fps=60) 
ani.save('TestsViscek_noise_%0.1f.mp4' % eta, writer=writervideo) 

plt.close() 


