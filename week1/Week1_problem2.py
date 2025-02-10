
import numpy as np
import matplotlib.pyplot as plt

T=1
g = 1
G = 1
Gr = 1
Dr = g / (Gr**2)



num_particles = 10000
num_steps = 10000
tmax = 10
times = np.linspace(0,tmax,num_steps)
dt = times[1]-times[0]



sigma_x = np.sqrt(2*g*dt / G**2)
sigma_phi = np.sqrt(2*g*dt / Gr**2)


xs = np.zeros([num_steps, num_particles])
phis = np.zeros_like(xs)

# phis[0] = np.pi/4 *np.ones(num_particles)
drift = T / (G*Dr) * np.cos(phis[0,0])
print(drift)

for i in range(1, num_steps):
    mean_x =  T/(G*Dr)*np.cos(phis[i-1])*(1-np.exp(-Dr*dt))
    alpha = np.random.normal(loc = mean_x, scale = sigma_x, size=[num_particles])
    beta = np.random.normal(scale = sigma_phi, size=[num_particles])
    xs[i] = xs[i-1] +  alpha 
    phis[i] = phis[i-1] + beta









for i in range(5):
    plt.plot(times, xs[:,i], alpha = 0.5)
plt.plot(times, np.mean(xs, axis=1), color='k', label='<x(t)>')
plt.plot(times, drift*np.ones_like(times), color='grey', linestyle=':', label=r'$ \frac {T}{\Gamma D_r} cos(\phi_0) $')
plt.legend()
plt.xlabel('t')
plt.ylabel('x(t)')
plt.savefig('xvst.png', dpi=600)