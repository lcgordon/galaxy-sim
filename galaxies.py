# -*- coding: utf-8 -*-
"""
Created on Sun Apr 19 15:45:05 2020

@author: conta
"""

#%matplotlib inline
# imports -- do your imports here
import numpy as np
import matplotlib.pyplot as plt
import scipy
import scipy.stats
import scipy.optimize
from scipy import signal
import astropy
from astropy import units as u
from scipy import integrate
from vpython import *

#https://github.com/vpython/vpython-jupyter/blob/master/Demos/Stars.ipynb
scene = canvas() # This is needed in Jupyter notebook and lab to make programs easily rerunnable
scene.width = scene.height = 600
scene.title = "Galaxy" # Display text below the 3D graphics:
scene.caption = """Right button drag or Ctrl-drag to rotate "camera" to view scene.
To zoom, drag with middle button or Alt/Option depressed, or use scroll wheel.
  On a two-button mouse, middle is left + right.
Touch screen: pinch/extend to zoom, swipe or two-finger rotate."""


#parameters
Nstars = 200  # change this to have more or fewer stars TOTAL
N1 = Nstars/2
N2 = Nstars - N1
G = 4.3e-30 #grav constant in units of pc/Msun * (pc/s)^2
Msun = 2e30 #kg
Rsun = 1e3 #pc (this is huge but oh well)
L = 1e5 #pc, this is the size scale
gal_rad = 0.5 * L #pc
gal_mass = 10e11 #Msun
star_mass = gal_mass/Nstars #units of Msun
#vsun = 0.8*sqrt(G*Msun/Rsun)

scene.range = 2*L
#box is 2L by 2L
scene.forward = vec(0,0,-L)

xaxis = curve(color=color.gray(0.5), radius=1e3)
xaxis.append(vec(0,0,0))
xaxis.append(vec(L,0,0))
yaxis = curve(color=color.gray(0.5), radius=1e3)
yaxis.append(vec(0,0,0))
yaxis.append(vec(0,L,0))
#zaxis = curve(color=color.gray(0.5), radius=1e8)
#zaxis.append(vec(0,0,0))
#zaxis.append(vec(0,0,L))

Stars = []
star_colors = [color.cyan, color.magenta]


#First Galaxy
psum = vec(0,0,0)
r = np.random.randint(0,gal_rad,Nstars)
phi_degrees = np.random.randint(0,360,Nstars)
phi = np.radians(phi_degrees)

x = r*np.sin(phi)
for i in range(Nstars): 
    if i % 2 == 0:
        x[i] = x[i] - L
    elif i%2 == 1:
        x[i] = x[i] + L
        
y = r*np.cos(phi)

mom_phi = np.radians(90-phi_degrees)

for i in range(Nstars):

    star = sphere(pos=vector(x[i], y[i],0), make_trail=True, retain=50, trail_radius=100)
    star.radius = Rsun
    star.mass = star_mass
    star.momentum = vector(0,0,0)
    #star.momentum = vector(mom_r*np.sin(mom_phi[0]), mom_r*np.cos(mom_phi[0]), 0) * star.mass
    star.color = star.trail_color = star_colors[i % 2]
    Stars.append( star )
    
    
# calculating initial velocities
mass_internal = []
for sk in range(Nstars):
    m = 0
    rad = np.abs(r[sk])
    for sl in range(Nstars):
        if sk==sl: continue
        rad2 = np.abs(r[sl])
        if rad > rad2: 
            m = m + 1
    mass_internal.append(m*star_mass)

v_i = np.zeros(Nstars)
for i in range(Nstars):
    v_i[i] = np.sqrt(G * mass_internal[i]/np.abs(r[i])) #units of km/s


for i in range(Nstars):
    velocity_v = vector(v_i[i]*np.sin(mom_phi[i]), -1* v_i[i]*np.cos(mom_phi[i]), 0)
    #if i % 2 ==0:
     #   velocity_v += vector(5e-12, 0, 0)
    #elif i%2==1:
     #   velocity_v += vector(-5e-12,0,0)
    Stars[i].momentum = velocity_v * Stars[i].mass
    psum = psum + Stars[i].momentum

#make total initial momentum equal zero
for i in range(Nstars):
    Stars[i].momentum = Stars[i].momentum - psum/Nstars



dt = 1e13 #1e12 or 1e13 is best value here
hitlist = []

def computeForces():
    global hitlist, Stars
    hitlist = []
    N = len(Stars)
    for i in range(N):
        si = Stars[i]
        if si is None: continue
        F = vec(0,0,0)
        pos1 = si.pos
        m1 = si.mass
        radius1 = si.radius
        for j in range(N):
            if i == j: continue
            sj = Stars[j]
            if sj is None: continue
            pos2 = sj.pos
            r = sj.pos - pos1
            rmag2 = mag2(r)
            radius2 = sj.radius
            if rmag2 <= (radius1 + radius2)**2:
                r = vec(0,0,0)
             
            F = F + (G*m1*sj.mass/(rmag2**1.5))*r
        si.momentum = si.momentum + F*dt
      

t = 0
v_left = vector(1e-11, 0, 0)
v_right = vector(-1e-11, 0,0)
while t < 1e16:
    #rate(100)
    
    # Compute all forces on all stars
    computeForces() #REMOVE THIS COMMENT IN ORDER TO RUN SIMULATION
    
    # Having updated all momenta, now update all positions
    for i in range(len(Stars)):
        if Stars[i] is None: continue
        
        if i % 2 == 0:
            Stars[i].pos = Stars[i].pos + Stars[i].momentum*(dt/star_mass) + (v_left * dt)
        elif i% 2 == 1:
            Stars[i].pos = Stars[i].pos + Stars[i].momentum*(dt/star_mass) + (v_right * dt)
    #for star in Stars:
     #   if star is None: continue
      #  star.pos = star.pos + star.momentum*(dt/star_mass)
            

    t += dt
    print(t)

#the z axis is pointing UP TOWARDS YOU