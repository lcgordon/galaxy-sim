# -*- coding: utf-8 -*-
"""
Created on Sun Apr 19 15:45:05 2020

@author: Lindsey Gordon

Last Updated: 4/25/2020

Loosely based on https://github.com/vpython/vpython-jupyter/blob/master/Demos/Stars.ipynb
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

#establishes VPython canvas
scene = canvas() 

#parameters (don't play with these)
Nstars = 400  # change this to have more or fewer stars TOTAL
N1 = Nstars/2
N2 = Nstars - N1
G = 4.3e-30 #grav constant in units of pc/Msun * (pc/s)^2
Msun = 2e30 #kg


#parameters (can play with these)
Rsun = 250 #pc (this is the size of the spheres)
L = 3e4 #pc, this is the size scale of the disk of the milky way
gal_rad = 0.5 * L #pc
gal_mass = 1 * 10**11 #Msun, VISIBLE mass
gal_mass_dark = 14 * 10**11 #dark matter mass in galaxy
star_mass = gal_mass/Nstars #units of Msun


#plotting purposes
scene.range = 2*L
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

scene.width = scene.height = 600
scene.title = "Galaxy" # Display text below the 3D graphics:
scene.caption = """Right button drag or Ctrl-drag to rotate "camera" to view scene.
To zoom, drag with middle button or Alt/Option depressed, or use scroll wheel.
  On a two-button mouse, middle is left + right.
Touch screen: pinch/extend to zoom, swipe or two-finger rotate."""


#Producing all the mass elements. Yes they are labelled stars. No I am not going to fix that.
 
Stars = []
star_colors = [color.cyan, color.magenta] #each galaxy is a different color. 
psum = vec(0,0,0)

r = np.zeros(Nstars)
#produce density distro
for i in range(Nstars):
    if i <= int(0.40 * Nstars):
        r[i] = np.random.randint(0, 0.25 * gal_rad)
    elif int(0.40 * Nstars) < i <= int(0.70 * Nstars):
        r[i] = np.random.randint(0.25*gal_rad, 0.5*gal_rad)
    elif int(0.70*Nstars) < i <= int(0.9*Nstars):
        r[i] = np.random.randint(0.5*gal_rad, 0.75*gal_rad)
    else:
        r[i] = np.random.randint(0.75*gal_rad, gal_rad)


#r = np.random.randint(0,gal_rad,Nstars) #N radial distances


phi_degrees = np.random.randint(0,360,Nstars) #N positions around axis
phi = np.radians(phi_degrees)
mom_phi = np.radians(90-phi_degrees)

x = r*np.sin(phi)
y = r*np.cos(phi)
for i in range(Nstars): #this separates the two galaxies out
    if i % 2 == 0:
        x[i] = x[i] - (0.75*L)
    elif i%2 == 1:
        x[i] = x[i] + (0.75*L)
        


for i in range(Nstars):
    star = sphere(pos=vector(x[i], y[i],0), make_trail=True, retain=20, trail_radius=50)
    star.radius = Rsun
    star.mass = star_mass
    star.momentum = vector(0,0,0) #initialized as zero
    star.color = star.trail_color = star_colors[i % 2]
    Stars.append(star)
    
    
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
v_left = vector(6e-12, 0, 0)
v_right = vector(-6e-12, 0,0)
for i in range(Nstars):
    v_i[i] = np.sqrt(G * mass_internal[i]/np.abs(r[i])) #units of pc/s


for i in range(Nstars):
    velocity_v = vector(v_i[i]*np.sin(mom_phi[i]), -1* v_i[i]*np.cos(mom_phi[i]), 0)
    if i % 2 == 0:
        velocity_v += v_left
    elif i%2==1:
        velocity_v += v_right
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
            if rmag2 <= (radius1 + radius2)**2: #so currently they ignore each other if they get too close together.
                r = vec(0,0,0) #this is not a perfect solution but it's what i've got right now.
             
            F = F + (G*m1*sj.mass/(rmag2**1.5))*r
        si.momentum = si.momentum + F*dt
      

t = 0 #initializes time for loop

while t < 1e16:

    # Compute all forces on all stars
    computeForces() #REMOVE THIS COMMENT IN ORDER TO RUN SIMULATION
    
    # Having updated all momenta, now update all positions
    for i in range(len(Stars)):
        if Stars[i] is None: continue
        Stars[i].pos = Stars[i].pos + Stars[i].momentum*(dt/star_mass) 

    t += dt
    print(t)

#the z axis is pointing UP TOWARDS YOU