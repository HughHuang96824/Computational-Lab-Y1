#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 7 14:24:33 2016

@author: Hugh
"""




import scipy as sp
import numpy as np
import pylab as pl
import scipy.integrate as spi
import matplotlib.pyplot as plt

G=6.67*(10**-11)
M=8.4*(10**23)
R=3.4*10**6
m=260
Vx_mars=-700
Vy_mars=300
restitution=0.7

def f(r,t):

    xx=r[0]
    vx=r[1]
    yy=r[2]
    vy=r[3]
    ax=-G*M*(r[0]-Vx_mars*t)/((r[0]-Vx_mars*t)**2 + (r[2]-Vy_mars*t)**2)**1.5
    ay=-G*M*(r[2]-Vy_mars*t)/((r[0]-Vx_mars*t)**2 + (r[2]-Vy_mars*t)**2)**1.5
#sets of ODE that shows displacements, velocities and accelerations in x and y-axis
    return [vx,ax,vy,ay]

t=sp.linspace(0.,500000.,500000)
vi=-200
vj=-vi
x0=8*R
y0=3*R
xx0=[x0,vi,y0,vj]
#escape velocity
V_escape=(2*G*M/(x0**2+y0**2)**0.5)**0.5

soln=spi.odeint(f,xx0,t)
x=soln[:,0]
v1=soln[:,1]
y=soln[:,2]
v2=soln[:,3]

#Equations of kinetic, potential and total energy
KE=0.5*m*(v1**2.+v2**2.)
PE=-G*m*M/((x-Vx_mars*t)**2+(y-Vy_mars*t)**2)**0.5
E=PE+KE

#Calculating the distance between Mars(moving) and the satellite
r=((x-Vx_mars*t)**2+(y-Vy_mars*t)**2)**0.5
   
#the position that is nearest to Mars 
index1=np.where(r == r.min())#http://stackoverflow.com/questions/19546863/find-the-index-of-minimum-values-in-given-array-in-python
if V_escape < ((vi-Vx_mars)**2+(vj-Vy_mars)**2)**0.5:#Checking whether it escapes
    if v1[-1]>0:
        angle = sp.arctan(v2[-1]/v1[-1])+1.25*sp.pi#Angular deviation
    else:
        angle = sp.arctan(v2[-1]/v1[-1])+0.25*sp.pi#Angular deviation  
    #print "The satellite is not captured; the angular deviation is ",angle,"."
    #print "At position (",x[index1][0],",",y[index1][0],"), the satellite is closest to the Mars; the distance is ",min(r),"."
else: 
    #print "The satellite is captured."
 if min(r) < R:#Check whether it collides with Mars and find the index for collison
    i=0
    while r[i]-R > 0:
        i+=1
    r_=r[:i+1]
    rr=r_-R
    index_collision=np.where(rr == rr.min())
    index_collision=index_collision[0]-1
    
    #print "The satellite collides with Mars at (",x[index_collision][0],",",y[index_collision][0],")."
    
 ## Method of Plotting Angular Deviation V.S. Speed   
Vi=np.arange(-(2*(Vy_mars-Vx_mars)+np.sqrt(4*(Vx_mars-Vy_mars)**2-8*(Vy_mars**2+Vx_mars**2-V_escape**2)))/4,-100000,-10)
def VelocityToTheta(v_i):

   t=sp.linspace(0.,10000.,20)
   v_j=-v_i
   xx0=[x0,v_i,y0,v_j]
   soln=spi.odeint(f,xx0,t)
    
   x=soln[:,0]
   v1=soln[:,1]
   y=soln[:,2]
   v2=soln[:,3]
   if v1[-1]>0:
       AllAngles = sp.arctan(v2[-1]/v1[-1])+1.25*sp.pi
   else:
       AllAngles = sp.arctan(v2[-1]/v1[-1])+0.25*sp.pi
   return AllAngles

NewTheta = []
for v_i in Vi:
    NewTheta.append(VelocityToTheta(v_i))
Vi=abs(Vi)*2**0.5




if min(r) < R:
    x=x[:index_collision[0]+1]#This code makes the satellite stop at the surface of the Mars if it collides
    y=y[:index_collision[0]+1]
    
    #Find the reflected velocity after collision 
    Positionx_collision=x[index_collision][0]-t[index_collision]*Vx_mars
    Positiony_collision=y[index_collision][0]-t[index_collision]*Vy_mars
    Calculation=((v1[index_collision]*Positionx_collision)+(v2[index_collision]*Positiony_collision))/(Positionx_collision**2+Positiony_collision**2)
    Vx_normal=Calculation*Positionx_collision
    Vy_normal=Calculation*Positiony_collision
    Vx_tangent=v1[index_collision]-Vx_normal
    Vy_tangent=v2[index_collision]-Vy_normal
    Vx_aftercollision=Vx_tangent-Vx_normal*restitution
    Vy_aftercollision=Vy_tangent-Vy_normal*restitution
    #Making a new trajectory after the collision
    xx0_new=[x[index_collision][0],Vx_aftercollision,y[index_collision][0],Vy_aftercollision]
    t_new=t+t[index_collision][0]
    soln_new=spi.odeint(f,xx0_new,t_new)
    
    x_new=soln_new[:,0]
    vx_new=soln_new[:,1]
    y_new=soln_new[:,2]
    vy_new=soln_new[:,3]

    r_new=((x_new-Vx_mars*t_new)**2+(y_new-Vy_mars*t_new)**2)**0.5
    
    #Checking whether it collides with Mars again
    if min(r_new)<R:
        i=0
        while r_new[i]-R >= 0:
            i+=1
        x_new=x_new[:i+1]
        y_new=y_new[:i+1]
        #Combining trajectories
        x=np.concatenate((x,x_new),axis=0)
        y=np.concatenate((y,y_new),axis=0)

  
pl.figure(1)
pl.plot(x,y)
circle3=pl.Circle((0,0),radius=R,color='r')#The origin location of the Mars
plt.gca().add_patch(circle3)
if min(r)<R:
    circle1=pl.Circle((t[index_collision]*Vx_mars,t[index_collision]*Vy_mars),radius=R,color='orange')#Mars initial position
    plt.gca().add_patch(circle1)#This is the position of the Mars when the satellite hits it the first time
    if min(r_new)<R:
        circle2=pl.Circle((t_new[i]*Vx_mars,t_new[i]*Vy_mars),radius=R,color='yellow')
        plt.gca().add_patch(circle2)#This is the position of the Mars when the satellite hits it the second time
pl.axis('equal')#http://stackoverflow.com/questions/17990845/how-to-equalize-the-scales-of-x-axis-and-y-axis-in-python-matplotlib
pl.title("Y Distance v.s X Distance")
pl.xlabel("X Distance (m)")
pl.ylabel("Y Distance (m)")

pl.figure(2)
pl.title("Angular Deviation If Not Captured v.s Speed")
pl.xlabel("Speed (m/s)")
pl.ylabel("Angular Deviation If Not Captured (rad)")
pl.plot(Vi,NewTheta,"r-")

pl.figure(3)
pl.plot(t,KE)
pl.title("Kinetic Engrgy v.s Time")
pl.xlabel("Time")
pl.ylabel("Kinetic Engrgy")

pl.figure(4)
pl.plot(t,PE)
pl.title("Potential Engrgy v.s. Time")
pl.xlabel("Time")
pl.ylabel("Potential Engrgy")

pl.figure(5)
pl.plot(t,E)
pl.title("Total Engrgy v.s. Time")
pl.xlabel("Time")
pl.ylabel("Total Engrgy")

pl.show()
