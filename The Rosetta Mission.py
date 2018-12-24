#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 20 13:33:31 2016

@author: Hugh
"""

import scipy as sp
import numpy as np
import pylab as pl
import scipy.integrate as spi
import matplotlib.pyplot as plt

G=6.67*(10**-11)
M=9.98*10**12
M_small=4.22*10**12
M_large=5.76*10**12
R_small=1100
R_large=1220
m=260
restitution=0.7

def f(r,t):

    xx=r[0]
    vx=r[1]
    yy=r[2]
    vy=r[3]
    ax=-((G*M_small*(r[0]-1659)/((r[0]-1659)**2 + r[2]**2)**1.5)+(G*M_large*(r[0]+660)/((r[0]+660)**2 + r[2]**2)**1.5))
    ay=-((G*M_small*r[2]/((r[0]-1659)**2 + r[2]**2)**1.5)+(G*M_large*r[2]/((r[0]+660)**2 + r[2]**2)**1.5))
    return [vx,ax,vy,ay]

t=sp.linspace(0.,500000.,50000)
vi=-0.599
vj=-0.03
x0=22473
y0=1100
xx0=[x0,vi,y0,vj]


soln=spi.odeint(f,xx0,t)
x=soln[:,0]
v1=soln[:,1]
y=soln[:,2]
v2=soln[:,3]


r_small=((x-1659)**2+y**2)**0.5
r_large=((x+660)**2+y**2)**0.5

def cut(m,n):
    i=0
    while m[i]-n > 0:
        i+=1
    i=i-1
    xa=x[:i+1]
    ya=y[:i+1]
    return [xa,ya,i]

x_1collision=cut(r_small,R_small)[0][-1]
y_1collision=cut(r_small,R_small)[1][-1]
    
if min(r_small) < R_small:      
    i=cut(r_small,R_small)[2]
    Positionx_collision=x[i]-1659
    Positiony_collision=y[i]
    Calculation=((v1[i]*Positionx_collision)+(v2[i]*Positiony_collision))/(Positionx_collision**2+Positiony_collision**2)
    Vx_normal=Calculation*Positionx_collision
    Vy_normal=Calculation*Positiony_collision
    Vx_tangent=v1[i]-Vx_normal
    Vy_tangent=v2[i]-Vy_normal
    Vx_aftercollision=Vx_tangent-Vx_normal*restitution
    Vy_aftercollision=Vy_tangent-Vy_normal*restitution
    xx0_new=[x[i],Vx_aftercollision,y[i],Vy_aftercollision]
    
    t_new=t+t[i]
    soln_new=spi.odeint(f,xx0_new,t_new)
    
    x_new=soln_new[:,0]
    vx_new=soln_new[:,1]
    y_new=soln_new[:,2]
    vy_new=soln_new[:,3]
    r_small_new=((x_new-1659)**2+y_new**2)**0.5
    
  
            
    
    if min(r_small_new)<R_small:
        i=0
        while r_small_new[i]-R_small > 0:
            i+=1
        i=i-1
        x_new=x_new[:i+1]
        y_new=y_new[:i+1]
        
        
        
        x=np.concatenate((cut(r_small,R_small)[0],x_new),axis=0)
        y=np.concatenate((cut(r_small,R_small)[1],y_new),axis=0)
        
        
        Positionx_collision=x_new[i]-1659
        Positiony_collision=y_new[i]
        Calculation=((vx_new[i]*Positionx_collision)+(vy_new[i]*Positiony_collision))/(Positionx_collision**2+Positiony_collision**2)
        Vx_normal=Calculation*Positionx_collision
        Vy_normal=Calculation*Positiony_collision
        Vx_tangent=vx_new[i]-Vx_normal
        Vy_tangent=vy_new[i]-Vy_normal
        Vx_aftercollision=Vx_tangent-Vx_normal*restitution
        Vy_aftercollision=Vy_tangent-Vy_normal*restitution
        xx0_new=[x_new[i],Vx_aftercollision,y_new[i],Vy_aftercollision]
        
        
        t_new=t+t_new[i]
        soln_new=spi.odeint(f,xx0_new,t_new)
        x_new=soln_new[:,0]
        vx_new=soln_new[:,1]
        y_new=soln_new[:,2]
        vy_new=soln_new[:,3]
        r_small_new=((x_new-1659)**2+y_new**2)**0.5
        
        if min(r_small_new)<R_small:
            i=0
            while r_small_new[i]-R_small >= 0:
                i+=1
           
            x_new=x_new[:i+1]
            y_new=y_new[:i+1]
            
            
            x=np.concatenate((x,x_new),axis=0)
            y=np.concatenate((y,y_new),axis=0)
            print (2*sp.pi*R_small*sp.arcsin(0.5*((x_1collision-x_new[-1])**2+(y_1collision-y_new[-1])**2)**0.5/R_small)/sp.pi)
      
            

pl.figure(1)
pl.plot(x,y)
circle1=pl.Circle((1659,0),radius=R_small,color='r')
circle2=pl.Circle((-660,0),radius=R_large,color='r')
plt.gca().add_patch(circle1)
plt.gca().add_patch(circle2)
pl.axis('equal')
pl.title("Y Distance v.s X Distance")
pl.xlabel("X Distance (m)")
pl.ylabel("Y Distance (m)")




