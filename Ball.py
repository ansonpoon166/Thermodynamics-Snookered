# -*- coding: utf-8 -*-
"""
Created on Tue Jan 21 14:47:14 2020

@author: awp18
"""


import numpy as np
import pylab as pl

class Ball:
    
    def __init__(self, m = 1.0, R = 0.005, r = [0, 0], v = [0, 0]):
        self.m = m
        self.R = R
        self.__r = np.array(r)
        self.__v = np.array(v)
        self.mom_change = 0
        if R > 0:
            self.patch =  pl.Circle(self.pos(), R, fc='r')
        else:
            self.patch = pl.Circle(r, -R, ec='b', fill=False, ls='solid') 
            #the container is given a negative radius.
            
        """
        This initailise the ball, with attributes including the radius, mass, position and velocity.
        The if function creates patches in the simulation, which creates a ball\
        for a gas particle and a hollow circle for the conatiner.
        """
        
    def pos(self):
        return self.__r

    def vel(self):
        return self.__v
    
    def setr(self, new_r):
        self.__r = np.array(new_r)
        self.patch.center = new_r
        
    def move(self, dt):

        g = self.vel()*dt
        r1 = np.add(self.pos(), g)
        self.setr(r1)
        
    def setv(self, new_v):
        self.__v = np.array(new_v)
        
    def time_to_collision(self, other):

        dr = self.pos() - other.pos()
        dv = self.vel() - other.vel()
        dR = self.R + other.R
        
        dis = (4*np.dot(dr,dv)**2-4*(np.dot(dv,dv))*((np.dot(dr,dr))-(dR**2)))
        # dis is the discriminant
        if dis>=0:
            
            dt1 = (-2*np.dot(dr,dv) + (dis)**0.5)/(2*np.dot(dv,dv)) 
            dt2 = (-2*np.dot(dr,dv) - (dis)**0.5)/(2*np.dot(dv,dv))
        
            if dt1>1e-9 and dt2>1e-9 :
                dt = min(dt1,dt2)
                return dt
            elif dt1 > 1e-9:
                dt = dt1
                return dt
            elif dt2 > 1e-9:
                dt = dt2
                return dt
            else:
                dt = 1e9
                return dt
        else:
            return 1e9
        
        
        

    def collide(self, other):
        dr = self.pos() - other.pos()
        dr_unit = dr/(np.linalg.norm(dr))
        v1_parallel_before = np.dot(self.vel(), dr_unit)
        v1_perpendicular = self.vel() - v1_parallel_before*dr_unit
        v2_parallel_before = np.dot(other.vel(), dr_unit)
        v2_perpendicular = other.vel() - v2_parallel_before*dr_unit
        
        dv = self.vel() - other.vel()
        v1_parallel_after = ((self.m - other.m)*v1_parallel_before)/(self.m + other.m) + ((2*other.m*v2_parallel_before)/(self.m + other.m))
        v2_parallel_after = ((2*self.m*v1_parallel_before)/(self.m + other.m)) - (((self.m - other.m)*other.vel())/(self.m + other.m))
        v1_after = v1_parallel_after*dr_unit + v1_perpendicular
        v2_after = v2_parallel_after*dr_unit + v2_perpendicular
        

        KE_before = 0.5*self.m*np.dot(self.vel(), self.vel()) + 0.5*other.m*np.dot(other.vel(), other.vel())
        KE_after = 0.5*self.m*np.dot(v1_after,v1_after) + 0.5*other.m*np.dot(v2_after, v2_after)
        if self.R>0 and other.R < 0: #this allows the program to consider only the collsion with the container
        
            self.mom_change=self.m*(np.linalg.norm(v1_after-self.vel()))
        elif self.R<0 and other.R>0:
            self.mom_change=other.m*(np.linalg.norm(v2_after-other.vel()))
        else:
            self.mom_change=0 #This finds the momentum change after the collision
        
        self.setv(v1_after)
        other.setv(v2_after) #This sets the new velocities

        
    
    
    