# -*- coding: utf-8 -*-
"""
Created on Tue Jan 21 16:36:18 2020

@author: awp18
"""
import numpy as np
import scipy as sp
import Ball as b
import pylab as pl
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


class Simulation:
    def __init__(self, n = 2, m = 1, ball_R = 1, r = [0, 0], v = [0, 0], con_R = 10, M = 10**100, v_factor = 1):
        self.ball_R=ball_R
        self.con_R=con_R
        self.m = m
        self.v_factor=v_factor #the scaling factor of the velocities
        self.container = b.Ball(M, -con_R, [0,0], [0,0])
        self.ball_list = [] # the list of balls inclduing the container
        self.mom_change = [] #The list of all impulses with the container
        self.distance = []
        self.sim_time = 0 #The total of the simulation
        self.speed = []
        self.each_time_for_collision = [] #The list of the time to collision
        self.KE_list = [] #The list of the sum f KE in each frame
        self.testing_KE = [] #The list of the KE of all individual particles in all frames
        self.temperature= [] # The list of tempeature in each frame
        self.n = n #square root of the number of balls

        
        for x in np.linspace(-con_R*0.6, con_R*0.6, num=n, endpoint =True):
            for y in np.linspace(-con_R*0.6, con_R*0.6, num=n, endpoint = True):
                ball=b.Ball(m, ball_R, r = [x, y], v = np.random.normal(0, scale = v_factor, size = 2))
                self.ball_list.append(ball) #This creates the grid

        self.container = b.Ball(M, -self.con_R, [0,0], [0,0])
        self.ball_list.append(self.container)

        

        
    def next_collision(self):
        col = []  
        Kin_E =[]
        sim_time = 0
        Kin_E_each_frame = 0
        
        big_dt = 10**5
        for i in self.ball_list:
            for j in self.ball_list:
                if self.ball_list.index(i) < self.ball_list.index(j):

                    t_ball = i.time_to_collision(j)

               
                    if t_ball < big_dt:
                        big_dt=t_ball
                        col = [t_ball, i, j] #This is to find the smallest time to collision
                        
        for i in self.ball_list:
            i.move(col[0])
        self.sim_time = self.sim_time + col[0]
        self.each_time_for_collision.append(self.sim_time)

        col[1].collide(col[2])

        self.mom_change.append(col[1].mom_change)
     
        for i in self.ball_list:
            KE = 0.5*i.m*np.dot(i.vel(),i.vel())
            Kin_E.append(KE)
            self.testing_KE.append(KE)# This is to calculate the kinetic energy of each particle
            
        Kin_E_each_frame = round(np.sum(Kin_E), 5)
        self.KE_list.append(Kin_E_each_frame)#This is to find the total KE of each frame
        
    def run(self, num_frames, animate=False, graphs = False, Max_Boltzmann =False, Pressure = False):
        if animate == True:
            f = pl.figure()
            ax = pl.axes(xlim=(self.container.R, -self.container.R), 
                         ylim=(self.container.R, -self.container.R))
            
            ax.add_patch(self.container.patch)
            for i in self.ball_list:
                ax.add_patch(i.patch)
        
        for frame in range(num_frames):
            self.next_collision()
            print("the frame is", frame)
            tot_mom_change = np.sum(self.mom_change)    

            if frame > num_frames - 10: #This is to only take data from the last 10 frames.
                for i in self.ball_list:
                    if self.ball_list.index(i) < len(self.ball_list)-2:
                        self.speed.append(np.linalg.norm(i.vel()))# This is to find the speed

                    if self.ball_list.index(i) < len(self.ball_list)-1:
                        absolute = np.linalg.norm(i.pos())
                        self.distance.append(absolute) #This is to find the distance from centre

            if animate:
                pl.pause(0.05)
        if graphs:
            fig, axs = plt.subplots(2,2)
            axs[0, 0].hist(self.distance, 50)
            axs[0, 0].set_title("distacne from centre")
            axs[0, 1].plot(self.each_time_for_collision, self.KE_list)
            axs[0, 1].set_title("kinetic energy")
            axs[1,0].hist(self.speed, 20, density = True)
            axs[1,0].set_title("velocity distribution")
            #This is to creates a collection of plot
        
        if Max_Boltzmann:
            bins_mid =[]
            n, bins, plot=plt.hist(self.speed, 15)
            for i in range(len(bins)):
                if i<len(bins)-1:
                    bins_mid.append((bins[i]+bins[i+1])/2)
                
            K_b = 1
            def distribution(v_2, temp, constant):
                return constant*(v_2)*np.exp((-0.5*self.m*v_2**2)/(K_b*temp))
            
            Guesses = [1, 1000]
            
            popt, pcov = sp.optimize.curve_fit(distribution, bins_mid, n, Guesses)
            print ("The params are", popt)
            
            temp = self.KE_list[0]/(1.38e-23*self.n**2)
            print("The temperature is", temp)
           
            x_func = np.linspace(0, 6, num=1000, endpoint=True)
            plt.plot(x_func, distribution(x_func, *popt), 'r')
            plt.xlabel("veloicty(m/s)")
            plt.ylabel("number of particles")
            plt.title("velocity distribution")
            plt.show()
            
        if Pressure:
            pressure = tot_mom_change/(self.sim_time*2*np.pi*self.con_R)
            temp = self.KE_list[0]/(1.38e-23*self.n**2)
            
            return pressure, temp
        

        


        

