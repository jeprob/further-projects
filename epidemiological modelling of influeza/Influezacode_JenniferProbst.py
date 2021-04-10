#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Name:      Probst, Jennifer
Email:     jennifer.probst@uzh.ch
Date:      11/03/2019
Semester:  FS19
Topic:     Epidemiological modeling of influenza with cellular automata
"""

import random
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.animation import FFMpegFileWriter
import seaborn
import math
import matplotlib.patches as mpatches
plt.rcParams['animation.ffmpeg_path'] = 'C:/Users/probs/Desktop/FFmpeg/bin/ffmpeg.exe'

class Society(): 
    ''' set up societies with their parameters'''
    
    def __init__(self, name, population, c, y, mS, mI, mR):
        #set variables
        self.name=name
        self.population=population
        self.c=c #contact rate
        self.y=y #recovery rate
        self.mS=mS #migration rate suceptible
        self.mR=mR #migration rate recovered
        self.mI=mI #migration rate infected
        
class Pandemic(): 
    '''set up pandemics with their parameters'''
    
    def __init__(self, name, b):
        #set variables
        self.name=name
        self.b=b #force of infection

class CellularAutomata(): 
    
    def __init__(self, society, pandemic, times, n):
        #define variables
        self.society=society
        self.pandemic=pandemic
        self.n=n
        self.times=times
        self.timeepidemic=0
        self.maxI=1
        self.fig = plt.figure(figsize=(self.n, self.n), frameon=False)
        self.colormap = seaborn.diverging_palette(220, 10, as_cmap = True)

        #set up arrays and total lists
        self.fractions= np.zeros((self.n,self.n, 3, self.times))
        self.people=np.zeros((self.n,self.n,3, self.times))
        self.plotting=np.zeros((self.n,self.n, 3, self.times))

        #distribute people to the grid
        self.fractions[:,:,:,0]=[1,0,0]
        self.setup_people()
        
        #choose where infection starts
        self.setup_infection_start()
        
    def setup_people(self):
        '''sets people in every grid with normal distribution around n'''
        a=[int(np.random.normal(self.society.population, 0.01*self.n, 1)) for i in range(self.n)]
        self.people[:,:,0,0]=a
        
    def setup_infection_start(self):
        '''starts the infections with 1 sick individual'''
        self.start_infection = math.floor(self.n / 2)
        a=self.people[self.start_infection, self.start_infection, :, 0]
        self.fractions[self.start_infection, self.start_infection, :, 0] = [(a[0]-1)/a[0], 1/a[0], 0]
        self.people[self.start_infection, self.start_infection, :, 0] = [a[0]-1, 1, 0]

    def loop(self, t):
        '''main loop for animation'''
        self.t=t
        self.simmulate1(self.fractions[:,:,:,t], self.people[:,:,:,t], t)
        self.update()

    def simmulate_whole(self):
        ''' simmulates n times and returns: total number of infected, time of pandemic/epidemic, maximum number of sick people'''
        
        #simmulate n times
        for t in range(self.times-1):
            self.simmulate1(self.fractions[:,:,:,t], self.people[:,:,:,t], t)
            Istep=(sum(sum(sum(self.people[:,:,1,:]))))
            if Istep>self.maxI:
                self.maxI=Istep
            if self.Istep>0.001 or self.recover>0.01: 
                self.timeepidemic+=1
        self.tI=sum(sum(sum(self.people[:,:,2,:])))
        
        #add totals
        self.total_S=[(sum(sum(self.people[:,:,0,t]))) for t in range(self.times)]
        self.total_I=[(sum(sum(self.people[:,:,1,t]))) for t in range(self.times)]
        self.total_R=[(sum(sum(self.people[:,:,2,t]))) for t in range(self.times)]
        
        self.plot()
        
        return(self.tI, self.timeepidemic, self.maxI)
        
    def simmulate1(self, fractions, people, t):
        ''' simmulates one round of updates & gives back the updated fractions and total numbers
            logic:  
            1. for all i update fractions in all grids with diff equ and total numbers
            2. calculate migration people and 
            3. store migration in migration matrix
            4. add total numbers and mirgation matrix
            5. update fractions'''
        
        #update fractions with diff equations
        self.recover=0
        for i in range(len(fractions)):
            for j in range(len(fractions)):
                fractions[i][j]=self.update_grid(fractions[i][j])
        
        #update total numbers
        for i in range(len(fractions)): 
            for j in range(len(fractions)):
                total_peops=sum(people[i][j])
                frac=fractions[i][j]
                people[i][j]=[frac[0]*total_peops,frac[1]*total_peops,frac[2]*total_peops]
                
        #calculate migration people, store migration in matrix (the ones leaving and coming)
        mm=np.zeros((self.n,self.n,3))
        for i in range(len(people)): 
            for j in range(len(people)):
                mm+=self.simmulate_migration(i, j, people[i][j])
                
        #add total number and migration matrix
        people+=mm
        
        #update frequences      
        for i in range(len(fractions)):
            for j in range(len(fractions)):
                total_peops=sum(people[i][j])
                peops=people[i][j]
                fractions[i][j]=[peops[0]/total_peops, peops[1]/total_peops, peops[2]/total_peops]
    
        #update values
        self.plotting=self.fractions[:, :, 1, t]
        self.fractions[:, :, :, t + 1] = fractions
        self.people[:, :, :, t + 1] = people
        self.Istep=(sum(sum(fractions[:,:,1])))
    
    def update_grid(self, grid):
        '''gets fractions in one grid as imput and returns updated fractions'''
        return ([self.update_S(grid), self.update_I(grid), self.update_R(grid)])
    
    def update_S(self, grid):
        '''given fractions in one grid, returns updated S'''
        a=grid[0]-self.society.c*grid[0]*grid[1]*self.pandemic.b
        return a
    
    def update_I(self, grid):
        '''given fractions in one grid, returns updated I'''
        b=grid[1]-self.society.y*grid[1] + self.society.c*grid[0]*grid[1]*self.pandemic.b
        return b
    
    def update_R(self, grid):
        '''given fractions in one grid, returns updated R'''
        self.recover=self.society.y*grid[1]
        c=grid[2]+self.recover
        return c

    def select_neighbor(self, i, j):
        ''' outputs neighbor coordinates as a tuple that randomly gets migrants
            change in coordinates:
            a= change in i
            b= change in j
            no migrants can get out of the grid
            ''' 
        while True: 
            a=random.randint(-1, 1)
            b=random.randint(-1, 1)
            if a==0 and b==0: 
                continue
            if i==(self.n-1) and a==1 and j==(self.n-1) and b==1:
                return (-1,-1) 
            if i==(self.n-1) and a==1:
                return (-1,b)
            if j==(self.n-1) and b==1: 
                return (a,-1)
            else:
                return (i+a,j+b)
                
    def simmulate_migration(self, i, j, grid):
        '''given absolute numbers of people in one grid it calculates the people to migrate 
            and where to and stores them at their destination grid in an otherwise empty 
            matrix which is returned.'''
        #specify temporary migration matrix
        mat=np.zeros((self.n,self.n,3))
        #calculate people to migrate
        migrants_S=grid[0]*self.society.mS
        migrants_I=grid[1]*self.society.mI
        migrants_R=grid[2]*self.society.mR
        #where to migrate
        a,b=self.select_neighbor(i,j)
        mat[a][b]=[migrants_S, migrants_I, migrants_R]
        mat[i][j]+=[-migrants_S, -migrants_I, -migrants_R]
        return mat
    
    def plot(self):   
        ''' generate plots '''
        #timeplot of totals
        timeline=list(range(self.times))
        fig1, ax1 = plt.subplots(1, 1, figsize=(10,10))
        ax1.plot(timeline, self.total_S, color='green', label='susceptible')
        ax1.plot(timeline, self.total_I, color='red', label='infected')
        ax1.plot(timeline, self.total_R, color='blue', label='recovered')
        plt.legend()
        ax1.set_xlabel('time', fontsize=18)
        ax1.set_ylabel('total number of individuals', fontsize=18)
        ax1.set_title('{} in {}'.format(self.pandemic.name, self.society.name), fontsize=18)
    
        fig1.savefig('{} in {} on {}x{} grid.pdf'.format(self.pandemic.name, self.society.name, self.n, self.n), format='pdf')
        plt.show()

    def update(self):
        '''plot for the animation'''
        self.fig.clear()
        ax = plt.subplot()
        plt.imshow(self.plotting, cmap = self.colormap)
        ax.grid(False)
        ax.set_xlim([-1/2, self.n-1/2])
        ax.set_ylim([-1/2, self.n-1/2])
        ax.set_title('t = {}'.format(self.t+1))
        
    def animate(self):
        '''animation'''
        mywriter=FFMpegFileWriter(fps=30, codec='libx264')
        cartoon = animation.FuncAnimation(self.fig, self.loop, frames = (self.times-1), interval = 25, repeat = False)
        cartoon.save('animation30x30.mp4', writer=mywriter)
        
class Szenarios(): 
    
    def __init__(self, times, n):
        '''define different flus, forces and societies'''
        self.setup()
        self.flus=[self.spanishflu, self.hongkongflu, self.asianflu, self.swineflu, self.seasonalflu]
        self.forces=[0.25, 0.18, 0.165, 0.146, 0.128]
        self.societies=[self.paris, self.dhaka, self.seoul, self.monacco, self.zurich, self.raipur, self.bergen, self.menzingen]
        self.times=times
        self.n=n
        
    def simmulate(self):
        self.lists()
        self.plotforce()
        self.plotflu()
    
    def setup(self):
        '''setup pandemics and societies'''
        #create pandemics
        self.spanishflu=Pandemic('Spanish Flu', 0.25)
        self.hongkongflu=Pandemic('Hong Kong Flu', 0.18)
        self.asianflu=Pandemic('Asian Flu', 0.165)
        self.swineflu=Pandemic('Swineflu', 0.146)
        self.seasonalflu=Pandemic('Seasonal Flu', 0.128)
        
        #create societies: population, c, y, mS, mI, mR
        self.paris=Society('Paris', 26000, 1, 0.1, 0.1, 0.01, 0.1)
        self.dhaka=Society('Dhaka', 29105, 1, 0.1, 0.01, 0.001, 0.01)
        self.seoul=Society('Seoul', 16456, 0.95, 0.1, 0.1, 0.01, 0.1)
        self.monacco=Society('Monacco', 18713, 1, 0.1, 0.1, 0.01, 0.1)
        self.zurich=Society('Zurich', 4454, 1, 0.1, 0.1, 0.01, 0.1)
        self.raipur=Society('Raipur', 4500, 0.95, 0.1, 0.01, 0.001, 0.01)
        self.bergen=Society('Bergen', 604, 1, 0.1, 0.1, 0.01, 0.1)
        self.menzingen=Society('Menzingen', 164, 1, 0.1, 0.01, 0.001, 0.01)
    
    def lists(self):
        #set up experiments and store results in list:
        self.loc_list=[]
        for society in self.societies:
            temp=[]
            for flu in self.flus:
                Simmulation=CellularAutomata(society, flu, self.times, self.n)
                
#                Simmulation.animate()
    
                sick, time, maxsick = Simmulation.simmulate_whole()
                sick=int(sick/(society.population*self.n**2))
                maxsick=maxsick/(society.population*self.n**2)
                temp.append([society.name, flu.name, flu.b, sick, time, maxsick])
            self.loc_list.append(temp)
            
        print(self.loc_list)
        
    def plotforce(self):
        #plot for force of infection vs people ill or time of pandemic or max of infected
        fig1, ax1 = plt.subplots(1, 1, figsize=(10,10))
        fig2, ax2 = plt.subplots(1, 1, figsize=(10,10))
        fig3, ax3 = plt.subplots(1, 1, figsize=(10,10))
        for i in range(len(self.societies)): 
            time=[]
            peop=[]
            maxsick=[]
            for j in range(len(self.flus)):
                time.append(self.loc_list[i][j][4])
                peop.append(self.loc_list[i][j][3])
                maxsick.append(self.loc_list[i][j][5])
            ax1.scatter(self.forces, peop, label='{}'.format(self.societies[i].name))
            ax1.set_xlabel('Force of infection', fontsize=18)
            ax1.set_ylabel('Infected people', fontsize=18)
            ax1.set_title('Infected people vs. force of infection', fontsize=18)
            ax2.scatter(self.forces, time, label='{}'.format(self.societies[i].name))
            ax2.set_xlabel('Force of infection', fontsize=18)
            ax2.set_ylabel('Time of pandemic', fontsize=18)
            ax2.set_title('Time of pandemic vs. force of infection', fontsize=18)
            ax3.scatter(self.forces, maxsick, label='{}'.format(self.societies[i].name))
            ax3.set_xlabel('Force of infection', fontsize=18)
            ax3.set_ylabel('Max number of infected people', fontsize=18)
            ax3.set_title('Max number of infected vs. force of infection', fontsize=18)
            
        fig1.legend(loc=0) 
        fig1.show()
        fig2.legend(loc=0)
        fig2.show()
        fig3.legend(loc=0)
        fig3.show()
        
        fig1.savefig('Infected people on {}x{} grid.pdf'.format(self.n, self.n), format='pdf')
        fig2.savefig('Time of pandemic for {}x{} grid.pdf'.format(self.n, self.n), format='pdf')
        fig3.savefig('Max number of infected for {}x{} grid.pdf'.format(self.n, self.n), format='pdf')
        
    def plotflu(self):
        #plot of time of pandemic vs people for each flu
        fig = plt.figure(figsize=(40,20))
        ax1 = fig.add_subplot(231)
        ax2 = fig.add_subplot(232)
        ax3 = fig.add_subplot(233)
        ax4 = fig.add_subplot(234)
        ax5 = fig.add_subplot(235)
        fig.subplots_adjust(hspace=0.2, wspace=0.2)
        axes=[ax1, ax2, ax3, ax4, ax5]
        
        times=np.zeros((len(self.flus),len(self.societies)))
        people=np.zeros((len(self.flus),len(self.societies)))
        colors=['blue', 'orange', 'green', 'red', 'purple', 'brown', 'pink', 'grey']
        
        for i in range(len(self.societies)): 
            for j in range(len(self.flus)):
                times[j][i]= self.loc_list[i][j][4]
                people[j][i]=self.loc_list[i][j][3]
        for j in range(len(self.flus)):
            axes[j].set_title('{}'.format(self.flus[j].name), fontsize=18)
            axes[j].set_ylabel('Infected people', fontsize=18)
            axes[j].set_xlabel('Time of pandemic', fontsize=18)
            axes[j].scatter(times[j], people[j], color=colors, ) 
            b = mpatches.Patch(color='blue', label='Paris')
            o = mpatches.Patch(color='orange', label='Dhaka')
            g = mpatches.Patch(color='green', label='Seoul')
            r = mpatches.Patch(color='red', label='Monacco')
            p = mpatches.Patch(color='purple', label='Zurich')
            br = mpatches.Patch(color='brown', label='Raipur')
            pi = mpatches.Patch(color='pink', label='Bergen')
            gr = mpatches.Patch(color='grey', label='Menzingen')
            plt.legend(handles=[b,o,g,r,p,br,pi,gr], bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., fontsize=18)
        axes[4].set_xlabel('Time of epidemic', fontsize=18)
        plt.show()
        
        fig.savefig('Flus on {}x{} grid.pdf'.format(self.n, self.n), format='pdf')


if __name__== "__main__":
    ''' different gridsizes simulated '''
    
#    grid5=Szenarios(800, 5)
#    grid5.simmulate()
#    grid15=Szenarios(1700, 15)
#    grid15.simmulate( )
#    grid30=Szenarios(2500, 30)
#    grid30.simmulate()
#    grid50=Szenarios(3500, 50)
#    grid50.simmulate()