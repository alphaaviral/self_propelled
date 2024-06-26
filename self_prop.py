import random
import numpy as np
import matplotlib.pyplot as plt
import math
from itertools import combinations
import copy
import csv
from sys import getsizeof
import pandas as pd
import time

class Particle:
    def __init__(self, pos, linvel):
        self.pos=pos
        self.linvel=linvel
        # self.isFromWedge=isFromWedge

class Chain:
    def __init__(self, chainSize=None, pos=None, angle=None, isWedge=False):    
        if pos is None:
            self.pos=np.array([random.uniform(boxX[0], boxX[1]), random.uniform(boxY[0], boxY[1])], dtype=np.float64)
        else:
            self.pos=pos

        if angle is None:
            self.angle=random.uniform(0,2*np.pi)
        else:
            self.angle=angle
        self.orient= np.array([np.cos(self.angle), np.sin(self.angle)], dtype=np.float64)      #unit vector, equal to distance between COM of consecutive particles
        self.linvel=np.array([0,0], dtype=np.float64)
        self.angvel=0                                                                          #vector-wise this is in z
        if chainSize is None:
            chainSize=globalChainSize
        self.chainSize=chainSize
        self.force=Fo*self.orient
        self.moment=0
        self.particles=[]
        self.isWedge=isWedge
        if self.chainSize%2!=0:
            for i in range(-1*int(self.chainSize/2),int(self.chainSize/2)+1):
                particlePos=self.pos+(d*i*self.orient)
                particleVel=self.linvel+zxycross(self.angvel, particlePos-self.pos)
                self.particles.append(Particle(particlePos, particleVel))
        else:
            for i in range(-1*(int(self.chainSize/2)-1),int(self.chainSize/2)+1):
                particlePos=self.pos+(d*i*self.orient)-(d/2)*self.orient
                particleVel=self.linvel+zxycross(self.angvel, particlePos-self.pos)
                self.particles.append(Particle(particlePos, particleVel))
                
    def updateParticles(self):
        for i in range(len(self.particles)):
            self.particles[i].pos=self.pos+((i-int(len(self.particles)/2))*d*self.orient)
            self.particles[i].linvel=self.linvel+zxycross(self.angvel, self.particles[i].pos-self.pos)

# cross of a vector in x-y plane with a vector perpendicular to it
def zxycross(z:float, xy:np.ndarray):
    return np.array([-1*xy[1]*z, xy[0]*z])

# force on a due to b
def getForce(a:Particle, b:Particle)->np.ndarray:
    # xdist = a.pos[0]-b.pos[0]
    # ydist = a.pos[1]-b.pos[1]
    # if abs(xdist)>boxSize[0]-abs(xdist):
    #     xdist*=-1
    # if abs(ydist)>boxSize[1]-abs(ydist):
    #     ydist*=-1
    # xdist = min(xdist, boxSize[0]-xdist)
    # ydist = min(ydist, boxSize[1]-ydist)
    # rmod=np.linalg.norm((xdist, ydist))
    apos = copy.deepcopy(a.pos)
    bpos= copy.deepcopy(b.pos)
    if abs(apos[0]-bpos[0])>(boxSize[0])/2:
        if apos[0]>bpos[0]:
            bpos[0]+=boxSize[0]
        else:
            apos[0]+=boxSize[0]

    if abs(apos[1]-bpos[1])>(boxSize[1])/2:
        if apos[1]>bpos[1]:
            bpos[1]+=boxSize[1]
        else:
            apos[1]+=boxSize[1]

    rmod=np.linalg.norm(apos-bpos)
    rhat=(apos-bpos)/rmod
    f = ((Uo*math.exp(-1*rmod/lamb))/rmod)*((1/rmod)+(1/lamb))*rhat
    return f
    

def updateChainForces(a:Chain, b:Chain):
    # if not isinstance(a, Chain) or not isinstance(b, Chain):
    #     raise TypeError("force cannot be calculated for non chains")
    if (not a.isWedge) and (not b.isWedge):
        xdist = abs(a.pos[0]-b.pos[0])
        ydist = abs(a.pos[1]-b.pos[1])
        # if boxSize[0]-xdist <=0:
        #     raise ValueError("box size too small")
        if min(xdist, boxSize[0]-xdist)>2*l or min(ydist, boxSize[1]-ydist)>2*l:
            # print("skipped")
            return

    for i in range(len(a.particles)):
        for j in range(len(b.particles)):
            force=getForce(a.particles[i], b.particles[j])
            a.force+=force
            b.force-=force
            a.moment+=np.cross(a.particles[i].pos-a.pos, force)
            b.moment+=np.cross(b.particles[j].pos-b.pos, -1*force)

def updateChainPositions(chain:Chain):
    # if not isinstance(chain, Chain):
    #     raise TypeError("force cannot be calculated for non chains")
    chain.pos+=chain.linvel*timeStep
    chain.pos[0]=((chain.pos[0]-boxX[0])%(boxSize[0]))+boxX[0]
    chain.pos[1]=((chain.pos[1]-boxY[0])%(boxSize[1]))+boxY[0]
    chain.angle+=chain.angvel*timeStep
    if chain.angle>2*np.pi:
        chain.angle-=2*np.pi
    if chain.angle<0:
        chain.angle+=2*np.pi
    chain.orient=np.array([np.cos(chain.angle), np.sin(chain.angle)])

def updateChainVelocities(chain:Chain):
    # if not isinstance(chain, Chain):
    #     raise TypeError("force cannot be calculated for non chains")
    chain.linvel=((np.dot(chain.force, chain.orient)/ft1)*chain.orient)+((chain.force-np.dot(chain.force, chain.orient)*chain.orient)/ft2)
    chain.angvel=chain.moment/fr

def getFrictionCoeff():
    a=aRatio
    fParallel = 2*np.pi/(np.log(a)-0.207+(0.98/a)-(0.133/(a**2)))
    fPerpendicular = 2*np.pi/(np.log(a)+0.839+(0.185/a)+(0.233/(a**2)))
    fRot = np.pi*(a**2)/(3*(np.log(a)-0.662+(0.917/a)-(0.05/(a**2))))
    return fo*fParallel, fo*fPerpendicular, fo*fRot

def addWedge(chainArray):
    # if wedgeSize%2==0:
    #     raise ValueError("length of wedge must be odd")
    upperChain = Chain(chainSize=wedgeSize, pos=np.array(([(wedgeSize-1)*d*np.cos(wedgeAngle/2)/2,(wedgeSize-1)*d*np.sin(wedgeAngle/2)/2]), dtype=np.float64), 
                       angle=wedgeAngle/2, isWedge=True)
    
    lowerChain = Chain(chainSize=wedgeSize, pos=np.array(([(wedgeSize+1)*d*np.cos(-1*wedgeAngle/2)/2,(wedgeSize+1)*d*np.sin(-1*wedgeAngle/2)/2]), dtype=np.float64), 
                       angle=-1*wedgeAngle/2, isWedge=True)
    
    chainArray.append(upperChain)
    chainArray.append(lowerChain)
    # return upperChain, lowerChain



if __name__=="__main__":
    Uo=1
    lamb=1
    l=10
    aRatio=l/lamb
    globalChainSize=int(round(9*aRatio/8))
    d=l/(math.sqrt((globalChainSize+1)*(globalChainSize-1)))
    timeStep=0.05
    Fo=2
    chainNos=30
    fo=1
    iterations=50
    wedgeSize=globalChainSize*4
    wedgeAngle=np.pi/2
    ft1, ft2, fr = getFrictionCoeff()
    boxX=(-5*l,5*l)
    boxY=(-5*l,5*l)
    boxSize = (boxX[1]-boxX[0], boxY[1]-boxY[0])
    totalChainNos=chainNos+2 #including wedges
    
    chains=[]
    addWedge(chains)
    for i in range(chainNos):
        chains.append(Chain())

    dataArr=np.zeros([iterations, totalChainNos, 6])
    for t in range(0, iterations):
        # timeData=[[t]]
        print(t)
        t_init=time.time()
        comb = combinations(chains, 2)

        for i in list(comb):
            updateChainForces(i[0], i[1])
        for i in range(len(chains)):
            # print(i)
            chain=chains[i]
            if not chain.isWedge:
                updateChainVelocities(chain)
                updateChainPositions(chain)
                chain.updateParticles()
            chain.force=Fo*chain.orient
            chain.moment=0

            dataArr[t][i][0]=chain.pos[0]
            dataArr[t][i][1]=chain.pos[1]
            dataArr[t][i][2]=chain.angle
            dataArr[t][i][3]=chain.linvel[0]
            dataArr[t][i][4]=chain.linvel[1]
            dataArr[t][i][5]=chain.angvel
            
            for j in range(len(chain.particles)):
                particle=chain.particles[j]
                if chain.isWedge:
                    plt.plot(particle.pos[0], particle.pos[1], 'b', marker=".", markersize=6)
                else:
                # for j in range(len(chain.particles)):
                    plt.plot(particle.pos[0], particle.pos[1], 'r', marker=".", markersize=6)

        plt.xlim(boxX[0], boxX[1])
        plt.ylim(boxY[0], boxY[1])

        plt.savefig("./out/test"+str(t)+".png")
        plt.clf()
        print(time.time()-t_init)
    dataArr=dataArr.reshape([iterations, totalChainNos*6])
    df=pd.DataFrame(dataArr)
    df.to_csv("./out.csv",header=False)