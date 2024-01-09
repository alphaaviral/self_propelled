import random
import numpy as np
import matplotlib.pyplot as plt
import math
from itertools import combinations
import copy

class Particle:
    def __init__(self, pos, linvel, isFromWedge=False):
        self.pos=pos
        self.linvel=linvel
        self.isFromWedge=isFromWedge

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
    if not isinstance(a, Particle) or not isinstance(b, Particle):
        raise TypeError("force cannot be calculated for non particles") 
    apos = copy.deepcopy(a.pos)
    bpos= copy.deepcopy(b.pos)
    if abs(apos[0]-bpos[0])>(boxX[1]-boxX[0])/2:
        if apos[0]>bpos[0]:
            bpos[0]+=boxX[1]-boxX[0]
        else:
            apos[0]+=boxX[1]-boxX[0]

    if abs(apos[1]-bpos[1])>(boxY[1]-boxY[0])/2:
        if apos[1]>bpos[1]:
            bpos[1]+=boxY[1]-boxY[0]
        else:
            apos[1]+=boxY[1]-boxY[0]

    rmod=np.linalg.norm(apos-bpos)
    rhat=(apos-bpos)/rmod
    f = ((Uo*math.exp(-1*rmod/lamb))/rmod)*((1/rmod)+(1/lamb))*rhat
    if not (a.isFromWedge or b.isFromWedge):
        f/=(globalChainSize**2)
    return f
    

def updateChainForces(a:Chain, b:Chain):
    if not isinstance(a, Chain) or not isinstance(b, Chain):
        raise TypeError("force cannot be calculated for non chains")
    
    for i in range(len(a.particles)):
        for j in range(len(b.particles)):
            force=getForce(a.particles[i], b.particles[j])
            a.force+=force
            b.force-=force
            a.moment+=np.cross(a.particles[i].pos-a.pos, force)
            b.moment+=np.cross(b.particles[j].pos-b.pos, -1*force)

def updateChainPositions(chain:Chain):
    if not isinstance(chain, Chain):
        raise TypeError("force cannot be calculated for non chains")
    chain.pos+=chain.linvel*timeStep
    chain.pos[0]=((chain.pos[0]-boxX[0])%(boxX[1]-boxX[0]))+boxX[0]
    chain.pos[1]=((chain.pos[1]-boxY[0])%(boxY[1]-boxY[0]))+boxY[0]
    chain.angle+=chain.angvel*timeStep
    if chain.angle>2*np.pi:
        chain.angle-=2*np.pi
    if chain.angle<0:
        chain.angle+=2*np.pi
    chain.orient=np.array([np.cos(chain.angle), np.sin(chain.angle)])


    # if chain.pos[0]>boxX[1] or chain.pos[0]<boxX[0]:
    #     chain.orient[0]*=-1
    #     chain.angle=np.arctan2(chain.orient[1], chain.orient[0])
    # if chain.pos[1]>boxY[1] or chain.pos[1]<boxY[0]:
    #     chain.orient[1]*=-1
    #     chain.angle=np.arctan2(chain.orient[1], chain.orient[0])    

def updateChainVelocities(chain:Chain):
    if not isinstance(chain, Chain):
        raise TypeError("force cannot be calculated for non chains")
    chain.linvel=((np.dot(chain.force, chain.orient)/ft1)*chain.orient)+((chain.force-np.dot(chain.force, chain.orient)*chain.orient)/ft2)
    chain.angvel=chain.moment/fr

# def getAspectRatio():
#     return ((2*particleRadius)+((globalChainSize-1)*1))/(2*particleRadius)

def getFrictionCoeff():
    a=aRatio
    fParallel = 2*np.pi/(np.log(a)-0.207+(0.98/a)-(0.133/(a**2)))
    fPerpendicular = 2*np.pi/(np.log(a)+0.839+(0.185/a)+(0.233/(a**2)))
    fRot = np.pi*(a**2)/(3*(np.log(a)-0.662+(0.917/a)-(0.05/(a**2))))
    return fo*fParallel, fo*fPerpendicular, fo*fRot

def makeWedge():
    # if wedgeSize%2==0:
    #     raise ValueError("length of wedge must be odd")
    upperChain = Chain(chainSize=wedgeSize, pos=np.array(([(wedgeSize-1)*d*np.cos(wedgeAngle/2)/2,(wedgeSize-1)*d*np.sin(wedgeAngle/2)/2]), dtype=np.float64), 
                       angle=wedgeAngle/2, isWedge=True)
    
    lowerChain = Chain(chainSize=wedgeSize, pos=np.array(([(wedgeSize+1)*d*np.cos(-1*wedgeAngle/2)/2,(wedgeSize+1)*d*np.sin(-1*wedgeAngle/2)/2]), dtype=np.float64), 
                       angle=-1*wedgeAngle/2, isWedge=True)
    
    return upperChain, lowerChain

if __name__=="__main__":
    Uo=10
    lamb=1
    l=3
    aRatio=l/lamb
    globalChainSize=int(round(9*aRatio/8))
    d=l/(math.sqrt((globalChainSize+1)*(globalChainSize-1)))
    timeStep=1
    Fo=2
    chainNos=20
    # particleRadius=0.55
    fo=1
    wedgeSize=globalChainSize*3
    wedgeAngle=np.pi/2
    ft1, ft2, fr = getFrictionCoeff()
    boxX=(-5*l,5*l)
    boxY=(-5*l,5*l)
    # if globalChainSize%2==0:
    #     raise ValueError("length of chain must be odd")
    
    chains=[]
    uC, lC = makeWedge()
    chains.append(uC)
    chains.append(lC)
    # exit()
    # Chain(chainSize=)

    for i in range(chainNos):
        chains.append(Chain())

    # chains.append(Chain(pos=np.array([0,0], dtype=np.float64), angle=np.pi/3))
    # chains.append(Chain(pos=np.array([0,10], dtype=np.float64), angle=5*np.pi/3))
    
    for t in range(0, 500):
        print(t)
        comb = combinations(chains, 2)

        for i in list(comb):
            updateChainForces(i[0], i[1])
        for chain in chains:
            # addSelfForce(i)
            if not chain.isWedge:
                updateChainVelocities(chain)
                updateChainPositions(chain)
                chain.updateParticles()
            chain.force=Fo*chain.orient
            chain.moment=0
            if chain.isWedge:
                for particle in chain.particles:
                    plt.plot(particle.pos[0], particle.pos[1], 'b', marker=".", markersize=6)
            else:
                for particle in chain.particles:
                    plt.plot(particle.pos[0], particle.pos[1], 'r', marker=".", markersize=6)
        plt.xlim(boxX[0], boxX[1])
        plt.ylim(boxY[0], boxY[1])

        plt.savefig("./out/test"+str(t)+".png")
        plt.clf()


