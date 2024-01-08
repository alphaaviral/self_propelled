import random
import numpy as np
import matplotlib.pyplot as plt
import math
from itertools import combinations

class Particle:
    def __init__(self, pos, linvel):
        self.pos=pos
        self.linvel=linvel

class Chain:
    def __init__(self, pos=None, angle=None):    
        if pos is None:
            self.pos=np.array([random.uniform(-100,100), random.uniform(-100,100)], dtype=np.float64)
        else:
            self.pos=pos

        if angle is None:
            self.angle=random.uniform(0,2*np.pi)
        else:
            self.angle=angle
        # self.angle=random.uniform(0,2*np.pi)
        self.orient= np.array([np.cos(self.angle), np.sin(self.angle)], dtype=np.float64)                    #unit vector, equal to diameter of one particle
        self.linvel=np.array([0,0], dtype=np.float64)
        self.angvel=0                    #vector-wise this is in z
        self.length=chainLength
        self.mass=chainMass
        self.moi=5  #update this
        self.force=Fo*self.orient
        self.moment=0
        self.particles=[]
        for i in range(-1*int(chainLength/2),int(chainLength/2)+1):
            particlePos=self.pos+(i*self.orient)
            particleVel=self.linvel+zxycross(self.angvel, particlePos-self.pos)
            self.particles.append(Particle(particlePos, particleVel))
    
    def updateParticles(self):
        for i in range(len(self.particles)):
            self.particles[i].pos=self.pos+((i-int(len(self.particles)/2))*self.orient)
            self.particles[i].linvel=self.linvel+zxycross(self.angvel, self.particles[i].pos-self.pos)

# cross of a vector in x-y plane with a vector perpendicular to it
def zxycross(z:float, xy:np.ndarray):
    return np.array([-1*xy[1]*z, xy[0]*z])

# force on a due to b
def getForce(a:Particle, b:Particle)->np.ndarray:
    if not isinstance(a, Particle) or not isinstance(b, Particle):
        raise TypeError("force cannot be calculated for non particles")
    rmod=np.linalg.norm(a.pos-b.pos)
    rhat=(a.pos-b.pos)/rmod
    f = ((Uo*math.exp(-1*rmod/lamb))/rmod)*((1/rmod)+(1/lamb))*rhat
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

# def addSelfForce(chain:Chain):
#     if not isinstance(chain, Chain):
#         raise TypeError("force cannot be calculated for non chains")
#     chain.force+=Fo*chain.orient

def updateChainPositions(chain:Chain):
    if not isinstance(chain, Chain):
        raise TypeError("force cannot be calculated for non chains")
    chain.pos+=chain.linvel*timeStep
    # angle = np.arctan2(a.orient[1], a.orient[0])
    chain.angle+=chain.angvel*timeStep
    if chain.angle>2*np.pi:
        chain.angle-=2*np.pi
    if chain.angle<0:
        chain.angle+=2*np.pi
    chain.orient=np.array([np.cos(chain.angle), np.sin(chain.angle)])
    if chain.pos[0]>100 or chain.pos[0]<-100:
        chain.orient[0]*=-1
        chain.angle=np.arctan2(chain.orient[1], chain.orient[0])
    if chain.pos[1]>100 or chain.pos[1]<-100:
        chain.orient[1]*=-1
        chain.angle=np.arctan2(chain.orient[1], chain.orient[0])    

def updateChainVelocities(chain:Chain):
    if not isinstance(chain, Chain):
        raise TypeError("force cannot be calculated for non chains")
    chain.linvel=((np.dot(chain.force, chain.orient)/ft1)*chain.orient)+((chain.force-np.dot(chain.force, chain.orient)*chain.orient)/ft2)
    chain.angvel=chain.moment/fr
    # chain.force/=chain.mass
    # chain.moment/=chain.moi
    # chain.linvel+=chain.force*timeStep
    # chain.angvel+=chain.moment*timeStep

if __name__=="__main__":
    Uo=100
    lamb=1
    chainLength=5
    chainMass=5
    timeStep=0.001
    Fo=200
    chainNos=50
    ft1=1
    ft2=1
    fr=1        #modify these ft and fr values
    # t=0
    if chainLength%2==0:
        raise ValueError("length of chain must be odd")
    chains=[]
    # for i in range(chainNos):
    #     chains.append(Chain())

    chains.append(Chain(pos=np.array([0,0], dtype=np.float64), angle=np.pi/3))
    chains.append(Chain(pos=np.array([0,10], dtype=np.float64), angle=5*np.pi/3))
    
    # comb = combinations(chains, 2)

    for t in range(0, 1000):
        print(t)
        comb = combinations(chains, 2)

        for i in list(comb):
            updateChainForces(i[0], i[1])
        for chain in chains:
            # addSelfForce(i)
            updateChainVelocities(chain)
            updateChainPositions(chain)
            chain.updateParticles()
            chain.force=Fo*chain.orient
            chain.moment=0
            for particle in chain.particles:
                plt.plot(particle.pos[0], particle.pos[1], 'r', marker=".", markersize=5)
        plt.xlim(-10, 50)
        plt.ylim(-10,50)
        

        # for particle in chains[0].particles:
        #     plt.plot(particle.pos[0], particle.pos[1], 'ro')
        # for particle in chains[1].particles:
        #     plt.plot(particle.pos[0], particle.pos[1], 'bo')

            # plt.plot(chain.)
        # plt.plot(chains[0].pos[0], chains[0].pos[1], 'ro')
        # plt.plot(chains[1].pos[0], chains[1].pos[1], 'bo')
        # plt.show()
        # if t%10==0:
        plt.savefig("./out/test"+str(t)+".png")
        plt.clf()


