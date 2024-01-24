#include <cmath>
#include <iostream>
#include <fstream>
using namespace std;

const float FLPI = M_PI;
const float Uo = 1;
const float lamb = 1;
const float l = 10;
const float aRatio = l/lamb;
const int globalChainSize = int(round(9*aRatio/8));
const float d = l/(sqrt((globalChainSize+1)*(globalChainSize-1)));
const float timeStep = 0.05;
const float Fo = 2;
const int chainNos = 100;
const float fo = 1;
const int wedgeSize=globalChainSize*10;
const float wedgeAngle=FLPI/2;
const float chain_chain_cutoff=l+lamb;
const float chain_wedge_cutoff=l/2+5*l+lamb;
const int saveChunkSize=100;
const string csvName="out.csv";
std::ofstream csvFile;

float* zxycross(float z, float* xy){
        float* tmp = new float[2];
        tmp[0] = -1*xy[1]*z, tmp[1] = xy[0]*z;
        return tmp;
}

float dot(float* x, float* y){
        float tmp = x[0]*y[0]+x[1]*y[1];
        return tmp;
}

float newmod(float a, float b){
        float mod_val = fmod(a, b);
        return mod_val >= 0 ? mod_val : b + mod_val;
}


float* getFrictionCoeff(){
        float* friction = new float[3];
        friction[0] =  2*FLPI/(log(aRatio)-0.207+(0.98/aRatio)-(0.133/(aRatio*aRatio)));
        friction[1] =  2*FLPI/(log(aRatio)+0.839+(0.185/aRatio)+(0.233/(aRatio*aRatio)));
        friction[2] =  FLPI*(aRatio*aRatio)/(3*(log(aRatio)-0.662+(0.917/aRatio)-(0.05/(aRatio*aRatio))));
        return friction;
}

const float* friction = getFrictionCoeff();
const float boxX[2] = {-15*l, 15*l};
const float boxY[2] = {-15*l, 15*l};

class Particle{
        public:
                float pos[2];
                float linvel[2];
};
class SPR{
        public:
                float pos[2];
                float linvel[2];
                float orient[2];
                float force[2];
                float angle;
                float angvel;
                float moment;
                Particle particles[globalChainSize];
                void update_particle(){
                        if (globalChainSize%2!=0){
                                for (int i=-1*int(globalChainSize/2); i<=int(globalChainSize/2); i++){
                                        int j=i;
                                        if (j<0){j+=globalChainSize;}
                                        particles[j].pos[0] = pos[0]+d*i*orient[0];
                                        particles[j].pos[1] = pos[1]+d*i*orient[1];
                                        float* rCOM= new float[2];
                                        rCOM[0] = particles[j].pos[0]-pos[0], rCOM[1] =  particles[j].pos[1]-pos[1];
                                        float* cross=zxycross(angvel, rCOM);
                                        particles[j].linvel[0] = linvel[0]+cross[0];
                                        particles[j].linvel[1] = linvel[1]+cross[1];
                                        delete[] cross; delete[] rCOM;
                                }
                        }
                        else{
                                for (int i=-1*(int(globalChainSize/2)-1); i<=int(globalChainSize/2); i++){
                                        int j=i;
                                        if (j<0){j+=globalChainSize;}
                                        particles[j].pos[0] = pos[0]+d*i*orient[0]-d*orient[0]/2;
                                        particles[j].pos[1] = pos[1]+d*i*orient[1]-d*orient[1]/2;
                                        float* rCOM= new float[2];
                                        rCOM[0] = particles[j].pos[0]-pos[0], rCOM[1] =  particles[j].pos[1]-pos[1];
                                        float* cross=zxycross(angvel, rCOM);
                                        particles[j].linvel[0] = linvel[0]+cross[0];
                                        particles[j].linvel[1] = linvel[1]+cross[1];
                                        delete[] cross; delete[] rCOM;
                                }
                        }
                }

                //Method to update chain parameters given forces are updated already. Call at every timestep
                void update_chain(){
                        linvel[0] = dot(force, orient)*orient[0]/friction[0]+(force[0]-dot(force, orient)*orient[0])/friction[1];
                        linvel[1] = dot(force, orient)*orient[1]/friction[0]+(force[1]-dot(force, orient)*orient[1])/friction[1];
                        angvel = moment/friction[2];
                        pos[0]+=linvel[0]*timeStep;
                        pos[1]+=linvel[1]*timeStep;
                        pos[0]=newmod((pos[0]-boxX[0]), (boxX[1]-boxX[0]))+boxX[0];
                        pos[1]=newmod((pos[1]-boxY[0]), (boxY[1]-boxY[0]))+boxY[0];
                        angle+=angvel*timeStep;
                        if (angle>2*FLPI){
                                angle-=2*FLPI;
                        }
                        if (angle<0){
                                angle+=2*FLPI;
                        }
                        orient[0] = cos(angle), orient[1] = sin(angle);
                }

};
class Wedge{
        public:
                float pos[2];
                float orient[2];
                Particle particles[wedgeSize];
                void update_Wedge(){
                        if (wedgeSize%2!=0){
                                for (int i=-1*int(wedgeSize/2); i<=int(wedgeSize/2); i++){
                                        int j=i;
                                        if (j<0){j+=wedgeSize;}
                                        particles[j].pos[0] = pos[0]+d*i*orient[0];
                                        particles[j].pos[1] = pos[1]+d*i*orient[1];
                                        particles[j].linvel[0] = 0;
                                        particles[j].linvel[1] = 0;
                                }
                        }
                        else{
                                for (int i=-1*(int(wedgeSize/2)-1); i<=int(wedgeSize/2); i++){
                                        int j=i;
                                        if (j<0){j+=wedgeSize;}
                                        particles[j].pos[0] = pos[0]+d*i*orient[0]-d*orient[0]/2;
                                        particles[j].pos[1] = pos[1]+d*i*orient[1]-d*orient[1]/2;
                                        particles[j].linvel[0] = 0;
                                        particles[j].linvel[1] = 0;
                                }
                        }
                }
};

float* getForce (Particle a, Particle b){
        float* tmp = new float[2];
        if (abs(a.pos[0]-b.pos[0])>(boxX[1]-boxX[0])/2){
                (a.pos[0]>b.pos[0])?(b.pos[0]+=boxX[1]-boxX[0]):(a.pos[0]+=boxX[1]-boxX[0]);
        }
        if (abs(a.pos[1]-b.pos[1])>(boxY[1]-boxY[0])/2){
                (a.pos[1]>b.pos[1])?(b.pos[1]+=boxY[1]-boxY[0]):(a.pos[1]+=boxY[1]-boxY[0]);
        }
        float rRel[2] = {a.pos[0]-b.pos[0], a.pos[1]-b.pos[1]};
        float rmod = sqrt(dot(rRel, rRel));
        tmp[0] = ((Uo*exp(-1*rmod/lamb))/rmod)*((1/rmod)+(1/lamb))*rRel[0]/rmod;
        tmp[1] = ((Uo*exp(-1*rmod/lamb))/rmod)*((1/rmod)+(1/lamb))*rRel[1]/rmod;
        return tmp;
}
void SPR_SPR(SPR& a, SPR& b){
	if (abs(a.pos[0]-b.pos[0])>chain_chain_cutoff or abs(a.pos[1]-b.pos[1])>chain_chain_cutoff){return;}
        for (int i=0; i<globalChainSize; i++){
                for (int j=0; j<globalChainSize; j++){
                        float* frc;
                        frc = getForce(a.particles[i], b.particles[j]);
                        a.force[0]+=frc[0]; a.force[1]+=frc[1];
                        b.force[0]-=frc[0]; b.force[1]-=frc[1];
                        a.moment+=(a.particles[i].pos[0]-a.pos[0])*frc[1]-(a.particles[i].pos[1]-a.pos[1])*frc[0];
                        b.moment+=-(b.particles[i].pos[0]-b.pos[0])*frc[1]+(b.particles[i].pos[1]-b.pos[1])*frc[0];
                        delete[] frc;
                }
        }
}
void SPR_Wedge(SPR& a, Wedge& b){
	if (abs(a.pos[0]-b.pos[0])>chain_wedge_cutoff or abs(a.pos[1]-b.pos[1])>chain_wedge_cutoff){return;}
        for (int i=0; i<globalChainSize; i++){
                for (int j=0; j<wedgeSize; j++){
                        float* frc;
                        frc = getForce(a.particles[i], b.particles[j]);
                        a.force[0]+=frc[0]; a.force[1]+=frc[1];
                        a.moment+=(a.particles[i].pos[0]-a.pos[0])*frc[1]-(a.particles[i].pos[1]-a.pos[1])*frc[0];
                        delete[] frc;
                }
        }
}

void addData(float arr[], SPR& a, int t, int i){
    int tOff = t*6*chainNos;
    int chainOff=i*6;
    int totalOff=tOff+chainOff;
    arr[totalOff]=a.pos[0];
    arr[totalOff+1]=a.pos[1];
    arr[totalOff+2]=a.angle;
    arr[totalOff+3]=a.linvel[0];
    arr[totalOff+4]=a.linvel[1];
    arr[totalOff+5]=a.angvel;
}

void saveData(float arr[]){
    int tOff=0, chainOff=0, totalOff=0;
    csvFile.open(csvName, ios::app);
    for (int t=0; t<saveChunkSize; t++){
        tOff = t*6*chainNos;
        for (int i=0; i<chainNos; i++){
            chainOff=i*6;
            totalOff=tOff+chainOff;
            csvFile << arr[totalOff] << ", " << arr[totalOff+1] << ", " << arr[totalOff+2] << ", " << arr[totalOff+3] << ", " << arr[totalOff+4] << ", " << arr[totalOff+5] << ", ";              
        }
        csvFile << endl;
    }
    csvFile.close();

}

int main(){
        SPR Rods[chainNos];
        float csvData[chainNos*6*saveChunkSize];
        for (int i=0; i<chainNos; i++){
                Rods[i].pos[0]=boxX[0]+(i%15)*12, Rods[i].pos[1]=boxY[0]+(i/15)*12;
                Rods[i].linvel[0]=0, Rods[i].linvel[1]=0;
                Rods[i].angle=i;
                Rods[i].angvel=0;
                Rods[i].orient[0]=cos(Rods[i].angle); Rods[i].orient[1]=sin(Rods[i].angle);
                Rods[i].force[0]=Fo*Rods[i].orient[0]; Rods[i].force[1]=Fo*Rods[i].orient[1];
                Rods[i].moment=0;
                Rods[i].update_particle();
        }
	Wedge lowerWedge, upperWedge;
        lowerWedge.orient[0]=cos(-1*wedgeAngle/2); lowerWedge.orient[1]=sin(-1*wedgeAngle/2);
        upperWedge.orient[0]=cos(wedgeAngle/2); upperWedge.orient[1]=sin(wedgeAngle/2);
        lowerWedge.pos[0]=(wedgeSize+1)*d*lowerWedge.orient[0]; lowerWedge.pos[1]=(wedgeSize+1)*d*lowerWedge.orient[1];
        upperWedge.pos[0]=(wedgeSize-1)*d*upperWedge.orient[0]; upperWedge.pos[1]=(wedgeSize-1)*d*upperWedge.orient[1];
        lowerWedge.update_Wedge(); upperWedge.update_Wedge();
        for (int t=0; t<10000; t++){
                for (int i=0; i<chainNos; i++){
                        for (int j=i+1; j<chainNos; j++){
                                SPR_SPR(Rods[i], Rods[j]);
                        }
			SPR_Wedge(Rods[i], lowerWedge);
                        SPR_Wedge(Rods[i], upperWedge);
                        Rods[i].update_chain();
                        Rods[i].update_particle();
                        Rods[i].force[0]=Fo*Rods[i].orient[0]; Rods[i].force[1]=Fo*Rods[i].orient[1];
                        Rods[i].moment=0;
                        addData(csvData, Rods[i], t%saveChunkSize, i);
                }
                if (t%saveChunkSize==saveChunkSize-1){
                    
                    saveData(csvData);
                }
                cout<<t<<endl;
                // cout << Rods[0].pos[0]; cout << ", "; cout << Rods[0].linvel[0]; cout <<", "; cout << Rods[0].angle;
                // cout<<"\n";
        }
        return 0;
}
