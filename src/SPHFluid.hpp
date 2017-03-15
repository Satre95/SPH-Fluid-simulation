#pragma once

#include "ofMain.h"
#include "SpatialHashTable.hpp"
#include "SPHParticle.hpp"

class SPHFluid {
public:
    //--------------------------------------------
    //MARK: - Public Methods -
    //--------------------------------------------

	SPHFluid(int numParticles = 1000);
    int numParticles;
    const float binSize = SPHParticle::smoothingRadius;
    const int numBins = 10000;
    float h = SPHParticle::smoothingRadius;
    float hRaise3 = pow(h, 3);
    float hRaise4 = pow(h, 4);
    const ofVec3f cubeDims;
    
    
    //--------------------------------------------
    //MARK: Draw fns.
    void drawParticles();
    

    //--------------------------------------------
    //MARK: Update fns
    void update();

private:
    //--------------------------------------------
    //MARK: - Private Methods -
    //--------------------------------------------

    //--------------------------------------------
    //MARK: Private Update fns.
    void updateSHT();
    void updateVBO();
    void updateParticleDensities();
    void computeForces();
    void applyForces();
    
    
    //--------------------------------------------
    //MARK: - Math fns.
    ///Del of kernel function.
    ofVec3f gradientOfKernelFn(SPHParticle & p1, SPHParticle & p2);
    ///Kernel function
    float kernelFn(SPHParticle & p1, SPHParticle & p2);
    ///Helper method to caclulate f(q) for kernel function.
    float helperKernelFn(float dist);
    ///Helper method to calculate f'(q) for kernel function.
    float helperKernelFnDerivative(float dist);

    //--------------------------------------------
    //MARK: - Private Variables
	SpatialHashTable<SPHParticle*> sht;
    std::vector<SPHParticle> particles;
    ofVbo particlesVbo;
    
 
};
