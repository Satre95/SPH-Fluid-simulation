#pragma once

#include "ofMain.h"
#include "SpatialHashTable.hpp"
#include "SPHParticle.hpp"

typedef float pressure;
typedef float density;
typedef float volume;

class SPHFluid {
public:
    //--------------------------------------------
    //MARK: - Public Methods -
    //--------------------------------------------

	SPHFluid(int numParticles = 10000);
    const int numParticles;
    const int binSize = SPHParticle::smoothingRadius;
    const int numBins = 10000;
    
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
	SpatialHashTable<SPHParticle*> buckets;
    std::vector<SPHParticle> particles;
    ofVbo particlesVbo;
    
 
};
