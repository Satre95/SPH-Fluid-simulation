#pragma once

#include "ofMain.h"
#include "SpatialHashTable.hpp"
#include "SPHParticle.hpp"

typedef float pressure;
typedef float density;
typedef float volume;

class SPHFluid {
public:
    
	SPHFluid(int numParticles = 1000);
    
    
private:
    //Kernel function
    ofVec3f gradientOfKernelFn(SPHParticle & p1, SPHParticle & p2);
    float kernelFn(SPHParticle & p1, SPHParticle & p2);
    //Helper method to caclulate f(q) for kernel function.
    float helperKernelFn(float dist);
    float helperKernelFnDerivative(float dist);

    
	SpatialHashTable<SPHParticle*> buckets;
    std::vector<SPHParticle> particles;
    
    
 
};
