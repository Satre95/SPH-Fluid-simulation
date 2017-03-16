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
    float hRaise2 = pow(h, 2);
    float hRaise3 = pow(h, 3);
    float hRaise4 = pow(h, 4);
    ofVec3f cubeDims;
    density restDensity = SPHParticle::mass / hRaise3;
    float stiffnessConstant = 5.0f;
    static ofVec3f gravity;
    float viscosityConstant;
    
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
    void updateParticlePressure();
    void computeForces();
    void applyForces();
    
    
    //--------------------------------------------
    //MARK: Math fns.
    
    ///∇ of kernel function.
    ofVec3f gradientOfKernelFn(SPHParticle & p1, SPHParticle & p2);
    ///Kernel function
    float kernelFn(SPHParticle & p1, SPHParticle & p2);
    ///Helper method to caclulate f(q) for kernel function.
    float helperKernelFn(float dist);
    ///Helper method to calculate f'(q) for kernel function.
    float helperKernelFnDerivative(float dist);
    ///Calculates Del of given quantity (a_i) using kernel function
    ofVec3f gradientOfQuantityHelperFn(float a_i, float a_j, SPHParticle & p_i, SPHParticle & p_j);
    ///Calculates Del^2 * A_i using kernel function
    float gradientSquaredOfQuantityHelperFn(float a_i_j, SPHParticle & p_i, SPHParticle & p_j);
    ///Calculates Del^2 * A_i using kernel function
    ofVec3f gradientSquaredOfQuantityHelperFn(ofVec3f a_i_j, SPHParticle & p_i, SPHParticle & p_j);
    

    //--------------------------------------------
    //MARK: Private Variables
	SpatialHashTable<SPHParticle*> sht;
    std::vector<SPHParticle> particles;
    ofVbo particlesVbo;
    
 
};
