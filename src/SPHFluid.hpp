#pragma once

#include "ofMain.h"
#include "SpatialHashTable.hpp"
#include "SPHParticle.hpp"
#include "ofxGui.h"

class SPHFluid {
public:
    //--------------------------------------------
    //MARK: - Public Methods -
    //--------------------------------------------

	SPHFluid(int num = 512);
    int numParticles;
    int startNum;
    const float binSize = SPHParticle::smoothingRadius;
    const int numBins = 10000;
    float h = SPHParticle::smoothingRadius;
    float hRaise2 = pow(h, 2);
    float hRaise3 = pow(h, 3);
    float hRaise4 = pow(h, 4);
    ofVec3f cubeDims;
    density restDensity = SPHParticle::mass / hRaise3;
    ofParameter<float> stiffnessConstant;
    ofParameter<float> fps;
    static ofVec3f gravity;
    ofParameter<float> viscosityConstant;
    float timeStep = 0.01f;
    ofVec3f boundingBox;
    
    //--------------------------------------------
    //MARK: Draw fns.
    void drawParticles();
    void drawGui();
    

    //--------------------------------------------
    //MARK: Update fns
    void update();
    void setPosAsColor(bool allow);
    void resetParticles();
    void addParticles();

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
    void detectCollisions();
    void computeForces();
    void applyForces();
    
    
    //--------------------------------------------
    //MARK: Math fns.
    
    //Random point in a sphere of given radius
    ofVec3f randomPointInSphere(float radius);
    
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
    
    ofxPanel panel;
};
