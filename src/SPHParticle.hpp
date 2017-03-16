#pragma once

#include "ofMain.h"
typedef float pressure;
typedef float density;

struct SPHParticle {
    SPHParticle() {
        vel = pos = lastPos = ofVec3f(0);
        localDensity = localPressure = 0;
    }
    
	ofVec3f pos;
    ofVec3f lastPos;
	ofVec3f vel;
    ofVec3f force;
    
	density localDensity;
	pressure localPressure;
    
    static float supportRadius;
    static float smoothingRadius;
    static float mass;
};
