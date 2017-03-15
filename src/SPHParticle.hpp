#pragma once

#include "ofMain.h"

struct SPHParticle {
    SPHParticle() {
        vel = pos = lastPos = ofVec3f(0);
        mass = localDensity = localVolume = localPressure = 0;
    }
    
	ofVec3f pos;
    ofVec3f lastPos;
	float mass;
	ofVec3f vel;
	float localDensity;
	float localVolume;
	float localPressure;
    static float supportRadius;
    static float smoothingRadius;
};
