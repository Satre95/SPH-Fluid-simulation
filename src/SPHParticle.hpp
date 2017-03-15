#pragma once

#include "ofMain.h"
typedef float pressure;
typedef float density;
typedef float volume;

struct SPHParticle {
    SPHParticle() {
        vel = pos = lastPos = ofVec3f(0);
        localDensity = localVolume = localPressure = 0;
    }
    
	ofVec3f pos;
    ofVec3f lastPos;
	ofVec3f vel;
	density localDensity;
	volume localVolume;
	pressure localPressure;
    static float supportRadius;
    static float smoothingRadius;
    static float mass;
};
