#pragma once

#include "ofMain.h"

struct SPHParticle {
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
