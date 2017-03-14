#pragma once

#include "ofMain.h"

struct SPHParticle {
	ofVec3f pos;
	float mass;
	ofVec3f vel;
	float supportRadius;
	float localDensity;
	float localVolume;
	float localPressure;
};
