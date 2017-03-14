#pragma once

#include "ofMain.h"
#include "SpatialHashTable.hpp"
#include "SPHParticle.hpp"

class SPHFluid {
public:
	SPHFluid();
private:
	SpatialHashTable<SPHParticle*> buckets;
};