#include "SPHFluid.hpp"

SPHFluid::SPHFluid() : buckets(SpatialHashTable<SPHParticle*>(5, 1000))
{}