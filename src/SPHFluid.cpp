#include "SPHFluid.hpp"

float SPHParticle::smoothingRadius = 0.025f;
float SPHParticle::supportRadius = 2.0f * smoothingRadius;

SPHFluid::SPHFluid(int numParticles) : buckets(SpatialHashTable<SPHParticle*>(5, 1000))
{
    particles.reserve(numParticles);
}

float SPHFluid::kernelFn(SPHParticle & p1, SPHParticle & p2) {
    ofVec3f pos1 = p1.pos;
    ofVec3f pos2 = p2.pos;
    return (1.0f / pow(SPHParticle::smoothingRadius, 3)) * helperKernelFn(pos1.distance(pos2));
}

ofVec3f SPHFluid::gradientOfKernelFn(SPHParticle & p1, SPHParticle & p2) {
    ofVec3f pos1 = p1.pos;
    ofVec3f pos2 = p2.pos;
    float q = pos1.distance(pos2);
    
    float term1 = 1 / pow(SPHParticle::smoothingRadius, 3+1);
    float term2 = helperKernelFnDerivative(q) / q;
    return term1 * term2 * (pos1 - pos2);
    
}

float SPHFluid::helperKernelFn(float q) {
    
    if(q >= 0.0f && q < 1.0f) {
        return (3.0f / (2.0f * PI)) * (2.0f / 3.0f - pow(q, 2) + 0.5f * pow(q, 3));
    }
    else if(q >= 1.0f && q < 2.0f) {
        return (3.0f / (2.0f * PI)) * (1.0f / 6.0f) * pow(2.0f - q, 3);
    }
    else
    {
        return 0.0f;
    }
}

float SPHFluid::helperKernelFnDerivative(float q) {
    if(q >=0.0f && q < 1.0f )
        return (3.0f / (2.0f * PI)) * (-2.0f * q + 1.5f * pow(q, 2));
    else if(q >= 1.0f && q < 2.0f)
        return (3.0f / (2.0f * PI)) * (-0.5f * pow((2 - q), 2));
    else
        return 0.0f;
}
