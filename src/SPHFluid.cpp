#include "SPHFluid.hpp"

float SPHParticle::smoothingRadius = 0.025f;
float SPHParticle::supportRadius = 2.0f * smoothingRadius;

SPHFluid::SPHFluid(int num) : sht(SpatialHashTable<SPHParticle*>(binSize, numBins)), numParticles(num)
{
    int particlesPerDim = std::floor(pow(numParticles, 1.0f/3.0f));
    
    //Allocate the particles in a 2x2x2 cube that is at height 4
    for (int z = 0; z < particlesPerDim; z++) {
        for(int y = 0; y < particlesPerDim; y++) {
            for( int x = 0; x < particlesPerDim; x++) {
                SPHParticle particle;
                float xDim = 2.0f * (float)x / (float)particlesPerDim;
                float yDim = 2.0f * (float)y / (float)particlesPerDim;
                float zDim = 2.0f * (float)z / (float)particlesPerDim;
                particle.pos = ofVec3f(xDim, yDim, zDim);
                particle.lastPos = particle.pos;
                particles.push_back(particle);
                //insert into the hash table.
                SPHParticle * particlePtr = &particles.back();
                sht.insert(particles.back().pos, particlePtr);
            }
        }
    }
}


//--------------------------------------------------
//MARK: - Math Functions
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

//--------------------------------------------------
//MARK: - Update Functions
void SPHFluid::update() {
    //1. Update the position data in the hashTable.
    updateSHT();
    //2. Update the particle densities
    updateParticleDensities();
    //3. Compute and sum pressure, viscosity, and gravity forces
    computeForces();
    //4. Apply forces as acceleration and velocity.
    applyForces();
    //5. Update VBO for rendering.
    updateVBO();
}

void SPHFluid::updateVBO() {
    //Extract positions
    vector<ofVec3f> positions(particles.size());
    //Use colors as positions
    vector<ofFloatColor> colors(particles.size());
    
    for(int i = 0; i < particles.size(); i++) {
        SPHParticle & p = particles.at(i);
        positions.at(i) = p.pos;
        colors.at(i) = ofFloatColor(p.pos.x / 2.0f,
                                    p.pos.y / 2.0f,
                                    p.pos.z / 2.0f,
                                    0.8f);
    }
    

    
    particlesVbo.setVertexData(positions.data(), (int)positions.size(), GL_DYNAMIC_DRAW);
    particlesVbo.setColorData(colors.data(), (int)colors.size(), GL_DYNAMIC_DRAW);
}

void SPHFluid::updateSHT() {
    for(SPHParticle & aParticle: particles) {
        //If the hashing the particles current pos does not hash
        //to same key as the last pos, then the particles needs to be
        //moved to a different bin.
        if(!sht.compareKeyHashes(aParticle.lastPos, aParticle.pos)) {
            sht.remove(aParticle.lastPos, &aParticle);
            sht.insert(aParticle.pos, &aParticle);
        }
    }
}

void SPHFluid::updateParticleDensities() {
    
}

void SPHFluid::computeForces() {
    
}

void SPHFluid::applyForces() {
    
}

//--------------------------------------------------
//MARK: - Draw Functions
void SPHFluid::drawParticles() {
    glPointSize(2.0f);
    particlesVbo.draw(GL_POINTS, 0, (int)particles.size());
}
