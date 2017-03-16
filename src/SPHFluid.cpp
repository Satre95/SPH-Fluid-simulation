#include "SPHFluid.hpp"
#include <list>
#include <math.h>

float SPHParticle::smoothingRadius = 0.025f;
float SPHParticle::supportRadius = 2.0f * SPHParticle::smoothingRadius;
float SPHParticle::mass = 0.1f;
ofVec3f SPHFluid::gravity(0.0f, -9.8f, 0.0f);

SPHFluid::SPHFluid(int num) : sht(SpatialHashTable<SPHParticle*>(binSize, numBins))
{
    int particlesPerDim = std::floor(pow(num, 1.0f/3.0f));
    numParticles = pow(particlesPerDim, 3);
    cubeDims = ofVec3f(particlesPerDim * SPHParticle::smoothingRadius);
    
    //Allocate the particles in a 2x2x2 cube that is at height 4
    for (int z = 0; z < particlesPerDim; z++) {
        for(int y = 0; y < particlesPerDim; y++) {
            for( int x = 0; x < particlesPerDim; x++) {
                SPHParticle particle;
                float xDim = cubeDims.x * (float)x / (float)particlesPerDim;
                float yDim = cubeDims.y * (float)y / (float)particlesPerDim;
                float zDim = cubeDims.z * (float)z / (float)particlesPerDim;
                particle.pos = ofVec3f(xDim, yDim, zDim);
                particle.lastPos = particle.pos;
                particles.push_back(particle);
                //insert into the hash table.
                SPHParticle * particlePtr = &particles.back();
                sht.insert(particlePtr->pos, particlePtr);
            }
        }
    }
}


//--------------------------------------------------
//MARK: - Math Functions
float SPHFluid::kernelFn(SPHParticle & p1, SPHParticle & p2) {
    ofVec3f pos1 = p1.pos;
    ofVec3f pos2 = p2.pos;
    return (1.0f / hRaise3) * helperKernelFn(pos1.distance(pos2) / h);
}

ofVec3f SPHFluid::gradientOfKernelFn(SPHParticle & p1, SPHParticle & p2) {
    ofVec3f pos1 = p1.pos;
    ofVec3f pos2 = p2.pos;
    float q = pos1.distance(pos2) / h;
    
    float term1 = 1 / hRaise4;
    float term2 = helperKernelFnDerivative(q) / q;
    return term1 * term2 * (pos1 - pos2);
}

float SPHFluid::helperKernelFn(float q) {
    if(q <= 0.00001f) return 0;
    
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
    if(q <= 0.00001f) return 0;
    
    if(q >=0.0f && q < 1.0f )
        return (3.0f / (2.0f * PI)) * (-2.0f * q + 1.5f * pow(q, 2));
    else if(q >= 1.0f && q < 2.0f)
        return (3.0f / (2.0f * PI)) * (-0.5f * pow((2 - q), 2));
    else
        return 0.0f;
}

ofVec3f SPHFluid::gradientOfQuantityHelperFn(float a_i, float a_j, SPHParticle & p_i, SPHParticle & p_j) {
    if(isnan(p_i.localDensity) || p_i.localDensity < 0.0001f ||
       isnan(p_j.localDensity) || p_j.localDensity < 0.0001f)
       return ofVec3f(0);
    
    float term2 = (a_i / pow(p_i.localDensity, 2)) + (a_j / pow(p_j.localDensity, 2));
    ofVec3f term3 = gradientOfKernelFn(p_i, p_j);;
    return p_j.mass * term2 * term3;
}

float SPHFluid::gradientSquaredOfQuantityHelperFn(float a_i_j, SPHParticle & p_i, SPHParticle & p_j) {
    float m_j = p_j.mass;
    density rho_j = p_j.localDensity;
    ofVec3f x_i_j(p_i.pos - p_j.pos);
    
    float term1 = m_j / rho_j;
    //term 2 is a_i_j
    float term3 = x_i_j.dot(gradientOfKernelFn(p_i, p_j));
    float term4 = x_i_j.dot(x_i_j) + 0.01f * hRaise2;
    
    return term1 * a_i_j * term3 / term4;
}

ofVec3f SPHFluid::gradientSquaredOfQuantityHelperFn(ofVec3f a_i_j, SPHParticle & p_i, SPHParticle & p_j) {
    float m_j = p_j.mass;
    density rho_j = p_j.localDensity;
    ofVec3f x_i_j(p_i.pos - p_j.pos);
    
    float term1 = m_j / rho_j;
    //term 2 is a_i_j
    float term3 = x_i_j.dot(gradientOfKernelFn(p_i, p_j));
    float term4 = x_i_j.dot(x_i_j) + 0.01f * hRaise2;
    
    return term1 * a_i_j * term3 / term4;
}

//--------------------------------------------------
//MARK: - Update Functions
void SPHFluid::update() {
    //1. Update the position data in the hashTable.
    updateSHT();
    //2. Update the particle densities
    updateParticleDensities();
    //3. Update the particle pressures
    updateParticlePressure();
    //4. Compute and sum pressure, viscosity, and gravity forces
    computeForces();
    //5. Apply forces as acceleration and velocity.
    applyForces();
    //6. Update VBO for rendering.
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
        colors.at(i) = ofFloatColor(p.pos.x / cubeDims.x,
                                    p.pos.y / cubeDims.y,
                                    p.pos.z / cubeDims.z,
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

/** 
 This function updates the densities of all particles based on their neighboring particles
 * 1. Compute density from neighboring particles in the same bin
 * 2. add in density from particles in neighboring bins
 */
void SPHFluid::updateParticleDensities() {
    for(SPHParticle & p_i: particles) {
        float rho = 0;
        //Get the neighboring buckets. note that this includes our bucket
        auto neighborBuckets = sht.getNeighboringBuckets(p_i.pos);
        
        for(auto aBucket: neighborBuckets) {
            //For some neighboring bucket, iterate over the particles
            for(auto aNeighbor: aBucket.get()) {
                if(aNeighbor == &p_i) continue;
                
                if(p_i.pos.distance(aNeighbor->pos) < 2.0f * h)
                    rho += aNeighbor->mass * kernelFn(p_i, *aNeighbor);
            }
        }
        p_i.localDensity = rho;
    }
}

void SPHFluid::updateParticlePressure() {
    for(auto & p_i:particles) {
        if(abs(p_i.localDensity) < 0.00001f || isnan(p_i.localDensity)) continue;
        
        p_i.localPressure = stiffnessConstant * (pow(p_i.localDensity / restDensity, 7) - 1.0f);
    }
}

/**
 * This function updates the forces on each particle by calculating the net force based 
 * off the pressure, viscosity, and gravity forces
 */
void SPHFluid::computeForces() {
    for(auto & p_i: particles) {
        auto neighborBuckets = sht.getNeighboringBuckets(p_i.pos);
        ofVec3f pressureForce(0);
        ofVec3f viscosityForce(0);
        
        for(auto aBucket: neighborBuckets) {
            for(auto p_j: aBucket.get()) {
                if(p_j == &p_i) continue;
                
                //1. Compute pressure force
                pressureForce += gradientOfQuantityHelperFn(p_i.localPressure, p_j->localPressure, p_i, *p_j);
                
                //2. Compute viscosity force
                
            }
        }
        
        //Need to multiply in mass to finish calculating pressure force
        pressureForce *= -p_i.mass;
        ofVec3f gravityForce(p_i.mass * gravity);
        
        //Sum up to a net force
        ofVec3f netForce = pressureForce + viscosityForce + gravityForce;
        p_i.force = netForce;
    }
}

//TODO: Implement fn.
void SPHFluid::applyForces() {
    
}

//--------------------------------------------------
//MARK: - Draw Functions
void SPHFluid::drawParticles() {
    glPointSize(4.0f);
    particlesVbo.draw(GL_POINTS, 0, (int)particles.size());
}
