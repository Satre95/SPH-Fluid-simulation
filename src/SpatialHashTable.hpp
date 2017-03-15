#pragma once

#include "ofMain.h"
#include <unordered_map>
#include <list>
#include <vector>


template <typename T>
class SpatialHashTable
{
public:
	SpatialHashTable();
	SpatialHashTable(float binSize, int numBins);

	void insert(ofVec3f position, T && value);
    void insert(ofVec3f position, T & value);
	void remove(ofVec3f position, T & value);
    void remove(ofVec3f position, T && value);
    bool exists(ofVec3f position);
    
    ///Ascertains whether the two positions hash to the same bin.
    bool compareKeyHashes(ofVec3f pos1, ofVec3f pos2);
    
    ///Function to iterate over all bins and process the data inside them.
//    void processDataInBuckets(const std::function< void( std::list<T> )> f);
    
    ///Function to return the bins represent the spaces near the given position, if any
    std::vector<std::reference_wrapper<std::list<T>>> getNeighboringBuckets(ofVec3f pos);
private:
	typedef long HashKey;
	SpatialHashTable<T>::HashKey hashPosition(ofVec3f position);
    std::list<T> * find(ofVec3f pos);

	const ofVec3f PRIMES = ofVec3f(73856093, 19349663, 83492791);
	float binSize;
	int numBins;
    float diagDist;
    float diagDist2d;
	std::unordered_map<HashKey, std::list<T>> bins;
};

template<typename T>
SpatialHashTable<T>::SpatialHashTable()
{
	this->binSize = 5; this->numBins = 500;
    diagDist = sqrt(pow(binSize, 2) * 3.0f);
    diagDist2d = sqrt(pow(binSize, 2) * 2.0f);
}

template<typename T>
SpatialHashTable<T>::SpatialHashTable(float binSize, int numBins)
{
	this->binSize = binSize;
	this->numBins = numBins;
    diagDist = sqrt(pow(binSize, 2) * 3.0f);
    diagDist2d = sqrt(pow(binSize, 2) * 2.0f);
}

template<typename T>
void SpatialHashTable<T>::insert(ofVec3f position, T && value) {
	//Get the hash key for the given position
	HashKey key = hashPosition(position);

	//if the bin that this position hashes to has not been used before, need to create it first.
	if (bins.find(key) == bins.end())
		bins.insert(std::make_pair(key, std::list<T>()));

	//Insert into the list at the given key
	std::list<T> & list = bins.at(key);
	list.push_back(value);
}

template<typename T>
void SpatialHashTable<T>::insert(ofVec3f position, T & value){
    //Get the hash key for the given position
    HashKey key = hashPosition(position);
    
    //if the bin that this position hashes to has not been used before, need to create it first.
    if (bins.find(key) == bins.end())
        bins.insert(std::make_pair(key, std::list<T>()));
    
    //Insert into the list at the given key
    std::list<T> & list = bins.at(key);
    list.push_back(value);
}

template <typename T>
void SpatialHashTable<T>::remove(ofVec3f position, T & value) {
	//Get the hash key for the given position
	HashKey key = hashPosition(position);

	//If the given position does not map to a slot in the map, return
	if (bins.find(key) == bins.end()) {
		std::cerr << "ERROR: Trying to remove item for position that is not stored in hash table" << std::endl;
		std::cerr << "Offending Position: " << position << std::endl;
		return;
	}

	std::list<T> & items = bins.at(key);
	items.remove(value);
}

template <typename T>
void SpatialHashTable<T>::remove(ofVec3f position, T && value) {
    //Get the hash key for the given position
    HashKey key = hashPosition(position);
    
    //If the given position does not map to a slot in the map, return
    if (bins.find(key) == bins.end()) {
        std::cerr << "ERROR: Trying to remove item for position that is not stored in hash table" << std::endl;
        std::cerr << "Offending Position: " << position << std::endl;
        return;
    }
    
    std::list<T> & items = bins.at(key);
    items.remove(value);
}

template<typename T>
typename SpatialHashTable<T>::HashKey SpatialHashTable<T>::hashPosition(ofVec3f position) {
	HashKey part1 = std::floor(position.x / binSize) * PRIMES.x;
	HashKey part2 = std::floor(position.y / binSize) * PRIMES.y;
	HashKey part3 = std::floor(position.z / binSize) * PRIMES.z;

	return (part1 ^ part2 ^ part3) % numBins;
}

template<typename T>
bool SpatialHashTable<T>::compareKeyHashes(ofVec3f pos1, ofVec3f pos2) {
    return hashPosition(pos1) == hashPosition(pos2);
}

//template <typename T>
//void SpatialHashTable<T>::processDataInBuckets(const std::function< void( std::list<T> )> f) {
//    for(auto itr = bins.begin(); itr != bins.end(); itr++) {
//        f(itr->second);
//    }
//}

template <typename T>
bool SpatialHashTable<T>::exists(ofVec3f position)
{ return bins.find(hashPosition(position)) != bins.end(); }

template <typename T>
std::list<T> * SpatialHashTable<T>::find(ofVec3f pos) {
    if(bins.find(hashPosition(pos)) == bins.end())
        return nullptr;
    else
        return &(*(bins.find(hashPosition(pos)))).second;
}

template <typename T>
std::vector<std::reference_wrapper<std::list<T>>> SpatialHashTable<T>::getNeighboringBuckets(ofVec3f pos) {
    vector<std::reference_wrapper<std::list<T>>> lists;
    
    //1. go over single dimensions
    //  1.a go over x dimension
    if(exists(pos + ofVec3f(binSize, 0, 0))) lists.push_back(*find(pos + ofVec3f(binSize, 0, 0)));
    if(exists(pos - ofVec3f(binSize, 0, 0))) lists.push_back(*find(pos - ofVec3f(binSize, 0, 0)));
    //  1.b go over y dimension
    if(exists(pos + ofVec3f(0, binSize, 0))) lists.push_back(*find(pos + ofVec3f(0, binSize, 0)));
    if(exists(pos - ofVec3f(0, binSize, 0))) lists.push_back(*find(pos - ofVec3f(0, binSize, 0)));
    //  1.c go over z dimension
    if(exists(pos + ofVec3f(0, 0, binSize))) lists.push_back(*find(pos + ofVec3f(0, 0, binSize)));
    if(exists(pos - ofVec3f(0, 0, binSize))) lists.push_back(*find(pos - ofVec3f(0, 0, binSize)));
    
    //2. go over diagonal cells
    //  2.a right top front  and left bottom back
    ofVec3f diag(diagDist);
    if(exists(pos + diag)) lists.push_back(*find(pos + diag));
    if(exists(pos - diag)) lists.push_back(*find(pos - diag));
    //  2.b right top back  and left bottom front
    diag = ofVec3f(diagDist, diagDist, -diagDist);
    if(exists(pos + diag)) lists.push_back(*find(pos + diag));
    if(exists(pos - diag)) lists.push_back(*find(pos - diag));
    //  2.c left top back and right bottom front
    diag = ofVec3f(-diagDist, diagDist, -diagDist);
    if(exists(pos + diag)) lists.push_back(*find(pos + diag));
    if(exists(pos - diag)) lists.push_back(*find(pos - diag));
    //  2.d left top front and right bottom back
    diag = ofVec3f(-diagDist, diagDist, diagDist);
    if(exists(pos + diag)) lists.push_back(*find(pos + diag));
    if(exists(pos - diag)) lists.push_back(*find(pos - diag));
    
    //3. go over adjacent cells (diagonal cells in same plane)
    //  3.a XY plane
    ofVec3f planeDiag(diagDist2d, diagDist2d, 0);
    if(exists(pos + planeDiag)) lists.push_back(*find(pos + planeDiag));
    if(exists(pos - planeDiag)) lists.push_back(*find(pos - planeDiag));
    planeDiag.x *= -1.0f;
    if(exists(pos + planeDiag)) lists.push_back(*find(pos + planeDiag));
    if(exists(pos - planeDiag)) lists.push_back(*find(pos - planeDiag));
    //  3.b YZ Plane
    planeDiag = ofVec3f(0, diagDist2d, diagDist2d);
    if(exists(pos + planeDiag)) lists.push_back(*find(pos + planeDiag));
    if(exists(pos - planeDiag)) lists.push_back(*find(pos - planeDiag));
    planeDiag.y *= -1.0f;
    if(exists(pos + planeDiag)) lists.push_back(*find(pos + planeDiag));
    if(exists(pos - planeDiag)) lists.push_back(*find(pos - planeDiag));
    //  3.c XZ Plane
    planeDiag = ofVec3f(diagDist2d, 0, diagDist2d);
    if(exists(pos + planeDiag)) lists.push_back(*find(pos + planeDiag));
    if(exists(pos - planeDiag)) lists.push_back(*find(pos - planeDiag));
    planeDiag.z *= -1.0f;
    if(exists(pos + planeDiag)) lists.push_back(*find(pos + planeDiag));
    if(exists(pos - planeDiag)) lists.push_back(*find(pos - planeDiag));
    
    //4. Finally, add our bucket
    if(exists(pos)) lists.push_back(*find(pos));
    
    return lists;
}
