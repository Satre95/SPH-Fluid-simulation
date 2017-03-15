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
    ///Ascertains whether the two positions hash to the same bin.
    bool compareKeyHashes(ofVec3f pos1, ofVec3f pos2);

	void testHash();
private:
	typedef long HashKey;
	SpatialHashTable<T>::HashKey hashPosition(ofVec3f position);

	const ofVec3f PRIMES = ofVec3f(73856093, 19349663, 83492791);
	float binSize;
	int numBins;
	std::unordered_map<HashKey, std::list<T>> bins;
};

template<typename T>
SpatialHashTable<T>::SpatialHashTable()
{
	this->binSize = 5; this->numBins = 500;
}

template<typename T>
SpatialHashTable<T>::SpatialHashTable(float binSize, int numBins)
{
	this->binSize = binSize;
	this->numBins = numBins;
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
