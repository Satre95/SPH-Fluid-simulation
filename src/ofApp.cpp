#include "ofApp.h"

#include "SpatialHashTable.hpp"
//--------------------------------------------------------------
void ofApp::setup(){
    ofSetFrameRate(60);
    camera.setDistance(5.0f);
    camera.setNearClip(0.01f);
    camera.setFarClip(1000.0f);
}

//--------------------------------------------------------------
void ofApp::update(){
    fluid.update();
}

//--------------------------------------------------------------
void ofApp::draw(){
    std::cerr << "Framerate: " << ofGetFrameRate() << std::endl;
    
    camera.begin();
    fluid.drawParticles();
    camera.end();
}

//--------------------------------------------------------------
void ofApp::keyPressed(int key){

}

//--------------------------------------------------------------
void ofApp::keyReleased(int key){

}

//--------------------------------------------------------------
void ofApp::mouseMoved(int x, int y ){

}

//--------------------------------------------------------------
void ofApp::mouseDragged(int x, int y, int button){

}

//--------------------------------------------------------------
void ofApp::mousePressed(int x, int y, int button){

}

//--------------------------------------------------------------
void ofApp::mouseReleased(int x, int y, int button){

}

//--------------------------------------------------------------
void ofApp::mouseEntered(int x, int y){

}

//--------------------------------------------------------------
void ofApp::mouseExited(int x, int y){

}

//--------------------------------------------------------------
void ofApp::windowResized(int w, int h){

}

//--------------------------------------------------------------
void ofApp::gotMessage(ofMessage msg){

}

//--------------------------------------------------------------
void ofApp::dragEvent(ofDragInfo dragInfo){ 

}
