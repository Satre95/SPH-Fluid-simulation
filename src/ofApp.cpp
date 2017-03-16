#include "ofApp.h"

#include "SpatialHashTable.hpp"
//--------------------------------------------------------------
void ofApp::setup(){
    ofSetFrameRate(60);
    camera.setDistance(2.0f);
    camera.setNearClip(0.01f);
    camera.setFarClip(1000.0f);
//    paused = true;
}

//--------------------------------------------------------------
void ofApp::update(){
    if(paused) return;
    
    fluid.update();
}

//--------------------------------------------------------------
void ofApp::draw(){
    ofBackground(ofColor::black);
    
    stringstream ss;
    ss << "Frame render time: ";
    ss << ofGetLastFrameTime();
    ofDrawBitmapString(ss.str(), 10, 10);
    
    camera.begin();
    fluid.drawParticles();
    camera.end();
}

//--------------------------------------------------------------
void ofApp::keyPressed(int key){

}

//--------------------------------------------------------------
void ofApp::keyReleased(int key){
    if(key == 'p')
        paused = !paused;
    else if(key == 'c') {
        posAsColor = !posAsColor;
        fluid.setPosAsColor(posAsColor);
    }
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
