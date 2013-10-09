#ifndef THREEDCLUSTERS_H
#define THREEDCLUSTERS_H

// Lattice dimensions
int Ns=40;
int leng1=Ns, leng2=Ns, leng3=Ns, nclusters=0;
int nconfig=1, selconfig=0;

string f3dlistname;

bool showboxes=false;

// angle for rotation for the camera direction
// float anglex = 0.0f;float angley = 0.0f;
float anglex = 0.0f;float angley = 45.0f;

// actual vector representing the camera's direction
float lx=0.0f,lz=-1.0f, ly = 0.0f;

// Initial XZ position of the camera
float x=0.0f, z=(double)Ns+25, y = 0.00f;

// FPS counter and frame limiter stuff
unsigned frameCount = 0;
double usedTime=0;
double pTime=0;
// const double FREQ=10; // Hz
const double FREQ=60; // Hz
const double TIMERMSECS=1000*1/(double)FREQ;
const double dt = TIMERMSECS/(double)1000;

float deltaAnglex = 0.0f;float deltaAngley = 0.0f;
float deltaMove = 0;
int xOrigin = -1; int yOrigin = -1; 

// Stuff which define sphere properties
double sphereradius=0.35; // Sphere radius
const int sphereSlices=6; // Sphere slices around Z axis
const int sphereStacks=6; //Sphere stacks/slices along the z axis

// double pointsize=8.5; // Sphere radius
double pointsize=4.0; // Sphere radius

// double alpha=0.1; // alpha value used for the sphere/points
double alpha=0.45; // alpha value used for the sphere/points

// Size of the scene/box around the scene 
double minx=-Ns - 0.5, maxx=Ns + 0.5,
      	miny=-Ns - 0.5,maxy=Ns + 0.5,
	minz=-Ns - 0.5,maxz=Ns + 0.5;

// The currently highlighted cluster
int cnt=0;

void cluster3input(int config);

// Will be set to true if we draw only the selected cluster
bool onecluster=false;

#endif //THREEDCLUSTERS_H
