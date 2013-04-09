#ifdef __APPLE__
	#include <GLUT/glut.h>
#else
	#include <GL/glut.h>
#endif

#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <vector>
#include <iostream>
#include <fstream>

using namespace std;

// angle for rotation for the camerca direction
float angle = 0.0f;

// actual vector representing the camera's direction
float lx=0.0f,lz=-1.0f, ly = 0.0f;

// XZ position of the camera
float x=0.0f, z=20.0f, y = 0.00f;

unsigned frameCount = 0;

double vtemp=0;

float deltaAngle = 0.0f;
float deltaMove = 0;
int xOrigin = -1; 

const double FREQ=60; // Hz
const double TIMERMSECS=1000*1/(double)FREQ;

const double dt = TIMERMSECS/1000;;

const int dim=3;

// Stuff for the spheres
const double sphereradius=0.5; // Sphere radius
const int sphereSlices=12; // Sphere slices around Z axis
const int sphereStacks=12; //Sphere stacks/slices along the z axis

double minx=-8.5, maxx=8.5,
      	miny=-8.5,maxy=8.5,
	minz=-8.5,maxz=8.5;

int mainWindow;

vector<vector<vector<int > > > lpoints;

vector<double> red;
vector<double> green;
vector<double> blue;

int Ns=8, Nt=4;
int leng1=Ns, leng2=Ns, leng3=Ns, leng4=Nt, nclusters=0;

void cluster3input(string f3dname){
	ifstream f3d;
	int is, cleng1, cleng2, cleng3, cleng4;
	f3d.open(f3dname.c_str());
	if(f3d.is_open()){
		f3d >> cleng1 >> cleng2 >> cleng3 >> cleng4;
		f3d >> nclusters;
		
		if(cleng1 != leng1 || cleng2 != leng2 || cleng3 != leng3 || cleng4 != leng4)
			cout << "WARNING: Lattice size missmatch!" << endl;
		
		int ii1, ii2, ii3, cluster;
		for(int i1=0;i1<leng1;i1++)
		for(int i2=0;i2<leng2;i2++)
		for(int i3=0;i3<leng3;i3++){
			f3d >> ii1 >> ii2 >> ii3 >> cluster;
			lpoints.at(ii1).at(ii2).at(ii3) = cluster;
		}
		f3d.close();
	}else{
		cout << "WARNING: Could not open 3dcluster file!" << endl;
	}
}

void mouseMove(int x, int y) {
        // this will only be true when the left button is down
        if (xOrigin >= 0) {

                // update deltaAngle
                deltaAngle = (x - xOrigin) * 0.001f;

                // update camera's direction
                lx = sin(angle + deltaAngle);
                lz = -cos(angle + deltaAngle);

                glutSetWindow(mainWindow);
                glutPostRedisplay();
        }   
}

void mouseButton(int button, int state, int x, int y) {
        // only start motion if the left button is pressed
        if (button == GLUT_LEFT_BUTTON) {
                // when the button is released
                if (state == GLUT_UP) {
                        angle += deltaAngle;
                        deltaAngle = 0.0f;
                        xOrigin = -1;
                }
                else  {// state = GLUT_DOWN
                        xOrigin = x;
                }
        }
}

void processNormalKeys(unsigned char key, int x, int y){
	if(key == 27)
		exit(0);
/*	else if (key == 'r'){
		int mod = glutGetModifiers();
		if (mod == GLUT_ACTIVE_ALT)
			red = 0.0;
		else
			red = 1.0;
	} */
}

void pressKey(int key, int x, int y){
        switch (key) {
                case GLUT_KEY_UP : deltaMove = 1.5f; break;
                case GLUT_KEY_DOWN : deltaMove = -1.5f; break;
        }   
        glutSetWindow(mainWindow);
        glutPostRedisplay();

}

void releaseKey(int key, int x, int y) {
	switch (key) {
		case GLUT_KEY_UP :
		case GLUT_KEY_DOWN : deltaMove = 0;break;
	}
}

void computePos(float deltaMove) {
	        x += deltaMove * lx * 0.1f;
	        z += deltaMove * lz * 0.1f;
}

void calcSphereColorTest(double &red, double &green, double &blue, int is){
	double s=1, v=1, h;
	
	h = is/(double)(leng1*leng2*leng3)*(double)360;

	h = h/(double)60;			// sector 0 to 5
	int i = floor( h );
	double f = h - i;			// factorial part of h
	double p = v * ( 1 - s );
	double q = v * ( 1 - s * f );
	double t = v * ( 1 - s * ( 1 - f ) );
	
	red=0; green=0; blue=0;
	
	switch( i ) {
		case 0:
			red=v; green=t; blue=p;
			break;
	
		case 1:
			red=q; green=v; blue=p;
			break;
	
		case 2:
			red=p; green=v; blue=t;
			break;
	
		case 3:
			red=p; green=q; blue=v;
			break;
	
		case 4:
			red=t; green=p; blue=v;
			break;
	
		default:
			red=v; green=p; blue=q;
			break;
	}
	
	cout << red << " " << green << " " << blue << endl;
}

void calcSphereColor(double &red, double &green, double &blue, int cluster){
	double s=1, v=1, h;
	
	h = cluster/(double)nclusters*(double)360;

	h = h/(double)60;			// sector 0 to 5
	int i = floor( h );
	double f = h - i;			// factorial part of h
	double p = v * ( 1 - s );
	double q = v * ( 1 - s * f );
	double t = v * ( 1 - s * ( 1 - f ) );
	
	red=0; green=0; blue=0;
	
	switch( i ) {
		case 0:
			red=v; green=t; blue=p;
			break;
	
		case 1:
			red=q; green=v; blue=p;
			break;
	
		case 2:
			red=p; green=v; blue=t;
			break;
	
		case 3:
			red=p; green=q; blue=v;
			break;
	
		case 4:
			red=t; green=p; blue=v;
			break;
	
		default:
			red=v; green=p; blue=q;
			break;
	}
	
	// cout << red << " " << green << " " << blue << endl;
}

void drawSquare(){
	glColor3f(1.0, 1.0, 1.0);
	glutWireCube(maxx-minx);
//	glutSolidCube(maxx-minx);
//	glBegin(GL_LINE_LOOP);
//		glColor3f(1.0, 0.0, 0.0);
//		glVertex3f(maxx,maxy,0);
//		glVertex3f(maxx,miny,0);
//		glVertex3f(minx,miny,0);
//		glVertex3f(minx,maxy,0);
//	glEnd();
}

void drawSphere(double x, double y, double z){
	glPushMatrix();
		glTranslatef(x-(double)Ns/2.0+0.5,y-(double)Ns/2.0+0.5,z-(double)Ns/2.0+0.5);
		// glRotatef(90.0,0.0,1.0,0.0);
		glutSolidSphere(sphereradius,sphereSlices,sphereStacks);
	glPopMatrix();
}

void drawBouncingPoint() {
	double x=0, crandx;
	double y=0, crandy;
	double z=0, crandz;

	drawSquare();

       	for(int i1=0;i1<leng1;i1++)
       	for(int i2=0;i2<leng2;i2++)
       	for(int i3=0;i3<leng3;i3++){
       		int is = i1 + i2*leng1 + i3*leng1*leng2;
		glColor4f(red[lpoints[i1][i2][i3]], green[lpoints[i1][i2][i3]], blue[lpoints[i1][i2][i3]], 0.6); // red, green, blue
		drawSphere(i1,i2,i3);
	}
}

double usedTime=0;
double pTime=0;

// Simple render function
void renderScene(int value){
	double currentTime = glutGet(GLUT_ELAPSED_TIME);
	double uTime=0; //used to modify time function call
	double fps;

	// Clear color and depth buffers
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	// Reset transformations
	glLoadIdentity();

	// Set the camera
        gluLookAt(x, y, z,
                 x + lx,y + ly,z + lz, 
                 0.0f,1.0f,0.0f);

        if(deltaMove) {
	                computePos(deltaMove);
	                glutSetWindow(mainWindow);
	                glutPostRedisplay();
        }

	glRotatef(angle, 0.0f, 1.0f, 0.0f);

	drawBouncingPoint();

	angle+=0.4f;

	usedTime += currentTime-pTime;
	pTime=currentTime;

	if(usedTime > 1000){
		fps=frameCount / (usedTime/(double)1000 );
		char* TempString = (char*)
		    malloc(512);
	 
		sprintf(
		    TempString,
		    "%f Frames Per Second",
		    fps 
		);  
	 
		glutSetWindowTitle(TempString);
		free(TempString);
		frameCount=0;
		usedTime=0;
	}else{
		frameCount++;
	}
	
	glutSwapBuffers();
	
	uTime = glutGet(GLUT_ELAPSED_TIME) - currentTime;
	if(uTime>TIMERMSECS){
		uTime=TIMERMSECS;
	}
	
	glutTimerFunc(TIMERMSECS-uTime, renderScene, 1);
}

void changeSize(int w, int h){
	if(h == 0)
		h = 1;

	float ratio = 1.0*w/h;

	// glPointSize( 6.0*ratio );
	
	// Projection matrix
	glMatrixMode(GL_PROJECTION);

	// Reset matrix
	glLoadIdentity();

	// Set viewpoint to be the entire window
	glViewport(0, 0, w, h);

	// Set the correct perspective
	// gluPerspective(field of view angle in yz plane, 
	// 			ratio, near, far clipping planes)
	gluPerspective(60,ratio,1,100);

	// Get back to the Modelview matrix
	glMatrixMode(GL_MODELVIEW);
}

int init(){;
	// OpenGL stuff
	glEnable(GL_DEPTH_TEST);
	glEnable( GL_POINT_SMOOTH );
	glEnable( GL_BLEND );
	glBlendFunc( GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA );
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glColorMaterial ( GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE );
	glEnable ( GL_COLOR_MATERIAL );

	glutIgnoreKeyRepeat(1);
	// Register callbacks
 	// glutDisplayFunc(renderScene);
	glutReshapeFunc(changeSize);
	// glutIdleFunc(renderScene);
	glutTimerFunc(0, renderScene, 0);
	glutKeyboardFunc(processNormalKeys);
        glutSpecialFunc(pressKey);
        glutSpecialUpFunc(releaseKey);
	glutMouseFunc(mouseButton);
	glutMotionFunc(mouseMove);
	
	// Data
	lpoints.resize(leng1);
	for(int i1=0;i1<leng1;i1++){
		lpoints[i1].resize(leng2);
		for(int i2=0;i2<leng2;i2++){
			lpoints[i1][i2].resize(leng3);
		}
	}
	
	string f3dname("clusters.data");
	cluster3input(f3dname);
	
	minx=-(double)Ns/2.0-0.5; maxx=(double)Ns/2.0+0.5;
      	miny=-(double)Ns/2.0-0.5; maxy=(double)Ns/2.0+0.5;
	minz=-(double)Ns/2.0-0.5; maxz=(double)Ns/2.0+0.5;
	
	red.resize(nclusters);
	green.resize(nclusters);
	blue.resize(nclusters);		
	
	for(int c=0;c<nclusters;c++)
		calcSphereColor(red[c], green[c], blue[c], c);
}

int main(int argc, char **argv){
	// Init GLUT and create window
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGBA);
	glutInitWindowPosition(100,100);
	glutInitWindowSize(640,640);
	mainWindow = glutCreateWindow("Tutorial");
	
	init();

	// Enter GLUT event processing cycle
	glutMainLoop();

	return 1;
}

