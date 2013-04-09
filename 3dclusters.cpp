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

#include "3dclusters_init.hpp"

int Ns=40, Nt=8;

// angle for rotation for the camerca direction
float anglex = 0.0f;float angley = 0.0f;

// actual vector representing the camera's direction
float lx=0.0f,lz=-1.0f, ly = 0.0f;

// XZ position of the camera
float x=0.0f, z=(double)Ns+20, y = 0.00f;

unsigned frameCount = 0;

double vtemp=0;

float deltaAnglex = 0.0f;float deltaAngley = 0.0f;
float deltaMove = 0;
int xOrigin = -1; int yOrigin = -1; 

const double FREQ=60; // Hz
const double TIMERMSECS=1000*1/(double)FREQ;

const double dt = TIMERMSECS/1000;;

const int dim=3;

// Stuff for the spheres
double sphereradius=0.35; // Sphere radius
const int sphereSlices=8; // Sphere slices around Z axis
const int sphereStacks=8; //Sphere stacks/slices along the z axis

double pointsize=10; // Sphere radius

double alpha=1;


double minx=-Ns - 0.5, maxx=Ns + 0.5,
      	miny=-Ns - 0.5,maxy=Ns + 0.5,
	minz=-Ns - 0.5,maxz=Ns + 0.5;

int mainWindow;

int nconfig=0, selconfig=0;
vector<string> filenames;

vector<vector<vector<int > > > lpoints;

vector<vector<vector<int > > > pointsdisabled;

vector<double> red, sred;
vector<double> green, sgreen;
vector<double> blue, sblue;

int leng1=Ns, leng2=Ns, leng3=Ns, leng4=Nt, nclusters=0;

int cnt=-1;

void calcSphereColor(double &red, double &green, double &blue, int cluster);

void cluster3input(int config){
	ifstream f3d;
	int is, cleng1, cleng2, cleng3, cleng4;
	f3d.open(filenames[config].c_str());
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
	
	red.resize(nclusters);
	green.resize(nclusters);
	blue.resize(nclusters);	
	
	for(int c=0;c<nclusters;c++)
		calcSphereColor(red[c], green[c], blue[c], c);
		
	sred=red;
	sgreen=green;
	sblue=blue;
	
	cnt=-1;
}

void mouseMove(int x, int y) {
        // this will only be true when the left button is down
        if (xOrigin >= 0 || yOrigin >= 0) {

                // update deltaAngle
                deltaAngley = (x - xOrigin) * 0.01f;
		deltaAnglex = (y - yOrigin) * 0.01f;                

                // update camera's direction

                glutSetWindow(mainWindow);
                glutPostRedisplay();
        }   
}

void mouseButton(int button, int state, int x, int y) {
        // only start motion if the left button is pressed
        if (button == GLUT_LEFT_BUTTON) {
                // when the button is released
                if (state == GLUT_UP) {
                        anglex += deltaAnglex;
                        angley += deltaAngley;
                        deltaAnglex = 0.0f;
                        deltaAngley = 0.0f;
                        xOrigin = -1;
                        yOrigin = -1;
                }
                else  {// state = GLUT_DOWN
                        xOrigin = x;
                        yOrigin = y;
                }
        }
}


void processNormalKeys(unsigned char key, int x, int y){
	if(key == 27){
		exit(0);
	}else if(key == '+'){
		int mod = glutGetModifiers();
		if (mod == GLUT_ACTIVE_ALT){
			 alpha += 0.1;
			 if(alpha > 1.0)
			 	alpha=1.0;
		}
		else{
			sphereradius += 0.1;
			pointsize += 1.0;
			if(sphereradius > 1.0)
				sphereradius = 1.0;
		}
		cout << "Radius = " << pointsize << " Alpha = "  << alpha << endl;
	}else if(key == '-'){
		int mod = glutGetModifiers();
		if (mod == GLUT_ACTIVE_ALT){
			 alpha -= 0.1;
			 if(alpha < 0.0)
			 	alpha=0;
		}
		else{
			sphereradius -= 0.1;
			pointsize -= 1.0;
			if(pointsize < 0.0)
				pointsize = 1.0;
			if(sphereradius < 0.0)
				sphereradius = 0.0;
		}
		cout << "Radius = " << pointsize << " Alpha = "  << alpha << endl;
	}else if(key == 'c' || key == 'C' ){
		int mod = glutGetModifiers();
		if (mod == GLUT_ACTIVE_SHIFT){
			cnt--;
			
			if(cnt < 0)
				cnt = nclusters-1;
			red=sred;
			green=sgreen;
			blue=sblue;
			red[cnt]=1; green[cnt]=1; blue[cnt]=1;
			cout << "Cluster " << cnt << " selected!" << endl;
		}else{
			cnt++;
			
			if(cnt > nclusters-1)
				cnt = 0;
			red=sred;
			green=sgreen;
			blue=sblue;
			red[cnt]=1; green[cnt]=1; blue[cnt]=1;
			cout << "Cluster " << cnt << " selected!" << endl;
		}
	}else if(key == 's' || key == 'S' ){
		int mod = glutGetModifiers();
		if (mod == GLUT_ACTIVE_SHIFT){
			cnt--;
			
			if(cnt < 0)
				cnt = nclusters-1;
			for(int i=0;i<nclusters;i++){
				red[i]=0;
				green[i]=0;
				blue[i]=0;
			}
			red[cnt]=1; green[cnt]=1; blue[cnt]=1;
			cout << "Cluster " << cnt << " selected!" << endl;
		}else{
			cnt++;
			
			if(cnt > nclusters-1)
				cnt = 0;
			for(int i=0;i<nclusters;i++){
				red[i]=0;
				green[i]=0;
				blue[i]=0;
			}
			red[cnt]=1; green[cnt]=1; blue[cnt]=1;
			cout << "Cluster " << cnt << " selected!" << endl;
		}	
	}else if(key == 'r'){
		cnt=-1;
		red=sred;
		green=sgreen;
		blue=sblue;
		cout << "No Cluster selected!" << endl;
	}else if(key == 'n' || key == 'N' ){
		int mod = glutGetModifiers();
		if (mod == GLUT_ACTIVE_SHIFT){
			selconfig--;
			if(selconfig==-1)
				selconfig=nconfig-1;
			cout << "Loading configuration " << selconfig << " ..." << endl;
			cluster3input(selconfig);
		}else{
			selconfig++;
			if(selconfig==nconfig)
				selconfig=0;
			cout << "Loading configuration " << selconfig << " ..." << endl;
			cluster3input(selconfig);
		}
	}	
}

void pressKey(int key, int x, int y){
        switch (key) {
                case GLUT_KEY_UP : deltaMove = 1.5f; break;
                case GLUT_KEY_DOWN : deltaMove = -1.5f; break;
                case GLUT_KEY_LEFT : deltaAngley = 1.5f; break;
                case GLUT_KEY_RIGHT : deltaAngley = -1.5f; break;                
        }   
        glutSetWindow(mainWindow);
        glutPostRedisplay();

}

void releaseKey(int key, int x, int y) {
	switch (key) {
		case GLUT_KEY_UP :
		case GLUT_KEY_DOWN : deltaMove = 0;break;
		case GLUT_KEY_LEFT :
		case GLUT_KEY_RIGHT : deltaAngley = 0;break;
	}
}

void computePos(float deltaMove) {
	        x += deltaMove * lx * 0.1f;
	        z += deltaMove * lz * 0.1f;
	        anglex += deltaAnglex;
	        angley += deltaAngley;
	        if(anglex>360) anglex-=360;
	        if(angley>360) angley-=360;
	        if(anglex<0) anglex+=360;
	        if(angley<0) angley+=360;
	        
	        cout << "anglex = " << anglex << " angley = " << angley << endl;
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
}

void drawSphere(double x, double y, double z){
	glPushMatrix();
		glTranslatef(x-(double)Ns/2.0+0.5,y-(double)Ns/2.0+0.5,z-(double)Ns/2.0+0.5);
		// glRotatef(90.0,0.0,1.0,0.0);
		glutSolidSphere(sphereradius,sphereSlices,sphereStacks);
		// glutSolidCube(sphereradius);
	glPopMatrix();
}

void drawBouncingPoint() {
	double x=0, crandx;
	double y=0, crandy;
	double z=0, crandz;
	
	glDepthMask(GL_TRUE);
	glEnable( GL_BLEND );
	
	drawSquare();

	//glPointSize(sphereradius);
	//glEnable(GL_POINT_SMOOTH);
	
	glPointSize(pointsize);

	if(alpha<1.0){
		glDepthMask(GL_FALSE);
	//	glDisable( GL_BLEND );
	}
	int i1=0, i2=0, i3=0;
	glBegin(GL_POINTS);
       	for(int ri1=0;ri1<leng1;ri1++)
       	for(int ri2=0;ri2<leng2;ri2++)
       	for(int ri3=0;ri3<leng3;ri3++){
       		if( angley >= 0 && angley <= 45){
       			i1=ri3; i2=ri2; i3=ri1;
       		}else if(angley >= 45 && angley <= 135){
       			i1=-ri1+leng1-1; i2=ri2; i3=-ri3+leng3-1;
       		}else if(angley >= 135 && angley <= 225){
       			i1=-ri3+leng3-1; i2=ri2; i3=-ri1+leng1-1;
       		}else if(angley >= 225 && angley <= 315){
       			i1=ri1; i2=ri2; i3=ri3;   
      		}else if(angley >= 315 && angley <= 360){
      				i1=ri3; i2=ri2; i3=ri1;
      		}else{
       			i1=ri1; i2=ri2; i3=ri3;
       		}
       		int is = i1 + i2*leng1 + i3*leng1*leng2;
       		if(red[lpoints[i1][i2][i3]]>0 || green[lpoints[i1][i2][i3]]>0 || blue[lpoints[i1][i2][i3]]>0){
       			if(red[lpoints[i1][i2][i3]]==1 && green[lpoints[i1][i2][i3]]==1 && blue[lpoints[i1][i2][i3]]==1){
       				glDepthMask(GL_TRUE);
				//glEnable( GL_BLEND );
				glColor4f(red[lpoints[i1][i2][i3]], green[lpoints[i1][i2][i3]], blue[lpoints[i1][i2][i3]], 1);
			}else{
				if(alpha<1.0){
					glDepthMask(GL_FALSE);
				//	glDisable( GL_BLEND );
				}
				glColor4f(red[lpoints[i1][i2][i3]], green[lpoints[i1][i2][i3]], blue[lpoints[i1][i2][i3]], alpha);
			}
			// drawSphere(i1,i2,i3);
			glVertex3f(i1-(double)Ns/2.0+0.5, i2-(double)Ns/2.0+0.5, i3-(double)Ns/2.0+0.5);
		}
	}
	glEnd();
	
	glDepthMask(GL_TRUE);
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

        if(deltaMove || deltaAnglex || deltaAngley) {
	                computePos(deltaMove);
	                glutSetWindow(mainWindow);
	                glutPostRedisplay();
        }

	glRotatef(anglex, 1.0f, 0.0f, 0.0f);
	glRotatef(angley, 0.0f, 1.0f, 0.0f);

	drawBouncingPoint();

	// angle+=0.4f;

	usedTime += currentTime-pTime;
	pTime=currentTime;

	if(usedTime > 1000){
		fps=frameCount / (usedTime/(double)1000 );
		char* TempString = (char*)
		    malloc(512);
	 
	 	if(cnt>-1){
		sprintf(
		    TempString,
		    "%f Frames Per Second, Configuration %i selected, Cluster %i selected",
		    fps, selconfig, cnt );
		}else{
		sprintf(
		    TempString,
		    "%f Frames Per Second, Configuration %i selected",
		    fps, selconfig );		
		}
	 
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

int getFilelist(string f3dlistname){
	// Read finname file and fill fevname
	ifstream fin;
	fin.open(f3dlistname.c_str());
	if(fin.is_open()!=true){
        	cout  << "ERROR: File " << f3dlistname <<  " to read configuration filename list could not be opened!" << endl;
		throw 1;
	}
	
	string strtmp;
	int n=0;
	while(n<nconfig && getline(fin, filenames.at(n)) ){
		n++;
	}
	if(n<nconfig){
		cout << "ERROR: Only found " << n << " names in " << f3dlistname << " !" << endl;
		return 1;
	}
}

int init(){;
	// OpenGL stuff
	glEnable(GL_DEPTH_TEST);
	glDepthMask(GL_TRUE);
	glEnable( GL_POINT_SMOOTH );
	glEnable( GL_BLEND );
	glBlendFunc( GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA );
//	glEnable(GL_LIGHTING);
//	glEnable(GL_LIGHT0);
// 	glEnable(GL_LIGHT1);
//	glColorMaterial ( GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE );
//	glEnable ( GL_COLOR_MATERIAL );
	
//	glEnable(GL_SMOOTH);
//	glShadeModel(GL_SMOOTH);

	// glutIgnoreKeyRepeat(1);
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
	
	
	nconfig=10; selconfig=0;
	filenames.resize(nconfig);
	string f3dlistname("cluster.list");
	getFilelist(f3dlistname);
	cluster3input(selconfig);
	
	minx=-(double)Ns/2.0-0.5; maxx=(double)Ns/2.0+0.5;
      	miny=-(double)Ns/2.0-0.5; maxy=(double)Ns/2.0+0.5;
	minz=-(double)Ns/2.0-0.5; maxz=(double)Ns/2.0+0.5;
}

int main(int argc, char **argv){
	// Init GLUT and create window
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGBA);
	glutInitWindowPosition(100,100);
	glutInitWindowSize(640,640);
	mainWindow = glutCreateWindow("3d clusters");
	
	parameterInit(argc, argv);
	init();

	// Enter GLUT event processing cycle
	glutMainLoop();

	return 1;
}

