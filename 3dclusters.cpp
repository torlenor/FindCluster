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

bool usespheres=false;

#include "3dclusters.h"
#include "3dclusters_init.hpp"

int mainWindow;

// Config names and data vectors
vector<string> filenames;
vector<vector<vector<int > > > lpoints;
vector<int> isinsector;

int box=-1;

vector<int> boxsize;
vector<int> boxes;

// One color for every point (we store it two times)
vector<double> red, sred;
vector<double> green, sgreen;
vector<double> blue, sblue;

#include "3dclusters_keys.hpp"

void calcSphereColor(double &red, double &green, double &blue, int cluster);

void cluster3input(int config){
	ifstream f3d;
	int is, cleng1, cleng2, cleng3, cleng4;
	f3d.open(filenames[config].c_str());
	if(f3d.is_open()){
		f3d >> cleng1 >> cleng2 >> cleng3 >> cleng4;
		f3d >> nclusters;
		
		if(cleng1 != leng1 || cleng2 != leng2 || cleng3 != leng3)
			cout << "WARNING: Lattice size missmatch!" << endl;
		
		int ii1, ii2, ii3, cluster, sector, is;
		for(int i1=0;i1<leng1;i1++)
		for(int i2=0;i2<leng2;i2++)
		for(int i3=0;i3<leng3;i3++){
			f3d >> ii1 >> ii2 >> ii3 >> cluster >> sector;
			lpoints.at(ii1).at(ii2).at(ii3) = cluster;
			is = ii1 + ii2*leng1 + ii3*leng1*leng2;
			isinsector.at(is) = sector;
		}
		f3d.close();
	}else{
		cout << "WARNING: Could not open 3dcluster file!" << endl;
	}
	
	red.resize(leng1*leng2*leng3);
	green.resize(leng1*leng2*leng3);
	blue.resize(leng1*leng2*leng3);	
	
	for(int i1=0;i1<leng1;i1++)
	for(int i2=0;i2<leng2;i2++)
	for(int i3=0;i3<leng3;i3++){
		int is = i1 + i2*leng1 + i3*leng1*leng2;
		calcSphereColor(red[is], green[is], blue[is], lpoints.at(i1).at(i2).at(i3));
	}
		
	sred=red;
	sgreen=green;
	sblue=blue;
	
	if(cnt>-1 && cnt<nclusters){
		if(onecluster==true){
			setOnlyCluster(cnt);		
		}else{
			setHightlightCluster(cnt);
		}
	}else{
		cnt=-1;	
	}
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

void pressKey(int key, int x, int y){
        switch (key) {
                case GLUT_KEY_UP : deltaMove = 5.0f; break;
                case GLUT_KEY_DOWN : deltaMove = -5.0f; break;
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
}

void drawSquare(){
	glColor3f(1.0, 1.0, 1.0);
	glutWireCube(maxx-minx);
}

void drawHalfBox(float x, float y, float z, float size, float red, float green, float blue){
	float ux, uy, uz; // Ursprung
	float ax, ay, az;
	float bx, by, bz;
	float cx, cy, cz;

	ux = x - size/2.0;
	uy = y - size/2.0;
	uz = z - size/2.0;
	
	ax = ux + size;
	ay = uy;
	az = uz;
	
	bx = ux;
	by = uy + size;
	bz = uz;
	
	cx = ux;
	cy = uy;
	cz = uz + size;

	glColor3f(red, green, blue);
	glBegin(GL_LINES);
		glVertex3i( ux, uy, uz);
		glVertex3i( ax, ay, az);

		glVertex3i( ux, uy, uz);
		glVertex3i( bx, by, bz);

		glVertex3i( ux, uy, uz);
		glVertex3i( cx, cy, cz);
	glEnd(); 
}

void drawBoxes(){
	double m1, m2, m3;
	bool cinbox=false;
	for(int box1=0; box1<boxes[box]; box1++){
	m1=box1*boxsize[box]+boxsize[box]/2.0-0.5;
	for(int box2=0; box2<boxes[box]; box2++){
	m2=box2*boxsize[box]+boxsize[box]/2.0-0.5;
	for(int box3=0; box3<boxes[box]; box3++){
		m3=box3*boxsize[box]+boxsize[box]/2.0-0.5;
		glPushMatrix();
			glTranslatef(m1-(double)Ns/2.0+0.5,m2-(double)Ns/2.0+0.5,m3-(double)Ns/2.0+0.5);
			cinbox=false;
			for(int i1=box1*boxsize[box];i1<(box1+1)*boxsize[box];i1++)
			for(int i2=box2*boxsize[box];i2<(box2+1)*boxsize[box];i2++)
			for(int i3=box3*boxsize[box];i3<(box3+1)*boxsize[box];i3++){
				if(lpoints.at(i1).at(i2).at(i3)==cnt){
					cinbox=true;
					break;
				}
			}
			//if(lpoints.at(m1).at(m2).at(m3)==cnt){
			if(cinbox){
			 	glColor3f(1.0, 1.0, 1.0);
			 	glutWireCube(boxsize[box]);
			}
		glPopMatrix();
	}
	}
	}
	
	/*
	float i1, i2, i3, ii1, ii2, ii3, m1, m2, m3;
	for(int box1=0; box1<boxes[box]; box1++){
	ii1=box1*boxsize[box]+boxsize[box]/2.0-0.5;
	for(int box2=0; box2<boxes[box]; box2++){
	ii2=box2*boxsize[box]+boxsize[box]/2.0-0.5;
	for(int box3=0; box3<boxes[box]; box3++){
	ii3=box3*boxsize[box]+boxsize[box]/2.0-0.5;
		i1=ii1-(double)Ns/2.0+0.5;
		i2=ii2-(double)Ns/2.0+0.5;
		i3=ii3-(double)Ns/2.0+0.5;
		drawHalfBox((float)i1, (float)i2, (float)i3, boxsize[box], 1.0, 1.0, 1.0);
	}
	}
	} */
	
}

void drawSphere(double x, double y, double z){
	glPushMatrix();
		glTranslatef(x-(double)Ns/2.0+0.5,y-(double)Ns/2.0+0.5,z-(double)Ns/2.0+0.5);
		// glRotatef(90.0,0.0,1.0,0.0);
		glutSolidSphere(sphereradius,sphereSlices,sphereStacks);
		// glutSolidCube(sphereradius);
	glPopMatrix();
}

void drawLattice() {
	double x=0, crandx;
	double y=0, crandy;
	double z=0, crandz;
	
	// Depth buffer modification for solid objects
	glDepthMask(GL_TRUE);
	glDisable( GL_BLEND );
	glDisable(GL_LIGHT0);
	glDisable(GL_LIGHTING);
	
	// Draw wireframe box around the scene
	drawSquare();
	if(showboxes==true){
		drawBoxes();
	}

	glEnable( GL_BLEND );
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);

	glPointSize(pointsize);

	int i1=0, i2=0, i3=0, is;
	if(usespheres==false)
		glBegin(GL_POINTS);
       	for(int ri1=0;ri1<leng1;ri1++)
       	for(int ri2=0;ri2<leng2;ri2++)
       	for(int ri3=0;ri3<leng3;ri3++){
       		// Render the points from back to front, regardless of camera angle
       		// TODO: Has to be implemented for rotations around x axis, too!
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
       		
       		is = i1 + i2*leng1 + i3*leng1*leng2;
       		if(red[is]>0 || green[is]>0 || blue[is]>0){
       		if(isinsector[is]<2){
       			if(red[is]==1 && green[is]==1 && blue[is]==1){
       				// The white points are always solid
       				glDepthMask(GL_TRUE);
				glColor4f(red[is], green[is], blue[is], 1);
			}else{
				if(alpha<1.0){
					// Deactivate depth buffer modification for transparent objects
					glDepthMask(GL_FALSE);
				}
				glColor4f(red[is], green[is], blue[is], alpha);
			}
			if(usespheres==true){
				drawSphere(i1,i2,i3);
			}else{
				glVertex3f(i1-(double)Ns/2.0+0.5, i2-(double)Ns/2.0+0.5, i3-(double)Ns/2.0+0.5);
			}
		}
		}
	}
	if(usespheres==false)
		glEnd(); // GL_POINTS
}

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

	// Compute movement and rotation of the camera
        if(deltaMove || deltaAnglex || deltaAngley) {
	                computePos(deltaMove);
	                glutSetWindow(mainWindow);
	                glutPostRedisplay();
        }

	// Rotate the camera
	glRotatef(anglex, 1.0f, 0.0f, 0.0f);
	glRotatef(angley, 0.0f, 1.0f, 0.0f);

	// Render function
	drawLattice();

	// FPS stuff
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
	// Keeps the aspect ratio right if one changes the window size
	if(h == 0)
		h = 1;

	float ratio = 1.0*w/h;

	// Projection matrix
	glMatrixMode(GL_PROJECTION);

	// Reset matrix
	glLoadIdentity();

	// Set viewpoint to be the entire window
	glViewport(0, 0, w, h);

	// Set the correct perspective
	// gluPerspective(field of view angle in yz plane, 
	// 			ratio, near, far clipping planes)
	gluPerspective(60,ratio,1,200);

	// Get back to the Modelview matrix
	glMatrixMode(GL_MODELVIEW);
}

int getFilelist(string f3dlistname){
	// Read f3dlistname file and fill filenames
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

int openglInit(){
	// OpenGL stuff
	glEnable(GL_DEPTH_TEST);
	glDepthMask(GL_TRUE);
	glEnable( GL_POINT_SMOOTH );
	glEnable( GL_BLEND );
	glBlendFunc( GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA );
	if(usespheres==true){
		glEnable(GL_LIGHTING);
		glEnable(GL_LIGHT0);
		glColorMaterial ( GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE );
		glEnable ( GL_COLOR_MATERIAL );
	}
// 	glEnable(GL_LIGHT1);
	
	glEnable(GL_SMOOTH);
	glShadeModel(GL_SMOOTH);

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
}

int init(){;
	// Vectors to hold the data
	lpoints.resize(leng1);
	for(int i1=0;i1<leng1;i1++){
		lpoints[i1].resize(leng2);
		for(int i2=0;i2<leng2;i2++){
			lpoints[i1][i2].resize(leng3);
		}
	}
	isinsector.resize(leng1*leng2*leng3);
	
	// Get the file list
	selconfig=0;
	filenames.resize(nconfig);
	getFilelist(f3dlistname);
	cluster3input(selconfig);
	
	minx=-(double)Ns/2.0-0.5; maxx=(double)Ns/2.0+0.5;
      	miny=-(double)Ns/2.0-0.5; maxy=(double)Ns/2.0+0.5;
	minz=-(double)Ns/2.0-0.5; maxz=(double)Ns/2.0+0.5;

	// Boxes
	for(int s=0;s<Ns;s++){
		if( Ns % (s+1) == 0)
			boxsize.push_back(s+1);
	}   

	for(unsigned int i=0;i<boxsize.size();i++){
		boxes.push_back(Ns/boxsize[i]);
	}   
}

int main(int argc, char **argv){
	// Init GLUT and create window
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGBA);
	glutInitWindowPosition(100,100);
	glutInitWindowSize(640,640);
	mainWindow = glutCreateWindow("3d clusters");
	
	if(parameterInit(argc, argv)!=0){
		cout << "Error in parameterInit()!" << endl;
		return 1;
	}
	
	openglInit();
	init();

	// Enter GLUT event processing cycle
	glutMainLoop();

	return 1;
}

