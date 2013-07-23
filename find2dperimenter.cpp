#include <iostream>
#include <vector>

using namespace std;

int Ns=8;
int dim=2;

void buildcluster(vector<vector<int> > &l){
	// lattice[x][y]=1;
	l[0][0]=0;                       l[3][0]=0; l[4][0]=0; l[5][0]=0; l[6][0]=0; l[7][0]=0;
	l[0][1]=0;                       l[3][1]=1; l[4][1]=1; l[5][1]=1;            l[7][1]=0;
	l[0][2]=0; l[1][2]=1; l[2][2]=1; l[3][2]=1; l[4][2]=1; l[5][2]=1;            l[7][2]=0;
	l[0][3]=0; l[1][3]=1; l[2][3]=1;                       l[5][3]=1; l[6][3]=1; l[7][3]=0;
	           l[1][4]=1; l[2][4]=1;                       l[5][4]=1; l[6][4]=1;
	l[0][5]=0; l[1][5]=1; l[2][5]=1; l[3][5]=1; l[4][5]=1; l[5][5]=1; l[6][5]=1; l[7][5]=0;
	                      l[2][6]=1; l[3][6]=1; l[4][6]=1; l[5][6]=1; l[6][6]=1;
	                                 l[3][7]=0; l[4][7]=0; l[5][7]=1;
}

void printcluster(vector<vector<int> > &lattice){
	cout << endl << "Lattice:" << endl << endl;
	for(int i=0;i<Ns;i++){
		for(int j=0;j<Ns;j++){
			if(lattice[j][i]==1)
				cout << "1" << " ";
			else
				cout << " " << " ";
		}
		cout << endl;
	}
	cout << endl;
}

void printclusterperimeter(vector<vector<int> > &lattice, vector<vector<int> > &isperimeter){
	cout << endl << "Lattice\t\t\tPerimeter:" << endl << endl;
	for(int i=0;i<Ns;i++){
		for(int j=0;j<Ns;j++){
			if(lattice[j][i]==1)
				cout << "1" << " ";
			else
				cout << " " << " ";
		}
		cout << "\t";
		for(int j=0;j<Ns;j++){
			if(isperimeter[j][i]==1)
				cout << "*" << " ";
			else
				cout << " " << " ";
		}
		cout << endl;
	}
	cout << endl;
}

void findperimeterbruteforce(vector<vector<int > > &isperimeter, vector<vector<int> > &lattice){
	for(int i=0;i<Ns;i++){
		for(int j=0;j<Ns;j++){
			isperimeter[j][i]=0;
			if(lattice[j][i]==1){
				if(lattice[j+1][i]==0 || lattice[j-1][i]==0 || lattice[j][i+1]==0 || lattice[j][i-1]==0)
					isperimeter[j][i]=1;
			}
		}
	}
}

void findperimeter(vector<vector<int > > &isperimeter, vector<vector<int> > &lattice, int startx, int starty){
	/* --------------------------------- LEFT ----------------------------------------- */
	cout << "Finding perimeter (going left from (" << startx << "," << starty << ")... " << endl;

	// Biased random walk
	// d is direction: d=(-x,-y,x.y) and goes from 0 to 3
	int d=1; // we start with direction upwards

	// Start at a point which belongs to the cluster
	int xp, xpt;
	int yp, ypt;

	// Needed to get the directions right
	d=3;

	xp=startx;
	yp=starty;

	do{
		// We try to go to the nearest neighbors starting in direction left and trying it clockwise 
		// as long as we find a point which belongs to the cluster.
		do{
			d=d+1;
			if(d>3)	d=d-4;

			// Do the try step
			if(d==3){ xpt=xp; ypt=yp+1; }
			if(d==2){ xpt=xp+1; ypt=yp; }
			if(d==1){ xpt=xp; ypt=yp-1; }
			if(d==0){ xpt=xp-1; ypt=yp; }

			if(xpt < 0)    xpt = xpt + Ns;
			if(xpt > Ns-1) xpt = xpt - Ns;
			if(ypt < 0)    ypt = ypt + Ns;
			if(ypt > Ns-1) ypt = ypt - Ns;

		}while(lattice[xpt][ypt] != 1);

		xp=xpt;
		yp=ypt;
		isperimeter[xpt][ypt]=1;

		d = d - 2;
		if(d<0)	d=d+4;
	}while(xp != startx || yp != starty); // Perform the walk until we hit (startx,starty) again
	
	
	/* --------------------------------- RIGHT ---------------------------------------- */
	cout << "Finding perimeter (going right from (" << startx << "," << starty << ")... " << endl;

	// Biased random walk
	// d is direction: d=(x,y,-x.-y) and goes from 0 to 3
	d=1; // we start with direction upwards

	// Start at a point which belongs to the cluster
	xp=0; xpt=0;
	yp=0; ypt=0;

	// Needed to get the directions right
	d=3;

	xp=startx;
	yp=starty;

	do{
		// We try to go to the nearest neighbors starting in direction right and trying it counter-clockwise 
		// as long as we find a point which belongs to the cluster.
		do{
			d=d+1;
			if(d>3)	d=d-4;

			// Do the try step
			if(d==3){ xpt=xp; ypt=yp-1; }
			if(d==2){ xpt=xp-1; ypt=yp; }
			if(d==1){ xpt=xp; ypt=yp+1; }
			if(d==0){ xpt=xp+1; ypt=yp; }

			if(xpt < 0)    xpt = xpt + Ns;
			if(xpt > Ns-1) xpt = xpt - Ns;
			if(ypt < 0)    ypt = ypt + Ns;
			if(ypt > Ns-1) ypt = ypt - Ns;

		}while(lattice[xpt][ypt] != 1);

		xp=xpt;
		yp=ypt;
		isperimeter[xpt][ypt]=1;

		d = d - 2;
		if(d<0)	d=d+4;
	}while(xp != startx || yp != starty); // Perform the walk until we hit (startx,starty) again
}

int main(){
	vector<vector<int > > lattice;
	lattice.resize(Ns);
	for(int i=0;i<Ns;i++)
		lattice[i].resize(Ns);
	for(int i=0;i<Ns;i++)
		for(int j=0;j<Ns;j++)
			lattice[i][j]=0;
	
	vector<vector<int > > isperimeter;
	isperimeter.resize(Ns);
	for(int i=0;i<Ns;i++)
		isperimeter[i].resize(Ns);

	
	buildcluster(lattice);

	// printcluster(lattice);

	// First, search for a starting point
	int startx=0;
	int starty=0;

	for(int y=0;y<Ns;y++)
	for(int x=0;x<Ns;x++){
		if(lattice[x][y]==1){
			startx=x;
			starty=y;
		}
	}
	findperimeter(isperimeter, lattice, startx, starty);
	
	// Do it again but start at a different corner
	for(int y=Ns-1;y>=0;y--)
	for(int x=Ns-1;x>=0;x--){
		if(lattice[x][y]==1){
			startx=x;
			starty=y;
		}
	}
	findperimeter(isperimeter, lattice, startx, starty);

	printclusterperimeter(lattice, isperimeter);
	
	cout << endl << "---------------------------------------------------------------------------" << endl
		<< " Brute force method (finds also internal perimeters " << endl << endl;
	findperimeterbruteforce(isperimeter, lattice);

	printclusterperimeter(lattice, isperimeter);

	cout << "Have a nice day!" << endl;

}
