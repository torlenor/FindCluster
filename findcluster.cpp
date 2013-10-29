/* findcluster.cpp
   Finds cluster and performs calculations with it.
   Uses Polyakov loop eigenvalues as input.

   v0.3.0 - 2013-05-03 Hans-Peter Schadler
*/

#include <iostream>
#include <fstream>
#include <complex>
#include <vector>
#include <cstdlib>
#include <string>
#include <algorithm>
#include <iomanip>
#include <limits>

using namespace std;

#include "version.h"
#include "include/generic.hpp"

#include "include/jackknife.h"

#include "findcluster.h"

int Ns=4, Nt=4;
int matrixdim=3, leng1=4, leng2=4, leng3=4, leng4=4, Nspace=4*4*4*4;

int nmeas=1;

bool usealternativesectors=false;

bool do3d=false;

bool detail=false; // Controlls if we want detailed information for every configuration
bool doboxes=false; // Controlls if we want box counting calculations
bool dodistance=false; // Controlls if we want box counting calculations
bool doradius=true; // Controlls if we want radius calculation

bool memorysaver=false; // Controlls if we drop the clusterdata arrays after observable calculations

// Filename string vector
vector<string> fevname;

// Vector for Polyakov loop eigenvalue data
vector<vector<complex<double> > > pollev;

// Neib vector
vector<vector<int> > neib;

// Box size vector
vector<int> boxsize;
vector<int> boxes;

Clusterstruct *clusterdata;
Observablestruct *obs;
Resultstruct results;

double fraction = 0.0;
double delta0 = M_PI/3.0;
double delta = delta0*fraction;

double r = 0;

#include "findcluster_init.hpp"
#include "findcluster_helper.hpp"
#include "findcluster_cluster.hpp"
#include "findcluster_radius.hpp"
#include "findcluster_path.hpp"
#include "findcluster_obs.hpp"
#include "findcluster_write.hpp"

int main(int argc, char *argv[]){
	// Handle command line information and initialize arrays and other stuff...
	if(init(argc, argv) != 0){
		cout << "Error in init()!" << endl;
		return 1;
	}

	cout << endl;
	printsettings();
	cout << endl;

	// Main part
	// Loop over all nmeas configurations
	cout << "------------------------------------------------------------------------------" << endl;
	for(int n=0;n<nmeas;n++){
		// Allocate memory for faster access
                (&clusterdata[n])->isinsector.resize(Nspace);
                (&clusterdata[n])->isincluster.resize(Nspace);

		if(detail){
			cout << endl << "------------------------------------------------------------------------------" << endl;
			cout << fevname[n] << "..." << endl;
		}else{
			cout << "\r" <<  fevname[n] << "..." << flush;
		}
		// Read Polyakov loop eigenvalues file
		if(readPollEvBinary(leng1, leng2, leng3, leng4, matrixdim, pollev, fevname[n]) != 0){
			cout << "ERROR: Problems with writePollEvBinary !" << endl;
			return 1;
		}
		
		// Check Polyakov loop eigenvalues
		checkPollEv(leng1, leng2, leng3, leng4, matrixdim, pollev);
		
		if(usealternativesectors==true){
			fillSectorsAlt(clusterdata[n], r); // Categorize lattice points by sectors using alternative prescription
		}else{
			fillSectors(clusterdata[n], delta); // Categorize lattice points by sectors
		}
		
		findClusters(clusterdata[n]);	// Identify clusters
		checkClusters(clusterdata[n]); // Check clusters
		findPercolatingCluster(clusterdata[n]); // Find percolating clusters

		sortClusterSize(clusterdata[n]); // Sort clusters per number of members

		calcObservables(obs[n], clusterdata[n]); // Calculate observables
		
		if(detail){
			cout << endl << "Details for " << fevname[n] << ":" << endl;
			writeConfigResultsstdout(obs[n], clusterdata[n]);
		}
		
		if(do3d){
			// Write data for 3dclusters program
			stringstream f3dclustername;
			f3dclustername << "3dcluster_" << leng1 << "x" << leng4 << "_m" << setprecision(0) << fixed << n << ".data";
			cluster3doutput(clusterdata[n], f3dclustername.str());
		}

		if(memorysaver){
			freeMem(clusterdata[n]);
		}
	}

	if(do3d){
		// Write clusterfilenamelist for 3dcluster program
		ofstream f3dclusterlist;
		f3dclusterlist.open("3dcluster.list");
		for(int n=0;n<nmeas;n++){
			stringstream f3dclustername;
			f3dclustername << "3dcluster_" << leng1 << "x" << leng4 << "_m" << setprecision(0) << fixed << n << ".data";
			f3dclusterlist << f3dclustername.str() << endl;
		}
		f3dclusterlist.close();
	}

	cout << endl;
	cout << "------------------------------------------------------------------------------" << endl;

	// Calculate the expectation values and write the results to file and stdout
	calcExp();
	writeresultsstdout();
	writeresults();

	delete [] clusterdata; clusterdata=0;
	delete [] obs; obs=0;
	
	return 0;
}

