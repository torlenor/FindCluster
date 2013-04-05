/* findcluster.cpp
   Find cluster and performs calculations with it.
   Uses Polyakov loop eigenvalue's as input.

   v0.0.0 - 2013-03-05 Hans-Peter Schadler
*/

#include <iostream>
#include <fstream>
#include <complex>
#include <vector>
#include <cstdlib>
#include <string>

using namespace std;

#include "include/generic.hpp"

int Ns=0, Nt=0;
int matrixdim=3, leng1=0, leng2=0, leng3=0, leng4=0, Nspace=0;
string fevname;

// Vector for Polyakov loop eigenvalue data
vector<vector<complex<double> > > pollev;

vector<vector<int> > neib;

// Vectors to store cluster information
vector<int> isinsector;  // isinsector[is] stores the sector of the lattice point is
vector<int> clustersector; // clustersector[c] stores the sector of the cluster c
vector<int> isincluster;  // isincluster[is] stores the cluster of the lattice point is
vector<vector<int> > clustermembers; // clustermembers[c][i] stores the members/lattice points of the cluster c (0 < i < N_c)
vector<int> numberofcmembers; // numberofcmembers[c] stores the number of members of the cluster c

void fillneib();
int init(int &argc, char *argv[]);

void fillsectors(double fraction);

void findclusters();
void checkclusters();

int latmap(int i1, int i2, int i3);

int main(int argc, char *argv[]){
	double fraction;
	// Handle command line information and initialize arrays and other stuff...
	if(init(argc, argv) != 0){
		cout << "Error in init()!" << endl;
	}

	cout << "Reading 3d lattice with Polyakov loop ev's... " << flush;
	if(readPollEvBinary(leng1, leng2, leng3, leng4, matrixdim, pollev, fevname) != 0){
		cout << "ERROR: Problems with writePollEvBinary !" << endl;
	}
	cout << "done!" << endl;
	checkPollEv(leng1, leng2, leng3, leng4, matrixdim, pollev);
	
	fraction = 0.9;
	
	// Main part
	// Categorize lattice points by sectors
	fillsectors(fraction);
	// Identify clusters
	cout << "Finding clusters... " << flush;
	findclusters();
	cout << "done!" << endl;
	cout << "Checking clusters... " << flush;
	checkclusters();
	cout << "done!" << endl;

	return 0;
}

int init(int &argc, char *argv[]){
	if(argc<4){
		cout << "./findcluster.x Ns Nt binevconfin" << endl;
		return 1;
	}

	Ns=atoi(argv[1]);
	Nt=atoi(argv[2]);
	fevname = argv[3];
	
	leng1=Ns; leng2=Ns; leng3=Ns; leng4=Nt;
	
	Nspace = Ns*Ns*Ns;

	neib.resize(Nspace);
	for(int is=0;is<Nspace;is++)
		neib[is].resize(6);

	fillneib();

	cout << "Allocating 3d (Ns^3 * 3) lattice for Polyakov loop ev's... " << flush;
	pollev.resize(Nspace);
	for(int is=0;is<Nspace;is++)
		pollev[is].resize(matrixdim);
	cout << "done!" << endl;
	
	cout << "Creating cluster data arrays..." << flush;
	isinsector.resize(Nspace);
	// clustersector will not be allocated here, but on the fly with push_back
	isincluster.resize(Nspace);
	// clustermembers will not be allocated here, but on the fly with push_back
	cout << "done!" << endl;
	
	return 0;
}

void fillsectors(double fraction){
	int csectp1=0;
	int csectm1=0;
	int csect0=0;
	
	double tracephase=0;
	double delta0 = M_PI/3.0;
	double delta = delta0*fraction;
	
	#ifdef DEBUG
	cout << "Delta = " << delta << endl;
	#endif
	
	for(int is=0;is<Nspace;is++)
		isinsector[is] = 2;
	
	for(int is=0;is<Nspace;is++){
		tracephase = arg(pollev[is][0] + pollev[is][1] + pollev[is][2]);
		
		if(abs(tracephase - 2.0*M_PI/3.0) < delta){
			isinsector[is]=1;
			csectp1++;
		}
		
		if(abs(tracephase) < delta){
			isinsector[is]=0;
			csect0++;
		}
		
		if(abs(tracephase + 2.0*M_PI/3.0) < delta){
			isinsector[is]=-1;
			csectm1++;
		}
	}
	
	#ifdef DEBUG
	cout << "Sector 1 # : " << csectp1 << " , Sector 0 # : " << csect0 << " , Sector -1 # : " << csectm1 << endl;
	#endif
}

void findclusters(){
	int sector = 0, cluster = 0, memberis = 0, isneib = 0;

	vector<bool> isfiled;
	isfiled.resize(Nspace);
	for(int is=0;is<Nspace;is++){
		isfiled[is] = false;
		isincluster[is] = -1;
	}
	
	int totalmembers=0;
	
	cluster = 0;
	for(int is=0;is<Nspace;is++){
		if(isfiled[is] == false){
			sector = isinsector[is];
			cluster++; // current cluster
			clustersector.push_back(sector); // Current cluster is in sector sector
			clustermembers.push_back(vector<int>()); // Add a new cluster
			
			clustermembers[cluster-1].push_back(is); // Add is to this new cluster
			isincluster[is] = cluster; // Set the cluster id for the point  is
			isfiled[is] = true; // Lattice point is now filed
			
			// Iterate over all found points in the actual cluster (loop does not have fixed length)
			int member = 0;
			while((unsigned int)member < clustermembers[cluster-1].size()){
				memberis = clustermembers[cluster-1][member];
			
				for(int mu=0;mu<6;mu++){
					isneib = neib[memberis][mu];
					if(isfiled[isneib] == false && isinsector[isneib] == sector){
						clustermembers[cluster-1].push_back(isneib);
						isincluster[isneib] = cluster;
						isfiled[isneib] = true;
					}
				}
				
				#ifdef DEBUG
				cout << "Cluster " << cluster << " member " << member + 1 << " of " << clustermembers[cluster-1].size() << " members." << endl;
				#endif
				
				member++;
				totalmembers++;
			}
		}
	}
	
	#ifdef DEBUG
	cout << "" << totalmembers << " of " << Nspace << " lattice points visited!" << endl;
	#endif
	
	for(int is=0;is<Nspace;is++){
		if(isfiled[is] == false)
			cout << "WARNING: Some points did not get filled!" << endl;
	}
}

void checkclusters(){
	for(int is=0;is<Nspace;is++){
		if(isincluster[is] < 0)
			cout << "WARNING: Some points do not have a cluster ID!" << endl;
	}
	
	int sum=0;
	for(unsigned int c=0;c<clustermembers.size();c++)
		sum+=clustermembers[c].size();
	if(sum != Nspace){
		cout << "WARNING: Number of lattice points differ from number of clustermember entries!" << endl;
		cout << "Number of entries = " << sum << " Number of lattice points = " << Nspace << endl;
	}
	
}


int latmap(int i1, int i2, int i3){
	return i1 + i2*leng1 + i3*leng1*leng2;
}

void fillneib(){
	// Fills the neib array
        int i1p,i2p,i3p,i1m,i2m,i3m,is,isp1,isp2,isp3,ism1,ism2,ism3;
        for(int i1 = 0;i1<leng1;i1++){
                i1p = i1 + 1;
                i1m = i1 - 1;
                if (i1p == leng1) i1p = 0;
                if (i1m == -1) i1m = leng1-1;

                for(int i2 = 0;i2<leng2;i2++){
                        i2p = i2 + 1;
                        i2m = i2 - 1;
                        if (i2p == leng2) i2p = 0;
                        if (i2m == -1) i2m = leng2-1;

                        for(int i3 = 0;i3<leng3;i3++){
                                i3p = i3 + 1;
                                i3m = i3 - 1;
                                if (i3p == leng3) i3p = 0;
                                if (i3m == -1) i3m = leng3-1;
                                
                                // Compute the site address and the addresses of the sites shifted
                                // by one unit in each direction
                                is = i1 + i2*leng1 + i3*leng1*leng2;

                                isp1 = i1p + i2*leng1 + i3*leng1*leng2;
                                isp2 = i1 + i2p*leng1 + i3*leng1*leng2;
                                isp3 = i1 + i2*leng1 + i3p*leng1*leng2;

                                ism1 = i1m + i2*leng1 + i3*leng1*leng2;
                                ism2 = i1 + i2m*leng1 + i3*leng1*leng2;
                                ism3 = i1 + i2*leng1 + i3m*leng1*leng2;

                                // Fill the neib array
                                neib[is][0] = isp1;
                                neib[is][1] = isp2;
                                neib[is][2] = isp3;

                                neib[is][3] = ism1;
                                neib[is][4] = ism2;
                                neib[is][5] = ism3;
                        }
                }
	}
}


