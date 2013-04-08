/* findcluster.cpp
   Finds cluster and performs calculations with it.
   Uses Polyakov loop eigenvalues as input.

   v0.0.0 - 2013-03-05 Hans-Peter Schadler
*/

#include <iostream>
#include <fstream>
#include <complex>
#include <vector>
#include <cstdlib>
#include <string>

using namespace std;

#include "version.h"
#include "include/generic.hpp"

int Ns=4, Nt=4;
int matrixdim=3, leng1=4, leng2=4, leng3=4, leng4=4, Nspace=4*4*4*4;

int nmeas=1;

// Filename string vector
vector<string> fevname;

// Vector for Polyakov loop eigenvalue data
vector<vector<complex<double> > > pollev;

// Neib vector
vector<vector<int> > neib;

// Vectors to store cluster information
struct Clusterstruct{
	vector<int> isinsector;  // lclusterdata.isinsector[is] stores the sector of the lattice point is
	vector<int> clustersector; // lclusterdata.clustersector[c] stores the sector of the cluster c
	vector<int> isincluster;  // lclusterdata.isincluster[is] stores the cluster of the lattice point is
	vector<vector<int> > clustermembers; // lclusterdata.clustermembers[c][i] stores the members/lattice points of the cluster c (0 < i < N_c)

	vector<int> percolatingclusters; // Percolating clusters
	vector<int> percclusterdirection;
};

Clusterstruct *clusterdata;

struct Observablestruct{
	int maxclustersize;
	int maxclusterid;
	
	double avgclustersize;

	double cut;
};

Observablestruct *obs;

void fillNeib();
int init(int &argc, char *argv[]);

// Stuff to find and categorize sectors/clusters for one configuration
void fillSectors(Clusterstruct &lclusterdata, double delta);
void findClusters(Clusterstruct &lclusterdata);
void checkClusters(Clusterstruct &lclusterdata);
void findPercolatingCluster(Clusterstruct &lclusterdata);
void calcObservables(Observablestruct &lobs, Clusterstruct &lclusterdata);

void writeOneConfigResultsstdout(Observablestruct &lobs, Clusterstruct &lclusterdata);

void calcExp();

int latmap(int i1, int i2, int i3);

double fraction = 1.0;
double delta0 = M_PI/3.0;
double delta = delta0*fraction;

#include "findcluster_init.hpp"

int main(int argc, char *argv[]){
	// Handle command line information and initialize arrays and other stuff...
	if(init(argc, argv) != 0){
		cout << "Error in init()!" << endl;
		return 1;
	}

	cout << endl;
	
	// Main part
	// Loop over all nmeas configurations
	for(int n=0;n<nmeas;n++){
		cout << "------------------------------------------------------------------------------" << endl;
		cout << "Working on configuration " << fevname[n] << "..." << endl;
		if(readPollEvBinary(leng1, leng2, leng3, leng4, matrixdim, pollev, fevname[n]) != 0){
			cout << "ERROR: Problems with writePollEvBinary !" << endl;
			return 1;
		}
	
		checkPollEv(leng1, leng2, leng3, leng4, matrixdim, pollev);
	
		fillSectors(clusterdata[n], delta); // Categorize lattice points by sectors
		findClusters(clusterdata[n]);	// Identify clusters
		checkClusters(clusterdata[n]); // Check clusters
		findPercolatingCluster(clusterdata[n]); // Find percolating clusters
	
		calcObservables(obs[n], clusterdata[n]);
		// writeOneConfigResultsstdout(obs[n], clusterdata[n]);
	}

	calcExp();

	delete [] clusterdata; clusterdata=0;
	delete [] obs; obs=0;
	
	return 0;
}

void calcExp(){
	double maxclustersize=0, maxclustersizeerr=0;
	double avgclustersize=0, avgclustersizeerr=0;
	double cut=0, cuterr=0;
	// Plain mean value
	for(int n=0; n<nmeas; n++){
		maxclustersize += (&obs[n])->maxclustersize;
		avgclustersize += (&obs[n])->avgclustersize;
		cut += (&obs[n])->cut;
	}

	maxclustersize = maxclustersize/(double)nmeas;
	avgclustersize = avgclustersize/(double)nmeas;
	cut = cut/(double)nmeas;

	// Plain standard deviation
	for(int n=0;n<nmeas;n++){
		maxclustersizeerr=pow((&obs[n])->maxclustersize-maxclustersize,2);
		avgclustersizeerr=pow((&obs[n])->avgclustersize-avgclustersize,2);
		cuterr=pow((&obs[n])->cut-cut,2);
	}
	maxclustersizeerr=sqrt(maxclustersizeerr)/(double)nmeas;
	avgclustersizeerr=sqrt(avgclustersizeerr)/(double)nmeas;
	cuterr=sqrt(cuterr)/(double)nmeas;



	cout << "Expectation values: " << endl;
	cout << "Average cluster size = " << avgclustersize << ", Maximum cluster size = " << maxclustersize << endl;
	cout << "Average cluster err  = " << avgclustersizeerr << ", Maximum cluster err  = " << maxclustersizeerr << endl;
	cout << "Cut = " << cut << " Cut err = " << cuterr << endl;

}

void calcObservables(Observablestruct &lobs, Clusterstruct &lclusterdata){
	cout << "Calculating observables... " << flush;
	// Find largest cluster
	int size=-1, largestcluster=-1;
	for(unsigned int c=0; c<lclusterdata.clustermembers.size(); c++){
		if( (int)lclusterdata.clustermembers[c].size() > size){
			size = lclusterdata.clustermembers[c].size();
			largestcluster = c;
		}
	}
	
	lobs.maxclustersize = size;
	lobs.maxclusterid = largestcluster;
	
	// Calculate average cluster size
	double avgclustersize=0;
	for(unsigned int c=0; c<lclusterdata.clustermembers.size(); c++){
		avgclustersize += lclusterdata.clustermembers[c].size();
	}
	
	lobs.avgclustersize = avgclustersize/(double)lclusterdata.clustermembers.size();

	// Calculate cut
	double cut=0;
	for(unsigned int is=0; is<lclusterdata.isinsector.size(); is++){
		if(lclusterdata.isinsector[is] == 2)
			cut = cut + 1;
	}

	lobs.cut = cut/(double)Nspace;

	cout << "done!" << endl;
}

void writeOneConfigResultsstdout(Observablestruct &lobs, Clusterstruct &lclusterdata){
	cout << endl << lclusterdata.percolatingclusters.size() << " of " << lclusterdata.clustermembers.size() << " clusters are percolating!" << endl;
	for(unsigned int c=0;c<lclusterdata.percolatingclusters.size();c++){
		cout << "Cluster " << lclusterdata.percolatingclusters[c] << " is percolating!" << endl;
		cout << "The cluster has " << lclusterdata.clustermembers[lclusterdata.percolatingclusters[c]].size() << " members and is in sector " << lclusterdata.clustersector[lclusterdata.percolatingclusters[c]] << " !" << endl;
	}
	
	cout << endl;
	cout << "Average cluster size = " << lobs.avgclustersize << endl;
	cout << "Largest cluster is cluster " << lobs.maxclusterid << " with " << lobs.maxclustersize << " members." << endl << endl;
}

void fillSectors(Clusterstruct &lclusterdata, double delta){
	cout << "Categorizing lattice points by sector (delta = " << delta << " )... " << flush;

	int csectp1=0;
	int csectm1=0;
	int csect0=0;
	
	double tracephase=0;
	
	#ifdef DEBUG
	cout << "Delta = " << delta << endl;
	#endif
	
	for(int is=0;is<Nspace;is++)
		lclusterdata.isinsector[is] = 2;
	
	for(int is=0;is<Nspace;is++){
		tracephase = arg(pollev[is][0] + pollev[is][1] + pollev[is][2]);
		
		if(abs(tracephase - 2.0*M_PI/3.0) < delta){
			lclusterdata.isinsector[is]=1;
			csectp1++;
		}
		
		if(abs(tracephase) < delta){
			lclusterdata.isinsector[is]=0;
			csect0++;
		}
		
		if(abs(tracephase + 2.0*M_PI/3.0) < delta){
			lclusterdata.isinsector[is]=-1;
			csectm1++;
		}
	}
	
	#ifdef DEBUG
	cout << "Sector 1 # : " << csectp1 << " , Sector 0 # : " << csect0 << " , Sector -1 # : " << csectm1 << endl;
	#endif
	
	cout << "done!" << endl;
}

void findPercolatingCluster(Clusterstruct &lclusterdata){
	// Idea from Fortran code
	cout << "Searching for percolating clusters... " << flush;
	
	int is=0;
	
	int cluster=0;
	
	vector<vector<vector<int> > > members;
	members.resize(lclusterdata.clustermembers.size());
	for(unsigned int i=0;i<members.size();i++){
		members[i].resize(3);
		members[i][0].resize(leng1);
		for(int j=0;j<leng1;j++)
			members[i][0][j]=0;
		members[i][1].resize(leng2);
		for(int j=0;j<leng2;j++)
			members[i][1][j]=0;
		members[i][2].resize(leng3);
		for(int j=0;j<leng3;j++)
			members[i][2][j]=0;
	}
	
	for(int i1=0;i1<leng1;i1++)
	for(int i2=0;i2<leng2;i2++)
	for(int i3=0;i3<leng3;i3++){
		is = latmap(i1, i2, i3);
		cluster = lclusterdata.isincluster[is];
		if(lclusterdata.clustersector[cluster] < 2){
			members[cluster][0][i1] = 1;
			members[cluster][1][i2] = 1;
			members[cluster][2][i3] = 1;
		}
	}
	
	int sum1, sum2, sum3;
	for(unsigned int c=0;c<lclusterdata.clustermembers.size();c++){
		sum1=0; sum2=0; sum3=0;
		for(int i1=0;i1<leng1;i1++)
			sum1 += members[c][0][i1];
		for(int i2=0;i2<leng2;i2++)
			sum2 += members[c][1][i2];
		for(int i3=0;i3<leng3;i3++)
			sum3 += members[c][2][i3];
			
		if(sum1 == leng1 || sum2 == leng2 || sum3 == leng3){
			lclusterdata.percolatingclusters.push_back(c);
			#ifdef DEBUG
			cout << "sum1 = " << sum1 << " sum2 = " << sum2 << " sum3 = " << sum3 << endl;
			cout << "Cluster c = " << c << " added!" << endl;
			#endif
		}
	}
	
	cout << "done!" << endl;
}

void findClusters(Clusterstruct &lclusterdata){
	cout << "Finding clusters... " << flush;
	int sector = 0, cluster = 0, memberis = 0, isneib = 0;

	vector<bool> isfiled;
	isfiled.resize(Nspace);
	for(int is=0;is<Nspace;is++){
		isfiled[is] = false;
		lclusterdata.isincluster[is] = -1;
	}
	
	int totalmembers=0;
	
	cluster = -1; // Note: We begin counting with ID 0!
	for(int is=0;is<Nspace;is++){
		if(isfiled[is] == false){
			sector = lclusterdata.isinsector[is];
			cluster++; // current cluster id (Note: We begin counting with ID 0!)
			lclusterdata.clustersector.push_back(sector); // Current cluster is in sector sector
			lclusterdata.clustermembers.push_back(vector<int>()); // Add a new cluster
			
			lclusterdata.clustermembers[cluster].push_back(is); // Add is to this new cluster
			lclusterdata.isincluster[is] = cluster; // Set the cluster id for the point  is
			isfiled[is] = true; // Lattice point is now filed
			
			// Iterate over all found points in the actual cluster (loop does not have fixed length)
			int member = 0;
			while((unsigned int)member < lclusterdata.clustermembers[cluster].size()){
				memberis = lclusterdata.clustermembers[cluster][member];
			
				for(int mu=0;mu<6;mu++){
					isneib = neib[memberis][mu];
					if(isfiled[isneib] == false && lclusterdata.isinsector[isneib] == sector){
						lclusterdata.clustermembers[cluster].push_back(isneib);
						lclusterdata.isincluster[isneib] = cluster;
						isfiled[isneib] = true;
					}
				}
				
				#ifdef DEBUG
				cout << endl <<  "Cluster " << cluster << " member " << member + 1 << " of " << lclusterdata.clustermembers[cluster].size() << " members." << endl;
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
	
	cout << "done!" << endl;
}

void checkClusters(Clusterstruct &lclusterdata){
	cout << "Checking clusters... " << flush;

	for(int is=0;is<Nspace;is++){
		if(lclusterdata.isincluster[is] < 0)
			cout << "WARNING: Some points do not have a valid cluster ID!" << endl;
	}
	
	int sum=0;
	for(unsigned int c=0;c<lclusterdata.clustermembers.size();c++)
		sum+=lclusterdata.clustermembers[c].size();
	if(sum != Nspace){
		cout << "WARNING: Number of lattice points differ from number of clustermember entries!" << endl;
		cout << "Number of entries = " << sum << " Number of lattice points = " << Nspace << endl;
	}
	
	// Check for conistency of the cluster arrays lclusterdata.isinsector and lclusterdata.clustersector
	for(unsigned int c=0;c<lclusterdata.clustermembers.size();c++){
		for(unsigned int member=0; member<lclusterdata.clustermembers[c].size(); member++){
			if(lclusterdata.isinsector[lclusterdata.clustermembers[c][member]] != lclusterdata.clustersector[c]){
				cout << "WARNING: Problem with lclusterdata.clustermembers and lclusterdata.clustersector consistency!" << endl;
			}
		}
	}
	
	// Check if all members in clustermmembers have a correct ID
	// Here we use vector.at(iter) so that also a range check is performed and an
	// exeption is thrown if it is outside vector range
	int mult=1;
	vector<int> membercheck;
	membercheck.resize(Nspace);
	for(unsigned int c=0;c<lclusterdata.clustermembers.size();c++){
		for(unsigned int member=0; member<lclusterdata.clustermembers[c].size(); member++){
			membercheck.at(lclusterdata.clustermembers[c][member]) = 1;
		}
	}
	for(int is=0;is<Nspace;is++){
		mult = mult*membercheck[is];
	}
	if(mult != 1){
		cout << "WARNING: Problem with clustermember IDs!" << endl;
	}
	
	cout << "done!" << endl;
}

int latmap(int i1, int i2, int i3){
	return i1 + i2*leng1 + i3*leng1*leng2;
}

void fillNeib(){
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


