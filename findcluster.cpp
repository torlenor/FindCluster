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
#include "findcluster_radius.hpp"
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

	
	calcExp();
	writeresultsstdout();
	writeresults();

	delete [] clusterdata; clusterdata=0;
	delete [] obs; obs=0;
	
	return 0;
}

void sortClusterSize(Clusterstruct &lclusterdata){
	lclusterdata.sortedcluster.resize(lclusterdata.clustermembers.size());
	lclusterdata.isinsortedcluster.resize(Nspace);
	
	for(unsigned int c=0; c<lclusterdata.clustermembers.size(); c++){
		lclusterdata.sortedcluster[c]=c;
	}
	
	int tmp1;
	bool notsorted=true;
	while(notsorted){
		notsorted=false;
		for(unsigned int c=0;c<lclusterdata.sortedcluster.size()-1;c++){
			if(lclusterdata.clustermembers[lclusterdata.sortedcluster[c+1]].size()<lclusterdata.clustermembers[lclusterdata.sortedcluster[c]].size()){
				tmp1=lclusterdata.sortedcluster[c+1];
				lclusterdata.sortedcluster[c+1]=lclusterdata.sortedcluster[c];
				lclusterdata.sortedcluster[c]=tmp1;
				notsorted=true;
			}
		}
	}
	
	reverse(lclusterdata.sortedcluster.begin(),lclusterdata.sortedcluster.end());
	
	for(unsigned c=0; c<lclusterdata.sortedcluster.size(); c++){
		for(unsigned member=0; member < lclusterdata.clustermembers[lclusterdata.sortedcluster[c]].size(); member++){
			lclusterdata.isinsortedcluster[lclusterdata.clustermembers[lclusterdata.sortedcluster[c]][member]] = c;
		}
	}

	for(unsigned c=0; c<lclusterdata.sortedcluster.size(); c++){
		if(lclusterdata.isinsector[lclusterdata.clustermembers[lclusterdata.sortedcluster[c]][0]] < 2)
			lclusterdata.sortedrealcluster.push_back(lclusterdata.sortedcluster[c]);
	}
}

void fillSectorsAlt(Clusterstruct &lclusterdata, double r){
	#ifdef DEBUG
	cout << "Categorizing lattice points by sector ( radius = " << r << " )... " << flush;
	#endif

	int csectp1=0;
	int csectm1=0;
	int csect0=0;

	delta = M_PI/3.0;
	
	double tracephase=0;
	double radius=0;
	
	#ifdef DEBUG
	cout << "Delta = " << delta << endl;
	cout << "r = " << r << endl;
	#endif

	lclusterdata.poll.resize(Nspace);
	
	for(int is=0;is<Nspace;is++)
		lclusterdata.isinsector[is] = 2;
	
	for(int is=0;is<Nspace;is++){
		tracephase = arg(pollev[is][0] + pollev[is][1] + pollev[is][2]);
		radius = abs(pollev[is][0] + pollev[is][1] + pollev[is][2]);
		
		if(abs(tracephase - 2.0*M_PI/3.0) <= delta){
			if(radius >= r){
				lclusterdata.isinsector[is]=1;
				csectp1++;
			}
		}
		
		if(abs(tracephase) < delta){
			if(radius >= r){
				lclusterdata.isinsector[is]=0;
				csect0++;
			}
		}
		
		if(abs(tracephase + 2.0*M_PI/3.0) <= delta){
			if(radius >= r){
				lclusterdata.isinsector[is]=-1;
				csectm1++;
			}
		}
		lclusterdata.poll[is]=pollev[is][0] + pollev[is][1] + pollev[is][2];
	}
	
	lclusterdata.nsectm1=csectm1; lclusterdata.nsect0=csect0; lclusterdata.nsectp1=csectp1;
	#ifdef DEBUG
	cout << "done!" << endl;
	#endif
}

void fillSectors(Clusterstruct &lclusterdata, double delta){
	#ifdef DEBUG
	cout << "Categorizing lattice points by sector (delta = " << delta << " )... " << flush;
	#endif

	int csectp1=0;
	int csectm1=0;
	int csect0=0;
	
	double tracephase=0;
	
	#ifdef DEBUG
	cout << "Delta = " << delta << endl;
	#endif

	lclusterdata.poll.resize(Nspace);
	
	for(int is=0;is<Nspace;is++)
		lclusterdata.isinsector[is] = 2;
	
	for(int is=0;is<Nspace;is++){
		tracephase = arg(pollev[is][0] + pollev[is][1] + pollev[is][2]);
		
		if(abs(tracephase - 2.0*M_PI/3.0) <= delta){
			lclusterdata.isinsector[is]=1;
			csectp1++;
		}
		
		if(abs(tracephase) < delta){
			lclusterdata.isinsector[is]=0;
			csect0++;
		}
		
		if(abs(tracephase + 2.0*M_PI/3.0) <= delta){
			lclusterdata.isinsector[is]=-1;
			csectm1++;
		}

		lclusterdata.poll[is]=pollev[is][0] + pollev[is][1] + pollev[is][2];
	}
	
	lclusterdata.nsectm1=csectm1; lclusterdata.nsect0=csect0; lclusterdata.nsectp1=csectp1;
	#ifdef DEBUG
	cout << "done!" << endl;
	#endif
}

void findPercolatingCluster(Clusterstruct &lclusterdata){
	#ifdef DBEUG
	cout << "Searching for percolating clusters... " << flush;
	#endif
	
	int is=0;
	
	int cluster=0;
	
	lclusterdata.clusterispercolating.resize(lclusterdata.clustermembers.size());
	for(unsigned int c=0;c<lclusterdata.clustermembers.size();c++){
		lclusterdata.clusterispercolating[c] = 0;
	}
	
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
	
	int sum1, sum2, sum3, percclustercnt=-1;
	for(unsigned int c=0;c<lclusterdata.clustermembers.size();c++){
		sum1=0; sum2=0; sum3=0;
		for(int i1=0;i1<leng1;i1++)
			sum1 += members[c][0][i1];
		for(int i2=0;i2<leng2;i2++)
			sum2 += members[c][1][i2];
		for(int i3=0;i3<leng3;i3++)
			sum3 += members[c][2][i3];
			
		if(sum1 == leng1 || sum2 == leng2 || sum3 == leng3){
			percclustercnt++;
			lclusterdata.percolatingclusters.push_back(c);
			lclusterdata.percolatingdirections.push_back(vector<int>());
			lclusterdata.percolatingdirections[percclustercnt].resize(3);
			lclusterdata.percolatingdirections[percclustercnt][0]=0;
			lclusterdata.percolatingdirections[percclustercnt][1]=0;
			lclusterdata.percolatingdirections[percclustercnt][2]=0;
			if(sum1 == leng1)
				lclusterdata.percolatingdirections[percclustercnt][0]=1;
			if(sum2 == leng2)
				lclusterdata.percolatingdirections[percclustercnt][1]=1;
			if(sum3 == leng3)
				lclusterdata.percolatingdirections[percclustercnt][2]=1;
			#ifdef DEBUG
			cout << "sum1 = " << sum1 << " sum2 = " << sum2 << " sum3 = " << sum3 << endl;
			cout << "Cluster c = " << c << " added!" << endl;
			cout << "Percolating in direction (" << lclusterdata.percolatingdirections[percclustercnt][0] << ","
				<< lclusterdata.percolatingdirections[percclustercnt][1] << "," 
				<< lclusterdata.percolatingdirections[percclustercnt][2] << ")." << endl;
			#endif
		
			if(lclusterdata.percolatingdirections[percclustercnt][0] == 1 && 
			  lclusterdata.percolatingdirections[percclustercnt][1] == 1 &&
			  lclusterdata.percolatingdirections[percclustercnt][2] == 1){
				lclusterdata.clusterispercolating[c]=1;
			}
		}
	}
	#ifdef DBEUG	
	cout << "done!" << endl;
	#endif
}

void findClusters(Clusterstruct &lclusterdata){
	#ifdef DEBUG
	cout << "Finding clusters... " << flush;
	#endif
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
			
			// Modification due to cluster radius with PBCs
			lclusterdata.clusterisperiodic.push_back(vector<int>()); // Add a new cluster
			lclusterdata.clusterisperiodic[cluster].resize(3);
			for(int i=0;i<3;i++)
				lclusterdata.clusterisperiodic[cluster][i]=0;
			int ci1=0, ci2=0, ci3=0;
		
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

						// Check if over periodic boundaries
						getCoords(memberis, ci1, ci2, ci3);
						switch(mu){
							case 0: if(ci1 + 1 == leng1)
									lclusterdata.clusterisperiodic[cluster][0]=1;
								break;

							case 1: if(ci2 + 1 == leng2)
									lclusterdata.clusterisperiodic[cluster][1]=1;
								break;
							
							case 2: if(ci3 + 1 == leng3)
									lclusterdata.clusterisperiodic[cluster][2]=1;
								break;
							
							case 3: if(ci1 - 1 == -1)
									lclusterdata.clusterisperiodic[cluster][0]=1;
								break;
							
							case 4: if(ci2 - 1 == -1)
									lclusterdata.clusterisperiodic[cluster][1]=1;
								break;
							
							case 5: if(ci3 - 1 == -1)
									lclusterdata.clusterisperiodic[cluster][2]=1;
								break;

							default: cout << "ERROR: Error in cluster switch!" << endl;
						}
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
	#ifdef DEBUG	
	cout << "done!" << endl;
	#endif
}

void checkClusters(Clusterstruct &lclusterdata){
	#ifdef DEBUG	
	cout << "Checking clusters... " << flush;
	#endif

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
	
	#ifdef DEBUG	
	cout << "done!" << endl;
	#endif
}

