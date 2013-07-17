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

double fraction = 0.0;
double delta0 = M_PI/3.0;
double delta = delta0*fraction;

double r = 0;

#include "findcluster_init.hpp"
#include "findcluster_radius.hpp"
#include "findcluster_obs.hpp"

void printsettings(){
	cout << "Settings:" << endl << endl;
	cout << "Lattice size = " << leng1 << "x" << leng2 << "x" << leng3 << "x" << leng4 << endl;
	cout << "Number of configurations = " << nmeas << endl;
	if(! usealternativesectors){
		cout << "Cut fraction = " << fraction << endl;
	}else{
		cout << "Alternative cut prescription, radius r = " << r << endl;
	}
	if(doboxes)
		cout << "Calculating 'box' observables." << endl;
	if(detail)
		cout << "Writing detailed results for every configuration." <<  endl;
	if(do3d)
		cout << "Writing 3dcluster visualization data files." <<  endl;
	if(memorysaver)
		cout << "Memory saver active." <<  endl;
}

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

			cout << endl << "The following clusters are connected over the PBCs: " << endl;
			for(unsigned int c=0;c<clusterdata[n].clusterisperiodic.size();c++){
				if(clusterdata[n].clusterisperiodic[clusterdata[n].sortedcluster[c]][0]==1 
				|| clusterdata[n].clusterisperiodic[clusterdata[n].sortedcluster[c]][1]==1 
				|| clusterdata[n].clusterisperiodic[clusterdata[n].sortedcluster[c]][2]==1){
					cout << "Cluster (size sorted) " << c << ", real id = " << clusterdata[n].sortedcluster[c] << endl;
				}
			}
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

	delete [] clusterdata; clusterdata=0;
	delete [] obs; obs=0;
	
	return 0;
}

void freeMem(Clusterstruct &lclusterdata){
	lclusterdata.isinsector.resize(0);
	lclusterdata.clustersector.resize(0);
	lclusterdata.isincluster.resize(0);

	for(unsigned c=0;c<lclusterdata.clustermembers.size();c++){
		lclusterdata.clustermembers[c].resize(0);
	}
	lclusterdata.clustermembers.resize(0);

	lclusterdata.percolatingclusters.resize(0);
	for(unsigned p=0; p<lclusterdata.percolatingdirections.size(); p++){
		lclusterdata.percolatingdirections[p].resize(0);
	}
	lclusterdata.percolatingdirections.resize(0);
    
	lclusterdata.sortedcluster.resize(0);
	lclusterdata.sortedrealcluster.resize(0);
	lclusterdata.isinsortedcluster.resize(0);
}

void getCoords(int is, int &i1, int &i2, int &i3){
	i1 = (is % (leng1*leng2) ) % leng1;
	i2 = (is % (leng1*leng2) ) / leng1;
	i3 = is / (leng1*leng2);

	if(is != i1 + i2*leng1 + i3*leng1*leng2)
		cout << "ERROR: Problem in getCoords!" << endl;
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

void writeClusterList(Clusterstruct &lclusterdata){
	vector<int> csize;
	for(unsigned c=0; c<lclusterdata.clustermembers.size(); c++){
		csize.push_back(lclusterdata.clustermembers[c].size());
	}

	sort(csize.begin(), csize.end());

	cout << "Cluster sizes" << endl;
	for(unsigned c=0; c<lclusterdata.clustermembers.size(); c++){
		cout << csize[c] << endl;
	}
}

void cluster3doutput(Clusterstruct &lclusterdata, string f3dname){
	ofstream f3d;
	int is;
	f3d.open(f3dname.c_str());
	if(f3d.is_open()){
		f3d << leng1 << " " << leng2 << " " << leng3 << " " << leng4 << endl;
		f3d << lclusterdata.clustermembers.size() << endl;
		for(int i1=0;i1<leng1;i1++)
		for(int i2=0;i2<leng2;i2++)
		for(int i3=0;i3<leng3;i3++){
			is = i1 + i2*leng1 + i3*leng1*leng2;

			f3d << i1 << " " << i2 << " " << i3 << " " << lclusterdata.isinsortedcluster[is] << " " << lclusterdata.isinsector[is] << endl;
		}
		f3d.close();
	}else{
		cout << "WARNING: Could not open 3dcluster file!" << endl;
	}
}

void calcExp(){
	double maxclustersize=0, maxclustersizeerr=0;

	double mlargestnpclustersize=0, mlargestnpclustersizeerr=0;

	double avgclustersize=0, avgclustersizeerr=0;
	double avgclustersizeF=0, avgclustersizeFerr=0;
	double avgclustersizenopercc=0, avgclustersizenoperccerr=0;
	double avgrootmeansquaredistanceR=0, avgrootmeansquaredistanceRerr=0;

	double avgpercc=0, avgperccerr=0;

	double cut=0, cuterr=0;

	double mlaserdim=0, mlaserdimerr=0;
	
	double ddata[nmeas];

	// Start of Jackknife
	// Maximal cluster size per volume expectation value
	for(int n=0;n<nmeas;n++){
		ddata[n]=0;
		for(int j=0;j<nmeas;j++){
			if(n!=j)
				ddata[n] += (&obs[j])->maxclustersize/(double)Nspace;
		}
		ddata[n] = ddata[n]/(double)(nmeas-1);
	}
	Jackknife(ddata, maxclustersize, maxclustersizeerr, nmeas);
	
	// Maximal non percolating cluster size per volume expectation value
	for(int n=0;n<nmeas;n++){
		ddata[n]=0;
		for(int j=0;j<nmeas;j++){
			if(n!=j)
				ddata[n] += (&obs[j])->largestnonpercclustersize/(double)Nspace;
		}
		ddata[n] = ddata[n]/(double)(nmeas-1);
	}
	Jackknife(ddata, mlargestnpclustersize, mlargestnpclustersizeerr, nmeas);
	
	// Average cluster size expectation value
	for(int n=0;n<nmeas;n++){
		ddata[n]=0;
		for(int j=0;j<nmeas;j++){
			if(n!=j)
				ddata[n] += (&obs[j])->avgclustersize;
		}
		ddata[n] = ddata[n]/(double)(nmeas-1);
	}
	Jackknife(ddata, avgclustersize, avgclustersizeerr, nmeas);
	
	// Average cluster size expectation value (Fortunato (1.7))
	for(int n=0;n<nmeas;n++){
		ddata[n]=0;
		for(int j=0;j<nmeas;j++){
			if(n!=j)
				ddata[n] += (&obs[j])->avgclustersizeF;
		}
		ddata[n] = ddata[n]/(double)(nmeas-1);
	}
	Jackknife(ddata, avgclustersizeF, avgclustersizeFerr, nmeas);
	
	// Average cluster size expectation value without percolating clusters
	for(int n=0;n<nmeas;n++){
		ddata[n]=0;
		for(int j=0;j<nmeas;j++){
			if(n!=j)
				ddata[n] += (&obs[j])->avgclustersizenopercc;
		}
		ddata[n] = ddata[n]/(double)(nmeas-1);
	}
	Jackknife(ddata, avgclustersizenopercc, avgclustersizenoperccerr, nmeas);
	
	// Average root mean square distance traveled R
	for(int n=0;n<nmeas;n++){
		ddata[n]=0;
		for(int j=0;j<nmeas;j++){
			if(n!=j)
				ddata[n] += (&obs[j])->rootmeansquaredistanceR;
		}
		ddata[n] = ddata[n]/(double)(nmeas-1);
	}
	Jackknife(ddata, avgrootmeansquaredistanceR, avgrootmeansquaredistanceRerr, nmeas);
	
	// Cut expectation value
	for(int n=0;n<nmeas;n++){
		ddata[n]=0;
		for(int j=0;j<nmeas;j++){
			if(n!=j)
				ddata[n] += (&obs[j])->cut;
		}
		ddata[n] = ddata[n]/(double)(nmeas-1);
	}
	Jackknife(ddata, cut, cuterr, nmeas);
	
	// Area expectation value
	for(int n=0;n<nmeas;n++){
		ddata[n]=0;
		for(int j=0;j<nmeas;j++){
			if(n!=j)
				ddata[n] += (&obs[j])->area;
		}
		ddata[n] = ddata[n]/(double)(nmeas-1);
	}
	Jackknife(ddata, mlaserdim, mlaserdimerr, nmeas);
	
	// Area largest non percolating cluster expectation value
	double marealargestnonpercc=0, marealargestnonperccerr=0;
	for(int n=0;n<nmeas;n++){
		ddata[n]=0;
		for(int j=0;j<nmeas;j++){
			if(n!=j)
				ddata[n] += (&obs[j])->arealargestnonperccluster;
		}
		ddata[n] = ddata[n]/(double)(nmeas-1);
	}
	Jackknife(ddata, marealargestnonpercc, marealargestnonperccerr, nmeas);
	
	// Radius of largest cluster expectation value
	double mlargestclusterradius=0, mlargestclusterradiuserr=0;
	for(int n=0;n<nmeas;n++){
		ddata[n]=0;
		for(int j=0;j<nmeas;j++){
			if(n!=j)
				ddata[n] += (&obs[j])->largestclusterradius;
		}
		ddata[n] = ddata[n]/(double)(nmeas-1);
	}
	Jackknife(ddata, mlargestclusterradius, mlargestclusterradiuserr, nmeas);

	// Polyakov loop expectation value only for points with sector < 2
	double mpoll=0, mpollerr=0;
	for(int n=0;n<nmeas;n++){
		ddata[n]=0;
		for(int j=0;j<nmeas;j++){
			if(n!=j)
				ddata[n] += (&obs[j])->poll;
		}
		ddata[n] = ddata[n]/(double)(nmeas-1);
	}
	Jackknife(ddata, mpoll, mpollerr, nmeas);

	cout << "Expectation values (single eliminitation jackknife): " << endl;
	cout << "Average cluster size = " << setprecision(14) << avgclustersize << ", Maximum cluster size / V = " << maxclustersize << endl;
	cout << "Average cluster err  = " << avgclustersizeerr << ", Maximum cluster / V err  = " << maxclustersizeerr << endl;
	cout << "Average cluster size Fortunato (1.7) = " << avgclustersizeF << ", Error  = " << avgclustersizeFerr << endl;
	cout << "Radius of largest cluster = " << mlargestclusterradius << ", Error = " << mlargestclusterradiuserr << endl;
	cout << "Root mean distance traveled R = " << avgrootmeansquaredistanceR << ", Error = " << avgrootmeansquaredistanceRerr << endl;
	cout << "Cut = " << cut << " Cut err = " << cuterr << endl;
	cout << "Laserdim = " << mlaserdim << " Laserdim err = " << mlaserdimerr << endl;
	cout << "Polyakov loop = " << mpoll << " Polyakov loop err = " << mpollerr << endl;

	cout << endl;

	stringstream fclustersizename;
	fclustersizename << "clustersize_" << Ns << "x" << Nt << "_f" << fraction << ".res";
	ofstream fclustersize;
	fclustersize.open(fclustersizename.str().c_str());
	fclustersize << "# Nt largestclusterweight largestclusterweighterr avgclusterweight avgclusterweighterr avgfortunatoclustersize avgfortunatoclustersizeerr largestnonpercclusterweight largestnonpercclusterweighterr largestclusterradius largestclusterradiuserr rootmeansquaredistanceR rootmeansquaredistanceR" << endl;
	fclustersize.flags (std::ios::scientific);
	fclustersize.precision(numeric_limits<double>::digits10 + 1);
	fclustersize << Nt << " " << maxclustersize << " " << maxclustersizeerr << " " << avgclustersize << " " << avgclustersizeerr << " " << avgclustersizeF << " " << avgclustersizeFerr << " " << mlargestnpclustersize << " " << mlargestnpclustersizeerr << " " << mlargestclusterradius << " " << mlargestclusterradiuserr << " " << avgrootmeansquaredistanceR << " " << avgrootmeansquaredistanceRerr << endl;
	fclustersize.close();
	
	stringstream favgclustersizename;
	favgclustersizename << "clustersize_avg_" << Ns << "x" << Nt << "_f" << fraction << ".res";
	ofstream favgclustersize;
	favgclustersize.open(favgclustersizename.str().c_str());
	favgclustersize << "# Nt avgclusterweight avgclusterweighterr avgfortunatoclustersize avgfortunatoclustersizeerr avgclusterweightnopercc avgclusterweightnoperccerr" << endl;
	favgclustersize.flags (std::ios::scientific);
	favgclustersize.precision(numeric_limits<double>::digits10 + 1);
	favgclustersize << Nt << " " << avgclustersize << " " << avgclustersizeerr << " " << avgclustersizeF << " " << avgclustersizeFerr << " " << avgclustersizenopercc << " " << avgclustersizenoperccerr << endl;
	favgclustersize.close();

	// Number of percolating clusters expectation value
	for(int n=0;n<nmeas;n++){
		ddata[n]=0;
		for(int j=0;j<nmeas;j++){
			if(n!=j)
				ddata[n] += (&obs[j])->percc;
		}
		ddata[n] = ddata[n]/(double)(nmeas-1);
	}
	Jackknife(ddata, avgpercc, avgperccerr, nmeas);

	stringstream fnperccname;
	fnperccname << "npercc_" << Ns << "x" << Nt << "_f" << fraction << ".res";
	ofstream fnpercc;
	fnpercc.open(fnperccname.str().c_str());
	fnpercc << "# Nt percclusters percclusterserr" << endl;
	fnpercc.flags (std::ios::scientific);
	fnpercc.precision(numeric_limits<double>::digits10 + 1);
	fnpercc << Nt << " " << avgpercc << " " << avgperccerr << endl;
	fnpercc.close();

	if(doboxes){
		// Box counts for largest cluster which is not in sector 2
		vector<double> avgboxcnt, avgboxcnterr;
		avgboxcnt.resize(boxsize.size());
		avgboxcnterr.resize(boxsize.size());
		
		for(unsigned int size=0; size<boxsize.size(); size++){
			for(int n=0;n<nmeas;n++){
				ddata[n]=0;

				for(int j=0;j<nmeas;j++){
					if(n!=j)
					ddata[n] += (&obs[j])->numberofboxes[(&obs[j])->largestclusterid][size];
				}
				ddata[n] = ddata[n]/(double)(nmeas-1);
			}
			Jackknife(ddata, avgboxcnt[size], avgboxcnterr[size], nmeas);
		}

		stringstream fboxcntname;
		fboxcntname << "boxcnt_" << Ns << "x" << Nt << "_f" << fraction << ".res";
		ofstream fboxcnt;
		fboxcnt.open(fboxcntname.str().c_str());
		fboxcnt << "# boxsize boxcnt boxcnterr" << endl;
		fboxcnt.flags (std::ios::scientific);
		fboxcnt.precision(numeric_limits<double>::digits10 + 1);
		for(unsigned int size=0; size<boxsize.size(); size++){
			fboxcnt << boxsize[size] << " " << avgboxcnt[size] << " " << avgboxcnterr[size] << endl;
		}
		fboxcnt.close();
		
		// Box counts for largest cluster which is not percolating
		for(unsigned int size=0; size<boxsize.size(); size++){
			avgboxcnt[size]=0;
			avgboxcnterr[size]=0;
			for(int n=0;n<nmeas;n++){
				ddata[n]=0;

				for(int j=0;j<nmeas;j++){
					if(n!=j)
					ddata[n] += (&obs[j])->numberofboxes[(&obs[j])->largestnonpercclusterid][size];
				}
				ddata[n] = ddata[n]/(double)(nmeas-1);
			}
			Jackknife(ddata, avgboxcnt[size], avgboxcnterr[size], nmeas);
		}

		fboxcntname.str("");
		fboxcntname << "boxcnt_nonpercc_" << Ns << "x" << Nt << "_f" << fraction << ".res";
		fboxcnt.open(fboxcntname.str().c_str());
		fboxcnt << "# boxsize boxcntnonpercc boxcntnonperccerr" << endl;
		fboxcnt.flags (std::ios::scientific);
		fboxcnt.precision(numeric_limits<double>::digits10 + 1);
		for(unsigned int size=0; size<boxsize.size(); size++){
			fboxcnt << boxsize[size] << " " << avgboxcnt[size] << " " << avgboxcnterr[size] << endl;
		}
		fboxcnt.close();
	}

	// Write surface area to file
	stringstream fareaname;
	fareaname << "area_" << Ns << "x" << Nt << "_f" << fraction << ".res";
	ofstream farea;
	farea.open(fareaname.str().c_str());
	farea << "# Nt area areaerr arealargestnonpercc arealargestnonperccerr largestclusterweight largestclusterweighterr largestnonpercclusterweight largestnonpercclusterweighterr" << endl;
	farea.flags (std::ios::scientific);
	farea.precision(numeric_limits<double>::digits10 + 1);
	farea << Nt << " " << mlaserdim << " " << mlaserdimerr << " " << marealargestnonpercc << " " << marealargestnonperccerr << " " << maxclustersize << " " << maxclustersizeerr << " " << mlargestnpclustersize << " " << mlargestnpclustersizeerr << endl;
	farea.close();

	
	// Write Polyakov loop to file
	stringstream fpollname;
	fpollname << "poll_" << Ns << "x" << Nt << "_f" << fraction << ".res";
	ofstream fpoll;
	fpoll.open(fpollname.str().c_str());
	fpoll << "# Nt poll pollerr" << endl;
	fpoll.flags (std::ios::scientific);
	fpoll.precision(numeric_limits<double>::digits10 + 1);
	fpoll << Nt << " " << mpoll << " " << mpollerr << endl;
	fpoll.close();
}

void writeConfigResultsstdout(Observablestruct &lobs, Clusterstruct &lclusterdata){
	cout << "Number of points in the different sectors:" << endl;
	cout << "Sector -1 # = " << lclusterdata.nsectm1 << " Sector 0 # = " << lclusterdata.nsect0 << " Sector 1 # = " << lclusterdata.nsectp1 << ", Sectors cut = " << Nspace - lclusterdata.nsectm1 - lclusterdata.nsect0 - lclusterdata.nsectp1 << endl << endl;

	cout  << lclusterdata.percolatingclusters.size() << " of " << lclusterdata.clustermembers.size() << " clusters are percolating!" << endl;
	for(unsigned int c=0;c<lclusterdata.percolatingclusters.size();c++){
		cout << "Cluster " << lclusterdata.percolatingclusters[c] << " is percolating in directions ("
			<< lclusterdata.percolatingdirections[c][0] << "," 
			<< lclusterdata.percolatingdirections[c][1] << ","
			<< lclusterdata.percolatingdirections[c][2] << ")!" << endl;
		cout << "The cluster has " << lclusterdata.clustermembers[lclusterdata.percolatingclusters[c]].size() << " members and is in sector " << lclusterdata.clustersector[lclusterdata.percolatingclusters[c]] << " !" << endl;
	}
	
	cout << endl;
	cout << "Average cluster size = " << lobs.avgclustersize << endl;
	cout << "Average cluster size Fortunato = " << lobs.avgclustersizeF << endl;
	cout << "Largest cluster is cluster " << lobs.maxclusterid << " with " << lobs.maxclustersize << " members." << endl;

	if(doboxes){
		cout << endl << "Boxes (Clusters sorted based on # of members):" << endl;
		for(unsigned int c=0; c<lclusterdata.sortedcluster.size(); c++){
			cout << "Cluster = " << lclusterdata.sortedcluster[c] << " with " << lclusterdata.clustermembers[lclusterdata.sortedcluster[c]].size() << " members:" << endl;
			for(unsigned int size=0; size<boxsize.size(); size++){
			cout << lobs.numberofboxes[lclusterdata.sortedcluster[c]][size] << " boxes of size " << boxsize[size] << " needed." << endl;
			}
		}
	}

//	cout << endl << "Cluster radius (Fortunato):" << endl;
//	for(unsigned int c=0; c<lclusterdata.sortedcluster.size(); c++){
//		if( lclusterdata.clustersector[lclusterdata.sortedcluster[c]] < 2){
//			cout << "Cluster ( " << lclusterdata.clustermembers[lclusterdata.sortedcluster[c]].size() << " members ) = " << lclusterdata.sortedcluster[c] << " radius = " << setprecision(5) << lobs.clusterradius[lclusterdata.sortedcluster[c]] << setprecision(2) << ", Center of Mass = (" << lobs.centerofmass[lclusterdata.sortedcluster[c]][0] << "," << lobs.centerofmass[lclusterdata.sortedcluster[c]][1] << "," << lobs.centerofmass[lclusterdata.sortedcluster[c]][2] << ")"  << endl;
//		}
//	}
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
