/*
 * findcluster_cluster.cpp - Cluster identification
 *
 * Copyright © 2014 H.-P. Schadler  <hanspeter.schadler@uni-graz.at>
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 */

#include "findcluster_cluster.h"

#include <algorithm>
#include <complex>
#include <iostream>
#include <vector>

#include "findcluster.h"

using std::abs;
using std::arg;

void sortClusterSize(Clusterstruct &lclusterdata, Options opt) {
	// Sorts clusters by their size and creates a sortedcluster vector
	// which contains the sorted cluster id's
	lclusterdata.sortedcluster.resize(lclusterdata.clustermembers.size());
	lclusterdata.isinsortedcluster.resize(opt.Nspace);
	
	for (unsigned int c=0; c<lclusterdata.clustermembers.size(); c++) {
		lclusterdata.sortedcluster[c]=c;
	}
	
	int tmp1;
	bool notsorted=true;
	while (notsorted) {
		notsorted=false;
		for (unsigned int c=0; c<lclusterdata.sortedcluster.size()-1; c++) {
			if (lclusterdata.clustermembers[lclusterdata.sortedcluster[c+1]].size()<lclusterdata.clustermembers[lclusterdata.sortedcluster[c]].size()) {
				tmp1=lclusterdata.sortedcluster[c+1];
				lclusterdata.sortedcluster[c+1]=lclusterdata.sortedcluster[c];
				lclusterdata.sortedcluster[c]=tmp1;
				notsorted=true;
			}
		}
	}
	
	reverse(lclusterdata.sortedcluster.begin(), lclusterdata.sortedcluster.end());
	
	for (unsigned c=0; c<lclusterdata.sortedcluster.size(); c++) {
		for (unsigned member=0; member < lclusterdata.clustermembers[lclusterdata.sortedcluster[c]].size(); member++) {
			lclusterdata.isinsortedcluster[lclusterdata.clustermembers[lclusterdata.sortedcluster[c]][member]] = c;
		}
	}

	for (unsigned c=0; c<lclusterdata.sortedcluster.size(); c++) {
		if (lclusterdata.isinsector[lclusterdata.clustermembers[lclusterdata.sortedcluster[c]][0]] < 2)
			lclusterdata.sortedrealcluster.push_back(lclusterdata.sortedcluster[c]);
	}
}

void fillSectorsAlt(Clusterstruct &lclusterdata, std::vector<std::vector<std::complex<double> > > &pollev, Options opt, double r) {
	// Categorizes the lattice points by sector using the alternative
	// prescription with the radius
	#ifdef DEBUG
		std::cout << "Categorizing lattice points by sector ( radius = " << r << " )... " << std::flush;
	#endif

	int csectp1=0;
	int csectm1=0;
	int csect0=0;

	opt.delta = M_PI/3.0;
	
  std::complex<double> trace=0;
	double tracephase=0;
	double radius=0;
	
	#ifdef DEBUG
		std::cout << "Delta = " << opt.delta << std::endl;
		std::cout << "r = " << r << std::endl;
	#endif

	lclusterdata.poll.resize(opt.Nspace);
	
	for (int is=0; is<opt.Nspace; is++)
		lclusterdata.isinsector[is] = 2;
	
	for (int is=0; is<opt.Nspace; is++) {
    trace=0;
    for(int i=0;i<opt.matrixdim;i++) {
      trace+=pollev.at(is).at(i);
    }
		tracephase = arg(trace);
		radius = abs(trace);
		
		if (abs(tracephase - 2.0*M_PI/3.0) <= opt.delta) {
			if (radius >= r) {
				lclusterdata.isinsector[is]=1;
				csectp1++;
			}
		}
		
		if (abs(tracephase) < opt.delta) {
			if (radius >= r) {
				lclusterdata.isinsector[is]=0;
				csect0++;
			}
		}
		
		if (abs(tracephase + 2.0*M_PI/3.0) <= opt.delta) {
			if (radius >= r) {
				lclusterdata.isinsector[is]=-1;
				csectm1++;
			}
		}
		lclusterdata.poll[is]=trace;
	}
	
	lclusterdata.nsectm1=csectm1; lclusterdata.nsect0=csect0; lclusterdata.nsectp1=csectp1;
	#ifdef DEBUG
		std::cout << "done!" << std::endl;
	#endif
}

void fillSectors(Clusterstruct &lclusterdata, std::vector<std::vector<std::complex<double> > > &pollev, Options opt, double delta) {
	// Categorizes the lattice points by sector
	#ifdef DEBUG
		std::cout << "Categorizing lattice points by sector (delta = " << delta << " )... " << std::flush;
	#endif

	int csectp1=0;
	int csectm1=0;
	int csect0=0;
	
  std::complex<double> trace;
	double tracephase=0;
	
	#ifdef DEBUG
		std::cout << "Delta = " << delta << std::endl;
	#endif

	lclusterdata.poll.resize(opt.Nspace);
	
	for (int is=0; is<opt.Nspace; is++)
		lclusterdata.isinsector[is] = 2;
	
	for (int is=0; is<opt.Nspace; is++) {
    trace=0;
    for(int i=0;i<opt.matrixdim;i++) {
      trace+=pollev.at(is).at(i);
    }
		tracephase = arg(trace);
		
		if (abs(tracephase - 2.0*M_PI/3.0) <= delta) {
			lclusterdata.isinsector[is]=1;
			csectp1++;
		}
		
		if (abs(tracephase) < delta) {
			lclusterdata.isinsector[is]=0;
			csect0++;
		}
		
		if (abs(tracephase + 2.0*M_PI/3.0) <= delta) {
			lclusterdata.isinsector[is]=-1;
			csectm1++;
		}

		lclusterdata.poll[is]=trace;
	}
	
	lclusterdata.nsectm1=csectm1; lclusterdata.nsect0=csect0; lclusterdata.nsectp1=csectp1;
	#ifdef DEBUG
		std::cout << "done!" << std::endl;
	#endif
}

void findPercolatingCluster(Clusterstruct &lclusterdata, Options opt) {
	// Finds the percolating clusters
	#ifdef DBEUG
		std::cout << "Searching for percolating clusters... " << flush;
	#endif
	
	int is=0;
	
	int cluster=0;
	
	lclusterdata.clusterispercolating.resize(lclusterdata.clustermembers.size());
	for (unsigned int c=0; c<lclusterdata.clustermembers.size(); c++) {
		lclusterdata.clusterispercolating[c] = 0;
	}
	
	std::vector<std::vector<std::vector<int> > > members;
	members.resize(lclusterdata.clustermembers.size());
	for (unsigned int i=0; i<members.size(); i++) {
		members[i].resize(3);
		members[i][0].resize(opt.leng1);
		for (int j=0; j<opt.leng1; j++)
			members[i][0][j]=0;
		members[i][1].resize(opt.leng2);
		for (int j=0; j<opt.leng2; j++)
			members[i][1][j]=0;
		members[i][2].resize(opt.leng3);
		for (int j=0; j<opt.leng3; j++)
			members[i][2][j]=0;
	}
	
	for (int i1=0; i1<opt.leng1; i1++)
	for (int i2=0; i2<opt.leng2; i2++)
	for (int i3=0; i3<opt.leng3; i3++) {
		is = latmap(i1, i2, i3);
		cluster = lclusterdata.isincluster[is];
		if (lclusterdata.clustersector[cluster] < 2) {
			members[cluster][0][i1] = 1;
			members[cluster][1][i2] = 1;
			members[cluster][2][i3] = 1;
		}
	}
	
	int sum1, sum2, sum3, percclustercnt=-1;
	for (unsigned int c=0; c<lclusterdata.clustermembers.size(); c++) {
		sum1=0; sum2=0; sum3=0;
		for (int i1=0; i1<opt.leng1; i1++)
			sum1 += members[c][0][i1];
		for (int i2=0; i2<opt.leng2; i2++)
			sum2 += members[c][1][i2];
		for (int i3=0; i3<opt.leng3; i3++)
			sum3 += members[c][2][i3];
			
		if (sum1 == opt.leng1 || sum2 == opt.leng2 || sum3 == opt.leng3) {
			percclustercnt++;
			lclusterdata.percolatingclusters.push_back(c);
			lclusterdata.percolatingdirections.push_back(std::vector<int>());
			lclusterdata.percolatingdirections[percclustercnt].resize(3);
			lclusterdata.percolatingdirections[percclustercnt][0]=0;
			lclusterdata.percolatingdirections[percclustercnt][1]=0;
			lclusterdata.percolatingdirections[percclustercnt][2]=0;
			if (sum1 == opt.leng1)
				lclusterdata.percolatingdirections[percclustercnt][0]=1;
			if (sum2 == opt.leng2)
				lclusterdata.percolatingdirections[percclustercnt][1]=1;
			if (sum3 == opt.leng3)
				lclusterdata.percolatingdirections[percclustercnt][2]=1;
			
			#ifdef DEBUG
				std::cout << "sum1 = " << sum1 << " sum2 = " << sum2 << " sum3 = " << sum3 << std::endl;
				std::cout << "Cluster c = " << c << " added!" << std::endl;
				std::cout << "Percolating in direction (" << lclusterdata.percolatingdirections[percclustercnt][0] << ","
					<< lclusterdata.percolatingdirections[percclustercnt][1] << "," 
					<< lclusterdata.percolatingdirections[percclustercnt][2] << ")." << std::endl;
			#endif
		
			if (lclusterdata.percolatingdirections[percclustercnt][0] == 1 && 
					lclusterdata.percolatingdirections[percclustercnt][1] == 1 &&
					lclusterdata.percolatingdirections[percclustercnt][2] == 1) {
				lclusterdata.clusterispercolating[c]=1;
			}
		}
	}
	#ifdef DBEUG	
		std::cout << "done!" << std::endl;
	#endif
}

void findClusters(Clusterstruct &lclusterdata, std::vector<std::vector<int> > &neib, Options opt) {
	// Find clusters
	#ifdef DEBUG
		std::cout << "Finding clusters... " << std::flush;
	#endif
	int sector = 0, cluster = 0, memberis = 0, isneib = 0;

	std::vector<bool> isfiled;
	isfiled.resize(opt.Nspace);
	for (int is=0; is<opt.Nspace; is++) {
		isfiled[is] = false;
		lclusterdata.isincluster[is] = -1;
	}
	
	int totalmembers=0;

	cluster = -1; // Note: We begin counting with ID 0!
	for (int is=0; is<opt.Nspace; is++) {
		if (isfiled[is] == false) {
			sector = lclusterdata.isinsector[is];
			cluster++; // current cluster id (Note: We begin counting with ID 0!)
			lclusterdata.clustersector.push_back(sector); // Current cluster is in sector sector
			
			lclusterdata.clustermembers.push_back(std::vector<int>()); // Add a new cluster
			lclusterdata.clustermembers[cluster].push_back(is); // Add is to this new cluster

			lclusterdata.isincluster[is] = cluster; // Set the cluster id for the point  is
			isfiled[is] = true; // Lattice point is now filed
			
			// Modification due to cluster radius with PBCs
			lclusterdata.clusterisperiodic.push_back(std::vector<int>()); // Add a new cluster
			lclusterdata.clusterisperiodic[cluster].resize(3);
			for (int i=0; i<3; i++)
				lclusterdata.clusterisperiodic[cluster][i]=0;
			int ci1=0, ci2=0, ci3=0;
		
			// Iterate over all found points in the actual cluster (loop does not have fixed length)
			int member = 0;
			while ((unsigned int)member < lclusterdata.clustermembers[cluster].size()) {
				memberis = lclusterdata.clustermembers[cluster][member];
			
				for (int mu=0; mu<6; mu++) {
					isneib = neib[memberis][mu];
					if (isfiled[isneib] == false && lclusterdata.isinsector[isneib] == sector) {
						lclusterdata.clustermembers[cluster].push_back(isneib);
						lclusterdata.isincluster[isneib] = cluster;
						isfiled[isneib] = true;

						// Check if over periodic boundaries
						getCoords(memberis, ci1, ci2, ci3);
						switch (mu) {
							case 0: if(ci1 + 1 == opt.leng1)
									lclusterdata.clusterisperiodic[cluster][0]=1;
								break;

							case 1: if(ci2 + 1 == opt.leng2)
									lclusterdata.clusterisperiodic[cluster][1]=1;
								break;
							
							case 2: if(ci3 + 1 == opt.leng3)
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

							default: std::cout << "ERROR: Error in cluster switch!" << std::endl;
						}
					}
				}
				
				#ifdef DEBUG
					std::cout << std::endl <<  "Cluster " << cluster << " member " << member + 1 << " of " << lclusterdata.clustermembers[cluster].size() << " members." << std::endl;
				#endif
				
				member++;
				totalmembers++;
			}
		}
	}
	
	#ifdef DEBUG
		std::cout << "" << totalmembers << " of " << opt.Nspace << " lattice points visited!" << std::endl;
	#endif
	
	for (int is=0; is<opt.Nspace; is++) {
		if (isfiled[is] == false)
			std::cout << "WARNING: Some points did not get filled!" << std::endl;
	}
	#ifdef DEBUG	
		std::cout << "done!" << std::endl;
	#endif
}

void checkClusters(Clusterstruct &lclusterdata, Options opt) {
	// Checks the clusters
	#ifdef DEBUG	
		std::cout << "Checking clusters... " << flush;
	#endif

	for (int is=0; is<opt.Nspace; is++) {
		if (lclusterdata.isincluster[is] < 0)
			std::cout << "WARNING: Some points do not have a valid cluster ID!" << std::endl;
	}
	
	int sum=0;
	for (unsigned int c=0; c<lclusterdata.clustermembers.size(); c++)
		sum+=lclusterdata.clustermembers[c].size();
	if (sum != opt.Nspace) {
		std::cout << "WARNING: Number of lattice points differ from number of clustermember entries!" << std::endl;
		std::cout << "Number of entries = " << sum << " Number of lattice points = " << opt.Nspace << std::endl;
	}
	
	// Check for conistency of the cluster arrays lclusterdata.isinsector and lclusterdata.clustersector
	for (unsigned int c=0; c<lclusterdata.clustermembers.size(); c++) {
		for (unsigned int member=0; member<lclusterdata.clustermembers[c].size(); member++) {
			if (lclusterdata.isinsector[lclusterdata.clustermembers[c][member]] != lclusterdata.clustersector[c]) {
				std::cout << "WARNING: Problem with lclusterdata.clustermembers and lclusterdata.clustersector consistency!" << std::endl;
			}
		}
	}
	
	// Check if all members in clustermmembers have a correct ID
	// Here we use vector.at(iter) so that also a range check is performed and an
	// exeption is thrown if it is outside vector range
	int mult=1;
	std::vector<int> membercheck;
	membercheck.resize(opt.Nspace);
	for (unsigned int c=0; c<lclusterdata.clustermembers.size(); c++) {
		for (unsigned int member=0; member<lclusterdata.clustermembers[c].size(); member++) {
			membercheck.at(lclusterdata.clustermembers[c][member]) = 1;
		}
	}
	for (int is=0; is<opt.Nspace; is++) {
		mult = mult*membercheck[is];
	}
	if (mult != 1) {
		std::cout << "WARNING: Problem with clustermember IDs!" << std::endl;
	}
	
	#ifdef DEBUG	
		std::cout << "done!" << std::endl;
	#endif
}
