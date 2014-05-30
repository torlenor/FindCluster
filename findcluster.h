/*
 * findcluster.h - Polyakov loop center cluster calculation - main header
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

#ifndef FINDCLUSTER_H
#define FINDCLUSTER_H

#include <complex>
#include <vector>
#include <string>

// options struct
struct Options{
  int Ns, Nt; // Spatial and temporal lattice extent
  int matrixdim, leng1, leng2, leng3, leng4, Nspace;

  bool wupperdata; // if true, marks that we are reading wuppertal data, if off we read the old polyakov loop ev structure

  int nmeas; // Number of measurements/configurations to use
  
  double fraction; // fraction to cut
  double delta; // delta = pi/3 *(1 - fraction), will be set in init

  bool usealternativesectors; // Controlls if we want the alternative sector identification
  double r; // radius for alternative sector identification

  bool do3d; // Controlls if we want to write out data for 3d visualization

  bool detail; // Controlls if we want detailed information for every configuration
  bool doboxes; // Controlls if we want box counting calculations
  bool dodistance; // Controlls if we want box counting calculations
  bool doradius; // Controlls if we want radius calculation
  bool domean; // Controlls if we want mean distance traveled calculations

  bool fastmode; // Fast mode: Only largest radius calc for f determination

  bool writemeas; // Controlls if we want to writeout all measurements

  std::vector<std::string> fevname; // Filenames of configurations
};

// Vectors to store cluster information
struct Clusterstruct{
	std::vector<int> isinsector;  // lclusterdata.isinsector[is] stores the sector of the lattice point is
	std::vector<int> clustersector; // lclusterdata.clustersector[c] stores the sector of the cluster c
	std::vector<int> isincluster;  // lclusterdata.isincluster[is] stores the cluster of the lattice point is
	std::vector<std::vector<int> > clustermembers; // lclusterdata.clustermembers[c][i] stores the members/lattice points of the cluster c (0 < i < N_c)

	std::vector<int> percolatingclusters; // Percolating clusters
	std::vector<std::vector<int> > percolatingdirections;
	std::vector<int> clusterispercolating;
	
	std::vector<int> sortedcluster;
	std::vector<int> sortedrealcluster;
	std::vector<int> isinsortedcluster;

	std::vector<std::vector<int> > clusterisperiodic;

	int nsectm1, nsect0, nsectp1;

	std::vector<std::complex<double> > poll;
};

struct Observablestruct{
	int maxclustersize;
	int maxclusterid;
	int maxclustersector;

  double largestnonpercclustersize;
	double largestnonpercclusterid;
	double largetsnonpercclustersector;
	
	double avgclustersize;
	double avgclustersizeF;

	double rootmeansquaredistanceR;
	
	double avgclustersizenp;
	double avgclustersizeFnp;

	double cut;

	std::vector<std::vector<int> > numberofboxes;
	
	// std::vector<std::vector<double> > centerofmass;
	std::vector<double> clusterradius;
	double largestclusterradius;
	double largestnpclusterradius;
	double avgclusterradius;
	double avgnpclusterradius;

	double percc;

	int largestclusterid;

	double area;
	double arealargestnonperccluster;
	double areaavgnonperccluster;

	double poll;
  
  std::vector<double> meanfreepathnew;
  double largestclustermeanfreepathnew;
  double largestnpclustermeanfreepathnew;
  double avgclustermeanfreepathnew;
  double avgnpclustermeanfreepathnew;
  double avgFclustermeanfreepathnew;
  double avgFnpclustermeanfreepathnew;
};

struct Resultstruct{
	double maxclustersize, maxclustersizeerr;
	double maxnonpercclustersize, maxnonpercclustersizeerr;
	double avgclustersize, avgclustersizeerr;
	double avgclusersizeFortunato, avgclusersizeFortunatoerr;
	double avgclustersizenp, avgclustersizenperr;
	double avgclustersizeFnp, avgclustersizeFnperr;
	double avgrootmeansquaredistance, avgrootmeansquaredistanceerr;
	double cut, cuterr;
	double totalperimeter, totalperimetererr;
	double largestnonpercperimeter, largestnonpercperimetererr;
	double avgnonpercperimeter, avgnonpercperimetererr;
	double largestclusterradius, largestclusterradiuserr;
	double largestnpclusterradius, largestnpclusterradiuserr;
	double avgclusterradius, avgclusterradiuserr;
	double avgnpclusterradius, avgnpclusterradiuserr;
	double polyakovloopaftercut, polyakovloopaftercuterr;
	double avgperccluster, avgpercclustererr;
	std::vector<double> largestclusterboxcount, largestclusterboxcounterr;
	std::vector<double> largestnonpercboxcount, largestnonpercboxcounterr;
  
  double largestclustermeanfreepathnew, largestclustermeanfreepathnewerr;
  double largestnpclustermeanfreepathnew, largestnpclustermeanfreepathnewerr;
  double avgclustermeanfreepathnew, avgclustermeanfreepathnewerr;
  double avgnpclustermeanfreepathnew, avgnpclustermeanfreepathnewerr;
  double avgFclustermeanfreepathnew, avgFclustermeanfreepathnewerr;
  double avgFnpclustermeanfreepathnew, avgFnpclustermeanfreepathnewerr;
};

void writeConfigResultsstdout(Observablestruct &lobs, Clusterstruct &lclusterdata);
void writeClusterList(Clusterstruct &lclusterdata);

void calcExp();

void sortClusterSize(Clusterstruct &lclusterdata);
void cluster3doutput(Clusterstruct &clusterdata, std::string f3dname);

int latmap(int i1, int i2, int i3);
void getCoords(int is, int &i1, int &i2, int &i3);

void freeMem(Clusterstruct &lclusterdata);

void fillNeib();

#endif // FINDCLUSTER_H
