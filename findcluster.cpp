/*
 * findcluster.cpp - Polyakov loop center cluster calculation - main file
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

#include "findcluster.h"

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <vector>

#include "findcluster_cluster.h"
#include "findcluster_obs.h"
#include "findcluster_write.h"
#include "findcluster_writemeasure.h"

#include "include/readowndata.hpp"
#include "include/readwupperdata.hpp"

#include "version.h"

Options opt; // Global variable declared extern in findcluster.h

// Vector for Polyakov loop eigenvalue data
std::vector<std::vector<std::complex<double> > > pollev;
// Neib vector
std::vector<std::vector<int> > neib;
// Box size vector
std::vector<int> boxsize;
std::vector<int> boxes;

#include "findcluster_init.hpp"

int main(int argc, char *argv[]) {
  std::vector<Observablestruct> obs;
  Resultstruct results;

	// Handle command line information and initialize arrays and other stuff...
	if (init(argc, argv, obs, results) != 0) {
		std::cout << "Error in init()!" << std::endl;
		return 1;
	}

	std::cout << std::endl;
	Printsettings(opt);
	std::cout << std::endl;

  if(opt.writemeas) {
    prepwriteMeasure(opt);
  }

	// Main part
	// Loop over all nmeas configurations
	std::cout << "------------------------------------------------------------------------------" << std::endl;
	for (int n=0; n<opt.nmeas; n++) {
		// Allocate memory for faster access
    Clusterstruct lclusterdata;
    lclusterdata.isinsector.resize(opt.Nspace);
    lclusterdata.isincluster.resize(opt.Nspace);

		if (opt.detail) {
			std::cout << std::endl << "------------------------------------------------------------------------------" << std::endl;
			std::cout << n+1 << " of "<< opt.nmeas << ": " << opt.fevname[n] << " ... " << std::endl;
		} else {
			std::cout << "\r" <<  "                                                                   " << std::flush;
			std::cout << "\r" << n+1 << " of "<< opt.nmeas << ": " << opt.fevname[n] << " ... " << std::flush;
		}
		// Read Polyakov loop eigenvalues file
    std::cout << "r" << std::flush;
		if (opt.wupperdata) {
      if (readWupperPollBinary(opt.leng1, opt.leng2, opt.leng3, opt.leng4, opt.matrixdim, pollev, opt.fevname[n]) != 0) {
        std::cout << "ERROR: Problems with readWupperPollBinary !" << std::endl;
        return 1;
      }
    } else {
      if (readPollEvBinary(opt.leng1, opt.leng2, opt.leng3, opt.leng4, opt.matrixdim, pollev, opt.fevname[n]) != 0) {
        std::cout << "ERROR: Problems with readPollEVBinary !" << std::endl;
        return 1;
      }
      // Check Polyakov loop eigenvalues
      checkPollEv(opt.leng1, opt.leng2, opt.leng3, opt.leng4, opt.matrixdim, pollev);
    }
		
    fillSectors(lclusterdata, pollev, opt, opt.delta); // Categorize lattice points by sectors
	
    std::cout << "c" << std::flush;
		findClusters(lclusterdata, neib, opt);	// Identify clusters
		checkClusters(lclusterdata, opt); // Check clusters
    std::cout << "p" << std::flush;
		findPercolatingCluster(lclusterdata, opt); // Find percolating clusters

    std::cout << "s" << std::flush;
		sortClusterSize(lclusterdata, opt); // Sort clusters per number of members

    std::cout << "o" << std::flush;
		CalcObservables(obs[n], lclusterdata, boxsize, boxes); // Calculate observables
		
		if (opt.detail) {
			std::cout << std::endl << "Details for " << opt.fevname[n] << ":" << std::endl;
			writeConfigResultsstdout(boxsize, boxes, obs[n], lclusterdata, opt);
		}
		
		if (opt.do3d) {
			// Write data for 3dclusters program
      std::stringstream f3dclustername;
			f3dclustername << "3dcluster_" << opt.leng1 << "x" << opt.leng4 << "_m" << std::setprecision(0) << std::fixed << n << ".data";
			cluster3doutput(lclusterdata, f3dclustername.str(), opt);
		}
	
    if (opt.writemeas) {
		  writeMeasures(obs[n], opt, n);
  	}
	}

	if (opt.do3d) {
		// Write clusterfilenamelist for 3dcluster program
    std::ofstream f3dclusterlist;
		f3dclusterlist.open("3dcluster.list");
		for (int n=0; n<opt.nmeas; n++) {
      std::stringstream f3dclustername;
			f3dclustername << "3dcluster_" << opt.leng1 << "x" << opt.leng4 << "_m" << std::setprecision(0) << std::fixed << n << ".data";
			f3dclusterlist << f3dclustername.str() << std::endl;
		}
		f3dclusterlist.close();
	}

	std::cout << std::endl;
	std::cout << "------------------------------------------------------------------------------" << std::endl;

	// Calculate the expectation values and write the results to file and stdout
	CalcExp(obs, results, opt, boxsize, boxes);
  if (opt.nmeas < 2) {
	  writeresultsstdout_singleconf(obs[0], opt);
  } else {
	  writeresultsstdout(results, opt);
  }
	writeresults(boxsize, boxes, results, opt);

	obs.resize(0);
	
	return 0;
}

