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
#include "findcluster_writemeasure.h"

#include "include/jackknife.h"
#include "include/readowndata.hpp"
#include "include/readwupperdata.hpp"
#include "version.h"

using namespace std;

Options opt;

// Vector for Polyakov loop eigenvalue data
vector<vector<complex<double> > > pollev;
// Neib vector
vector<vector<int> > neib;
// Box size vector
vector<int> boxsize;
vector<int> boxes;

Observablestruct *obs;
Resultstruct results;

#include "findcluster_helper.hpp"
#include "findcluster_init.hpp"
#include "findcluster_obs.hpp"
#include "findcluster_write.hpp"

int main(int argc, char *argv[]) {
	// Handle command line information and initialize arrays and other stuff...
	if (init(argc, argv) != 0) {
		cout << "Error in init()!" << endl;
		return 1;
	}

	cout << endl;
	printsettings();
	cout << endl;

  if(opt.writemeas) {
    prepwriteMeasure(opt);
  }

	// Main part
	// Loop over all nmeas configurations
	cout << "------------------------------------------------------------------------------" << endl;
	for (int n=0; n<opt.nmeas; n++) {
		// Allocate memory for faster access
    Clusterstruct lclusterdata;
    lclusterdata.isinsector.resize(opt.Nspace);
    lclusterdata.isincluster.resize(opt.Nspace);

		if (opt.detail) {
			cout << endl << "------------------------------------------------------------------------------" << endl;
			cout << n+1 << " of "<< opt.nmeas << ": " << opt.fevname[n] << " ... " << endl;
		} else {
			cout << "\r" <<  "                                                                   " << flush;
			cout << "\r" << n+1 << " of "<< opt.nmeas << ": " << opt.fevname[n] << " ... " << flush;
		}
		// Read Polyakov loop eigenvalues file
    cout << "r" << flush;
		if (opt.wupperdata) {
      if (readWupperPollBinary(opt.leng1, opt.leng2, opt.leng3, opt.leng4, opt.matrixdim, pollev, opt.fevname[n]) != 0) {
        cout << "ERROR: Problems with readWupperPollBinary !" << endl;
        return 1;
      }
    } else {
      if (readPollEvBinary(opt.leng1, opt.leng2, opt.leng3, opt.leng4, opt.matrixdim, pollev, opt.fevname[n]) != 0) {
        cout << "ERROR: Problems with readPollEVBinary !" << endl;
        return 1;
      }
      // Check Polyakov loop eigenvalues
      checkPollEv(opt.leng1, opt.leng2, opt.leng3, opt.leng4, opt.matrixdim, pollev);
    }
		
		if (opt.usealternativesectors==true) {
			fillSectorsAlt(lclusterdata, pollev, opt, opt.r); // Categorize lattice points by sectors using alternative prescription
		} else {
			fillSectors(lclusterdata, pollev, opt, opt.delta); // Categorize lattice points by sectors
		}
	
    cout << "c" << flush;
		findClusters(lclusterdata, neib, opt);	// Identify clusters
		checkClusters(lclusterdata, opt); // Check clusters
    cout << "p" << flush;
		findPercolatingCluster(lclusterdata, opt); // Find percolating clusters

    cout << "s" << flush;
		sortClusterSize(lclusterdata, opt); // Sort clusters per number of members

    cout << "o" << flush;
		calcObservables(obs[n], lclusterdata); // Calculate observables
		
		if (opt.detail) {
			cout << endl << "Details for " << opt.fevname[n] << ":" << endl;
			writeConfigResultsstdout(obs[n], lclusterdata);
		}
		
		if (opt.do3d) {
			// Write data for 3dclusters program
			stringstream f3dclustername;
			f3dclustername << "3dcluster_" << opt.leng1 << "x" << opt.leng4 << "_m" << setprecision(0) << fixed << n << ".data";
			cluster3doutput(lclusterdata, f3dclustername.str());
		}
	
    if (opt.writemeas) {
		  writeMeasures(obs[n], opt, n);
  	}
	}

	if (opt.do3d) {
		// Write clusterfilenamelist for 3dcluster program
		ofstream f3dclusterlist;
		f3dclusterlist.open("3dcluster.list");
		for (int n=0; n<opt.nmeas; n++) {
			stringstream f3dclustername;
			f3dclustername << "3dcluster_" << opt.leng1 << "x" << opt.leng4 << "_m" << setprecision(0) << fixed << n << ".data";
			f3dclusterlist << f3dclustername.str() << endl;
		}
		f3dclusterlist.close();
	}

	cout << endl;
	cout << "------------------------------------------------------------------------------" << endl;

	// Calculate the expectation values and write the results to file and stdout
	calcExp();
  if (opt.nmeas < 2) {
	  writeresultsstdout_singleconf(obs[0]);
  } else {
	  writeresultsstdout();
  }
	writeresults();

	// delete [] clusterdata; clusterdata=0;
	delete [] obs; obs=0;
	
	return 0;
}

