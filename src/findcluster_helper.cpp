/*
 * findcluster_helper.cpp - various functions
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

#include "findcluster_helper.h"

#include <iostream>
#include <vector>

#include "findcluster.h"

void fillNeib(const Options &opt, std::vector<std::vector<int> > &neib) {
	// Fills the neib array
  int i1p,i2p,i3p,i1m,i2m,i3m,is,isp1,isp2,isp3,ism1,ism2,ism3;
  for (int i1=0; i1<opt.leng1; i1++) {
    i1p = i1 + 1;
    i1m = i1 - 1;
    if (i1p == opt.leng1) i1p = 0;
    if (i1m == -1) i1m = opt.leng1-1;

    for (int i2=0; i2<opt.leng2; i2++) {
      i2p = i2 + 1;
      i2m = i2 - 1;
      if (i2p == opt.leng2) i2p = 0;
      if (i2m == -1) i2m = opt.leng2-1;

      for (int i3 = 0;i3<opt.leng3;i3++) {
        i3p = i3 + 1;
        i3m = i3 - 1;
        if (i3p == opt.leng3) i3p = 0;
        if (i3m == -1) i3m = opt.leng3-1;
        
        // Compute the site address and the addresses of the sites shifted
        // by one unit in each direction
        is = i1 + i2*opt.leng1 + i3*opt.leng1*opt.leng2;

        isp1 = i1p + i2*opt.leng1 + i3*opt.leng1*opt.leng2;
        isp2 = i1 + i2p*opt.leng1 + i3*opt.leng1*opt.leng2;
        isp3 = i1 + i2*opt.leng1 + i3p*opt.leng1*opt.leng2;

        ism1 = i1m + i2*opt.leng1 + i3*opt.leng1*opt.leng2;
        ism2 = i1 + i2m*opt.leng1 + i3*opt.leng1*opt.leng2;
        ism3 = i1 + i2*opt.leng1 + i3m*opt.leng1*opt.leng2;

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

int latmap(const int i1, const int i2, const int i3, const Options &opt) {
	return i1 + i2*opt.leng1 + i3*opt.leng1*opt.leng2;
}

void Printsettings(const Options &opt) {
	std::cout << "Settings:" << std::endl << std::endl;
	std::cout << "Lattice size = " << opt.leng1 << "x" << opt.leng2 << "x" << opt.leng3 << "x" << opt.leng4 << std::endl;
	if (opt.wupperdata)
		std::cout << "Reading Wuppertal Polakov loop data format." <<  std::endl;
	std::cout << "Number of configurations = " << opt.nmeas << std::endl;
	std::cout << "Cut fraction = " << opt.fraction << std::endl;
	if (opt.doboxes)
		std::cout << "Calculating 'box' observables." << std::endl;
	if (opt.detail)
		std::cout << "Writing detailed results for every configuration." <<  std::endl;
	if (opt.do3d)
		std::cout << "Writing 3dcluster visualization data files." <<  std::endl;
	if (opt.writemeas)
		std::cout << "Writing all measurements to file." <<  std::endl;
	if (opt.fastmode)
		std::cout << "FAST MODE: Fast mode activated! Calculating only largest cluster radius!" <<  std::endl;
}

void freeMem(Clusterstruct &lclusterdata) {
	lclusterdata.isinsector.resize(0);
	lclusterdata.clustersector.resize(0);
	lclusterdata.isincluster.resize(0);

	for (unsigned c=0;c<lclusterdata.clustermembers.size();c++) {
		lclusterdata.clustermembers[c].resize(0);
	}
	lclusterdata.clustermembers.resize(0);

	lclusterdata.percolatingclusters.resize(0);
	for (unsigned p=0; p<lclusterdata.percolatingdirections.size(); p++) {
		lclusterdata.percolatingdirections[p].resize(0);
	}
	lclusterdata.percolatingdirections.resize(0);
    
	lclusterdata.sortedcluster.resize(0);
	lclusterdata.sortedrealcluster.resize(0);
	lclusterdata.isinsortedcluster.resize(0);
}

void getCoords(const Options &opt, const int is, int &i1, int &i2, int &i3) {
	i1 = (is % (opt.leng1*opt.leng2) ) % opt.leng1;
	i2 = (is % (opt.leng1*opt.leng2) ) / opt.leng1;
	i3 = is / (opt.leng1*opt.leng2);

	if (is != i1 + i2*opt.leng1 + i3*opt.leng1*opt.leng2)
		std::cout << "ERROR: Problem in getCoords!" << std::endl;
}