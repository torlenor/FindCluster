/*
 * findcluster_radius.h - Cluster radius calculations
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

#include "findcluster_radius.h"

void getCoordsShift(int is, int &i1, int &i2, int &i3, int *shift, const Options &opt){
	i1 = (is % (opt.leng1*opt.leng2) ) % opt.leng1;
	i2 = (is % (opt.leng1*opt.leng2) ) / opt.leng1;
	i3 = is / (opt.leng1*opt.leng2);

	if(i1>shift[0])
	i1=i1-opt.leng1;
	if(i2>shift[1])
	i2=i2-opt.leng2;
	if(i3>shift[2])
	i3=i3-opt.leng3;

#ifdef DEBUG
	if (is != i1 + i2*opt.leng1 + i3*opt.leng1*opt.leng2)
		std::cout << "ERROR: Problem in getCoords!" << std::endl;
#endif
}

double ObsClusterRadius(const Clusterstruct &lclusterdata, const Options &opt, const int clusterid){
	// Calculation of the cluster radius for the cluster with the largest weight.
	
	double centerofmass[3], radiussquare, radiussquaremin;
	int i1, i2, i3;

	int shift[3];

	int c = clusterid;

	radiussquare=0;
	radiussquaremin=1E30;

	// Calculate center of mass
	for(int s1=0;s1<opt.leng1/2;s1++)
	for(int s2=0;s2<opt.leng2/2;s2++)
	for(int s3=0;s3<opt.leng3/2;s3++){
		shift[0]=s1 + 0.5;
		shift[1]=s2 + 0.5;
		shift[2]=s3 + 0.5;

		radiussquare=0;

		centerofmass[0]=0;
		centerofmass[1]=0;
		centerofmass[2]=0;
		for(unsigned int member=0; member<lclusterdata.clustermembers[c].size();member++){
			getCoordsShift(lclusterdata.clustermembers[c][member], i1, i2, i3, shift, opt);
			centerofmass[0] += i1;
			centerofmass[1] += i2; 
			centerofmass[2] += i3; 
		}

		centerofmass[0] = centerofmass[0]/(double)lclusterdata.clustermembers[c].size();
		centerofmass[1] = centerofmass[1]/(double)lclusterdata.clustermembers[c].size();
		centerofmass[2] = centerofmass[2]/(double)lclusterdata.clustermembers[c].size();

		for(unsigned int member=0; member<lclusterdata.clustermembers[c].size();member++){
			getCoordsShift(lclusterdata.clustermembers[c][member], i1, i2, i3, shift, opt);
			radiussquare += (pow(centerofmass[0] - i1, 2)
					+ pow(centerofmass[1] - i2, 2)
					+ pow(centerofmass[2] - i3, 2));
		}
		radiussquare = sqrt(radiussquare/(double)lclusterdata.clustermembers[c].size());
		if(radiussquare<radiussquaremin){
			radiussquaremin=radiussquare;
		}
	}
	
	return radiussquaremin;
}

void ObsClusterRadiusAll(Observablestruct &lobs, Clusterstruct &lclusterdata, const Options &opt){
	// Calculation of the cluster radius. We save the largest cluster (in terms
	// of the cluster radius).
	
	lobs.clusterradius.resize(lclusterdata.clustermembers.size());

	for(unsigned int c=0;c<lclusterdata.clustermembers.size();c++){
		if(lclusterdata.clustermembers[c].size()>1 && lclusterdata.clustersector[c] < 2){
			lobs.clusterradius[c] = ObsClusterRadius(lclusterdata, opt, c);
		}
	}
	lobs.largestclusterradius=lobs.clusterradius[lobs.maxclusterid];

	// Additional radius observables
	// Radius largest non-percolating cluster
	lobs.largestnpclusterradius=lobs.clusterradius[lobs.largestnonpercclusterid];

	// Average over all and over non-percolating clusters
	int cnt=0, cntnp=0;
	lobs.avgclusterradius=0;
	lobs.avgnpclusterradius=0;
	for(unsigned int c=0;c<lclusterdata.clustermembers.size();c++){
		if(lclusterdata.clustersector[c]<2){
			lobs.avgclusterradius += lobs.clusterradius[c];
			cnt++;
			if(lclusterdata.clusterispercolating[c] == 0){
				lobs.avgnpclusterradius += lobs.clusterradius[c];
				cntnp++;
			}
		}
	}

	lobs.avgclusterradius = lobs.avgclusterradius/(double)cnt;
	lobs.avgnpclusterradius = lobs.avgnpclusterradius/(double)cntnp;
}