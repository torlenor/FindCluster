/*
 * findcluster_box.cpp - Box counting calculations
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

#include <iostream>
#include <vector>

#include "findcluster.h"
#include "findcluster_helper.h"

void obsBoxesOnlyLargest(Observablestruct &lobs, Clusterstruct &lclusterdata, Options opt, std::vector<int> &boxsize, std::vector<int> &boxes) {
	lobs.numberofboxes.resize(lclusterdata.clustermembers.size());
	for(unsigned int c=0;c<lobs.numberofboxes.size();c++){
		lobs.numberofboxes[c].resize(boxsize.size());
	}

	int c=0;
	int boxcnt=0, is, i1, i2, i3;
	bool clusterinbox=false;

	// Calculation for largestcluster
	c=lobs.largestclusterid;
	// We should not calculate trivial box sizes!
	lobs.numberofboxes[c][0]=lclusterdata.clustermembers[c].size();
	lobs.numberofboxes[c][boxsize.size()-1]=1;
	for(unsigned int size=1;size<boxsize.size()-1;size++){
	// for(unsigned int size=0;size<boxsize.size();size++){
		boxcnt=0;
		// Loop over all boxes
		for(int box1=0;box1<boxes[size];box1++)
		for(int box2=0;box2<boxes[size];box2++)
		for(int box3=0;box3<boxes[size];box3++){
			clusterinbox=false;
			// Loop over all points in box
			for(int b1=0;b1<boxsize[size];b1++)
			for(int b2=0;b2<boxsize[size];b2++)
			for(int b3=0;b3<boxsize[size];b3++){
				i1 = b1 + box1*boxsize[size];
				i2 = b2 + box2*boxsize[size];
				i3 = b3 + box3*boxsize[size];
				
				is = latmap(i1, i2, i3, opt);

				if(lclusterdata.isincluster[is] == (int)c)
					clusterinbox=true;
			}

			if(clusterinbox==true)
				boxcnt++;

		} // Boxes in all directions
		lobs.numberofboxes[c][size]=boxcnt;
	} // Boxsize
	if(lobs.numberofboxes[c][0] != (int)lclusterdata.clustermembers[c].size()){
		std::cout << "ERROR: Number of boxes for boxsize = 1 has to be equal to number of cluster elements!" << std::endl;
	}
	if(lobs.numberofboxes[c][lobs.numberofboxes[c].size()-1] != 1){
		std::cout << "ERROR: Number of boxes for boxsize = Ns has to be equal 1!" << std::endl;
	}
	
	// Calculation for largest non percolating cluster
	c=lobs.largestnonpercclusterid;
	// We should not calculate trivial box sizes!
	lobs.numberofboxes[c][0]=lclusterdata.clustermembers[c].size();
	lobs.numberofboxes[c][boxsize.size()-1]=1;
	for(unsigned int size=1;size<boxsize.size()-1;size++){
	// for(unsigned int size=0;size<boxsize.size();size++){
		boxcnt=0;
		// Loop over all boxes
		for(int box1=0;box1<boxes[size];box1++)
		for(int box2=0;box2<boxes[size];box2++)
		for(int box3=0;box3<boxes[size];box3++){
			clusterinbox=false;
			// Loop over all points in box
			for(int b1=0;b1<boxsize[size];b1++)
			for(int b2=0;b2<boxsize[size];b2++)
			for(int b3=0;b3<boxsize[size];b3++){
				i1 = b1 + box1*boxsize[size];
				i2 = b2 + box2*boxsize[size];
				i3 = b3 + box3*boxsize[size];
				
				is = latmap(i1, i2, i3, opt);

				if(lclusterdata.isincluster[is] == (int)c)
					clusterinbox=true;
			}

			if(clusterinbox==true)
				boxcnt++;

		} // Boxes in all directions
		lobs.numberofboxes[c][size]=boxcnt;
	} // Boxsize
	if(lobs.numberofboxes[c][0] != (int)lclusterdata.clustermembers[c].size()){
		std::cout << "ERROR: Number of boxes for boxsize = 1 has to be equal to number of cluster elements!" << std::endl;
	}
	if(lobs.numberofboxes[c][lobs.numberofboxes[c].size()-1] != 1){
		std::cout << "ERROR: Number of boxes for boxsize = Ns has to be equal 1!" << std::endl;
	}
}

void obsBoxes(Observablestruct &lobs, Clusterstruct &lclusterdata, Options opt, std::vector<int> &boxsize, std::vector<int> &boxes) {
	lobs.numberofboxes.resize(lclusterdata.clustermembers.size());
	for(unsigned int c=0;c<lobs.numberofboxes.size();c++){
		lobs.numberofboxes[c].resize(boxsize.size());
	}

	int boxcnt=0, is, i1, i2, i3;
	bool clusterinbox=false;
	for(unsigned int c=0; c<lclusterdata.clustermembers.size();c++){
		// We should not calculate trivial box sizes!
		lobs.numberofboxes[c][0]=lclusterdata.clustermembers[c].size();
		lobs.numberofboxes[c][boxsize.size()-1]=1;
		for(unsigned int size=1;size<boxsize.size()-1;size++){
		// for(unsigned int size=0;size<boxsize.size();size++){
			boxcnt=0;
			// Loop over all boxes
			for(int box1=0;box1<boxes[size];box1++)
			for(int box2=0;box2<boxes[size];box2++)
			for(int box3=0;box3<boxes[size];box3++){
				clusterinbox=false;
				// Loop over all points in box
				for(int b1=0;b1<boxsize[size];b1++)
				for(int b2=0;b2<boxsize[size];b2++)
				for(int b3=0;b3<boxsize[size];b3++){
					i1 = b1 + box1*boxsize[size];
					i2 = b2 + box2*boxsize[size];
					i3 = b3 + box3*boxsize[size];
					
					is = latmap(i1, i2, i3, opt);

					if(lclusterdata.isincluster[is] == (int)c)
						clusterinbox=true;
				}

				if(clusterinbox==true)
					boxcnt++;

			} // Boxes in all directions
			lobs.numberofboxes[c][size]=boxcnt;
		} // Boxsize
		if(lobs.numberofboxes[c][0] != (int)lclusterdata.clustermembers[c].size()){
			std::cout << "ERROR: Number of boxes for boxsize = 1 has to be equal to number of cluster elements!" << std::endl;
		}
		if(lobs.numberofboxes[c][lobs.numberofboxes[c].size()-1] != 1){
			std::cout << "ERROR: Number of boxes for boxsize = Ns has to be equal 1!" << std::endl;
		}
	} // Cluster
}
