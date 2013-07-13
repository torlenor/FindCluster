#ifndef FINDCLUSTER_OBS_HPP
#define FINDCLUSTER_OBS_HPP

#include "findcluster_radius.hpp"

void obsLargestCluster(Observablestruct &lobs, Clusterstruct &lclusterdata){
	// Find largest cluster
	int size = -1, largestcluster = -1;
	for(unsigned int c=0; c<lclusterdata.clustermembers.size(); c++){
		if( (int)lclusterdata.clustermembers[c].size() > size && lclusterdata.clustersector[c] < 2){
			size = lclusterdata.clustermembers[c].size();
			largestcluster = c;
		}
	}
	
	lobs.maxclustersize = size;
	lobs.maxclusterid = largestcluster;
	lobs.maxclustersector = lclusterdata.clustersector[largestcluster];
	
	// Find largest non percollating cluster
	largestcluster = -1;
	size = -1;
	for(unsigned int i=0;i<lclusterdata.sortedrealcluster.size();i++){
		if(lclusterdata.clusterispercolating[lclusterdata.sortedrealcluster[i]] == 0){
			largestcluster = lclusterdata.sortedrealcluster[i];
			size = lclusterdata.clustermembers[largestcluster].size();
			break;
		}
	}
	
	lobs.largestnonpercclustersize = size;
	lobs.largestnonpercclusterid = largestcluster;
	lobs.largetsnonpercclustersector = lclusterdata.clustersector[largestcluster];
}

void obsAverageClusterSize(Observablestruct &lobs, Clusterstruct &lclusterdata){
	// Calculate average cluster size
	double avgclustersize=0;
	int cnt=0;
	for(unsigned int c=0; c<lclusterdata.clustermembers.size(); c++){
		if(lclusterdata.clustersector[c] < 2){
			avgclustersize += lclusterdata.clustermembers[c].size();
			cnt++;
		}
	}
	lobs.avgclustersize = avgclustersize/(double)cnt;
}

void obsAverageClusterSizeFortunato(Observablestruct &lobs, Clusterstruct &lclusterdata){
	// Calculate average cluster size from Fortunato (1.7)
	// We consider only non-percolating clusters.
	// First calculate the number of clusters of size s per lattice site
	vector<int> sizes;
	vector<double> sizedist;
	int curcsize;
	int knownsize=0;
	for(unsigned int c=0; c<lclusterdata.clustermembers.size(); c++){
		if(lclusterdata.clustersector[c] < 2){
			curcsize = lclusterdata.clustermembers[c].size();
			knownsize=0;
			for(unsigned int s=0; s<sizes.size(); s++){
				if(sizes[s] == curcsize){
					knownsize = s;
					break;
				}
			}
			if(knownsize>0){
				sizedist[knownsize] += 1.0/(double)Nspace;
			}else{
				sizes.push_back(curcsize);
				sizedist.push_back(1.0/(double)Nspace);
			}
		}
	}

	double norm=0;
	for(unsigned int size=0; size<sizedist.size(); size++){
		norm += sizedist[size]*sizes[size];
	}

	double avgclustersizeF=0;
	for(unsigned int size=0; size<sizedist.size(); size++){
		avgclustersizeF += sizedist[size]*pow(sizes[size],2)/norm;
	}

	lobs.avgclustersizeF = avgclustersizeF;
}

void obsAverageClusterSizeNoPercc(Observablestruct &lobs, Clusterstruct &lclusterdata){
	// Calculate average cluster size without percolating clusters
	// First calculate the number of clusters of size s per lattice site
	vector<int> sizes;
	vector<double> sizedist;
	int curcsize;
	int knownsize=0;
	for(unsigned int c=0; c<lclusterdata.clustermembers.size(); c++){
		if(lclusterdata.clustersector[c] < 2 && lclusterdata.clusterispercolating[c] == 0){
			curcsize = lclusterdata.clustermembers[c].size();
			knownsize=0;
			for(unsigned int s=0; s<sizes.size(); s++){
				if(sizes[s] == curcsize){
					knownsize = s;
					break;
				}
			}
			if(knownsize>0){
				sizedist[knownsize] += 1.0/(double)Nspace;
			}else{
				sizes.push_back(curcsize);
				sizedist.push_back(1.0/(double)Nspace);
			}
		}
	}

	double norm=0;
	for(unsigned int size=0; size<sizedist.size(); size++){
		norm += sizedist[size]*sizes[size];
	}

	double avgclustersizenopercc=0;
	for(unsigned int size=0; size<sizedist.size(); size++){
		avgclustersizenopercc += sizedist[size]*pow(sizes[size],2)/norm;
	}

	lobs.avgclustersizenopercc = avgclustersizenopercc;
}

void obsCutPercentage(Observablestruct &lobs, Clusterstruct &lclusterdata){
	// Calculate cut
	double cut=0;
	for(unsigned int is=0; is<lclusterdata.isinsector.size(); is++){
		if(lclusterdata.isinsector[is] == 2)
			cut = cut + 1;
	}

	lobs.cut = cut/(double)Nspace;
}

void obsArea(Observablestruct &lobs, Clusterstruct &lclusterdata){
	lobs.area=0;
	int i1=0, i2=0, i3=0, is=0;
	bool incluster=false;
	// First direction i1
	for(i2=0;i2<Ns;i2++)
	for(i3=0;i3<Ns;i3++){
		i1=0;
		is=latmap(i1, i2, i3);
		if(lclusterdata.isincluster[is]==lobs.largestclusterid)
			incluster=true;
		for(i1=1;i1<Ns;i1++){
			is=latmap(i1, i2, i3);
			if(incluster==true && lclusterdata.isincluster[is] != lobs.largestclusterid){
				incluster=false;
				lobs.area++;
			}else if(incluster==false && lclusterdata.isincluster[is] == lobs.largestclusterid){
				incluster=true;
				lobs.area++;
			}
		}
	}
	
	// Second direction i2
	incluster=false;
	for(i1=0;i1<Ns;i1++)
	for(i3=0;i3<Ns;i3++){
		i2=0;
		is=latmap(i1, i2, i3);
		if(lclusterdata.isincluster[is]==lobs.largestclusterid)
			incluster=true;
		for(i2=1;i2<Ns;i2++){
			is=latmap(i1, i2, i3);
			if(incluster==true && lclusterdata.isincluster[is] != lobs.largestclusterid){
				incluster=false;
				lobs.area++;
			}else if(incluster==false && lclusterdata.isincluster[is] == lobs.largestclusterid){
				incluster=true;
				lobs.area++;
			}
		}
	}
	
	// Third direction i3
	incluster=false;
	for(i1=0;i1<Ns;i1++)
	for(i2=0;i2<Ns;i2++){
		i3=0;
		is=latmap(i1, i2, i3);
		if(lclusterdata.isincluster[is]==lobs.largestclusterid)
			incluster=true;
		for(i3=1;i3<Ns;i3++){
			is=latmap(i1, i2, i3);
			if(incluster==true && lclusterdata.isincluster[is] != lobs.largestclusterid){
				incluster=false;
				lobs.area++;
			}else if(incluster==false && lclusterdata.isincluster[is] == lobs.largestclusterid){
				incluster=true;
				lobs.area++;
			}
		}
	}
	lobs.area = lobs.area/(double)(6*Ns*Ns); // Mean over all 3 directions
}

void obsAreaLargestNonPercCluster(Observablestruct &lobs, Clusterstruct &lclusterdata){
	lobs.arealargestnonperccluster=0;
	int i1=0, i2=0, i3=0, is=0;
	bool incluster=false;
	
	int cluster=0;
	// Find largest non percollating cluster
	for(unsigned int i=0;i<lclusterdata.sortedrealcluster.size();i++){
		if(lclusterdata.clusterispercolating[lclusterdata.sortedrealcluster[i]] == 0){
			cluster=lclusterdata.sortedrealcluster[i];
			break;
		}
	}
	
	// First direction i1
	for(i2=0;i2<Ns;i2++)
	for(i3=0;i3<Ns;i3++){
		i1=0;
		is=latmap(i1, i2, i3);
		if(lclusterdata.isincluster[is]==cluster)
			incluster=true;
		for(i1=1;i1<Ns;i1++){
			is=latmap(i1, i2, i3);
			if(incluster==true && lclusterdata.isincluster[is] != cluster){
				incluster=false;
				lobs.arealargestnonperccluster++;
			}else if(incluster==false && lclusterdata.isincluster[is] == cluster){
				incluster=true;
				lobs.arealargestnonperccluster++;
			}
		}
	}
	
	// Second direction i2
	incluster=false;
	for(i1=0;i1<Ns;i1++)
	for(i3=0;i3<Ns;i3++){
		i2=0;
		is=latmap(i1, i2, i3);
		if(lclusterdata.isincluster[is]==cluster)
			incluster=true;
		for(i2=1;i2<Ns;i2++){
			is=latmap(i1, i2, i3);
			if(incluster==true && lclusterdata.isincluster[is] != cluster){
				incluster=false;
				lobs.arealargestnonperccluster++;
			}else if(incluster==false && lclusterdata.isincluster[is] == cluster){
				incluster=true;
				lobs.arealargestnonperccluster++;
			}
		}
	}
	
	// Third direction i3
	incluster=false;
	for(i1=0;i1<Ns;i1++)
	for(i2=0;i2<Ns;i2++){
		i3=0;
		is=latmap(i1, i2, i3);
		if(lclusterdata.isincluster[is]==cluster)
			incluster=true;
		for(i3=1;i3<Ns;i3++){
			is=latmap(i1, i2, i3);
			if(incluster==true && lclusterdata.isincluster[is] != cluster){
				incluster=false;
				lobs.arealargestnonperccluster++;
			}else if(incluster==false && lclusterdata.isincluster[is] == cluster){
				incluster=true;
				lobs.arealargestnonperccluster++;
			}
		}
	}
	lobs.arealargestnonperccluster = lobs.arealargestnonperccluster/(double)(6*Ns*Ns); // Mean over all 3 directions
}

void obsNumberOfPercClusters(Observablestruct &lobs, Clusterstruct &lclusterdata){
	// Number of percolating clusters
	int pcount=0;
	for(unsigned int p=0; p<lclusterdata.percolatingclusters.size(); p++){
		if(lclusterdata.percolatingdirections[p][0] == 1 && lclusterdata.percolatingdirections[p][1] == 1 && lclusterdata.percolatingdirections[p][2] == 1){
			pcount++;
		}
	}
	lobs.percc=pcount;
}

void obsPollAfterCut(Observablestruct &lobs, Clusterstruct &lclusterdata){
	/* Calculation of Polyakov loop expectation value for points which 
	 * survived the cut. For that loop over all points which are not in 
	 * sector 2 and calculate the spatial average of |P(x)| */
	complex<double> pollsum=0;
	int pollcnt=0;
	for(int is=0;is<Nspace;is++){
		if(lclusterdata.isinsector[is] < 2){
			pollsum += lclusterdata.poll[is];
			pollcnt++;
		}
	}
	lobs.poll = abs(pollsum)/(double)pollcnt;
}

void obsBoxesOnlyLargest(Observablestruct &lobs, Clusterstruct &lclusterdata){
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
				
				is = latmap(i1, i2, i3);

				if(lclusterdata.isincluster[is] == (int)c)
					clusterinbox=true;
			}

			if(clusterinbox==true)
				boxcnt++;

		} // Boxes in all directions
		lobs.numberofboxes[c][size]=boxcnt;
	} // Boxsize
	if(lobs.numberofboxes[c][0] != (int)lclusterdata.clustermembers[c].size()){
		cout << "ERROR: Number of boxes for boxsize = 1 has to be equal to number of cluster elements!" << endl;
	}
	if(lobs.numberofboxes[c][lobs.numberofboxes[c].size()-1] != 1){
		cout << "ERROR: Number of boxes for boxsize = Ns has to be equal 1!" << endl;
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
				
				is = latmap(i1, i2, i3);

				if(lclusterdata.isincluster[is] == (int)c)
					clusterinbox=true;
			}

			if(clusterinbox==true)
				boxcnt++;

		} // Boxes in all directions
		lobs.numberofboxes[c][size]=boxcnt;
	} // Boxsize
	if(lobs.numberofboxes[c][0] != (int)lclusterdata.clustermembers[c].size()){
		cout << "ERROR: Number of boxes for boxsize = 1 has to be equal to number of cluster elements!" << endl;
	}
	if(lobs.numberofboxes[c][lobs.numberofboxes[c].size()-1] != 1){
		cout << "ERROR: Number of boxes for boxsize = Ns has to be equal 1!" << endl;
	}
}

void obsBoxes(Observablestruct &lobs, Clusterstruct &lclusterdata){
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
					
					is = latmap(i1, i2, i3);

					if(lclusterdata.isincluster[is] == (int)c)
						clusterinbox=true;
				}

				if(clusterinbox==true)
					boxcnt++;

			} // Boxes in all directions
			lobs.numberofboxes[c][size]=boxcnt;
		} // Boxsize
		if(lobs.numberofboxes[c][0] != (int)lclusterdata.clustermembers[c].size()){
			cout << "ERROR: Number of boxes for boxsize = 1 has to be equal to number of cluster elements!" << endl;
		}
		if(lobs.numberofboxes[c][lobs.numberofboxes[c].size()-1] != 1){
			cout << "ERROR: Number of boxes for boxsize = Ns has to be equal 1!" << endl;
		}
	} // Cluster
}

void calcObservables(Observablestruct &lobs, Clusterstruct &lclusterdata){
	#ifdef DEBUG	
	cout << "Calculating observables... " << flush;
	#endif

	// lobs.maxclustersize, lobs.maxclusterid, lobs.maxclustersector
	obsLargestCluster(lobs, lclusterdata);
	// lobs.avgclustersize
	obsAverageClusterSize(lobs, lclusterdata);
	// lobs.avgclustersizeF
	obsAverageClusterSizeFortunato(lobs, lclusterdata);
	// lobs.avgclustersizenopercc
	obsAverageClusterSizeNoPercc(lobs, lclusterdata);
	// lobs.cut
	obsCutPercentage(lobs, lclusterdata);
	// lobs.largestclusterid
	lobs.largestclusterid=lclusterdata.sortedrealcluster[0];
	// lobs.percc
	obsNumberOfPercClusters(lobs, lclusterdata);
	// lobs.area
	obsArea(lobs, lclusterdata);
	// lobs.arealargestnonperccluster
	obsAreaLargestNonPercCluster(lobs, lclusterdata);
	// lobs.poll
	obsPollAfterCut(lobs, lclusterdata);
	
	// lobs.numberofboxes
	// if(doboxes)
	// obsBoxes(lobs, lclusterdata);
	// lobs.numberofboxes
	if(doboxes)
		obsBoxesOnlyLargest(lobs, lclusterdata);
	
	if(doradius)		
		clusterRadius(lobs, lclusterdata);
	
	#ifdef DEBUG
	cout << "done!" << endl;
	#endif
}

#endif // FINDCLUSTER_OBS_HPP
