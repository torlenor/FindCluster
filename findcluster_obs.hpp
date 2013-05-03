#ifndef FINDCLUSTER_OBS_HPP
#define FINDCLUSTER_OBS_HPP

void obsLargestCluster(Observablestruct &lobs, Clusterstruct &lclusterdata){
	// Find largest cluster
	int size=-1, largestcluster=-1;
	for(unsigned int c=0; c<lclusterdata.clustermembers.size(); c++){
		if( (int)lclusterdata.clustermembers[c].size() > size && lclusterdata.clustersector[c] < 2){
			size = lclusterdata.clustermembers[c].size();
			largestcluster = c;
		}
	}
	
	lobs.maxclustersize = size;
	lobs.maxclusterid = largestcluster;
	lobs.maxclustersector = lclusterdata.clustersector[largestcluster];
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
	lobs.laserdim=0;
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
				lobs.laserdim++;
			}else if(incluster==false && lclusterdata.isincluster[is] == lobs.largestclusterid){
				incluster=true;
				lobs.laserdim++;
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
				lobs.laserdim++;
			}else if(incluster==false && lclusterdata.isincluster[is] == lobs.largestclusterid){
				incluster=true;
				lobs.laserdim++;
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
				lobs.laserdim++;
			}else if(incluster==false && lclusterdata.isincluster[is] == lobs.largestclusterid){
				incluster=true;
				lobs.laserdim++;
			}
		}
	}
	lobs.laserdim = lobs.laserdim/(double)(6*Ns*Ns); // Mean over all 3 directions
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

void calcObservables(Observablestruct &lobs, Clusterstruct &lclusterdata){
	#ifdef DEBUG	
	cout << "Calculating observables... " << flush;
	#endif

	obsLargestCluster(lobs, lclusterdata);
	obsAverageClusterSize(lobs, lclusterdata);
	obsAverageClusterSizeFortunato(lobs, lclusterdata);
	obsCutPercentage(lobs, lclusterdata);
	lobs.largestclusterid=lclusterdata.sortedrealcluster[0];
	obsNumberOfPercClusters(lobs, lclusterdata);
	obsArea(lobs, lclusterdata);
	obsPollAfterCut(lobs, lclusterdata);
	
	#ifdef DEBUG
	cout << "done!" << endl;
	#endif

}

#endif // FINDCLUSTER_OBS_HPP
