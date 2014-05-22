#ifndef FINDCLUSTER_OBS_HPP
#define FINDCLUSTER_OBS_HPP

#include "findcluster_path.hpp"
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
	
	// Calculate average cluster size non percolating
	avgclustersize=0;
	cnt=0;
	for(unsigned int c=0; c<lclusterdata.clustermembers.size(); c++){
		if(lclusterdata.clustersector[c] < 2 && lclusterdata.clusterispercolating[c] == 0){
			avgclustersize += lclusterdata.clustermembers[c].size();
			cnt++;
		}
	}
	lobs.avgclustersizenp = avgclustersize/(double)cnt;
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
				sizedist[knownsize] += 1.0/(double)opt.Nspace;
			}else{
				sizes.push_back(curcsize);
				sizedist.push_back(1.0/(double)opt.Nspace);
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
				sizedist[knownsize] += 1.0/(double)opt.Nspace;
			}else{
				sizes.push_back(curcsize);
				sizedist.push_back(1.0/(double)opt.Nspace);
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

	lobs.avgclustersizeFnp = avgclustersizenopercc;
}

void obsCutPercentage(Observablestruct &lobs, Clusterstruct &lclusterdata){
	// Calculate cut
	double cut=0;
	for(unsigned int is=0; is<lclusterdata.isinsector.size(); is++){
		if(lclusterdata.isinsector[is] == 2)
			cut = cut + 1;
	}

	lobs.cut = cut/(double)opt.Nspace;
}

void obsArea(Observablestruct &lobs, Clusterstruct &lclusterdata){
	lobs.area=0;
	int i1=0, i2=0, i3=0, is=0;
	bool incluster=false;
	// First direction i1
	for(i2=0;i2<opt.Ns;i2++)
	for(i3=0;i3<opt.Ns;i3++){
		i1=0;
		is=latmap(i1, i2, i3);
		if(lclusterdata.isincluster[is]==lobs.largestclusterid)
			incluster=true;
		for(i1=1;i1<opt.Ns;i1++){
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
	for(i1=0;i1<opt.Ns;i1++)
	for(i3=0;i3<opt.Ns;i3++){
		i2=0;
		is=latmap(i1, i2, i3);
		if(lclusterdata.isincluster[is]==lobs.largestclusterid)
			incluster=true;
		for(i2=1;i2<opt.Ns;i2++){
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
	for(i1=0;i1<opt.Ns;i1++)
	for(i2=0;i2<opt.Ns;i2++){
		i3=0;
		is=latmap(i1, i2, i3);
		if(lclusterdata.isincluster[is]==lobs.largestclusterid)
			incluster=true;
		for(i3=1;i3<opt.Ns;i3++){
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
	lobs.area = lobs.area/(double)(6*opt.Ns*opt.Ns); // Mean over all 3 directions
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
	for(i2=0;i2<opt.Ns;i2++)
	for(i3=0;i3<opt.Ns;i3++){
		i1=0;
		is=latmap(i1, i2, i3);
		if(lclusterdata.isincluster[is]==cluster)
			incluster=true;
		for(i1=1;i1<opt.Ns;i1++){
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
	for(i1=0;i1<opt.Ns;i1++)
	for(i3=0;i3<opt.Ns;i3++){
		i2=0;
		is=latmap(i1, i2, i3);
		if(lclusterdata.isincluster[is]==cluster)
			incluster=true;
		for(i2=1;i2<opt.Ns;i2++){
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
	for(i1=0;i1<opt.Ns;i1++)
	for(i2=0;i2<opt.Ns;i2++){
		i3=0;
		is=latmap(i1, i2, i3);
		if(lclusterdata.isincluster[is]==cluster)
			incluster=true;
		for(i3=1;i3<opt.Ns;i3++){
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
	lobs.arealargestnonperccluster = lobs.arealargestnonperccluster/(double)(6*opt.Ns*opt.Ns); // Mean over all 3 directions
}

void obsAreaAvgNonPercCluster(Observablestruct &lobs, Clusterstruct &lclusterdata){
	lobs.areaavgnonperccluster=0;
	int i1=0, i2=0, i3=0, is=0;
	bool incluster=false;
	
	int clusters=0;
	
  for(int cluster=0; cluster<(int)lclusterdata.clustermembers.size(); cluster++){
		if(lclusterdata.clustersector[cluster] < 2 && lclusterdata.clusterispercolating[cluster] == 0){
      clusters++;
      incluster=false;
    
      // First direction i1
      for(i2=0;i2<opt.Ns;i2++)
      for(i3=0;i3<opt.Ns;i3++){
        i1=0;
        is=latmap(i1, i2, i3);
        if(lclusterdata.isincluster[is]==cluster)
          incluster=true;
        for(i1=1;i1<opt.Ns;i1++){
          is=latmap(i1, i2, i3);
          if(incluster==true && lclusterdata.isincluster[is] != cluster){
            incluster=false;
            lobs.areaavgnonperccluster++;
          }else if(incluster==false && lclusterdata.isincluster[is] == cluster){
            incluster=true;
            lobs.areaavgnonperccluster++;
          }
        }
      }
      
      // Second direction i2
      incluster=false;
      for(i1=0;i1<opt.Ns;i1++)
      for(i3=0;i3<opt.Ns;i3++){
        i2=0;
        is=latmap(i1, i2, i3);
        if(lclusterdata.isincluster[is]==cluster)
          incluster=true;
        for(i2=1;i2<opt.Ns;i2++){
          is=latmap(i1, i2, i3);
          if(incluster==true && lclusterdata.isincluster[is] != cluster){
            incluster=false;
            lobs.areaavgnonperccluster++;
          }else if(incluster==false && lclusterdata.isincluster[is] == cluster){
            incluster=true;
            lobs.areaavgnonperccluster++;
          }
        }
      }
      
      // Third direction i3
      incluster=false;
      for(i1=0;i1<opt.Ns;i1++)
      for(i2=0;i2<opt.Ns;i2++){
        i3=0;
        is=latmap(i1, i2, i3);
        if(lclusterdata.isincluster[is]==cluster)
          incluster=true;
        for(i3=1;i3<opt.Ns;i3++){
          is=latmap(i1, i2, i3);
          if(incluster==true && lclusterdata.isincluster[is] != cluster){
            incluster=false;
            lobs.areaavgnonperccluster++;
          }else if(incluster==false && lclusterdata.isincluster[is] == cluster){
            incluster=true;
            lobs.areaavgnonperccluster++;
          }
        }
      }
    }
  }
    
  lobs.areaavgnonperccluster = lobs.areaavgnonperccluster/(double)(6*opt.Ns*opt.Ns*clusters); // Mean over all 3 directions
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
	for(int is=0;is<opt.Nspace;is++){
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

void obsRootMeanSquareDistance(Observablestruct &lobs, Clusterstruct &lclusterdata){
	// THIS NEEDS THE CLUSTER RADIUS OF ALL CLUSTERS!!!
	// Calculate root mean square distance traveled R
	// Stauffer pp.117 eq. 105
	
	// First calculate the number of clusters of size s per lattice site
	vector<int> sizes;
	vector<double> sizedist;
	vector<double> radiusavg, radiuscnt;
	double meansquaredistanceR=0;
	int curcsize;
	int knownsize=0;

	for(unsigned int c=0; c<lclusterdata.clustermembers.size(); c++){
		if(lclusterdata.clustersector[c] < 2 && lclusterdata.clustermembers[c].size() > 1){
			curcsize = lclusterdata.clustermembers[c].size();
			knownsize=0;
			for(unsigned int s=0; s<sizes.size(); s++){
				if(sizes[s] == curcsize){
					knownsize = s;
					break;
				}
			}
			if(knownsize>0){
				sizedist[knownsize] += 1.0/(double)opt.Nspace;

				radiusavg[knownsize] += lobs.clusterradius.at(c);
				radiuscnt[knownsize]++;
			}else{
				sizes.push_back(curcsize);
				sizedist.push_back(1.0/(double)opt.Nspace);
				
				radiusavg.push_back(lobs.clusterradius.at(c));
				radiuscnt.push_back(1);
			}
		}
	}

	// Norm of raidusavg
	for(unsigned int size=0; size<sizedist.size(); size++){
		radiusavg[size]=radiusavg[size]/(double)radiuscnt[size];
	}

	for(unsigned int size=0; size<sizedist.size(); size++){
		meansquaredistanceR += sizedist[size]*sizes[size]*pow(radiusavg.at(size),2.0);
	}

	lobs.rootmeansquaredistanceR = sqrt(meansquaredistanceR);
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
	// lobs.areaavgnonperccluster
	obsAreaAvgNonPercCluster(lobs, lclusterdata);
	// lobs.poll
	obsPollAfterCut(lobs, lclusterdata);

	
	// lobs.numberofboxes
	// if(doboxes)
	// obsBoxes(lobs, lclusterdata);
	// lobs.numberofboxes
	if(opt.doboxes) {
    cout << "b" << flush;
		obsBoxesOnlyLargest(lobs, lclusterdata);
  }
	
	if(opt.doradius){
		if(opt.dodistance){
			obsClusterRadius(lobs, lclusterdata);
		}else{
			obsClusterRadiusOnlyLargest(lobs, lclusterdata);
      obsClusterRadiusOnlyLargestNP(lobs, lclusterdata);
		}
	}

	if(opt.domean){
    obsClusterMeanFreePathLargest(lobs, lclusterdata);
    obsClusterMeanFreePath(lobs, lclusterdata);
    obsClusterMeanFreePathNew(lobs, lclusterdata);
    obsAverageMeanfreepathNew(lobs, lclusterdata);
  }
	
	// lobs.rootmeansquaredistanceR
	if(opt.dodistance){
		obsRootMeanSquareDistance(lobs, lclusterdata);
	}

	#ifdef DEBUG
	cout << "done!" << endl;
	#endif
}

void calcExp(){
	double ddata[opt.nmeas];

	// Start of Jackknife
	// Maximal cluster size per volume expectation value
	for(int n=0;n<opt.nmeas;n++){
		ddata[n]=0;
		for(int j=0;j<opt.nmeas;j++){
			if(n!=j)
				ddata[n] += (&obs[j])->maxclustersize/(double)opt.Nspace;
		}
		ddata[n] = ddata[n]/(double)(opt.nmeas-1);
	}
	Jackknife(ddata, results.maxclustersize, results.maxclustersizeerr, opt.nmeas);
	
	// Maximal non percolating cluster size per volume expectation value
	for(int n=0;n<opt.nmeas;n++){
		ddata[n]=0;
		for(int j=0;j<opt.nmeas;j++){
			if(n!=j)
				ddata[n] += (&obs[j])->largestnonpercclustersize/(double)opt.Nspace;
		}
		ddata[n] = ddata[n]/(double)(opt.nmeas-1);
	}
	Jackknife(ddata, results.maxnonpercclustersize, results.maxnonpercclustersizeerr, opt.nmeas);
	
	// Average cluster size expectation value
	for(int n=0;n<opt.nmeas;n++){
		ddata[n]=0;
		for(int j=0;j<opt.nmeas;j++){
			if(n!=j)
				ddata[n] += (&obs[j])->avgclustersize;
		}
		ddata[n] = ddata[n]/(double)(opt.nmeas-1);
	}
	Jackknife(ddata, results.avgclustersize, results.avgclustersizeerr, opt.nmeas);
	
	// Average cluster size expectation value (Fortunato (1.7))
	for(int n=0;n<opt.nmeas;n++){
		ddata[n]=0;
		for(int j=0;j<opt.nmeas;j++){
			if(n!=j)
				ddata[n] += (&obs[j])->avgclustersizeF;
		}
		ddata[n] = ddata[n]/(double)(opt.nmeas-1);
	}
	Jackknife(ddata, results.avgclusersizeFortunato, results.avgclusersizeFortunatoerr, opt.nmeas);
	
	// Average cluster size expectation value without percolating clusters Fortunato
	for(int n=0;n<opt.nmeas;n++){
		ddata[n]=0;
		for(int j=0;j<opt.nmeas;j++){
			if(n!=j)
				ddata[n] += (&obs[j])->avgclustersizeFnp;
		}
		ddata[n] = ddata[n]/(double)(opt.nmeas-1);
	}
	Jackknife(ddata, results.avgclustersizeFnp, results.avgclustersizeFnperr, opt.nmeas);
	
	// Average cluster size expectation value without percolating clusters
	for(int n=0;n<opt.nmeas;n++){
		ddata[n]=0;
		for(int j=0;j<opt.nmeas;j++){
			if(n!=j)
				ddata[n] += (&obs[j])->avgclustersizenp;
		}
		ddata[n] = ddata[n]/(double)(opt.nmeas-1);
	}
	Jackknife(ddata, results.avgclustersizenp, results.avgclustersizenperr, opt.nmeas);
	
	if(opt.dodistance){
		// Average root mean square distance traveled R
		for(int n=0;n<opt.nmeas;n++){
			ddata[n]=0;
			for(int j=0;j<opt.nmeas;j++){
				if(n!=j)
					ddata[n] += (&obs[j])->rootmeansquaredistanceR;
			}
			ddata[n] = ddata[n]/(double)(opt.nmeas-1);
		}
		Jackknife(ddata, results.avgrootmeansquaredistance, results.avgrootmeansquaredistanceerr, opt.nmeas);
	}
	
	// Cut expectation value
	for(int n=0;n<opt.nmeas;n++){
		ddata[n]=0;
		for(int j=0;j<opt.nmeas;j++){
			if(n!=j)
				ddata[n] += (&obs[j])->cut;
		}
		ddata[n] = ddata[n]/(double)(opt.nmeas-1);
	}
	Jackknife(ddata, results.cut, results.cuterr, opt.nmeas);
	
	// Perimeter expectation value
	for(int n=0;n<opt.nmeas;n++){
		ddata[n]=0;
		for(int j=0;j<opt.nmeas;j++){
			if(n!=j)
				ddata[n] += (&obs[j])->area;
		}
		ddata[n] = ddata[n]/(double)(opt.nmeas-1);
	}
	Jackknife(ddata, results.totalperimeter, results.totalperimetererr, opt.nmeas);
	
	// Perimeter largest non percolating cluster expectation value
	for(int n=0;n<opt.nmeas;n++){
		ddata[n]=0;
		for(int j=0;j<opt.nmeas;j++){
			if(n!=j)
				ddata[n] += (&obs[j])->arealargestnonperccluster;
		}
		ddata[n] = ddata[n]/(double)(opt.nmeas-1);
	}
	Jackknife(ddata, results.largestnonpercperimeter, results.largestnonpercperimetererr, opt.nmeas);
	
  // Perimeter avg non percolating cluster expectation value
	for(int n=0;n<opt.nmeas;n++){
		ddata[n]=0;
		for(int j=0;j<opt.nmeas;j++){
			if(n!=j)
				ddata[n] += (&obs[j])->areaavgnonperccluster;
		}
		ddata[n] = ddata[n]/(double)(opt.nmeas-1);
	}
	Jackknife(ddata, results.avgnonpercperimeter, results.avgnonpercperimetererr, opt.nmeas);
	
	// Radius of largest cluster expectation value
	for(int n=0;n<opt.nmeas;n++){
		ddata[n]=0;
		for(int j=0;j<opt.nmeas;j++){
			if(n!=j)
				ddata[n] += (&obs[j])->largestclusterradius;
		}
		ddata[n] = ddata[n]/(double)(opt.nmeas-1);
	}
	Jackknife(ddata, results.largestclusterradius, results.largestclusterradiuserr, opt.nmeas);
	
	// Radius of largest non-perc. cluster expectation value
	for(int n=0;n<opt.nmeas;n++){
		ddata[n]=0;
		for(int j=0;j<opt.nmeas;j++){
			if(n!=j)
				ddata[n] += (&obs[j])->largestnpclusterradius;
		}
		ddata[n] = ddata[n]/(double)(opt.nmeas-1);
	}
	Jackknife(ddata, results.largestnpclusterradius, results.largestnpclusterradiuserr, opt.nmeas);
	
	// Avg. radius of cluster
	for(int n=0;n<opt.nmeas;n++){
		ddata[n]=0;
		for(int j=0;j<opt.nmeas;j++){
			if(n!=j)
				ddata[n] += (&obs[j])->avgclusterradius;
		}
		ddata[n] = ddata[n]/(double)(opt.nmeas-1);
	}
	Jackknife(ddata, results.avgclusterradius, results.avgclusterradiuserr, opt.nmeas);
	
	// Avg. radius of non-perc. clusters
	for(int n=0;n<opt.nmeas;n++){
		ddata[n]=0;
		for(int j=0;j<opt.nmeas;j++){
			if(n!=j)
				ddata[n] += (&obs[j])->avgnpclusterradius;
		}
		ddata[n] = ddata[n]/(double)(opt.nmeas-1);
	}
	Jackknife(ddata, results.avgnpclusterradius, results.avgnpclusterradiuserr, opt.nmeas);

	// Polyakov loop expectation value only for points with sector < 2
	for(int n=0;n<opt.nmeas;n++){
		ddata[n]=0;
		for(int j=0;j<opt.nmeas;j++){
			if(n!=j)
				ddata[n] += (&obs[j])->poll;
		}
		ddata[n] = ddata[n]/(double)(opt.nmeas-1);
	}
	Jackknife(ddata, results.polyakovloopaftercut, results.polyakovloopaftercuterr, opt.nmeas);
	
	// Number of percolating clusters expectation value
	for(int n=0;n<opt.nmeas;n++){
		ddata[n]=0;
		for(int j=0;j<opt.nmeas;j++){
			if(n!=j)
				ddata[n] += (&obs[j])->percc;
		}
		ddata[n] = ddata[n]/(double)(opt.nmeas-1);
	}
	Jackknife(ddata, results.avgperccluster, results.avgpercclustererr, opt.nmeas);
	
	if(opt.domean){
    // Mean free path of largest cluster expectation value
    for(int n=0;n<opt.nmeas;n++){
      ddata[n]=0;
      for(int j=0;j<opt.nmeas;j++){
        if(n!=j)
          ddata[n] += (&obs[j])->largestclustermeanfreepath;
      }
      ddata[n] = ddata[n]/(double)(opt.nmeas-1);
    }
    Jackknife(ddata, results.largestclustermeanfreepath, results.largestclustermeanfreepatherr, opt.nmeas);
    
    // Mean free path of largest non-perc. cluster expectation value
    for(int n=0;n<opt.nmeas;n++){
      ddata[n]=0;
      for(int j=0;j<opt.nmeas;j++){
        if(n!=j)
          ddata[n] += (&obs[j])->largestnpclustermeanfreepath;
      }
      ddata[n] = ddata[n]/(double)(opt.nmeas-1);
    }
    Jackknife(ddata, results.largestnpclustermeanfreepath, results.largestnpclustermeanfreepatherr, opt.nmeas);
    
    // Mean free path of avg cluster expectation value
    for(int n=0;n<opt.nmeas;n++){
      ddata[n]=0;
      for(int j=0;j<opt.nmeas;j++){
        if(n!=j)
          ddata[n] += (&obs[j])->avgclustermeanfreepath;
      }
      ddata[n] = ddata[n]/(double)(opt.nmeas-1);
    }
    Jackknife(ddata, results.avgclustermeanfreepath, results.avgclustermeanfreepatherr, opt.nmeas);
    
    // Mean free path of avg non-perc. cluster expectation value
    for(int n=0;n<opt.nmeas;n++){
      ddata[n]=0;
      for(int j=0;j<opt.nmeas;j++){
        if(n!=j)
          ddata[n] += (&obs[j])->avgnpclustermeanfreepath;
      }
      ddata[n] = ddata[n]/(double)(opt.nmeas-1);
    }
    Jackknife(ddata, results.avgnpclustermeanfreepath, results.avgnpclustermeanfreepatherr, opt.nmeas);

    // NEW DEFINITION OF MEAN FREE PATH
    // Mean free path of largest cluster expectation value (NEW)
    for(int n=0;n<opt.nmeas;n++){
      ddata[n]=0;
      for(int j=0;j<opt.nmeas;j++){
        if(n!=j)
          ddata[n] += (&obs[j])->largestclustermeanfreepathnew;
      }
      ddata[n] = ddata[n]/(double)(opt.nmeas-1);
    }
    Jackknife(ddata, results.largestclustermeanfreepathnew, results.largestclustermeanfreepathnewerr, opt.nmeas);
    
    // Mean free path of largest non-perc. cluster expectation value (NEW)
    for(int n=0;n<opt.nmeas;n++){
      ddata[n]=0;
      for(int j=0;j<opt.nmeas;j++){
        if(n!=j)
          ddata[n] += (&obs[j])->largestnpclustermeanfreepathnew;
      }
      ddata[n] = ddata[n]/(double)(opt.nmeas-1);
    }
    Jackknife(ddata, results.largestnpclustermeanfreepathnew, results.largestnpclustermeanfreepathnewerr, opt.nmeas);
      
    // Mean free path of avg cluster expectation value
    for(int n=0;n<opt.nmeas;n++){
      ddata[n]=0;
      for(int j=0;j<opt.nmeas;j++){
        if(n!=j)
          ddata[n] += (&obs[j])->avgclustermeanfreepathnew;
      }
      ddata[n] = ddata[n]/(double)(opt.nmeas-1);
    }
    Jackknife(ddata, results.avgclustermeanfreepathnew, results.avgclustermeanfreepathnewerr, opt.nmeas);
    
    // Mean free path of avg non-perc. cluster expectation value
    for(int n=0;n<opt.nmeas;n++){
      ddata[n]=0;
      for(int j=0;j<opt.nmeas;j++){
        if(n!=j)
          ddata[n] += (&obs[j])->avgnpclustermeanfreepathnew;
      }
      ddata[n] = ddata[n]/(double)(opt.nmeas-1);
    }
    Jackknife(ddata, results.avgnpclustermeanfreepathnew, results.avgnpclustermeanfreepathnewerr, opt.nmeas);
    
    // Mean free path of avg cluster expectation value Fortunato/Stauffer
    for(int n=0;n<opt.nmeas;n++){
      ddata[n]=0;
      for(int j=0;j<opt.nmeas;j++){
        if(n!=j)
          ddata[n] += (&obs[j])->avgFclustermeanfreepathnew;
      }
      ddata[n] = ddata[n]/(double)(opt.nmeas-1);
    }
    Jackknife(ddata, results.avgFclustermeanfreepathnew, results.avgFclustermeanfreepathnewerr, opt.nmeas);
    
    // Mean free path of avg non-perc. cluster expectation value Fortunato/Stauffer
    for(int n=0;n<opt.nmeas;n++){
      ddata[n]=0;
      for(int j=0;j<opt.nmeas;j++){
        if(n!=j)
          ddata[n] += (&obs[j])->avgFnpclustermeanfreepathnew;
      }
      ddata[n] = ddata[n]/(double)(opt.nmeas-1);
    }
    Jackknife(ddata, results.avgFnpclustermeanfreepathnew, results.avgFnpclustermeanfreepathnewerr, opt.nmeas);
    // END NEW DEFINITION OF MEAN FREE PATH
  } // END IF(domean)

	if(opt.doboxes){
		// Box counts for largest cluster which is not in sector 2
		vector<double> avgboxcnt, avgboxcnterr;
		avgboxcnt.resize(boxsize.size());
		avgboxcnterr.resize(boxsize.size());
		
		results.largestclusterboxcount.resize(boxsize.size());
		results.largestclusterboxcounterr.resize(boxsize.size());
		
		results.largestnonpercboxcount.resize(boxsize.size());
		results.largestnonpercboxcounterr.resize(boxsize.size());
		
		for(unsigned int size=0; size<boxsize.size(); size++){
			for(int n=0;n<opt.nmeas;n++){
				ddata[n]=0;

				for(int j=0;j<opt.nmeas;j++){
					if(n!=j)
					ddata[n] += (&obs[j])->numberofboxes[(&obs[j])->largestclusterid][size];
				}
				ddata[n] = ddata[n]/(double)(opt.nmeas-1);
			}
			Jackknife(ddata, avgboxcnt[size], avgboxcnterr[size], opt.nmeas);
			results.largestclusterboxcount[size]=avgboxcnt[size];
			results.largestclusterboxcounterr[size]=avgboxcnterr[size];
		}
		
		// Box counts for largest cluster which is not percolating
		for(unsigned int size=0; size<boxsize.size(); size++){
			avgboxcnt[size]=0;
			avgboxcnterr[size]=0;
			for(int n=0;n<opt.nmeas;n++){
				ddata[n]=0;

				for(int j=0;j<opt.nmeas;j++){
					if(n!=j)
					ddata[n] += (&obs[j])->numberofboxes[(&obs[j])->largestnonpercclusterid][size];
				}
				ddata[n] = ddata[n]/(double)(opt.nmeas-1);
			}
			Jackknife(ddata, avgboxcnt[size], avgboxcnterr[size], opt.nmeas);
			results.largestnonpercboxcount[size]=avgboxcnt[size];
			results.largestnonpercboxcounterr[size]=avgboxcnterr[size];
		}
	}
}

#endif // FINDCLUSTER_OBS_HPP
