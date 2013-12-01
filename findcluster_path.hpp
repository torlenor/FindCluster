#ifndef FINDCLUSTER_PATH_HPP
#define FINDCLUSTER_PATH_HPP
void obsClusterMeanFreePathLargest(Observablestruct &lobs, Clusterstruct &lclusterdata){
  // Calculates the mean free path of the largest cluster
	lobs.meanfreepath.resize(lclusterdata.clustermembers.size());

	int c = lobs.maxclusterid;
  // c=lclusterdata.clustermembers.size()-1;

  double path1=0, path2=0, path3=0;
  int paths=0;

  int i1=0, i2=0, i3=0, startpoint=0, ii=0;

  vector<int> istagged;
  istagged.resize(leng3);

  // First direction
  paths=0;
  path1=0;
  for (i1=0; i1<leng1; i1++) {
  for (i2=0; i2<leng2; i2++) {
    for( int ij=0; ij<leng3; ij++) {
      istagged[ij] = 0;
    }
    for (i3=0; i3<leng3; i3++) {
      // go through i3 until we hit a point in the cluster or i3=leng1-1

      if (lclusterdata.isincluster[latmap(i1, i2, i3)] == c && istagged[i3] == 0) {
        paths++;

        // Go in +3 direction
        startpoint = i3;
        ii = startpoint;
        do {
          if (istagged[ii] == 1) {
            cout << "WTF!" << endl;
          }
          path1 = path1 + 1.0;
          istagged[ii] = 1;
          ii++;
          if(ii == leng3)
            ii = ii - leng3;
        } while ( lclusterdata.isincluster[latmap(i1, i2, ii)] == c && ii != startpoint);
      
        // Go in -3 direction
        ii = startpoint;
        do {
          if( istagged[ii] == 0) {
            path1 = path1 + 1.0;
            istagged[ii] = 1;
          }
          ii--;
          if(ii == -1)
            ii = ii + leng3;
        } while ( lclusterdata.isincluster.at(latmap(i1, i2, ii)) == c && ii != startpoint);
        
      }
    }
  }
  }
  path1=path1/(double)paths;
  // END 1st direction
  
  // 2nd direction
  paths=0;
  for (i1=0; i1<leng1; i1++) {
  for (i3=0; i3<leng3; i3++) {
    for( int ij=0; ij<leng2; ij++) {
      istagged[ij] = 0;
    }
    for (i2=0; i2<leng2; i2++) {
      // go through i2 until we hit a point in the cluster or i3=leng1-1

      if (lclusterdata.isincluster[latmap(i1, i2, i3)] == c && istagged[i2] == 0) {
        paths++;

        // Go in +2 direction
        startpoint = i2;
        ii = startpoint;
        do {
          if (istagged[ii] == 1) {
            cout << "WTF!" << endl;
          }
          path2 = path2 + 1.0;
          istagged[ii] = 1;
          ii++;
          if(ii == leng2)
            ii = ii - leng2;
        } while ( lclusterdata.isincluster[latmap(i1, ii, i3)] == c && ii != startpoint);
      
        // Go in -2 direction
        ii = startpoint;
        do {
          if( istagged[ii] == 0) {
            path2 = path2 + 1.0;
            istagged[ii] = 1;
          }
          ii--;
          if(ii == -1)
            ii = ii + leng2;
        } while ( lclusterdata.isincluster.at(latmap(i1, ii, i3)) == c && ii != startpoint);
       
      }
    }
  }
  }
  path2=path2/(double)paths;
  // END 2nd direction
  
  // 3rd direction
  paths=0;
  for (i2=0; i2<leng2; i2++) {
  for (i3=0; i3<leng3; i3++) {
    for( int ij=0; ij<leng1; ij++) {
      istagged[ij] = 0;
    }
    for (i1=0; i1<leng1; i1++) {
      // go through i1 until we hit a point in the cluster or i1=leng1-1

      if (lclusterdata.isincluster[latmap(i1, i2, i3)] == c && istagged[i1] == 0) {
        paths++;

        // Go in +1 direction
        startpoint = i1;
        ii = startpoint;
        do {
          if (istagged[ii] == 1) {
            cout << "WTF!" << endl;
          }
          path3 = path3 + 1.0;
          istagged[ii] = 1;
          ii++;
          if(ii == leng1)
            ii = ii - leng1;
        } while ( lclusterdata.isincluster[latmap(ii, i2, i3)] == c && ii != startpoint);
      
        // Go in -1 direction
        ii = startpoint;
        do {
          if( istagged[ii] == 0) {
            path3 = path3 + 1.0;
            istagged[ii] = 1;
          }
          ii--;
          if(ii == -1)
            ii = ii + leng1;
        } while ( lclusterdata.isincluster.at(latmap(ii, i2, i3)) == c && ii != startpoint);
        
      }
    }
  }
  }
  path3=path3/(double)paths;
  // END 3nd direction


  lobs.meanfreepath[c]=(path1+path2+path3)/(double)3.0;
  lobs.largestclustermeanfreepath=(path1+path2+path3)/(double)3.0;
}

void obsClusterMeanFreePath(Observablestruct &lobs, Clusterstruct &lclusterdata){
  // Calculates the mean free path of all clusters
	lobs.meanfreepath.resize(lclusterdata.clustermembers.size());

  int totalcount=0;

	for (int c=0; c<(int)lclusterdata.clustermembers.size(); c++) {
		if (lclusterdata.clustersector[c] < 2) {
      totalcount=0;
      double path1=0, path2=0, path3=0;
      int paths=0;

      int i1=0, i2=0, i3=0, startpoint=0, ii=0;

      vector<int> istagged;
      istagged.resize(leng3);

      // First direction
      paths=0;
      path1=0;
      for (i1=0; i1<leng1; i1++) {
      for (i2=0; i2<leng2; i2++) {
        for( int ij=0; ij<leng3; ij++) {
          istagged[ij] = 0;
        }
        for (i3=0; i3<leng3; i3++) {
          // go through i3 until we hit a point in the cluster or i3=leng3-1

          if (lclusterdata.isincluster[latmap(i1, i2, i3)] == c && istagged[i3] == 0) {
            paths++;

            // Go in +3 direction
            startpoint = i3;
            ii = startpoint;
            do {
              if (istagged[ii] == 1) {
                cout << "WTF!" << endl;
              }
              path1 = path1 + 1.0;
              istagged[ii] = 1;
              totalcount++;
              ii++;
              if(ii == leng3)
                ii = ii - leng3;
            } while ( lclusterdata.isincluster[latmap(i1, i2, ii)] == c && ii != startpoint);
          
            // Go in -3 direction
            ii = startpoint;
            do {
              if( istagged[ii] == 0) {
                path1 = path1 + 1.0;
                istagged[ii] = 1;
                totalcount++;
              }
              ii--;
              if(ii == -1)
                ii = ii + leng3;
            } while ( lclusterdata.isincluster.at(latmap(i1, i2, ii)) == c && ii != startpoint);
            
          }
        }
      }
      }
      path1=path1/(double)paths;
      // END 1st direction
      if ( (totalcount - (int)lclusterdata.clustermembers[c].size()) != 0) {
        cout << "ERROR: Something wrong in path calculation, direction 3. Number of counted cluster members != clustermembers" << endl;
      }

      // 2nd direction
      paths=0;
      totalcount=0;
      for (i1=0; i1<leng1; i1++) {
      for (i3=0; i3<leng3; i3++) {
        for( int ij=0; ij<leng2; ij++) {
          istagged[ij] = 0;
        }
        for (i2=0; i2<leng2; i2++) {
          // go through i2 until we hit a point in the cluster or i2=leng2-1

          if (lclusterdata.isincluster[latmap(i1, i2, i3)] == c && istagged[i2] == 0) {
            paths++;

            // Go in +2 direction
            startpoint = i2;
            ii = startpoint;
            do {
              if (istagged[ii] == 1) {
                cout << "WTF!" << endl;
              }
              path2 = path2 + 1.0;
              istagged[ii] = 1;
              totalcount++;
              ii++;
              if(ii == leng2)
                ii = ii - leng2;
            } while ( lclusterdata.isincluster[latmap(i1, ii, i3)] == c && ii != startpoint);
          
            // Go in -2 direction
            ii = startpoint;
            do {
              if( istagged[ii] == 0) {
                path2 = path2 + 1.0;
                totalcount++;
                istagged[ii] = 1;
              }
              ii--;
              if(ii == -1)
                ii = ii + leng2;
            } while ( lclusterdata.isincluster.at(latmap(i1, ii, i3)) == c && ii != startpoint);
            
          }
        }
      }
      }
      path2=path2/(double)paths;
      // END 2nd direction
      if ( (totalcount - (int)lclusterdata.clustermembers[c].size()) != 0) {
        cout << "ERROR: Something wrong in path calculation, direction 2. Number of counted cluster members != clustermembers" << endl;
      }

      // 3rd direction
      paths=0;
      totalcount=0;
      for (i2=0; i2<leng2; i2++) {
      for (i3=0; i3<leng3; i3++) {
        for( int ij=0; ij<leng1; ij++) {
          istagged[ij] = 0;
        }
        for (i1=0; i1<leng1; i1++) {
          // go through i1 until we hit a point in the cluster or i1=leng1-1

          if (lclusterdata.isincluster[latmap(i1, i2, i3)] == c && istagged[i1] == 0) {
            paths++;

            // Go in +1 direction
            startpoint = i1;
            ii = startpoint;
            do {
              if (istagged[ii] == 1) {
                cout << "WTF!" << endl;
              }
              path3 = path3 + 1.0;
              istagged[ii] = 1;
              totalcount++;
              ii++;
              if(ii == leng1)
                ii = ii - leng1;
            } while ( lclusterdata.isincluster[latmap(ii, i2, i3)] == c && ii != startpoint);
          
            // Go in -1 direction
            ii = startpoint;
            do {
              if( istagged[ii] == 0) {
                path3 = path3 + 1.0;
                totalcount++;
                istagged[ii] = 1;
              }
              ii--;
              if(ii == -1)
                ii = ii + leng1;
            } while ( lclusterdata.isincluster.at(latmap(ii, i2, i3)) == c && ii != startpoint);
            
          }
        }
      }
      }
      path3=path3/(double)paths;
      // END 3nd direction
      
      if ( (totalcount - (int)lclusterdata.clustermembers[c].size()) != 0) {
        cout << "ERROR: Something wrong in path calculation, direction 1. Number of counted cluster members != clustermembers" << endl;
      }

      lobs.meanfreepath[c]=(path1+path2+path3)/(double)3.0;
  } // cluster if
  } // cluster loop
  
  lobs.largestclustermeanfreepath=lobs.meanfreepath[lobs.maxclusterid];
  lobs.largestnpclustermeanfreepath=lobs.meanfreepath[lobs.largestnonpercclusterid];

  // Calculate average over all clusters
  double avgclustermeanfreepath=0;
  int cnt=0;
  for(unsigned int c=0; c<lclusterdata.clustermembers.size(); c++){
    if(lclusterdata.clustersector[c] < 2){
      avgclustermeanfreepath += lobs.meanfreepath[c];
      cnt++;
    }    
  }
  lobs.avgclustermeanfreepath = avgclustermeanfreepath/(double)cnt;
  
  // Calculate average cluster size non percolating
  avgclustermeanfreepath=0;
  cnt=0;
  for(unsigned int c=0; c<lclusterdata.clustermembers.size(); c++){
    if(lclusterdata.clustersector[c] < 2 && lclusterdata.clusterispercolating[c] == 0){
      avgclustermeanfreepath += lobs.meanfreepath[c];
      cnt++;
    }    
  }
  lobs.avgnpclustermeanfreepath = avgclustermeanfreepath/(double)cnt;

  // Write a list of clusterweight and meanfreepath to stdout
  /* for(unsigned int c=0; c<lclusterdata.sortedrealcluster.size(); c++){
    if(lclusterdata.clustersector[lclusterdata.sortedrealcluster[c]] < 2){
      cout << lclusterdata.clustermembers[lclusterdata.sortedrealcluster[c]].size() << " " << lobs.meanfreepath[lclusterdata.sortedrealcluster[c]] << endl;
    }
  } */

}

void obsClusterMeanFreePathNew(Observablestruct &lobs, Clusterstruct &lclusterdata){
  // Calculates the mean free path of all clusters using the new definition
	lobs.meanfreepathnew.resize(lclusterdata.clustermembers.size());

  int totalcount=0;

	for (int c=0; c<(int)lclusterdata.clustermembers.size(); c++) {
		if (lclusterdata.clustersector[c] < 2) {
      totalcount=0;
      double path1=0, path2=0, path3=0;
      int paths=0;
      double segment=0;

      int i1=0, i2=0, i3=0, startpoint=0, ii=0;

      vector<int> istagged;
      istagged.resize(leng3);

      // First direction
      paths=0;
      path1=0;
      segment=0;
      for (i1=0; i1<leng1; i1++) {
      for (i2=0; i2<leng2; i2++) {
        for( int ij=0; ij<leng3; ij++) {
          istagged[ij] = 0;
        }
        for (i3=0; i3<leng3; i3++) {
          // go through i3 until we hit a point in the cluster or i3=leng3-1

          if (lclusterdata.isincluster[latmap(i1, i2, i3)] == c && istagged[i3] == 0) {
            paths++;
            segment=0;

            // Go in +3 direction
            startpoint = i3;
            ii = startpoint;
            do {
              if (istagged[ii] == 1) {
                cout << "WTF!" << endl;
              }
              segment += 1.0;
              // path1 = path1 + 1.0;
              istagged[ii] = 1;
              totalcount++;
              ii++;
              if(ii == leng3)
                ii = ii - leng3;
            } while ( lclusterdata.isincluster[latmap(i1, i2, ii)] == c && ii != startpoint);
          
            // Go in -3 direction
            ii = startpoint;
            do {
              if( istagged[ii] == 0) {
                segment += 1.0;
                // path1 = path1 + 1.0;
                istagged[ii] = 1;
                totalcount++;
              }
              ii--;
              if(ii == -1)
                ii = ii + leng3;
            } while ( lclusterdata.isincluster.at(latmap(i1, i2, ii)) == c && ii != startpoint);
            path1 += segment*segment;
          } // new segment found if end
        }
      }
      }
      path1=path1/((double)lclusterdata.clustermembers[c].size()*(double)2.0);
      // END 1st direction
      if ( (totalcount - (int)lclusterdata.clustermembers[c].size()) != 0) {
        cout << "ERROR: Something wrong in path calculation, direction 3. Number of counted cluster members != clustermembers" << endl;
      }

      // 2nd direction
      paths=0;
      path2=0;
      totalcount=0;
      segment=0;
      for (i1=0; i1<leng1; i1++) {
      for (i3=0; i3<leng3; i3++) {
        for( int ij=0; ij<leng2; ij++) {
          istagged[ij] = 0;
        }
        for (i2=0; i2<leng2; i2++) {
          // go through i2 until we hit a point in the cluster or i2=leng2-1

          if (lclusterdata.isincluster[latmap(i1, i2, i3)] == c && istagged[i2] == 0) {
            paths++;
            segment=0;

            // Go in +2 direction
            startpoint = i2;
            ii = startpoint;
            do {
              if (istagged[ii] == 1) {
                cout << "WTF!" << endl;
              }
              segment += 1.0;
              // path2 = path2 + 1.0;
              istagged[ii] = 1;
              totalcount++;
              ii++;
              if(ii == leng2)
                ii = ii - leng2;
            } while ( lclusterdata.isincluster[latmap(i1, ii, i3)] == c && ii != startpoint);
          
            // Go in -2 direction
            ii = startpoint;
            do {
              if( istagged[ii] == 0) {
                segment += 1.0;
                // path2 = path2 + 1.0;
                totalcount++;
                istagged[ii] = 1;
              }
              ii--;
              if(ii == -1)
                ii = ii + leng2;
            } while ( lclusterdata.isincluster.at(latmap(i1, ii, i3)) == c && ii != startpoint);
            path2 += segment*segment;
          } // new segment found if end
        }
      }
      }
      path2=path2/((double)lclusterdata.clustermembers[c].size()*(double)2.0);
      // END 2nd direction
      if ( (totalcount - (int)lclusterdata.clustermembers[c].size()) != 0) {
        cout << "ERROR: Something wrong in path calculation, direction 2. Number of counted cluster members != clustermembers" << endl;
      }

      // 3rd direction
      paths=0;
      totalcount=0;
      segment=0;
      for (i2=0; i2<leng2; i2++) {
      for (i3=0; i3<leng3; i3++) {
        for( int ij=0; ij<leng1; ij++) {
          istagged[ij] = 0;
        }
        for (i1=0; i1<leng1; i1++) {
          // go through i1 until we hit a point in the cluster or i1=leng1-1

          if (lclusterdata.isincluster[latmap(i1, i2, i3)] == c && istagged[i1] == 0) {
            paths++;
            segment=0;

            // Go in +1 direction
            startpoint = i1;
            ii = startpoint;
            do {
              if (istagged[ii] == 1) {
                cout << "WTF!" << endl;
              }
              segment += 1.0;
              // path3 = path3 + 1.0;
              istagged[ii] = 1;
              totalcount++;
              ii++;
              if(ii == leng1)
                ii = ii - leng1;
            } while ( lclusterdata.isincluster[latmap(ii, i2, i3)] == c && ii != startpoint);
          
            // Go in -1 direction
            ii = startpoint;
            do {
              if( istagged[ii] == 0) {
                segment += 1.0;
                // path3 = path3 + 1.0;
                totalcount++;
                istagged[ii] = 1;
              }
              ii--;
              if(ii == -1)
                ii = ii + leng1;
            } while ( lclusterdata.isincluster.at(latmap(ii, i2, i3)) == c && ii != startpoint);
            path3 += segment*segment;
          } // new segment found if end
        }
      }
      }
      path3=path3/((double)lclusterdata.clustermembers[c].size()*(double)2.0);
      // END 3nd direction
      
      if ( (totalcount - (int)lclusterdata.clustermembers[c].size()) != 0) {
        cout << "ERROR: Something wrong in path calculation, direction 1. Number of counted cluster members != clustermembers" << endl;
      }

      lobs.meanfreepathnew[c]=(path1+path2+path3)/(double)3.0;
  } // cluster if
  } // cluster loop
  
  lobs.largestclustermeanfreepathnew=lobs.meanfreepathnew[lobs.maxclusterid];
  lobs.largestnpclustermeanfreepathnew=lobs.meanfreepathnew[lobs.largestnonpercclusterid];
  
    // Write a list of clusterweight and meanfreepath to stdout
  /* for(unsigned int c=0; c<lclusterdata.sortedrealcluster.size(); c++){
    if(lclusterdata.clustersector[lclusterdata.sortedrealcluster[c]] < 2){
      cout << lclusterdata.clustermembers[lclusterdata.sortedrealcluster[c]].size() << " " << lobs.meanfreepathnew[lclusterdata.sortedrealcluster[c]] << endl;
    }
  } */
}

void obsAverageMeanfreepathNew(Observablestruct &lobs, Clusterstruct &lclusterdata){

  // Naive averages
  // Calculate average over all clusters
  double avgclustermeanfreepathnew=0;
  int cnt=0;
  for(unsigned int c=0; c<lclusterdata.clustermembers.size(); c++){
    if(lclusterdata.clustersector[c] < 2){
      avgclustermeanfreepathnew += lobs.meanfreepathnew[c];
      cnt++;
    }    
  }
  lobs.avgclustermeanfreepathnew = avgclustermeanfreepathnew/(double)cnt;
  // Calculate average cluster size non percolating
  avgclustermeanfreepathnew=0;
  cnt=0;
  for(unsigned int c=0; c<lclusterdata.clustermembers.size(); c++){
    if(lclusterdata.clustersector[c] < 2 && lclusterdata.clusterispercolating[c] == 0){
      avgclustermeanfreepathnew += lobs.meanfreepathnew[c];
      cnt++;
    }    
  }
  lobs.avgnpclustermeanfreepathnew = avgclustermeanfreepathnew/(double)cnt;

  // --------------------------------------------------------------------- //
  // Fortunato/Stauffer average over all clusters
  // --------------------------------------------------------------------- //
	// First calculate the number of clusters of size s per lattice site
	vector<int> sizes;
	vector<double> sizedist;
	vector<double> pathavg;
	vector<int> pathcnt;
	double avgFclustermeanfreepathnew=0;
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
				sizedist[knownsize] += 1.0/(double)Nspace;

				pathavg[knownsize] += lobs.meanfreepathnew.at(c);
				pathcnt[knownsize]++;
			}else{
				sizes.push_back(curcsize);
				sizedist.push_back(1.0/(double)Nspace);
				
				pathavg.push_back(lobs.meanfreepathnew.at(c));
				pathcnt.push_back(1);
			}
		}
	}

	// Norm of pathavg
	for(unsigned int size=0; size<sizedist.size(); size++){
		pathavg[size]=pathavg[size]/(double)pathcnt[size];
	}
  
  // Norm	
  double norm=0;
	for(unsigned int size=0; size<sizedist.size(); size++){
		norm += sizedist[size]*sizes[size];
	}


	for(unsigned int size=0; size<sizedist.size(); size++){
		avgFclustermeanfreepathnew += sizedist[size]*sizes[size]*pathavg.at(size);
	}

	lobs.avgFclustermeanfreepathnew = avgFclustermeanfreepathnew/(double)norm;
	
	// --------------------------------------------------------------------- //
  // Fortunato/Stauffer average over all clusters (non-perc.)
  // --------------------------------------------------------------------- //
	// First calculate the number of clusters of size s per lattice site
	sizes.resize(0);
	sizedist.resize(0);
	pathavg.resize(0);
	pathcnt.resize(0);
	double avgFnpclustermeanfreepathnew=0;
	curcsize=0;
	knownsize=0;

	for(unsigned int c=0; c<lclusterdata.clustermembers.size(); c++){
		if(lclusterdata.clustersector[c] < 2 && lclusterdata.clustermembers[c].size() > 1  && lclusterdata.clusterispercolating[c] == 0){
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

				pathavg[knownsize] += lobs.meanfreepathnew.at(c);
				pathcnt[knownsize]++;
			}else{
				sizes.push_back(curcsize);
				sizedist.push_back(1.0/(double)Nspace);
				
				pathavg.push_back(lobs.meanfreepathnew.at(c));
				pathcnt.push_back(1);
			}
		}
	}

	// Norm of pathavg
	for(unsigned int size=0; size<sizedist.size(); size++){
		pathavg[size]=pathavg[size]/(double)pathcnt[size];
	}
	
	// Norm	
  norm=0;
	for(unsigned int size=0; size<sizedist.size(); size++){
		norm += sizedist[size]*sizes[size];
	}

	for(unsigned int size=0; size<sizedist.size(); size++){
		avgFnpclustermeanfreepathnew += sizedist[size]*sizes[size]*pathavg.at(size);
	}

	lobs.avgFnpclustermeanfreepathnew = avgFnpclustermeanfreepathnew/(double)norm;
}


#endif // FINDCLUSTER_PATH_HPP
