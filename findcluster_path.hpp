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
        /* ii = startpoint;
        do {
          if( istagged[ii] == 0) {
            path1 = path1 + 1.0;
            istagged[ii] = 1;
          }
          ii--;
          if(ii == -1)
            ii = ii + leng3;
        } while ( lclusterdata.isincluster.at(latmap(i1, i2, ii)) == c && ii != startpoint);
        */
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
       /* ii = startpoint;
        do {
          if( istagged[ii] == 0) {
            path2 = path2 + 1.0;
            istagged[ii] = 1;
          }
          ii--;
          if(ii == -1)
            ii = ii + leng2;
        } while ( lclusterdata.isincluster.at(latmap(i1, ii, i3)) == c && ii != startpoint);
        */
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
        /*ii = startpoint;
        do {
          if( istagged[ii] == 0) {
            path3 = path3 + 1.0;
            istagged[ii] = 1;
          }
          ii--;
          if(ii == -1)
            ii = ii + leng1;
        } while ( lclusterdata.isincluster.at(latmap(ii, i2, i3)) == c && ii != startpoint);
        */
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

	for (int c=0; c<(int)lclusterdata.clustermembers.size(); c++) {
		if (lclusterdata.clustersector[c] < 2) {
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
  } // cluster if
  } // cluster loop
  
  lobs.largestclustermeanfreepath=lobs.meanfreepath[lobs.maxclusterid];

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
}

#endif // FINDCLUSTER_PATH_HPP
