#ifndef FINDCLUSTER_PATH_HPP
#define FINDCLUSTER_PATH_HPP
void obsClusterMeanFreePathLargest(Observablestruct &lobs, Clusterstruct &lclusterdata){
  // Calculates the mean free path of the largest cluster
	lobs.meanfreepath.resize(lclusterdata.clustermembers.size());

	int c = lobs.maxclusterid;
  // c=lclusterdata.clustermembers.size()-1;

  double path1=0, path2=0, path3=0;
  int paths=0;

  int i1=0, i2=0, i3=0, startpoint=0;

  vector<int> istagged3;
  istagged3.resize(leng3);

  // First direction
  paths=0;
  for (i1=0; i1<leng1; i1++) {
  for (i2=0; i2<leng2; i2++) {
    i3=0;
    // go through i3 until we hit a point in the cluster or i3=leng1-1
    for( int ii=0; ii<leng3; ii++) {
      istagged3[ii] = 0;
    }

    while (i3<leng1 && lclusterdata.isincluster[latmap(i1, i2, i3)] != c) {
      i3++;
    }

    if(i3 >= leng3)
      continue;
    if(istagged3[i3] == 1) {
      continue;
    }

    // go through the cluster part and count the paths and path++ onces until we hit a point not
    // in the cluster anymore or we hit the start point (percolating over this direction)
    paths++;
    startpoint = i3;
    // Go in +3 direction
    do {
      path1 = path1 + 1.0;
      istagged3[i3] = 1;
      i3++;
      if(i3 == leng3)
        i3 = i3 - leng3;
    } while ( lclusterdata.isincluster[latmap(i1, i2, i3)] == c && i3 != startpoint);
    
    // Go in -3 direction
    i3 = startpoint;
    do {
      if( istagged3[i3] == 0) {
        path1 = path1 + 1.0;
        istagged3[i3] = 1;
      }
      i3--;
      if(i3 == -1)
        i3 = i3 + leng3;
    } while ( lclusterdata.isincluster.at(latmap(i1, i2, i3)) == c && i3 != startpoint);
    
    i3 = startpoint + 1;
  }
  }
  path1=path1/(double)paths;
  // END 1st direction
  
  // 2nd direction
  paths=0;
  for (i1=0; i1<leng1; i1++) {
  for (i3=0; i3<leng3; i3++) {
    i2=0;
    // go through i3 until we hit a point in the cluster or i3=leng1-1
    for( int ii=0; ii<leng2; ii++) {
      istagged3[ii] = 0;
    }

    while (i2<leng2 && lclusterdata.isincluster[latmap(i1, i2, i3)] != c) {
      i2++;
    }

    if(i2 >= leng2)
      continue;
    if(istagged3[i2] == 1) {
      continue;
    }

    // go through the cluster part and count the paths and path++ onces until we hit a point not
    // in the cluster anymore or we hit the start point (percolating over this direction)
    paths++;
    startpoint = i2;
    // Go in +3 direction
    do {
      path2 = path2 + 1.0;
      istagged3[i2] = 1;
      i2++;
      if(i2 == leng2)
        i2 = i2 - leng2;
    } while ( lclusterdata.isincluster[latmap(i1, i2, i3)] == c && i2 != startpoint);

    
    // Go in -3 direction
    i2 = startpoint;
    do {
      if( istagged3[i2] == 0) {
        path2 = path2 + 1.0;
        istagged3[i2] = 1;
      }
      i2--;
      if(i2 == -1)
        i2 = i2 + leng2;
    } while ( lclusterdata.isincluster.at(latmap(i1, i2, i3)) == c && i2 != startpoint);
    
    i2 = startpoint + 1;
  }
  }
  path2=path2/(double)paths;
  // END 2nd direction
  
  // 3rd direction
  paths=0;
  for (i2=0; i2<leng2; i2++) {
  for (i3=0; i3<leng3; i3++) {
    i1=0;
    // go through i3 until we hit a point in the cluster or i3=leng1-1
    for( int ii=0; ii<leng1; ii++) {
      istagged3[ii] = 0;
    }

    while (i1<leng1 && lclusterdata.isincluster[latmap(i1, i2, i3)] != c) {
      i1++;
    }

    if(i1 >= leng1)
      continue;
    if(istagged3[i1] == 1) {
      continue;
    }

    // go through the cluster part and count the paths and path++ onces until we hit a point not
    // in the cluster anymore or we hit the start point (percolating over this direction)
    paths++;
    startpoint = i1;
    // Go in +3 direction
    do {
      path3 = path3 + 1.0;
      istagged3[i1] = 1;
      i1++;
      if(i1 == leng1)
        i1 = i1 - leng1;
    } while ( lclusterdata.isincluster[latmap(i1, i2, i3)] == c && i1 != startpoint);
    
    // Go in -3 direction
    i1 = startpoint;
    do {
      if( istagged3[i1] == 0) {
        path3 = path3 + 1.0;
        istagged3[i1] = 1;
      }
      i1--;
      if(i1 == -1)
        i1 = i1 + leng1;
    } while ( lclusterdata.isincluster.at(latmap(i1, i2, i3)) == c && i1 != startpoint);
    
    i1 = startpoint + 1;
  }
  }
  path3=path3/(double)paths;
  // END 2nd direction


  lobs.meanfreepath[c]=(path1+path2+path3)/(double)3.0;
  lobs.largestclustermeanfreepath=(path1+path2+path3)/(double)3.0;
  cout << "Total path = " << lobs.largestclustermeanfreepath << endl;
}

#endif // FINDCLUSTER_PATH_HPP
