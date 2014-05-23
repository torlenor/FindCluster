#ifndef FINDCLUSTER_RADIUS_HPP
#define FINDCLUSTER_RADIUS_HPP

void getCoordsShift(int is, int &i1, int &i2, int &i3, int *shift){
	        i1 = (is % (opt.leng1*opt.leng2) ) % opt.leng1;
	        i2 = (is % (opt.leng1*opt.leng2) ) / opt.leng1;
	        i3 = is / (opt.leng1*opt.leng2);

		if(i1>shift[0])
			i1=i1-opt.leng1;
		if(i2>shift[1])
			i2=i2-opt.leng2;
		if(i3>shift[2])
			i3=i3-opt.leng3;

	        // if(is != i1 + i2*opt.leng1 + i3*opt.leng1*opt.leng2)
	        //        cout << "ERROR: Problem in getCoords!" << endl;
}

void obsClusterRadius(Observablestruct &lobs, Clusterstruct &lclusterdata){
	// Calculation of the cluster radius. We save the largest cluster (in terms
	// of the cluster radius).
	double centerofmass[3], radiussquare, radiussquaremin, centerofmassmin[3];
	int i1, i2, i3;

	int shift[3];

	lobs.centerofmass.resize(lclusterdata.clustermembers.size());
	lobs.clusterradius.resize(lclusterdata.clustermembers.size());

//	int c = lobs.maxclusterid;

//	double radiussquaremax=0;

	for(unsigned int c=0;c<lclusterdata.clustermembers.size();c++){
		if(lclusterdata.clustermembers[c].size()>1 && lclusterdata.clustersector[c] < 2){
			// Do it only for clusters with size > 1 in sectors < 2

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
					getCoordsShift(lclusterdata.clustermembers[c][member], i1, i2, i3, shift);
					centerofmass[0] += i1;
					centerofmass[1] += i2; 
					centerofmass[2] += i3; 
				}

				centerofmass[0] = centerofmass[0]/(double)lclusterdata.clustermembers[c].size();
				centerofmass[1] = centerofmass[1]/(double)lclusterdata.clustermembers[c].size();
				centerofmass[2] = centerofmass[2]/(double)lclusterdata.clustermembers[c].size();

				for(unsigned int member=0; member<lclusterdata.clustermembers[c].size();member++){
					getCoordsShift(lclusterdata.clustermembers[c][member], i1, i2, i3, shift);
					radiussquare += (pow(centerofmass[0] - i1, 2)
							+ pow(centerofmass[1] - i2, 2)
							+ pow(centerofmass[2] - i3, 2));
				}
				radiussquare = sqrt(radiussquare/(double)lclusterdata.clustermembers[c].size());
				if(radiussquare<radiussquaremin){
					radiussquaremin=radiussquare;
					centerofmassmin[0]=centerofmass[0];
					centerofmassmin[1]=centerofmass[1];
					centerofmassmin[2]=centerofmass[2];
				}
			}

			lobs.centerofmass[c].push_back(centerofmassmin[0]);
			lobs.centerofmass[c].push_back(centerofmassmin[1]);
			lobs.centerofmass[c].push_back(centerofmassmin[2]);
			
			lobs.clusterradius[c]=radiussquaremin;
//			if(radiussquaremin>radiussquaremax)
//				radiussquaremax=radiussquaremin;
		}
	}
	lobs.largestclusterradius=lobs.clusterradius[lobs.maxclusterid];
	// lobs.largestclusterradius=radiussquaremax;
	// lobs.largestclusterradius=radiussquaremin;

	// Additional radius observables
	// Radius largest non-percolating cluster
	lobs.largestnpclusterradius=lobs.clusterradius[lobs.largestnonpercclusterid];

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

void obsClusterRadiusOnlyLargest(Observablestruct &lobs, Clusterstruct &lclusterdata){
	// Calculation of the cluster radius. We save the largest cluster (in terms
	// of the cluster weight).
	double centerofmass[3], radiussquare, radiussquaremin, centerofmassmin[3];
	int i1, i2, i3;

	int shift[3];

	lobs.centerofmass.resize(lclusterdata.clustermembers.size());
	lobs.clusterradius.resize(lclusterdata.clustermembers.size());

	int c = lobs.maxclusterid;

//	double radiussquaremax=0;

//	for(unsigned int c=0;c<lclusterdata.clustermembers.size();c++){
//		if(lclusterdata.clustermembers[c].size()>1 && lclusterdata.clustersector[c] < 2){
			// Do it only for clusters with size > 1 in sectors < 2

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
					getCoordsShift(lclusterdata.clustermembers[c][member], i1, i2, i3, shift);
					centerofmass[0] += i1;
					centerofmass[1] += i2; 
					centerofmass[2] += i3; 
				}

				centerofmass[0] = centerofmass[0]/(double)lclusterdata.clustermembers[c].size();
				centerofmass[1] = centerofmass[1]/(double)lclusterdata.clustermembers[c].size();
				centerofmass[2] = centerofmass[2]/(double)lclusterdata.clustermembers[c].size();

				for(unsigned int member=0; member<lclusterdata.clustermembers[c].size();member++){
					getCoordsShift(lclusterdata.clustermembers[c][member], i1, i2, i3, shift);
					radiussquare += (pow(centerofmass[0] - i1, 2)
							+ pow(centerofmass[1] - i2, 2)
							+ pow(centerofmass[2] - i3, 2));
				}
				radiussquare = sqrt(radiussquare/(double)lclusterdata.clustermembers[c].size());
				if(radiussquare<radiussquaremin){
					radiussquaremin=radiussquare;
					centerofmassmin[0]=centerofmass[0];
					centerofmassmin[1]=centerofmass[1];
					centerofmassmin[2]=centerofmass[2];
				}
			}

			lobs.centerofmass[c].push_back(centerofmassmin[0]);
			lobs.centerofmass[c].push_back(centerofmassmin[1]);
			lobs.centerofmass[c].push_back(centerofmassmin[2]);
			
//			lobs.clusterradius[c]=radiussquaremin;
//			if(radiussquaremin>radiussquaremax)
//				radiussquaremax=radiussquaremin;
//		}
//	}
	// lobs.largestclusterradius=lobs.clusterradius[lclusterdata.sortedcluster[0]];
	// lobs.largestclusterradius=radiussquaremax;
	lobs.largestclusterradius=radiussquaremin;
}

void obsClusterRadiusOnlyLargestNP(Observablestruct &lobs, Clusterstruct &lclusterdata){
	// Calculation of the cluster radius. We save the largest cluster (in terms
	// of the cluster weight).
	double centerofmass[3], radiussquare, radiussquaremin, centerofmassmin[3];
	int i1, i2, i3;

	int shift[3];

	int c = lobs.largestnonpercclusterid;

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
      getCoordsShift(lclusterdata.clustermembers[c][member], i1, i2, i3, shift);
      centerofmass[0] += i1;
      centerofmass[1] += i2; 
      centerofmass[2] += i3; 
    }

    centerofmass[0] = centerofmass[0]/(double)lclusterdata.clustermembers[c].size();
    centerofmass[1] = centerofmass[1]/(double)lclusterdata.clustermembers[c].size();
    centerofmass[2] = centerofmass[2]/(double)lclusterdata.clustermembers[c].size();

    for(unsigned int member=0; member<lclusterdata.clustermembers[c].size();member++){
      getCoordsShift(lclusterdata.clustermembers[c][member], i1, i2, i3, shift);
      radiussquare += (pow(centerofmass[0] - i1, 2)
          + pow(centerofmass[1] - i2, 2)
          + pow(centerofmass[2] - i3, 2));
    }

    radiussquare = sqrt(radiussquare/(double)lclusterdata.clustermembers[c].size());

    if(radiussquare<radiussquaremin){
      radiussquaremin=radiussquare;
      centerofmassmin[0]=centerofmass[0];
      centerofmassmin[1]=centerofmass[1];
      centerofmassmin[2]=centerofmass[2];
    }
  }

  lobs.centerofmass[c].push_back(centerofmassmin[0]);
  lobs.centerofmass[c].push_back(centerofmassmin[1]);
  lobs.centerofmass[c].push_back(centerofmassmin[2]);
  
  lobs.largestnpclusterradius=radiussquaremin;
}

#endif // FINDCLUSTER_RADIUS_HPP
