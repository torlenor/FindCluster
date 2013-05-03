#ifndef FINDCLUSTER_RADIUS_HPP
#define FINDCLUSTER_RADIUS_HPP

void findShift(Clusterstruct &lclusterdata, int c){
	int shifti1=0, shifti2=0, shifti3=0;
	bool isempty=false;
	int i1, i2, i3, is;
	// if true, perform shift calculation in i1
	if(lclusterdata.clusterisperiodic[c][0] == 1){
		shifti1=0;
		// Find empty i2/i3 plane
		for(i1=0;i1<Ns;i1++){
			isempty=true;
			for(i2=0;i2<Ns;i2++){
				for(i3=0;i3<Ns;i3++){
					is = latmap(i1, i2, i3);
					if(lclusterdata.isincluster[is] == c)
						isempty=false;
				}
			}
			if(isempty==true){
				shifti1=i1;
				break;
			}
		}
		cout << "Cluster " << c << " shifti1 = " << shifti1 << endl;
	}
	
	// if true, perform shift calculation in i2
	if(lclusterdata.clusterisperiodic[c][1] == 1){
		shifti2=0;
		// Find empty i1/i3 plane
		for(i2=0;i2<Ns;i2++){
			isempty=true;
			for(i1=0;i1<Ns;i1++){
				for(i3=0;i3<Ns;i3++){
					is = latmap(i1, i2, i3);
					if(lclusterdata.isincluster[is] == c)
						isempty=false;
				}
			}
			if(isempty==true){
				shifti2=i2;
				break;
			}
		}
		cout << "Cluster " << c << " shifti2 = " << shifti2 << endl;
	}
	
	// if true, perform shift calculation in i3
	if(lclusterdata.clusterisperiodic[c][2] == 1){
		shifti3=0;
		// Find empty i1/i2 plane
		for(i3=0;i3<Ns;i3++){
			isempty=true;
			for(i1=0;i1<Ns;i1++){
				for(i2=0;i2<Ns;i2++){
					is = latmap(i1, i2, i3);
					if(lclusterdata.isincluster[is] == c)
						isempty=false;
				}
			}
			if(isempty==true){
				shifti3=i3;
				break;
			}
		}
		cout << "Cluster " << c << " shifti3 = " << shifti3 << endl;
	}
}

void clusterRadius(Observablestruct &lobs, Clusterstruct &lclusterdata){
	double centerofmass[3], radiussquare;
	int i1, i2, i3;

	lobs.centerofmass.resize(lclusterdata.clustermembers.size());

	for(unsigned int c=0;c<lclusterdata.clustermembers.size();c++){
		radiussquare=0;
		centerofmass[0]=0;
		centerofmass[1]=0;
		centerofmass[2]=0;
		// cout << "mini1 = " << mini1 << " mini2 = " << mini2 << " mini3 = " << mini3 << endl;
		for(unsigned int member=0; member<lclusterdata.clustermembers[c].size();member++){
			getCoords(lclusterdata.clustermembers[c][member], i1, i2, i3);
			centerofmass[0] += i1;
			centerofmass[1] += i2; 
			centerofmass[2] += i3; 
		}

		centerofmass[0] = centerofmass[0]/(double)lclusterdata.clustermembers[c].size();
		centerofmass[1] = centerofmass[1]/(double)lclusterdata.clustermembers[c].size();
		centerofmass[2] = centerofmass[2]/(double)lclusterdata.clustermembers[c].size();

		lobs.centerofmass[c].push_back(centerofmass[0]);
		lobs.centerofmass[c].push_back(centerofmass[1]);
		lobs.centerofmass[c].push_back(centerofmass[2]);
		
		for(unsigned int member=0; member<lclusterdata.clustermembers[c].size();member++){
			getCoords(lclusterdata.clustermembers[c][member], i1, i2, i3);
			radiussquare += (pow(centerofmass[0] - i1, 2)
					+ pow(centerofmass[1] - i2, 2)
					+ pow(centerofmass[2] - i3, 2));
		}
		radiussquare = sqrt(radiussquare/(double)lclusterdata.clustermembers[c].size());
		lobs.clusterradius.push_back(radiussquare);
	}
}

#endif // FINDCLUSTER_RADIUS_HPP
