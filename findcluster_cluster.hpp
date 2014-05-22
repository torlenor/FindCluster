#ifndef FINDCLUSTER_CLUSTER_HPP
#define FINDCLUSTER_CLUSTER_HPP

void sortClusterSize(Clusterstruct &lclusterdata){
	// Sorts clusters by their size and creates a sortedcluster vector
	// which contains the sorted cluster id's
	lclusterdata.sortedcluster.resize(lclusterdata.clustermembers.size());
	lclusterdata.isinsortedcluster.resize(opt.Nspace);
	
	for(unsigned int c=0; c<lclusterdata.clustermembers.size(); c++){
		lclusterdata.sortedcluster[c]=c;
	}
	
	int tmp1;
	bool notsorted=true;
	while(notsorted){
		notsorted=false;
		for(unsigned int c=0;c<lclusterdata.sortedcluster.size()-1;c++){
			if(lclusterdata.clustermembers[lclusterdata.sortedcluster[c+1]].size()<lclusterdata.clustermembers[lclusterdata.sortedcluster[c]].size()){
				tmp1=lclusterdata.sortedcluster[c+1];
				lclusterdata.sortedcluster[c+1]=lclusterdata.sortedcluster[c];
				lclusterdata.sortedcluster[c]=tmp1;
				notsorted=true;
			}
		}
	}
	
	reverse(lclusterdata.sortedcluster.begin(),lclusterdata.sortedcluster.end());
	
	for(unsigned c=0; c<lclusterdata.sortedcluster.size(); c++){
		for(unsigned member=0; member < lclusterdata.clustermembers[lclusterdata.sortedcluster[c]].size(); member++){
			lclusterdata.isinsortedcluster[lclusterdata.clustermembers[lclusterdata.sortedcluster[c]][member]] = c;
		}
	}

	for(unsigned c=0; c<lclusterdata.sortedcluster.size(); c++){
		if(lclusterdata.isinsector[lclusterdata.clustermembers[lclusterdata.sortedcluster[c]][0]] < 2)
			lclusterdata.sortedrealcluster.push_back(lclusterdata.sortedcluster[c]);
	}
}

void fillSectorsAlt(Clusterstruct &lclusterdata, double r){
	// Categorizes the lattice points by sector using the alternative
	// prescription with the radius
	#ifdef DEBUG
	cout << "Categorizing lattice points by sector ( radius = " << r << " )... " << flush;
	#endif

	int csectp1=0;
	int csectm1=0;
	int csect0=0;

	opt.delta = M_PI/3.0;
	
	double tracephase=0;
	double radius=0;
	
	#ifdef DEBUG
	cout << "Delta = " << opt.delta << endl;
	cout << "r = " << r << endl;
	#endif

	lclusterdata.poll.resize(opt.Nspace);
	
	for(int is=0;is<opt.Nspace;is++)
		lclusterdata.isinsector[is] = 2;
	
	for(int is=0;is<opt.Nspace;is++){
		tracephase = arg(pollev[is][0] + pollev[is][1] + pollev[is][2]);
		radius = abs(pollev[is][0] + pollev[is][1] + pollev[is][2]);
		
		if(abs(tracephase - 2.0*M_PI/3.0) <= opt.delta){
			if(radius >= r){
				lclusterdata.isinsector[is]=1;
				csectp1++;
			}
		}
		
		if(abs(tracephase) < opt.delta){
			if(radius >= r){
				lclusterdata.isinsector[is]=0;
				csect0++;
			}
		}
		
		if(abs(tracephase + 2.0*M_PI/3.0) <= opt.delta){
			if(radius >= r){
				lclusterdata.isinsector[is]=-1;
				csectm1++;
			}
		}
		lclusterdata.poll[is]=pollev[is][0] + pollev[is][1] + pollev[is][2];
	}
	
	lclusterdata.nsectm1=csectm1; lclusterdata.nsect0=csect0; lclusterdata.nsectp1=csectp1;
	#ifdef DEBUG
	cout << "done!" << endl;
	#endif
}

void fillSectors(Clusterstruct &lclusterdata, double delta){
	// Categorizes the lattice points by sector
	#ifdef DEBUG
	cout << "Categorizing lattice points by sector (delta = " << delta << " )... " << flush;
	#endif

	int csectp1=0;
	int csectm1=0;
	int csect0=0;
	
	double tracephase=0;
	
	#ifdef DEBUG
	cout << "Delta = " << delta << endl;
	#endif

	lclusterdata.poll.resize(opt.Nspace);
	
	for(int is=0;is<opt.Nspace;is++)
		lclusterdata.isinsector[is] = 2;
	
	for(int is=0;is<opt.Nspace;is++){
		tracephase = arg(pollev[is][0] + pollev[is][1] + pollev[is][2]);
		
		if(abs(tracephase - 2.0*M_PI/3.0) <= delta){
			lclusterdata.isinsector[is]=1;
			csectp1++;
		}
		
		if(abs(tracephase) < delta){
			lclusterdata.isinsector[is]=0;
			csect0++;
		}
		
		if(abs(tracephase + 2.0*M_PI/3.0) <= delta){
			lclusterdata.isinsector[is]=-1;
			csectm1++;
		}

		lclusterdata.poll[is]=pollev[is][0] + pollev[is][1] + pollev[is][2];
	}
	
	lclusterdata.nsectm1=csectm1; lclusterdata.nsect0=csect0; lclusterdata.nsectp1=csectp1;
	#ifdef DEBUG
	cout << "done!" << endl;
	#endif
}

void findPercolatingCluster(Clusterstruct &lclusterdata){
	// Finds the percolating clusters
	#ifdef DBEUG
	cout << "Searching for percolating clusters... " << flush;
	#endif
	
	int is=0;
	
	int cluster=0;
	
	lclusterdata.clusterispercolating.resize(lclusterdata.clustermembers.size());
	for(unsigned int c=0;c<lclusterdata.clustermembers.size();c++){
		lclusterdata.clusterispercolating[c] = 0;
	}
	
	vector<vector<vector<int> > > members;
	members.resize(lclusterdata.clustermembers.size());
	for(unsigned int i=0;i<members.size();i++){
		members[i].resize(3);
		members[i][0].resize(opt.leng1);
		for(int j=0;j<opt.leng1;j++)
			members[i][0][j]=0;
		members[i][1].resize(opt.leng2);
		for(int j=0;j<opt.leng2;j++)
			members[i][1][j]=0;
		members[i][2].resize(opt.leng3);
		for(int j=0;j<opt.leng3;j++)
			members[i][2][j]=0;
	}
	
	for(int i1=0;i1<opt.leng1;i1++)
	for(int i2=0;i2<opt.leng2;i2++)
	for(int i3=0;i3<opt.leng3;i3++){
		is = latmap(i1, i2, i3);
		cluster = lclusterdata.isincluster[is];
		if(lclusterdata.clustersector[cluster] < 2){
			members[cluster][0][i1] = 1;
			members[cluster][1][i2] = 1;
			members[cluster][2][i3] = 1;
		}
	}
	
	int sum1, sum2, sum3, percclustercnt=-1;
	for(unsigned int c=0;c<lclusterdata.clustermembers.size();c++){
		sum1=0; sum2=0; sum3=0;
		for(int i1=0;i1<opt.leng1;i1++)
			sum1 += members[c][0][i1];
		for(int i2=0;i2<opt.leng2;i2++)
			sum2 += members[c][1][i2];
		for(int i3=0;i3<opt.leng3;i3++)
			sum3 += members[c][2][i3];
			
		if(sum1 == opt.leng1 || sum2 == opt.leng2 || sum3 == opt.leng3){
			percclustercnt++;
			lclusterdata.percolatingclusters.push_back(c);
			lclusterdata.percolatingdirections.push_back(vector<int>());
			lclusterdata.percolatingdirections[percclustercnt].resize(3);
			lclusterdata.percolatingdirections[percclustercnt][0]=0;
			lclusterdata.percolatingdirections[percclustercnt][1]=0;
			lclusterdata.percolatingdirections[percclustercnt][2]=0;
			if(sum1 == opt.leng1)
				lclusterdata.percolatingdirections[percclustercnt][0]=1;
			if(sum2 == opt.leng2)
				lclusterdata.percolatingdirections[percclustercnt][1]=1;
			if(sum3 == opt.leng3)
				lclusterdata.percolatingdirections[percclustercnt][2]=1;
			#ifdef DEBUG
			cout << "sum1 = " << sum1 << " sum2 = " << sum2 << " sum3 = " << sum3 << endl;
			cout << "Cluster c = " << c << " added!" << endl;
			cout << "Percolating in direction (" << lclusterdata.percolatingdirections[percclustercnt][0] << ","
				<< lclusterdata.percolatingdirections[percclustercnt][1] << "," 
				<< lclusterdata.percolatingdirections[percclustercnt][2] << ")." << endl;
			#endif
		
			if(lclusterdata.percolatingdirections[percclustercnt][0] == 1 && 
			  lclusterdata.percolatingdirections[percclustercnt][1] == 1 &&
			  lclusterdata.percolatingdirections[percclustercnt][2] == 1){
				lclusterdata.clusterispercolating[c]=1;
			}
		}
	}
	#ifdef DBEUG	
	cout << "done!" << endl;
	#endif
}

void findClusters(Clusterstruct &lclusterdata){
	// Find clusters
	#ifdef DEBUG
	cout << "Finding clusters... " << flush;
	#endif
	int sector = 0, cluster = 0, memberis = 0, isneib = 0;

	vector<bool> isfiled;
	isfiled.resize(opt.Nspace);
	for(int is=0;is<opt.Nspace;is++){
		isfiled[is] = false;
		lclusterdata.isincluster[is] = -1;
	}
	
	int totalmembers=0;

	cluster = -1; // Note: We begin counting with ID 0!
	for(int is=0;is<opt.Nspace;is++){
		if(isfiled[is] == false){
			sector = lclusterdata.isinsector[is];
			cluster++; // current cluster id (Note: We begin counting with ID 0!)
			lclusterdata.clustersector.push_back(sector); // Current cluster is in sector sector
			
			lclusterdata.clustermembers.push_back(vector<int>()); // Add a new cluster
			lclusterdata.clustermembers[cluster].push_back(is); // Add is to this new cluster

			lclusterdata.isincluster[is] = cluster; // Set the cluster id for the point  is
			isfiled[is] = true; // Lattice point is now filed
			
			// Modification due to cluster radius with PBCs
			lclusterdata.clusterisperiodic.push_back(vector<int>()); // Add a new cluster
			lclusterdata.clusterisperiodic[cluster].resize(3);
			for(int i=0;i<3;i++)
				lclusterdata.clusterisperiodic[cluster][i]=0;
			int ci1=0, ci2=0, ci3=0;
		
			// Iterate over all found points in the actual cluster (loop does not have fixed length)
			int member = 0;
			while((unsigned int)member < lclusterdata.clustermembers[cluster].size()){
				memberis = lclusterdata.clustermembers[cluster][member];
			
				for(int mu=0;mu<6;mu++){
					isneib = neib[memberis][mu];
					if(isfiled[isneib] == false && lclusterdata.isinsector[isneib] == sector){
						lclusterdata.clustermembers[cluster].push_back(isneib);
						lclusterdata.isincluster[isneib] = cluster;
						isfiled[isneib] = true;

						// Check if over periodic boundaries
						getCoords(memberis, ci1, ci2, ci3);
						switch(mu){
							case 0: if(ci1 + 1 == opt.leng1)
									lclusterdata.clusterisperiodic[cluster][0]=1;
								break;

							case 1: if(ci2 + 1 == opt.leng2)
									lclusterdata.clusterisperiodic[cluster][1]=1;
								break;
							
							case 2: if(ci3 + 1 == opt.leng3)
									lclusterdata.clusterisperiodic[cluster][2]=1;
								break;
							
							case 3: if(ci1 - 1 == -1)
									lclusterdata.clusterisperiodic[cluster][0]=1;
								break;
							
							case 4: if(ci2 - 1 == -1)
									lclusterdata.clusterisperiodic[cluster][1]=1;
								break;
							
							case 5: if(ci3 - 1 == -1)
									lclusterdata.clusterisperiodic[cluster][2]=1;
								break;

							default: cout << "ERROR: Error in cluster switch!" << endl;
						}
					}
				}
				
				#ifdef DEBUG
				cout << endl <<  "Cluster " << cluster << " member " << member + 1 << " of " << lclusterdata.clustermembers[cluster].size() << " members." << endl;
				#endif
				
				member++;
				totalmembers++;
			}
		}
	}
	
	#ifdef DEBUG
	cout << "" << totalmembers << " of " << opt.Nspace << " lattice points visited!" << endl;
	#endif
	
	for(int is=0;is<opt.Nspace;is++){
		if(isfiled[is] == false)
			cout << "WARNING: Some points did not get filled!" << endl;
	}
	#ifdef DEBUG	
	cout << "done!" << endl;
	#endif
}

void checkClusters(Clusterstruct &lclusterdata){
	// Checks the clusters
	#ifdef DEBUG	
	cout << "Checking clusters... " << flush;
	#endif

	for(int is=0;is<opt.Nspace;is++){
		if(lclusterdata.isincluster[is] < 0)
			cout << "WARNING: Some points do not have a valid cluster ID!" << endl;
	}
	
	int sum=0;
	for(unsigned int c=0;c<lclusterdata.clustermembers.size();c++)
		sum+=lclusterdata.clustermembers[c].size();
	if(sum != opt.Nspace){
		cout << "WARNING: Number of lattice points differ from number of clustermember entries!" << endl;
		cout << "Number of entries = " << sum << " Number of lattice points = " << opt.Nspace << endl;
	}
	
	// Check for conistency of the cluster arrays lclusterdata.isinsector and lclusterdata.clustersector
	for(unsigned int c=0;c<lclusterdata.clustermembers.size();c++){
		for(unsigned int member=0; member<lclusterdata.clustermembers[c].size(); member++){
			if(lclusterdata.isinsector[lclusterdata.clustermembers[c][member]] != lclusterdata.clustersector[c]){
				cout << "WARNING: Problem with lclusterdata.clustermembers and lclusterdata.clustersector consistency!" << endl;
			}
		}
	}
	
	// Check if all members in clustermmembers have a correct ID
	// Here we use vector.at(iter) so that also a range check is performed and an
	// exeption is thrown if it is outside vector range
	int mult=1;
	vector<int> membercheck;
	membercheck.resize(opt.Nspace);
	for(unsigned int c=0;c<lclusterdata.clustermembers.size();c++){
		for(unsigned int member=0; member<lclusterdata.clustermembers[c].size(); member++){
			membercheck.at(lclusterdata.clustermembers[c][member]) = 1;
		}
	}
	for(int is=0;is<opt.Nspace;is++){
		mult = mult*membercheck[is];
	}
	if(mult != 1){
		cout << "WARNING: Problem with clustermember IDs!" << endl;
	}
	
	#ifdef DEBUG	
	cout << "done!" << endl;
	#endif
}

#endif // FINDCLUSTER_CLUSTER_HPP
