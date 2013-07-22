#ifndef FINDCLUSTER_WRITE_HPP
#define FINDCLUSTER_WRITE_HPP

void writeresults(){
	stringstream fclustersizename;
	fclustersizename << "clustersize_" << Ns << "x" << Nt << "_f" << fraction << ".res";
	ofstream fclustersize;
	fclustersize.open(fclustersizename.str().c_str());
	fclustersize << "# Nt largestclusterweight largestclusterweighterr avgclusterweight avgclusterweighterr avgfortunatoclustersize avgfortunatoclustersizeerr largestnonpercclusterweight largestnonpercclusterweighterr largestclusterradius largestclusterradiuserr rootmeansquaredistanceR rootmeansquaredistanceR" << endl;
	fclustersize.flags (std::ios::scientific);
	fclustersize.precision(numeric_limits<double>::digits10 + 1);
	fclustersize << Nt << " " << results.maxclustersize << " " << results.maxclustersizeerr << " " << results.avgclustersize << " " << results.avgclustersizeerr << " " << results.avgclusersizeFortunato << " " << results.avgclusersizeFortunatoerr << " " << results.maxnonpercclustersize << " " << results.maxnonpercclustersizeerr << " " << results.largestclusterradius << " " << results.largestclusterradiuserr << " " << results.avgrootmeansquaredistance << " " << results.avgrootmeansquaredistanceerr << endl;
	fclustersize.close();
	
	stringstream favgclustersizename;
	favgclustersizename << "clustersize_avg_" << Ns << "x" << Nt << "_f" << fraction << ".res";
	ofstream favgclustersize;
	favgclustersize.open(favgclustersizename.str().c_str());
	favgclustersize << "# Nt avgclusterweight avgclusterweighterr avgfortunatoclustersize avgfortunatoclustersizeerr avgclusterweightnopercc avgclusterweightnoperccerr" << endl;
	favgclustersize.flags (std::ios::scientific);
	favgclustersize.precision(numeric_limits<double>::digits10 + 1);
	favgclustersize << Nt << " " << results.avgclustersize << " " << results.avgclustersizeerr << " " << results.avgclusersizeFortunato << " " << results.avgclusersizeFortunatoerr << " " << results.avgnonpercclustersize << " " << results.avgnonpercclustersizeerr << endl;
	favgclustersize.close();

	stringstream fnperccname;
	fnperccname << "npercc_" << Ns << "x" << Nt << "_f" << fraction << ".res";
	ofstream fnpercc;
	fnpercc.open(fnperccname.str().c_str());
	fnpercc << "# Nt percclusters percclusterserr" << endl;
	fnpercc.flags (std::ios::scientific);
	fnpercc.precision(numeric_limits<double>::digits10 + 1);
	fnpercc << Nt << " " << results.avgperccluster << " " << results.avgpercclustererr << endl;
	fnpercc.close();
	
	if(doboxes){
		stringstream fboxcntname;
		fboxcntname << "boxcnt_" << Ns << "x" << Nt << "_f" << fraction << ".res";
		ofstream fboxcnt;
		fboxcnt.open(fboxcntname.str().c_str());
		fboxcnt << "# boxsize boxcnt boxcnterr" << endl;
		fboxcnt.flags (std::ios::scientific);
		fboxcnt.precision(numeric_limits<double>::digits10 + 1);
		for(unsigned int size=0; size<boxsize.size(); size++){
			fboxcnt << boxsize[size] << " " << results.largestclusterboxcount[size] << " " << results.largestclusterboxcounterr[size] << endl;
		}
		fboxcnt.close();
		
		fboxcntname.str("");
		fboxcntname << "boxcnt_nonpercc_" << Ns << "x" << Nt << "_f" << fraction << ".res";
		fboxcnt.open(fboxcntname.str().c_str());
		fboxcnt << "# boxsize boxcntnonpercc boxcntnonperccerr" << endl;
		fboxcnt.flags (std::ios::scientific);
		fboxcnt.precision(numeric_limits<double>::digits10 + 1);
		for(unsigned int size=0; size<boxsize.size(); size++){
			fboxcnt << boxsize[size] << " " << results.largestnonpercboxcount[size] << " " << results.largestnonpercboxcounterr[size] << endl;
		}
		fboxcnt.close();
	}
	
	stringstream fareaname;
	fareaname << "area_" << Ns << "x" << Nt << "_f" << fraction << ".res";
	ofstream farea;
	farea.open(fareaname.str().c_str());
	farea << "# Nt area areaerr arealargestnonpercc arealargestnonperccerr largestclusterweight largestclusterweighterr largestnonpercclusterweight largestnonpercclusterweighterr" << endl;
	farea.flags (std::ios::scientific);
	farea.precision(numeric_limits<double>::digits10 + 1);
	farea << Nt << " " << results.totalperimeter << " " << results.totalperimetererr << " " << results.largestnonpercperimeter << " " << results.largestnonpercperimetererr << " " << results.maxclustersize << " " << results.maxclustersizeerr << " " << results.maxnonpercclustersize << " " << results.maxnonpercclustersizeerr << endl;
	farea.close();
	
	stringstream fpollname;
	fpollname << "poll_" << Ns << "x" << Nt << "_f" << fraction << ".res";
	ofstream fpoll;
	fpoll.open(fpollname.str().c_str());
	fpoll << "# Nt poll pollerr domainwallpoll(largest cluster) domainwallpollerr(largest cluster) avgdomainwallpoll avgdomainwallpollerr" << endl;
	fpoll.flags (std::ios::scientific);
	fpoll.precision(numeric_limits<double>::digits10 + 1);
	fpoll << Nt << " " << results.polyakovloopaftercut << " " << results.polyakovloopaftercuterr << " " << results.largestclusterdomainwallpoll << " " << results.largestclusterdomainwallpollerr << " " << results.avgdomainwallpoll << " " << results.avgdomainwallpollerr << endl;
	fpoll.close();
}

void writeresultsstdout(){
	cout << "Expectation values (single eliminitation jackknife): " << endl;
	cout << "Average cluster size = " << setprecision(14) << results.avgclustersize << ", Maximum cluster size / V = " << results.maxclustersize << endl;
	cout << "Average cluster err  = " << results.avgclustersizeerr << ", Maximum cluster / V err  = " << results.maxclustersizeerr << endl;
	cout << "Average cluster size Fortunato (1.7) = " << results.avgclusersizeFortunato << ", Error  = " << results.avgclusersizeFortunatoerr << endl;
	cout << "Radius of largest cluster = " << results.largestclusterradius << ", Error = " << results.largestclusterradiuserr << endl;
	if(dodistance){
		cout << "Root mean distance traveled R = " << results.avgrootmeansquaredistance << ", Error = " << results.avgrootmeansquaredistanceerr << endl;
	}
	cout << "Cut = " << results.cut << " Cut err = " << results.cuterr << endl;
	cout << "Laserdim = " << results.totalperimeter << " Laserdim err = " << results.totalperimetererr << endl;
	cout << "Polyakov loop = " << results.polyakovloopaftercut << " Polyakov loop err = " << results.polyakovloopaftercuterr << endl;
	cout << "Polyakov loop (Domain Wall largest cluster) = " << results.largestclusterdomainwallpoll << ", Error = " << results.largestclusterdomainwallpollerr << endl;
	cout << "Polyakov loop (Domain Wall average) = " << results.avgdomainwallpoll << ", Error = " << results.avgdomainwallpollerr << endl;

	cout << endl;
}

void writeConfigResultsstdout(Observablestruct &lobs, Clusterstruct &lclusterdata){
	cout << "Number of points in the different sectors:" << endl;
	cout << "Sector -1 # = " << lclusterdata.nsectm1 << " Sector 0 # = " << lclusterdata.nsect0 << " Sector 1 # = " << lclusterdata.nsectp1 << ", Sectors cut = " << Nspace - lclusterdata.nsectm1 - lclusterdata.nsect0 - lclusterdata.nsectp1 << endl << endl;

	cout  << lclusterdata.percolatingclusters.size() << " of " << lclusterdata.clustermembers.size() << " clusters are percolating!" << endl;
	for(unsigned int c=0;c<lclusterdata.percolatingclusters.size();c++){
		cout << "Cluster " << lclusterdata.percolatingclusters[c] << " is percolating in directions ("
			<< lclusterdata.percolatingdirections[c][0] << "," 
			<< lclusterdata.percolatingdirections[c][1] << ","
			<< lclusterdata.percolatingdirections[c][2] << ")!" << endl;
		cout << "The cluster has " << lclusterdata.clustermembers[lclusterdata.percolatingclusters[c]].size() << " members and is in sector " << lclusterdata.clustersector[lclusterdata.percolatingclusters[c]] << " !" << endl;
	}
	
	cout << endl;
	cout << "Average cluster size = " << lobs.avgclustersize << endl;
	cout << "Average cluster size Fortunato = " << lobs.avgclustersizeF << endl;
	cout << "Largest cluster is cluster " << lobs.maxclusterid << " with " << lobs.maxclustersize << " members." << endl;

	if(doboxes){
		cout << endl << "Boxes (Clusters sorted based on # of members):" << endl;
		for(unsigned int c=0; c<lclusterdata.sortedcluster.size(); c++){
			cout << "Cluster = " << lclusterdata.sortedcluster[c] << " with " << lclusterdata.clustermembers[lclusterdata.sortedcluster[c]].size() << " members:" << endl;
			for(unsigned int size=0; size<boxsize.size(); size++){
			cout << lobs.numberofboxes[lclusterdata.sortedcluster[c]][size] << " boxes of size " << boxsize[size] << " needed." << endl;
			}
		}
	}

	cout << endl << "The following clusters are connected over the PBCs: " << endl;
	for(unsigned int c=0;c<lclusterdata.clusterisperiodic.size();c++){
		if(lclusterdata.clusterisperiodic[lclusterdata.sortedcluster[c]][0]==1 
		|| lclusterdata.clusterisperiodic[lclusterdata.sortedcluster[c]][1]==1 
		|| lclusterdata.clusterisperiodic[lclusterdata.sortedcluster[c]][2]==1){
			cout << "Cluster (size sorted) " << c << ", real id = " << lclusterdata.sortedcluster[c] << endl;
		}
	}
}

void writeClusterList(Clusterstruct &lclusterdata){
	vector<int> csize;
	for(unsigned c=0; c<lclusterdata.clustermembers.size(); c++){
		csize.push_back(lclusterdata.clustermembers[c].size());
	}

	sort(csize.begin(), csize.end());

	cout << "Cluster sizes" << endl;
	for(unsigned c=0; c<lclusterdata.clustermembers.size(); c++){
		cout << csize[c] << endl;
	}
}

void cluster3doutput(Clusterstruct &lclusterdata, string f3dname){
	ofstream f3d;
	int is;
	f3d.open(f3dname.c_str());
	if(f3d.is_open()){
		f3d << leng1 << " " << leng2 << " " << leng3 << " " << leng4 << endl;
		f3d << lclusterdata.clustermembers.size() << endl;
		for(int i1=0;i1<leng1;i1++)
		for(int i2=0;i2<leng2;i2++)
		for(int i3=0;i3<leng3;i3++){
			is = i1 + i2*leng1 + i3*leng1*leng2;

			f3d << i1 << " " << i2 << " " << i3 << " " << lclusterdata.isinsortedcluster[is] << " " << lclusterdata.isinsector[is] << endl;
		}
		f3d.close();
	}else{
		cout << "WARNING: Could not open 3dcluster file!" << endl;
	}
}

#endif // FINDCLUSTER_WRITE_HPP
