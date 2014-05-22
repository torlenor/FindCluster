#ifndef FINDCLUSTER_WRITE_HPP
#define FINDCLUSTER_WRITE_HPP

void writeresults(){
	stringstream fclustersizename;
	fclustersizename << "clustersize_" << opt.Ns << "x" << opt.Nt << "_f" << opt.fraction << ".res";
	ofstream fclustersize;
	fclustersize.open(fclustersizename.str().c_str());
	fclustersize << "# Nt largestclusterweight largestclusterweighterr avgclusterweight avgclusterweighterr avgfortunatoclustersize avgfortunatoclustersizeerr largestnonpercclusterweight largestnonpercclusterweighterr largestclusterradius largestclusterradiuserr rootmeansquaredistanceR rootmeansquaredistanceR" << endl;
	fclustersize.flags (std::ios::scientific);
	fclustersize.precision(numeric_limits<double>::digits10 + 1);
	fclustersize << opt.Nt << " " 
	  << results.maxclustersize << " " << results.maxclustersizeerr << " " 
	  << results.avgclustersize << " " << results.avgclustersizeerr << " " 
	  << results.avgclusersizeFortunato << " " << results.avgclusersizeFortunatoerr << " " 
	  << results.maxnonpercclustersize << " " << results.maxnonpercclustersizeerr << " " 
	  << results.largestclusterradius << " " << results.largestclusterradiuserr << " " 
	  << results.avgrootmeansquaredistance << " " << results.avgrootmeansquaredistanceerr 
	<< endl;
	fclustersize.close();
	
	stringstream favgclustersizename;
	favgclustersizename << "clustersize_avg_" << opt.Ns << "x" << opt.Nt << "_f" << opt.fraction << ".res";
	ofstream favgclustersize;
	favgclustersize.open(favgclustersizename.str().c_str());
	favgclustersize << "# Nt avgclusterweight avgclusterweighterr avgfortunatoclustersize avgfortunatoclustersizeerr avgclusterweightFnp avgclusterweightFnperr avgclusterweightnp avgclusterweightnperr avgnpclusterradius avgnpclusterradiuserr avgnpclusterarea avgnpclusterareaerr" << endl;
	favgclustersize.flags (std::ios::scientific);
	favgclustersize.precision(numeric_limits<double>::digits10 + 1);
	favgclustersize << opt.Nt << " " 
    << results.avgclustersize << " " << results.avgclustersizeerr << " " 
    << results.avgclusersizeFortunato << " " << results.avgclusersizeFortunatoerr << " " 
    << results.avgclustersizenp << " " << results.avgclustersizenperr << " " 
    << results.avgclustersizeFnp << " " << results.avgclustersizeFnperr << " " 
    << results.avgnpclusterradius << " " << results.avgnpclusterradiuserr << " "
    << results.avgnonpercperimeter << " " << results.avgnonpercperimetererr <<  " "
	<< endl;
	favgclustersize.close();
	
	stringstream fclusterradiusname;
	fclusterradiusname << "clusterradius_" << opt.Ns << "x" << opt.Nt << "_f" << opt.fraction << ".res";
	ofstream fclusterradius;
	fclusterradius.open(fclusterradiusname.str().c_str());
	fclusterradius << "# Nt largestclusterradius largestclusterradiuserr largestnpclusterradius largestnpclusterradiuserr avgclusterradius avgclusterradiuserr avgnpcluterradius avgnpclusterradiuserr rootmeansquaredistanceR rootmeansquaredistanceRerr maxclustersize maxclustersizeerr maxnonpercclustersize maxnonpercclustersizeerr avgclustersize avgclustersizeerr avgclusersizeFortunato avgclusersizeFortunatoerr avgclustersizeFnp avgclustersizeFnperr" << endl;
	fclusterradius.flags (std::ios::scientific);
	fclusterradius.precision(numeric_limits<double>::digits10 + 1);
	fclusterradius << opt.Nt << " " << results.largestclusterradius << " " << results.largestclusterradiuserr << " "
    << results.largestnpclusterradius << " " << results.largestnpclusterradiuserr << " "
    << results.avgclusterradius << " " << results.avgclusterradiuserr << " " 
    << results.avgnpclusterradius << " " << results.avgnpclusterradiuserr << " "
    << results.avgrootmeansquaredistance << " " << results.avgrootmeansquaredistanceerr << " "
    << results.maxclustersize << " " << results.maxclustersizeerr << " "
    << results.maxnonpercclustersize << " " << results.maxnonpercclustersizeerr << " "
    << results.avgclustersize << " " << results.avgclustersizeerr << " "
    << results.avgclusersizeFortunato << " " << results.avgclusersizeFortunatoerr << " "
    << results.avgclustersizeFnp << " " << results.avgclustersizeFnperr
	<< endl;
	fclusterradius.close();

	stringstream fnperccname;
	fnperccname << "npercc_" << opt.Ns << "x" << opt.Nt << "_f" << opt.fraction << ".res";
	ofstream fnpercc;
	fnpercc.open(fnperccname.str().c_str());
	fnpercc << "# Nt percclusters percclusterserr cut cuterr" << endl;
	fnpercc.flags (std::ios::scientific);
	fnpercc.precision(numeric_limits<double>::digits10 + 1);
	fnpercc << opt.Nt << " " << results.avgperccluster << " " << results.avgpercclustererr << " " 
    << results.cut << " " << results.cuterr 
  << endl;
	fnpercc.close();
	
  stringstream ffreepathname;
	ffreepathname << "freepath_" << opt.Ns << "x" << opt.Nt << "_f" << opt.fraction << ".res";
	ofstream ffreepath;
	ffreepath.open(ffreepathname.str().c_str());
	ffreepath << "# Nt largestcluster largestclustererr largestnpcluster largestnpclustererr avgcluster avgclustererr avgnpcluster avgnpclustererr" << endl;
	ffreepath.flags (std::ios::scientific);
	ffreepath.precision(numeric_limits<double>::digits10 + 1);
	ffreepath << opt.Nt << " " 
    << results.largestclustermeanfreepath << " " << results.largestclustermeanfreepatherr << " " 
    << results.largestnpclustermeanfreepath << " " << results.largestnpclustermeanfreepatherr << " " 
    << results.avgclustermeanfreepath << " " << results.avgclustermeanfreepatherr << " " 
    << results.avgnpclustermeanfreepath << " " << results.avgnpclustermeanfreepatherr
  << endl;
	ffreepath.close();
	
	// NEW DEFINITION OF MEAN FREE PATH
	stringstream ffreepathnewname;
	ffreepathnewname << "meanfreepath_new_" << opt.Ns << "x" << opt.Nt << "_f" << opt.fraction << ".res";
	ofstream ffreepathnew;
	ffreepathnew.open(ffreepathnewname.str().c_str());
	ffreepathnew << "# 1Nt 2largestcluster 3largestclustererr 4largestnpcluster 5largestnpclustererr 6avgcluster 7avgclustererr 8avgnpcluster 9avgnpclustererr 10avgFcluster 11avgFclustererr 12avgFnpcluster 13avgFnpclustererr" << endl;
	ffreepathnew.flags (std::ios::scientific);
	ffreepathnew.precision(numeric_limits<double>::digits10 + 1);
	ffreepathnew << opt.Nt << " " 
    << results.largestclustermeanfreepathnew << " " << results.largestclustermeanfreepathnewerr << " " 
    << results.largestnpclustermeanfreepathnew << " " << results.largestnpclustermeanfreepathnewerr << " " 
    << results.avgclustermeanfreepathnew << " " << results.avgclustermeanfreepathnewerr << " " 
    << results.avgnpclustermeanfreepathnew << " " << results.avgnpclustermeanfreepathnewerr << " " 
    << results.avgFclustermeanfreepathnew << " " << results.avgFclustermeanfreepathnewerr << " " 
    << results.avgFnpclustermeanfreepathnew << " " << results.avgFnpclustermeanfreepathnewerr
  << endl;
	ffreepathnew.close();
	// END NEW DEFINITION OF MEAN FREE PATH
	
	if(opt.doboxes){
		stringstream fboxcntname;
		fboxcntname << "boxcnt_" << opt.Ns << "x" << opt.Nt << "_f" << opt.fraction << ".res";
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
		fboxcntname << "boxcnt_nonpercc_" << opt.Ns << "x" << opt.Nt << "_f" << opt.fraction << ".res";
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
	fareaname << "area_" << opt.Ns << "x" << opt.Nt << "_f" << opt.fraction << ".res";
	ofstream farea;
	farea.open(fareaname.str().c_str());
	farea << "# Nt area areaerr arealargestnonpercc arealargestnonperccerr largestclusterweight largestclusterweighterr largestnonpercclusterweight largestnonpercclusterweighterr areaavgnonpercc areaavgnonperccerr avgclusersizeFnp avgclusersizeFnperr" << endl;
	farea.flags (std::ios::scientific);
	farea.precision(numeric_limits<double>::digits10 + 1);
	farea << opt.Nt << " " 
    << results.totalperimeter << " " << results.totalperimetererr << " " 
    << results.largestnonpercperimeter << " " << results.largestnonpercperimetererr << " " 
    << results.maxclustersize << " " << results.maxclustersizeerr << " " 
    << results.maxnonpercclustersize << " " << results.maxnonpercclustersizeerr <<  " " 
    << results.avgnonpercperimeter << " " << results.avgnonpercperimetererr <<  " "
    << results.avgclustersizeFnp << " " << results.avgclustersizeFnperr << " "
	<< endl;
	farea.close();
	
	stringstream fpollname;
	fpollname << "poll_" << opt.Ns << "x" << opt.Nt << "_f" << opt.fraction << ".res";
	ofstream fpoll;
	fpoll.open(fpollname.str().c_str());
	fpoll << "# Nt poll pollerr domainwallpoll(largest cluster) domainwallpollerr(largest cluster) avgdomainwallpoll avgdomainwallpollerr" << endl;
	fpoll.flags (std::ios::scientific);
	fpoll.precision(numeric_limits<double>::digits10 + 1);
	fpoll << opt.Nt << " " << results.polyakovloopaftercut << " " << results.polyakovloopaftercuterr << endl;
	fpoll.close();
}

void writeresultsstdout(){
	cout << "Expectation values (single eliminitation jackknife): " << endl;
	cout << "Average cluster size = " << setprecision(14) << results.avgclustersize << ", Maximum cluster size / V = " << results.maxclustersize << endl;
	cout << "Average cluster err  = " << results.avgclustersizeerr << ", Maximum cluster / V err  = " << results.maxclustersizeerr << endl;
	cout << "Average cluster size Fortunato (1.7) = " << results.avgclusersizeFortunato << ", Error  = " << results.avgclusersizeFortunatoerr << endl;
	cout << "Radius of largest cluster = " << results.largestclusterradius << ", Error = " << results.largestclusterradiuserr << endl;
	if(opt.dodistance){
		cout << "Root mean distance traveled R = " << results.avgrootmeansquaredistance << ", Error = " << results.avgrootmeansquaredistanceerr << endl;
	}
	cout << "Cut = " << results.cut << " Cut err = " << results.cuterr << endl;
	cout << "Laserdim = " << results.totalperimeter << " Laserdim err = " << results.totalperimetererr << endl;
	cout << "Polyakov loop = " << results.polyakovloopaftercut << " Polyakov loop err = " << results.polyakovloopaftercuterr << endl;

	cout << endl;
}

void writeConfigResultsstdout(Observablestruct &lobs, Clusterstruct &lclusterdata){
	cout << "Number of points in the different sectors:" << endl;
	cout << "Sector -1 # = " << lclusterdata.nsectm1 << " Sector 0 # = " << lclusterdata.nsect0 << " Sector 1 # = " << lclusterdata.nsectp1 << ", Sectors cut = " << opt.Nspace - lclusterdata.nsectm1 - lclusterdata.nsect0 - lclusterdata.nsectp1 << endl << endl;

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

	if(opt.doboxes){
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
		f3d << opt.leng1 << " " << opt.leng2 << " " << opt.leng3 << " " << opt.leng4 << endl;
		f3d << lclusterdata.clustermembers.size() << endl;
		for(int i1=0;i1<opt.leng1;i1++)
		for(int i2=0;i2<opt.leng2;i2++)
		for(int i3=0;i3<opt.leng3;i3++){
			is = i1 + i2*opt.leng1 + i3*opt.leng1*opt.leng2;

			f3d << i1 << " " << i2 << " " << i3 << " " << lclusterdata.isinsortedcluster[is] << " " << lclusterdata.isinsector[is] << endl;
		}
		f3d.close();
	}else{
		cout << "WARNING: Could not open 3dcluster file!" << endl;
	}
}

#endif // FINDCLUSTER_WRITE_HPP
