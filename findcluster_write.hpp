#ifndef FINDCLUSTER_WRITE_HPP
#define FINDCLUSTER_WRITE_HPP

void writeclustersize() {
  std::stringstream fclustersizename;
	fclustersizename << "clustersize_" << opt.Ns << "x" << opt.Nt << "_f" << opt.fraction << ".res";

  std::ofstream fclustersize;
	fclustersize.open(fclustersizename.str().c_str());

	fclustersize << "# 1_Nt 2_largestclusterweight 3_largestclusterweighterr 4_avgclusterweight 5_avgclusterweighterr 6_avgfortunatoclustersize 7_avgfortunatoclustersizeerr 8_largestnonpercclusterweight 9_largestnonpercclusterweighterr 10_largestclusterradius 11_largestclusterradiuserr 12_rootmeansquaredistanceR 13_rootmeansquaredistanceR" << std::endl;
	fclustersize.flags (std::ios::scientific);
	fclustersize.precision(numeric_limits<double>::digits10 + 1);
	fclustersize << opt.Nt << " " 
	  << results.maxclustersize << " " << results.maxclustersizeerr << " " 
	  << results.avgclustersize << " " << results.avgclustersizeerr << " " 
	  << results.avgclusersizeFortunato << " " << results.avgclusersizeFortunatoerr << " " 
	  << results.maxnonpercclustersize << " " << results.maxnonpercclustersizeerr << " " 
	  << results.largestclusterradius << " " << results.largestclusterradiuserr << " " 
	  << results.avgrootmeansquaredistance << " " << results.avgrootmeansquaredistanceerr 
	<< std::endl;

	fclustersize.close();
}

void writeavgclustersize() {
  std::stringstream favgclustersizename;
	favgclustersizename << "clustersize_avg_" << opt.Ns << "x" << opt.Nt << "_f" << opt.fraction << ".res";

  std::ofstream favgclustersize;
	favgclustersize.open(favgclustersizename.str().c_str());

	favgclustersize << "# 1_Nt 2_avgclusterweight 3_avgclusterweighterr 4_avgfortunatoclustersize 5_avgfortunatoclustersizeerr 6_avgclusterweightFnp 7_avgclusterweightFnperr 8_avgclusterweightnp 9_avgclusterweightnperr 10_avgnpclusterradius 11_avgnpclusterradiuserr 12_avgnpclusterarea 13_avgnpclusterareaerr" << std::endl;
	favgclustersize.flags (std::ios::scientific);
	favgclustersize.precision(numeric_limits<double>::digits10 + 1);
	favgclustersize << opt.Nt << " " 
    << results.avgclustersize << " " << results.avgclustersizeerr << " " 
    << results.avgclusersizeFortunato << " " << results.avgclusersizeFortunatoerr << " " 
    << results.avgclustersizenp << " " << results.avgclustersizenperr << " " 
    << results.avgclustersizeFnp << " " << results.avgclustersizeFnperr << " " 
    << results.avgnpclusterradius << " " << results.avgnpclusterradiuserr << " "
    << results.avgnonpercperimeter << " " << results.avgnonpercperimetererr <<  " "
	<< std::endl;

	favgclustersize.close();
}

void writeclusterradius() {
  std::stringstream fclusterradiusname;
	fclusterradiusname << "clusterradius_" << opt.Ns << "x" << opt.Nt << "_f" << opt.fraction << ".res";

  std::ofstream fclusterradius;
	fclusterradius.open(fclusterradiusname.str().c_str());

	fclusterradius << "# 1_Nt 2_largestclusterradius 3_largestclusterradiuserr 4_largestnpclusterradius 5_largestnpclusterradiuserr 6_avgclusterradius 7_avgclusterradiuserr 8_avgnpcluterradius 9_avgnpclusterradiuserr 10_rootmeansquaredistanceR 11_rootmeansquaredistanceRerr 12_maxclustersize 13_maxclustersizeerr 14_maxnonpercclustersize 15_maxnonpercclustersizeerr 16_avgclustersize 17_avgclustersizeerr 18_avgclusersizeFortunato 19_avgclusersizeFortunatoerr 20_avgclustersizeFnp 21_avgclustersizeFnperr" << std::endl;
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
	<< std::endl;

	fclusterradius.close();
}

void writenpercc() {
  std::stringstream fnperccname;
	fnperccname << "npercc_" << opt.Ns << "x" << opt.Nt << "_f" << opt.fraction << ".res";

  std::ofstream fnpercc;
	fnpercc.open(fnperccname.str().c_str());

	fnpercc << "# 1_Nt 2_percclusters 3_percclusterserr 4_cut 5_cuterr" << std::endl;
	fnpercc.flags (std::ios::scientific);
	fnpercc.precision(numeric_limits<double>::digits10 + 1);
	fnpercc << opt.Nt << " " << results.avgperccluster << " " << results.avgpercclustererr << " " 
    << results.cut << " " << results.cuterr 
  << std::endl;

	fnpercc.close();
}

void writemeanfreepathnew() {
  std::stringstream ffreepathnewname;
	ffreepathnewname << "meanfreepath_new_" << opt.Ns << "x" << opt.Nt << "_f" << opt.fraction << ".res";

  std::ofstream ffreepathnew;
	ffreepathnew.open(ffreepathnewname.str().c_str());

	ffreepathnew << "# 1_Nt 2_largestcluster 3_largestclustererr 4_largestnpcluster 5_largestnpclustererr 6_avgcluster 7_avgclustererr 8_avgnpcluster 9_avgnpclustererr 10_avgFcluster 11_avgFclustererr 12_avgFnpcluster 13_avgFnpclustererr" << std::endl;
	ffreepathnew.flags (std::ios::scientific);
	ffreepathnew.precision(numeric_limits<double>::digits10 + 1);
	ffreepathnew << opt.Nt << " " 
    << results.largestclustermeanfreepathnew << " " << results.largestclustermeanfreepathnewerr << " " 
    << results.largestnpclustermeanfreepathnew << " " << results.largestnpclustermeanfreepathnewerr << " " 
    << results.avgclustermeanfreepathnew << " " << results.avgclustermeanfreepathnewerr << " " 
    << results.avgnpclustermeanfreepathnew << " " << results.avgnpclustermeanfreepathnewerr << " " 
    << results.avgFclustermeanfreepathnew << " " << results.avgFclustermeanfreepathnewerr << " " 
    << results.avgFnpclustermeanfreepathnew << " " << results.avgFnpclustermeanfreepathnewerr
  << std::endl;

	ffreepathnew.close();
}

void writebox() {
    std::stringstream fboxcntname;
		fboxcntname << "boxcnt_" << opt.Ns << "x" << opt.Nt << "_f" << opt.fraction << ".res";

    std::ofstream fboxcnt;
		fboxcnt.open(fboxcntname.str().c_str());

		fboxcnt << "# 1_boxsize 2_boxcnt 3_boxcnterr" << std::endl;
		fboxcnt.flags (std::ios::scientific);
		fboxcnt.precision(numeric_limits<double>::digits10 + 1);
		for (unsigned int size=0; size<boxsize.size(); size++) {
			fboxcnt << boxsize[size] << " " << results.largestclusterboxcount[size] << " " << results.largestclusterboxcounterr[size] << std::endl;
		}

		fboxcnt.close();
}

void writeboxnp() {
    std::stringstream fboxcntname;
		fboxcntname << "boxcnt_nonpercc_" << opt.Ns << "x" << opt.Nt << "_f" << opt.fraction << ".res";

    std::ofstream fboxcnt;
		fboxcnt.open(fboxcntname.str().c_str());

		fboxcnt << "# 1_boxsize 2_boxcntnonpercc 3_boxcntnonperccerr" << std::endl;
		fboxcnt.flags (std::ios::scientific);
		fboxcnt.precision(numeric_limits<double>::digits10 + 1);
		for(unsigned int size=0; size<boxsize.size(); size++){
			fboxcnt << boxsize[size] << " " << results.largestnonpercboxcount[size] << " " << results.largestnonpercboxcounterr[size] << std::endl;
		}

		fboxcnt.close();
}

void writearea() {
  std::stringstream fareaname;
	fareaname << "area_" << opt.Ns << "x" << opt.Nt << "_f" << opt.fraction << ".res";

  std::ofstream farea;
	farea.open(fareaname.str().c_str());

	farea << "# 1_Nt 2_area 3_areaerr 4_arealargestnonpercc 5_arealargestnonperccerr 6_largestclusterweight 7_largestclusterweighterr 8_largestnonpercclusterweight 9_largestnonpercclusterweighterr 10_areaavgnonpercc 11_areaavgnonperccerr 12_avgclusersizeFnp 13_avgclusersizeFnperr" << std::endl;
	farea.flags (std::ios::scientific);
	farea.precision(numeric_limits<double>::digits10 + 1);
	farea << opt.Nt << " " 
    << results.totalperimeter << " " << results.totalperimetererr << " " 
    << results.largestnonpercperimeter << " " << results.largestnonpercperimetererr << " " 
    << results.maxclustersize << " " << results.maxclustersizeerr << " " 
    << results.maxnonpercclustersize << " " << results.maxnonpercclustersizeerr <<  " " 
    << results.avgnonpercperimeter << " " << results.avgnonpercperimetererr <<  " "
    << results.avgclustersizeFnp << " " << results.avgclustersizeFnperr << " "
	<< std::endl;

	farea.close();
}

void writepoll() {
	stringstream fpollname;
	fpollname << "poll_" << opt.Ns << "x" << opt.Nt << "_f" << opt.fraction << ".res";

	ofstream fpoll;
	fpoll.open(fpollname.str().c_str());

	fpoll << "# 1_Nt 2_poll 3_pollerr" << std::endl;
	fpoll.flags (std::ios::scientific);
	fpoll.precision(numeric_limits<double>::digits10 + 1);
	fpoll << opt.Nt << " " 
    << results.polyakovloopaftercut << " " << results.polyakovloopaftercuterr 
  << std::endl;

	fpoll.close();
}

void writeresults() {
  writeclustersize();
  writeavgclustersize();
  writeclusterradius();

  writenpercc();

  writemeanfreepathnew();
	
	if(opt.doboxes){
    writebox();
    writeboxnp();
	}
	
  writearea();
  writepoll();
}

void writeresultsstdout() {
  std::cout << "Expectation values (single eliminitation jackknife): " << std::endl;
  std::cout << "Average cluster size = " << setprecision(14) << results.avgclustersize << ", Maximum cluster size / V = " << results.maxclustersize << std::endl;
	std::cout << "Average cluster err  = " << results.avgclustersizeerr << ", Maximum cluster / V err  = " << results.maxclustersizeerr << std::endl;
  std::cout << "Average cluster size Fortunato (1.7) = " << results.avgclusersizeFortunato << ", Error  = " << results.avgclusersizeFortunatoerr << std::endl;
  std::cout << "Radius of largest cluster = " << results.largestclusterradius << ", Error = " << results.largestclusterradiuserr << std::endl;
	if(opt.dodistance){
    std::cout << "Root mean distance traveled R = " << results.avgrootmeansquaredistance << ", Error = " << results.avgrootmeansquaredistanceerr << std::endl;
	}
  std::cout << "Cut = " << results.cut << " Cut err = " << results.cuterr << std::endl;
  std::cout << "Laserdim = " << results.totalperimeter << " Laserdim err = " << results.totalperimetererr << std::endl;
  std::cout << "Polyakov loop = " << results.polyakovloopaftercut << " Polyakov loop err = " << results.polyakovloopaftercuterr << endl << std::endl;
}

void writeConfigResultsstdout(Observablestruct &lobs, Clusterstruct &lclusterdata){
	std::cout << "Number of points in the different sectors:" << std::endl;
	std::cout << "Sector -1 # = " << lclusterdata.nsectm1 << " Sector 0 # = " << lclusterdata.nsect0 << " Sector 1 # = " << lclusterdata.nsectp1 << ", Sectors cut = " << opt.Nspace - lclusterdata.nsectm1 - lclusterdata.nsect0 - lclusterdata.nsectp1 << std::endl << std::endl;

	std::cout  << lclusterdata.percolatingclusters.size() << " of " << lclusterdata.clustermembers.size() << " clusters are percolating!" << std::endl;
	for(unsigned int c=0;c<lclusterdata.percolatingclusters.size();c++){
		std::cout << "Cluster " << lclusterdata.percolatingclusters[c] << " is percolating in directions ("
			<< lclusterdata.percolatingdirections[c][0] << "," 
			<< lclusterdata.percolatingdirections[c][1] << ","
			<< lclusterdata.percolatingdirections[c][2] << ")!" << std::endl;
		std::cout << "The cluster has " << lclusterdata.clustermembers[lclusterdata.percolatingclusters[c]].size() << " members and is in sector " << lclusterdata.clustersector[lclusterdata.percolatingclusters[c]] << " !" << std::endl;
	}
	
	std::cout << std::endl;
	std::cout << "Average cluster size = " << lobs.avgclustersize << std::endl;
	std::cout << "Average cluster size Fortunato = " << lobs.avgclustersizeF << std::endl;
	std::cout << "Largest cluster is cluster " << lobs.maxclusterid << " with " << lobs.maxclustersize << " members." << std::endl;

	if(opt.doboxes){
		std::cout << std::endl << "Boxes (Clusters sorted based on # of members):" << std::endl;
		for(unsigned int c=0; c<lclusterdata.sortedcluster.size(); c++){
			std::cout << "Cluster = " << lclusterdata.sortedcluster[c] << " with " << lclusterdata.clustermembers[lclusterdata.sortedcluster[c]].size() << " members:" << std::endl;
			for(unsigned int size=0; size<boxsize.size(); size++){
			std::cout << lobs.numberofboxes[lclusterdata.sortedcluster[c]][size] << " boxes of size " << boxsize[size] << " needed." << std::endl;
			}
		}
	}

	std::cout << std::endl << "The following clusters are connected over the PBCs: " << std::endl;
	for(unsigned int c=0;c<lclusterdata.clusterisperiodic.size();c++){
		if(lclusterdata.clusterisperiodic[lclusterdata.sortedcluster[c]][0]==1 
		|| lclusterdata.clusterisperiodic[lclusterdata.sortedcluster[c]][1]==1 
		|| lclusterdata.clusterisperiodic[lclusterdata.sortedcluster[c]][2]==1){
			std::cout << "Cluster (size sorted) " << c << ", real id = " << lclusterdata.sortedcluster[c] << std::endl;
		}
	}
}

void writeClusterList(Clusterstruct &lclusterdata) {
  std::vector<int> csize;
	for (unsigned c=0; c<lclusterdata.clustermembers.size(); c++) {
		csize.push_back(lclusterdata.clustermembers[c].size());
	}

	sort(csize.begin(), csize.end());

	std::cout << "Cluster sizes" << std::endl;
	for (unsigned c=0; c<lclusterdata.clustermembers.size(); c++) {
		std::cout << csize[c] << std::endl;
	}
}

void cluster3doutput(Clusterstruct &lclusterdata, string f3dname) {
	int is;

  std::ofstream f3d;
	f3d.open(f3dname.c_str());

	if (f3d.is_open()) {
		f3d << opt.leng1 << " " << opt.leng2 << " " << opt.leng3 << " " << opt.leng4 << std::endl;
		f3d << lclusterdata.clustermembers.size() << std::endl;
		for (int i1=0; i1<opt.leng1; i1++)
		for (int i2=0; i2<opt.leng2; i2++)
		for (int i3=0; i3<opt.leng3; i3++) {
			is = i1 + i2*opt.leng1 + i3*opt.leng1*opt.leng2;

			f3d << i1 << " " << i2 << " " << i3 << " " << lclusterdata.isinsortedcluster[is] << " " << lclusterdata.isinsector[is] << std::endl;
		}
		f3d.close();
	}else{
		std::cout << "WARNING: Could not open 3dcluster file!" << std::endl;
	}
}

#endif // FINDCLUSTER_WRITE_HPP
