/*
 * findcluster_measurewrite.cpp - Write out measurements
 *
 * Copyright © 2014 H.-P. Schadler  <hanspeter.schadler@uni-graz.at>
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 */

#include "findcluster_writemeasure.h"

#include <limits>
#include <fstream>

#include "findcluster.h"

void writemeasclustersize(Observablestruct *obs, Options &opt) {
  std::stringstream fclustersizename;
	fclustersizename << "meas_clustersize_" << opt.Ns << "x" << opt.Nt << "_f" << opt.fraction << ".res";

  std::ofstream fclustersize;
	fclustersize.open(fclustersizename.str().c_str());

	fclustersize << "# 0_meas 1_Nt 2_largestclusterweight 3_largestclusterweighterr 4_avgclusterweight 5_avgclusterweighterr 6_avgfortunatoclustersize 7_avgfortunatoclustersizeerr 8_largestnonpercclusterweight 9_largestnonpercclusterweighterr 10_largestclusterradius 11_largestclusterradiuserr 12_rootmeansquaredistanceR 13_rootmeansquaredistanceR" << std::endl;
	fclustersize.flags (std::ios::scientific);
	fclustersize.precision(std::numeric_limits<double>::digits10 + 1);
	
	// Loop over all measurements
	for (int m=0; m<opt.nmeas; m++) {
		fclustersize << m << " " << opt.Nt << " " 
			<< (&obs[m])->maxclustersize/(double)opt.Nspace << " " 
			<< (&obs[m])->avgclustersize << " " 
			<< (&obs[m])->avgclustersizeF << " " 
			<< (&obs[m])->largestnonpercclustersize/(double)opt.Nspace << " " 
			<< (&obs[m])->largestclusterradius << " " 
			<< (&obs[m])->rootmeansquaredistanceR << " " 
		<< "\n";
	} // End loop m
	
	fclustersize.close();
}

/*
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
*/
void writeMeasures(Observablestruct *obs, Options &opt) {
  writemeasclustersize(obs, opt);
  /* writemeasavgclustersize();
  writemeasclusterradius();

  writemeasnpercc();

  writemeasmeanfreepathnew();
	
	if(opt.doboxes){
    writemeasbox();
    writemeasboxnp();
	}
	
  writemeasarea();
  writemeaspoll(); */
}
