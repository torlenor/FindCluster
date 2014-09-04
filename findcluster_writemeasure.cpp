/*
 * findcluster_writemeasure.cpp - Write out measurements
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

void prepwriteMeasure(Options &opt) {
  // Prepare file for writemeascluster
  std::stringstream fclustersizename;
	fclustersizename << "meas_cluster_" << opt.Ns << "x" << opt.Nt << "_f" << opt.fraction << ".res";
  std::ofstream fclustersize;
	fclustersize.open(fclustersizename.str().c_str());
	fclustersize << "# 0_meas opt.maxclustersize opt.largestnonpercclustersize opt.largestnonpercclusterid opt.largetsnonpercclustersector opt.avgclustersize opt.avgclustersizeF opt.rootmeansquaredistanceR opt.avgclustersizenp opt.avgclustersizeFnp opt.cut opt.largestclusterradius opt.largestnpclusterradius opt.avgclusterradius opt.avgnpclusterradius opt.percc opt.area opt.arealargestnonperccluster opt.areaavgnonperccluster opt.poll opt.largestclustermeanfreepathnew opt.largestnpclustermeanfreepathnew opt.avgclustermeanfreepathnew opt.avgnpclustermeanfreepathnew opt.avgFclustermeanfreepathnew opt.avgFnpclustermeanfreepathnew" << std::endl;
	fclustersize.close();
}

void writemeascluster(Observablestruct &lobs, Options &opt, int m) {
  std::stringstream fclustersizename;
	fclustersizename << "meas_cluster_" << opt.Ns << "x" << opt.Nt << "_f" << opt.fraction << ".res";

  std::ofstream fclustersize;
	fclustersize.open(fclustersizename.str().c_str(), std::ofstream::out | std::ofstream::app);
	
  fclustersize.flags (std::ios::scientific);
	fclustersize.precision(std::numeric_limits<double>::digits10 + 1);
	
  fclustersize << m << " " <<  
    lobs.maxclustersize << " " << 
    lobs.largestnonpercclustersize << " " << 
    lobs.largestnonpercclusterid << " " << 
    lobs.largetsnonpercclustersector << " " << 
    lobs.avgclustersize << " " << 
    lobs.avgclustersizeF << " " << 
    lobs.rootmeansquaredistanceR << " " << 
    lobs.avgclustersizenp << " " << 
    lobs.avgclustersizeFnp << " " << 
    lobs.cut << " " << 
    lobs.largestclusterradius << " " << 
    lobs.largestnpclusterradius << " " << 
    lobs.avgclusterradius << " " << 
    lobs.avgnpclusterradius << " " << 
    lobs.percc << " " << 
    lobs.area << " " << 
    lobs.arealargestnonperccluster << " " << 
    lobs.areaavgnonperccluster << " " << 
    lobs.poll << " " << 
    lobs.largestclustermeanfreepathnew << " " << 
    lobs.largestnpclustermeanfreepathnew << " " << 
    lobs.avgclustermeanfreepathnew << " " << 
    lobs.avgnpclustermeanfreepathnew << " " << 
    lobs.avgFclustermeanfreepathnew << " " << 
    lobs.avgFnpclustermeanfreepathnew 
  << "\n";
	
	fclustersize.close();
}

void writeMeasures(Observablestruct &lobs, Options &opt, int m) {
  writemeascluster(lobs, opt, m);
}
