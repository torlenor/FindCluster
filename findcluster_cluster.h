/*
 * findcluster_cluster.h - Cluster identification - headers
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

#ifndef FINDCLUSTER_CLUSTER_H
#define FINDCLUSTER_CLUSTER_H

#include <complex>
#include <vector>

#include "findcluster.h"

void fillSectorsAlt(Clusterstruct &lclusterdata, std::vector<std::vector<std::complex<double> > > &pollev, Options opt, double r);
void fillSectors(Clusterstruct &lclusterdata, std::vector<std::vector<std::complex<double> > > &pollev, Options opt, double delta);

void findClusters(Clusterstruct &lclusterdata, std::vector<std::vector<int> > &neib, Options opt);
void checkClusters(Clusterstruct &lclusterdata, Options opt);

void findPercolatingCluster(Clusterstruct &lclusterdata, Options opt);

void sortClusterSize(Clusterstruct &lclusterdata, Options opt);

#endif // FINDCLUSTER_CLUSTER_H