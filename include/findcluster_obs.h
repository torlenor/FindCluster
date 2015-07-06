/*
 * findcluster_obs.h - observable calculations - headers
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

#ifndef FINDCLUSTER_OBS_H
#define FINDCLUSTER_OBS_H

#include "findcluster.h"
#include "findcluster_box.h"
#include "findcluster_path.h"
#include "findcluster_radius.h"

void ObsLargestCluster(Observablestruct &lobs, Clusterstruct &lclusterdata);

void ObsAverageClusterSize(Observablestruct &lobs, Clusterstruct &lclusterdata);

void ObsAverageClusterSizeFortunato(Observablestruct &lobs, Clusterstruct &lclusterdata);

void ObsAverageClusterSizeNoPercc(Observablestruct &lobs, Clusterstruct &lclusterdata);

void ObsCutPercentage(Observablestruct &lobs, Clusterstruct &lclusterdata);

void ObsArea(Observablestruct &lobs, Clusterstruct &lclusterdata);

void ObsAreaLargestNonPercCluster(Observablestruct &lobs, Clusterstruct &lclusterdata);

void ObsAreaAvgNonPercCluster(Observablestruct &lobs, Clusterstruct &lclusterdata);

void ObsNumberOfPercClusters(Observablestruct &lobs, Clusterstruct &lclusterdata);

void ObsPollAfterCut(Observablestruct &lobs, Clusterstruct &lclusterdata);

void ObsRootMeanSquareDistance(Observablestruct &lobs, Clusterstruct &lclusterdata);

void CalcObservables(Observablestruct &lobs, Clusterstruct &lclusterdata, std::vector<int> &boxsize, std::vector<int> &boxes);

void CalcExp(const std::vector<Observablestruct> &obs, Resultstruct &results, Options &opt, std::vector<int> &boxsize, std::vector<int> &boxes);

#endif // FINDCLUSTER_OBS_H
