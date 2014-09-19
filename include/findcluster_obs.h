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
