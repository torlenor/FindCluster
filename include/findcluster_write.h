#ifndef FINDCLUSTER_WRITE_H
#define FINDCLUSTER_WRITE_H

#include <iostream>

#include "findcluster.h"

void writeclustersize(const Resultstruct &results, const Options &opt);

void writeavgclustersize(const Resultstruct &results, const Options &opt);

void writeclusterradius(const Resultstruct &results, const Options &opt);

void writenpercc(const Resultstruct &results, const Options &opt);

void writemeanfreepathnew(const Resultstruct &results, const Options &opt);

void writearea(const Resultstruct &results, const Options &opt);

void writepoll(const Resultstruct &results, const Options &opt);

void writebox(const Resultstruct &result, std::vector<int> &boxsize, std::vector<int> &boxes, const Options &opt);

void writeboxnp(const Resultstruct &results, std::vector<int> &boxsize, std::vector<int> &boxes, const Options &opt);

void writeresults(std::vector<int> &boxsize, std::vector<int> &boxes, const Resultstruct &results, const Options &opt);

void writeresultsstdout(const Resultstruct &results, const Options &opt);

void writeresultsstdout_singleconf(Observablestruct &lobs, const Options &opt);

void writeConfigResultsstdout(std::vector<int> &boxsize, std::vector<int> &boxes, Observablestruct &lobs, Clusterstruct &lclusterdata, const Options &opt);

void writeClusterList(Clusterstruct &lclusterdata);

void cluster3doutput(Clusterstruct &lclusterdata, std::string f3dname, const Options &opt);

#endif // FINDCLUSTER_WRITE_H
