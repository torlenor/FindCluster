/*
 * findcluster_write.h - write results to file and stdout - headers
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
