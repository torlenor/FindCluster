/*
 * findcluster_path.h - Mean free path calculations - headers
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

#ifndef FINDCLUSTER_PATH_H
#define FINDCLUSTER_PATH_H

#include "findcluster.h"

// void obsClusterMeanFreePathLargest(Observablestruct &lobs, Clusterstruct &lclusterdata, Options opt);
// void obsClusterMeanFreePath(Observablestruct &lobs, Clusterstruct &lclusterdata, Options opt);
void obsClusterMeanFreePathNew(Observablestruct &lobs, Clusterstruct &lclusterdata, Options opt);
void obsAverageMeanfreepathNew(Observablestruct &lobs, Clusterstruct &lclusterdata, Options opt);

#endif // FINDCLUSTER_PATH_H
