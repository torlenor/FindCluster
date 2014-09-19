/*
 * findcluster_box.h - Box counting calculations - headers
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

#ifndef FINDCLUSTER_BOX_H
#define FINDCLUSTER_BOX_H

#include <vector>

#include "findcluster.h"

void ObsBoxesOnlyLargest(Observablestruct &lobs, Clusterstruct &lclusterdata, const Options &opt, std::vector<int> &boxsize, std::vector<int> &boxes);
void ObsBoxes(Observablestruct &lobs, Clusterstruct &lclusterdata, const Options &opt, std::vector<int> &boxsize, std::vector<int> &boxes);

#endif // FINDCLUSTER_BOX_H
