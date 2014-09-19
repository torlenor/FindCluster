/*
 * findcluster_helper.h - various functions - headers
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

#ifndef FINDCLUSTER_HELPER_H
#define FINDCLUSTER_HELPER_H

#include <vector>

#include "findcluster.h"

void fillNeib(const Options &opt, std::vector<std::vector<int> > &neib);

int latmap(const int i1, const int i2, const int i3, const Options &opt);

void Printsettings(const Options &opt);

void freeMem(Clusterstruct &lclusterdata);

void getCoords(const Options &opt, const int is, int &i1, int &i2, int &i3);

#endif // FINDCLUSTER_HELPER_H
