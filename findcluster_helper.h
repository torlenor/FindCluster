#ifndef FINDCLUSTER_HELPER_H
#define FINDCLUSTER_HELPER_H

#include <vector>

void fillNeib(const Options &opt, std::vector<std::vector<int> > &neib);

int latmap(const int i1, const int i2, const int i3, const Options &opt);

void Printsettings(const Options &opt);

void freeMem(Clusterstruct &lclusterdata);

void getCoords(const Options &opt, const int is, int &i1, int &i2, int &i3);

#endif // FINDCLUSTER_HELPER_H
