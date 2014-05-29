#ifndef FINDCLUSTER_HELPER_HPP
#define FINDCLUSTER_HELPER_HPP

void fillNeib() {
	// Fills the neib array
  int i1p,i2p,i3p,i1m,i2m,i3m,is,isp1,isp2,isp3,ism1,ism2,ism3;
  for (int i1=0; i1<opt.leng1; i1++) {
    i1p = i1 + 1;
    i1m = i1 - 1;
    if (i1p == opt.leng1) i1p = 0;
    if (i1m == -1) i1m = opt.leng1-1;

    for (int i2=0; i2<opt.leng2; i2++) {
      i2p = i2 + 1;
      i2m = i2 - 1;
      if (i2p == opt.leng2) i2p = 0;
      if (i2m == -1) i2m = opt.leng2-1;

      for (int i3 = 0;i3<opt.leng3;i3++) {
        i3p = i3 + 1;
        i3m = i3 - 1;
        if (i3p == opt.leng3) i3p = 0;
        if (i3m == -1) i3m = opt.leng3-1;
        
        // Compute the site address and the addresses of the sites shifted
        // by one unit in each direction
        is = i1 + i2*opt.leng1 + i3*opt.leng1*opt.leng2;

        isp1 = i1p + i2*opt.leng1 + i3*opt.leng1*opt.leng2;
        isp2 = i1 + i2p*opt.leng1 + i3*opt.leng1*opt.leng2;
        isp3 = i1 + i2*opt.leng1 + i3p*opt.leng1*opt.leng2;

        ism1 = i1m + i2*opt.leng1 + i3*opt.leng1*opt.leng2;
        ism2 = i1 + i2m*opt.leng1 + i3*opt.leng1*opt.leng2;
        ism3 = i1 + i2*opt.leng1 + i3m*opt.leng1*opt.leng2;

        // Fill the neib array
        neib[is][0] = isp1;
        neib[is][1] = isp2;
        neib[is][2] = isp3;

        neib[is][3] = ism1;
        neib[is][4] = ism2;
        neib[is][5] = ism3;
      }
    }
  }
}

int latmap(int i1, int i2, int i3) {
	return i1 + i2*opt.leng1 + i3*opt.leng1*opt.leng2;
}

void printsettings() {
	cout << "Settings:" << endl << endl;
	cout << "Lattice size = " << opt.leng1 << "x" << opt.leng2 << "x" << opt.leng3 << "x" << opt.leng4 << endl;
	if (opt.wupperdata)
		cout << "Reading Wuppertal Polakov loop data format." <<  endl;
	cout << "Number of configurations = " << opt.nmeas << endl;
	if (! opt.usealternativesectors) {
		cout << "Cut fraction = " << opt.fraction << endl;
	} else {
		cout << "Alternative cut prescription, radius r = " << opt.r << endl;
	}
	if (opt.doboxes)
		cout << "Calculating 'box' observables." << endl;
	if (opt.detail)
		cout << "Writing detailed results for every configuration." <<  endl;
	if (opt.do3d)
		cout << "Writing 3dcluster visualization data files." <<  endl;
	if (opt.writemeas)
		cout << "Writing all measurements to file." <<  endl;
}

void freeMem(Clusterstruct &lclusterdata) {
	lclusterdata.isinsector.resize(0);
	lclusterdata.clustersector.resize(0);
	lclusterdata.isincluster.resize(0);

	for (unsigned c=0;c<lclusterdata.clustermembers.size();c++) {
		lclusterdata.clustermembers[c].resize(0);
	}
	lclusterdata.clustermembers.resize(0);

	lclusterdata.percolatingclusters.resize(0);
	for (unsigned p=0; p<lclusterdata.percolatingdirections.size(); p++) {
		lclusterdata.percolatingdirections[p].resize(0);
	}
	lclusterdata.percolatingdirections.resize(0);
    
	lclusterdata.sortedcluster.resize(0);
	lclusterdata.sortedrealcluster.resize(0);
	lclusterdata.isinsortedcluster.resize(0);
}

void getCoords(int is, int &i1, int &i2, int &i3) {
	i1 = (is % (opt.leng1*opt.leng2) ) % opt.leng1;
	i2 = (is % (opt.leng1*opt.leng2) ) / opt.leng1;
	i3 = is / (opt.leng1*opt.leng2);

	if (is != i1 + i2*opt.leng1 + i3*opt.leng1*opt.leng2)
		cout << "ERROR: Problem in getCoords!" << endl;
}

#endif // FINDCLUSTER_HELPER_HPP
