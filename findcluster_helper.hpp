#ifndef FINDCLUSTER_HELPER_HPP
#define FINDCLUSTER_HELPER_HPP

void freeMem(Clusterstruct &lclusterdata){
	lclusterdata.isinsector.resize(0);
	lclusterdata.clustersector.resize(0);
	lclusterdata.isincluster.resize(0);

	for(unsigned c=0;c<lclusterdata.clustermembers.size();c++){
		lclusterdata.clustermembers[c].resize(0);
	}
	lclusterdata.clustermembers.resize(0);

	lclusterdata.percolatingclusters.resize(0);
	for(unsigned p=0; p<lclusterdata.percolatingdirections.size(); p++){
		lclusterdata.percolatingdirections[p].resize(0);
	}
	lclusterdata.percolatingdirections.resize(0);
    
	lclusterdata.sortedcluster.resize(0);
	lclusterdata.sortedrealcluster.resize(0);
	lclusterdata.isinsortedcluster.resize(0);
}

void getCoords(int is, int &i1, int &i2, int &i3){
	i1 = (is % (leng1*leng2) ) % leng1;
	i2 = (is % (leng1*leng2) ) / leng1;
	i3 = is / (leng1*leng2);

	if(is != i1 + i2*leng1 + i3*leng1*leng2)
		cout << "ERROR: Problem in getCoords!" << endl;
}

#endif // FINDCLUSTER_HELPER_HPP
