#ifndef LPOLLHANDLER_HPP
#define LPOLLHANDLER_HPP

// Local Polyakov loop write
int readBinaryLPoll(int leng1, int leng2, int leng3, int leng4, int matrixdim, string fconfigname, vector<vector<complex<double> > > &lpollarray){
	int elems=0;

       	FILE* pFile;

	int nindex=leng1*leng2*leng3*matrixdim*matrixdim;

	int cleng1, cleng2, cleng3, cleng4;
	complex<double> pollref, plaqref, poll, plaq;

        pFile = fopen(fconfigname.c_str(), "rb");
   	if (pFile == NULL) perror ("Error opening file");
	else{
		elems += fread(&cleng1, sizeof(int), 1, pFile);
		elems += fread(&cleng2, sizeof(int), 1, pFile);
		elems += fread(&cleng3, sizeof(int), 1, pFile);
		elems += fread(&cleng4, sizeof(int), 1, pFile);
		
		// Check lattice dimensions
		if(leng1!=cleng1 || leng2!=cleng2 || leng3!=cleng3 || leng4!=cleng4){
			cout << "ERROR: lattice dimensions are different in config file and settings !" << endl;
			return 1;
		}
		
		elems += fread(&pollref, sizeof(std::complex<double>),1, pFile);
		elems += fread(&plaqref, sizeof(std::complex<double>),1, pFile);

		int is=0;

		for(int i1=0;i1<leng1;i1++)
		for(int i2=0;i2<leng2;i2++)
		for(int i3=0;i3<leng3;i3++){
			is = i1 + i2*leng1 + i3*leng1*leng2;
			elems += fread(&lpollarray[is][0], sizeof(std::complex<double>), matrixdim*matrixdim, pFile);
		}

		fclose(pFile);
	}
        return nindex + 6 - elems;
}

complex<double> calcPoll(int leng1, int leng2, int leng3, int leng4, int matrixdim, vector<vector<complex<double> > > &lpollarray){
        // Calculates the Polyakov loop spatial average from local Polyakov loop lattice
        complex<double> poll(0,0);
        
	for(int is=0;is<leng1*leng2*leng3;is++){
                poll += lpollarray[is][0*3 + 0] + lpollarray[is][1*3 + 1] + lpollarray[is][2*3 + 2];
        }   

        poll = poll/((double)leng1*leng2*leng3);

        return poll;
}

#endif // LPOLLHANDLER_HPP
