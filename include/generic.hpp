#ifndef GENERIC_HPP
#define GENERIC_HPP

#include <iostream>
#include <fstream>
#include <complex>
#include <vector>
#include <string>

int writePollEv(string fevname);
int writePollEvBinary(string fevname);
void checkPollEv();

void checkPollEv(int leng1, int leng2, int leng3, int leng4, int matrixdim, vector<vector<complex<double> > > &pollev){
	cout << "Performing some checks on the evs... " << flush;

	complex<double> sum;
	double eps=1E-12;
	
	double absval=0;

	for(int is=0;is<leng1*leng2*leng3;is++)
		for(int i=0;i<3;i++){
			absval = abs(pollev[is][i]);
			sum += absval;
			if(abs(absval - 1.0) > eps){
				cout << endl << "WARNING: abs(abs(ev) - 1.0) = " << abs(absval - 1.0) << " > eps"  << endl;
			}
		}

	sum = sum/(double)(3*leng1*leng2*leng3);

	if(abs(abs(sum) - 1.0) > eps)
		cout << endl << "WARNING: Sum over all abs(ev)'s / 3Ns*Ns*Ns = " << abs(sum) << " > 1"  << endl;
	
	cout << "done!" << endl;
}

int readPollEv(int leng1, int leng2, int leng3, int leng4, int matrixdim, vector<vector<complex<double> > > &pollev, string fevname){
	ifstream fev;
	fev.open(fevname.c_str());
	
        if(fev.is_open()!=true){
                cout  << "ERROR: File " << fevname <<  " to write Polyakov loop eigenvalues could not be opened!" << endl;
		throw 1;
        }

	string strtmp;
	int lineskip=2;

	for(int j=0;j<lineskip;j++){
		getline(fev, strtmp);
	}

	int itmp;
	for(int is=0;is<leng1*leng2*leng3;is++){
		fev >> itmp >> real(pollev[is][0]) >> imag(pollev[is][0]) >> real(pollev[is][1]) >>
			imag(pollev[is][1]) >> real(pollev[is][2]) >> imag(pollev[is][2]);
	}

	fev.close();

	return 0;
}

int readPollEvBinary(int leng1, int leng2, int leng3, int leng4, int matrixdim, vector<vector<complex<double> > > &pollev, string fevname){
	cout << "Reading 3d lattice with Polyakov loop evs... " << flush;
	
	int elems=0;
       	FILE* pFile;
       	
       	int cleng1, cleng2, cleng3, cleng4;

	unsigned int nindex=leng1*leng2*leng3*matrixdim;

	int is;

	if( nindex != pollev.size() * pollev[pollev.size()-1].size() ){
		cout << "ERROR: Something wrong in readPollEvBinary!" << endl;
		cout << "Want to read nindex = " << nindex << " entries, but has space for only " << pollev.size() * pollev[pollev.size()-1].size() << " entries!" << endl;

	}
	
       	pFile = fopen(fevname.c_str(), "rb");
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
		
		for(int i1=0;i1<leng1;i1++)
		for(int i2=0;i2<leng2;i2++)
		for(int i3=0;i3<leng3;i3++){
			is = i1 + i2*leng1 + i3*leng1*leng2;
			elems += fread(&pollev[is][0], sizeof(std::complex<double>), matrixdim, pFile);
		}

		fclose(pFile);
	}
	
	cout << "done!" << endl;
	
        return nindex + 4 - elems;
}

#endif // GENERIC_HPP
