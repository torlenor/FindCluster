/* genpollev_from_lpoll.cpp
   Calculates Polyakov loop eigenvalues from a local Polyakov loop configuration

   v0.0.0 - 2013-03-05 Hans-Peter Schadler
*/

#include <iostream>
#include <fstream>
#include <complex>
#include <vector>

#include <gsl/gsl_histogram.h>
#include <gsl/gsl_histogram2d.h>
#include <gsl/gsl_errno.h>

#include "include/ConfigData.hpp"

using namespace std;

#include "include/matrixfunctions.hpp"
#include "include/lpollhandler.hpp"

#include "include/generic.hpp"

extern "C" {
	// Eigenvalues for a general complex matrix
	void zgeev_(const char *JOBVBL, const char* JOBVR, int *N, complex<double>* A, int* lda, 
			complex<double>* EVAL, complex<double>* VL, int *ldvl, complex<double>* VR, 
			int *ldvr, complex<double>* WORK, int *lwork, double *RWORK, int* INFO);
}

int Ns=0;
int Nt=0;

vector<vector<complex<double> > > pollev;

void localPollEv(vector<vector<complex<double> > > &lpollarray, complex<double> *polleval, int i1, int i2, int i3);
void fillPollEv(vector<vector<complex<double> > > &lpollarray);

int matrixdim=3, leng1=0, leng2=0, leng3=0, leng4=0;

void gslhisto(string fhistofilename, int nbinsx, int nbinsy, double *datax, double *datay, int Ndatax, int Ndatay);

int main(int argc, char *argv[]){
	if(argc<5){
		cout << "./genpollev_from_lpoll.x Ns Nt lpollconfigin pollevalout" << endl;
		return 1;
	}

	Ns=atoi(argv[1]);
	Nt=atoi(argv[2]);
	string fconfigname=argv[3];
	string fevname=argv[4];

	leng1=Ns; leng2=Ns; leng3=Ns; leng4=Nt;

	vector<vector<complex<double> > > lpollarray;
	
	lpollarray.resize(leng1*leng2*leng3);
	for(int is=0;is<leng1*leng2*leng3;is++)
		lpollarray[is].resize(matrixdim*matrixdim);

	if(readBinaryLPoll(leng1, leng2, leng3, leng4, matrixdim, fconfigname,lpollarray) != 0){
		cout << "ERROR: Error reading the file!" << endl;
		return 1;
	}

	complex<double> pollTr;

	pollTr = calcPoll(Ns, Ns, Ns, Nt, 3, lpollarray);
	cout << "Total Polyakov loop P = " << pollTr << endl;

	cout << "Allocating 3d (Ns^3 + 3) lattice for Polyakov loop ev's... " << flush;
	pollev.resize(Ns*Ns*Ns);
	cout << "done!" << endl;

	cout << "Filling 3d lattice with Polyakov loop ev's... " << flush;
	fillPollEv(lpollarray);
	cout << "done!" << endl;

	checkPollEv(leng1, leng2, leng3, leng4, matrixdim, pollev);

	complex<double> totpoll=0;
	for(int is=0;is<leng1*leng2*leng3;is++){
		totpoll += pollev[is][0] + pollev[is][1] + pollev[is][2];
	}
	totpoll = totpoll/(double)(leng1*leng2*leng3);
	cout << "Total Polyakov loop from pollev P = " << totpoll << endl;

	double delta = 1e-10;
	double err = abs(pollTr - totpoll);
	cout << "Error = " << err << endl;
	if(err > delta){
		cout << "WARNING: Error in Polyakov loop large!!!" << endl;
	}

	cout << "Writing 3d lattice to disk... " << flush;
	//writePollEv(fevname);
	if(writePollEvBinary(fevname) != 0){
		cout << "ERROR: Problems with writePollEvBinary !" << endl;
	}
	cout << "done!" << endl;

	return 0;
}

int writePollEv(string fevname){
	ofstream fev;
	fev.open(fevname.c_str());
	
	fev << "# is real(ev_1) imag(ev_1) real(ev_2) imag(ev_2) real(ev_3) imag(ev_3)" << endl;
	fev << "# is = i1 + i2*Ns + i3*Ns*Ns" << endl;
       	fev.flags (std::ios::scientific);
	fev.precision(numeric_limits<double>::digits10 + 1);

        if(fev.is_open()!=true){
                cout  << "ERROR: File " << fevname <<  " to write Polyakov loop eigenvalues could not be opened!" << endl;
		throw 1;
        }

	for(int is=0;is<Ns*Ns*Ns;is++){
		fev << is << " " << real(pollev[is][0]) << " " << imag(pollev[is][0]) << " "
				<< real(pollev[is][1]) << " " << imag(pollev[is][1]) << " "
				<< real(pollev[is][2]) << " " << imag(pollev[is][2]) << endl;
	}

	fev.close();

	return 0;
}

int writePollEvBinary(string fevname){
	int elems=0;
       	FILE* pFile;

	int nindex=Ns*Ns*Ns*matrixdim;

	int is;

	if( nindex != pollev.size() * pollev[pollev.size()-1].size() ){
		cout << "ERROR: Something wrong in writePollEvBinary!" << endl;
		cout << "Want to write nindex = " << nindex << " entries, but there are only " << pollev.size() * pollev[pollev.size()-1].size() << " !" << endl;

	}
	
       	pFile = fopen(fevname.c_str(), "wb");
   	if (pFile == NULL) perror ("Error opening file");
	else{
		elems += fwrite(&Ns, sizeof(int), 1, pFile);
		elems += fwrite(&Ns, sizeof(int), 1, pFile);
		elems += fwrite(&Ns, sizeof(int), 1, pFile);
		elems += fwrite(&Nt, sizeof(int), 1, pFile);
		
		for(int i1=0;i1<leng1;i1++)
		for(int i2=0;i2<leng2;i2++)
		for(int i3=0;i3<leng3;i3++){
			is = i1 + i2*leng1 + i3*leng1*leng2;
			elems += fwrite(&pollev[is][0], sizeof(std::complex<double>), matrixdim, pFile);
		}

		fclose(pFile);
	}
	
        return nindex + 4 - elems;
}

void fillPollEv(vector<vector<complex<double> > > &lpollarray){
	// Fills the 3d lattice with the Polyakov loop eigenvalues
	complex<double> polleval[3];

	int is=0;

	for(int i1=0;i1<Ns;i1++)
	for(int i2=0;i2<Ns;i2++)
	for(int i3=0;i3<Ns;i3++){
		localPollEv(lpollarray, polleval, i1, i2, i3);
		is = i1 + i2*Ns + i3*Ns*Ns;
		for(int i=0;i<3;i++){
			pollev[is].push_back(polleval[i]);
		}
	}
}

void localPollEv(vector<vector<complex<double> > > &lpollarray, complex<double> *polleval, int i1, int i2, int i3){
	complex<double> poll[3*3];

	int is=0;

	is = i1 + i2*Ns + i3*Ns*Ns;

	for(int i=0;i<9;i++)
		poll[i] = lpollarray[is][i];
	
	int dim=3;

	char NN='N';

        int LWORK = 4*dim;
        complex<double> *WORK = new complex<double>[LWORK];
        double *RWORK = new double[2*dim];
       	int INFO;

	complex<double> VL[dim], VR[dim];
	int one=1;
	zgeev_(&NN, &NN, &dim, &poll[0], &dim, polleval, VL, &dim, VR, 
			&dim, WORK, &LWORK, RWORK, &INFO);

	delete [] WORK;
	delete [] RWORK;
}
