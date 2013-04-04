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

extern "C" {
	// Eigenvalues for a general complex matrix
	void zgeev_(const char *JOBVBL, const char* JOBVR, int *N, complex<double>* A, int* lda, 
			complex<double>* EVAL, complex<double>* VL, int *ldvl, complex<double>* VR, 
			int *ldvr, complex<double>* WORK, int *lwork, double *RWORK, int* INFO);
}

int Ns=0;
int Nt=0;

ConfigData *config;
vector<vector<complex<double> > > pollev;

complex<double> localPollEv(complex<double> *polleval, int i1, int i2, int i3);
complex<double> totalPollTr();
void fillPollEv();
int writePollEv(string fevname);
void checkPollEv();

void gslhisto(string fhistofilename, int nbinsx, int nbinsy, double *datax, double *datay, int Ndatax, int Ndatay);

int main(int argc, char *argv[]){
	if(argc<5){
		cout << "./cluster.x Ns Nt confin outfile" << endl;
		return 1;
	}

	Ns=atoi(argv[1]);
	Nt=atoi(argv[2]);
	string fconfigname=argv[3];
	string foutname=argv[4];

	int Nspace = Ns*Ns*Ns;

	config = new ConfigData(Ns, Ns, Ns, Nt, 3);

	cout << "Reading the file..." << endl;
	if(config->readBinaryConfig2(fconfigname)){
		return 1;
	}

	complex<double> pollTr;

	pollTr = totalPollTr();
	cout << "Total Polyakov loop P = " << pollTr << endl;

	cout << "Performing a random Z_3 rotation... " << flush;
	config->z3rot();
	cout << "done!" << endl;

	cout << "Allocating 3d (Ns^3 + 3) lattice for Polyakov loop ev's... " << flush;
	pollev.resize(Nspace);
	cout << "done!" << endl;

	cout << "Filling 3d lattice with Polyakov loop ev's... " << flush;
	fillPollEv();
	cout << "done!" << endl;

	checkPollEv();

	cout << "Writing 3d lattice to disk... " << flush;
	string fevname("pollev.data");
	writePollEv(fevname);
	cout << "done!" << endl;

	return 0;
}

void checkPollEv(){
	cout << "Performing some checks on the ev's... " << flush;

	complex<double> sum;
	double eps=1E-15;

	for(int is=0;is<Ns*Ns*Ns;is++)
		for(int i=0;i<3;i++)
			sum += abs(pollev[is][i]);

	sum = sum/(double)(3*Ns*Ns*Ns);

	if(abs(sum) - 1.0 > eps)
		cout << endl << "WARNING: Sum over all abs(ev)'s / 3Ns*Ns*Ns = " << abs(sum) << " > 1"  << endl;

	cout << "done!" << endl;
}

int writePollEv(string fevname){
	ofstream fev;
	fev.open(fevname.c_str());
	
	fev << "# is real(ev_1) imag(ev_1) real(ev_2) imag(ev_2) real(ev_3) imag(ev_3)" << endl;
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

void fillPollEv(){
	// Fills the 3d lattice with the Polyakov loop eigenvalues
	complex<double> polleval[3];

	int is=0;

	for(int i1=0;i1<Ns;i1++)
	for(int i2=0;i2<Ns;i2++)
	for(int i3=0;i3<Ns;i3++){
		localPollEv(polleval, i1, i2, i3);
		is = i1 + i2*Ns + i3*Ns*Ns;
		// cout << "Pushing ev's... " << flush;
		for(int i=0;i<3;i++){
			pollev[is].push_back(polleval[i]);
			// cout << polleval[i] << endl;
		}
		// cout << "done!" << endl;
	}
}

complex<double> localPollEv(complex<double> *polleval, int i1, int i2, int i3){
	complex<double> poll[3][3];

	complex<double> up[3][3], uu[3][3], upaux[3][3];

	int t=0, i4=0, is=0;

	is = config->latmap(i1, i2, i3, i4);
	config->extract(*up,3,is);

	for(i4=1;i4<Nt-1;i4++){
		is = config->latmap(i1, i2, i3, i4);
		config->extract(*uu,3,is);
		axb(*upaux,*up,*uu, 3);
		aeb(*up,*upaux, 3);
	}    

	i4 = Nt-1;
	is = config->latmap(i1, i2, i3, i4);

	config->extract(*uu,3,is);
	axb(*poll,*up,*uu, 3);

//	cout << poll[0][0] << endl;
	
	int dim=3;

	char NN='N';

        int LWORK = 4*dim;
        complex<double> *WORK = new complex<double>[LWORK];
        double *RWORK = new double[2*dim];
       	int INFO;

	complex<double> VL[dim], VR[dim];
	int one=1;
	// cout << "Calculating eigenvalue... " << flush;
	zgeev_(&NN, &NN, &dim, &poll[0][0], &dim, polleval, VL, &dim, VR, 
			&dim, WORK, &LWORK, RWORK, &INFO);
	// cout << "done!" << endl;

	delete [] WORK;
	delete [] RWORK;
}

complex<double> totalPollTr(){
	complex<double> poll;

	complex<double> up[3][3], uu[3][3], upaux[3][3];

	int t=0, i4=0, is=0;
	for(int i1=0;i1<Ns;i1++)
	for(int i2=0;i2<Ns;i2++)
	for(int i3=0;i3<Ns;i3++){
		i4=0; // i4 = time direction
		is = config->latmap(i1, i2, i3, i4);
		config->extract(*up,3,is);

		for(i4=1;i4<Nt-1;i4++){
			is = config->latmap(i1, i2, i3, i4);
			config->extract(*uu,3,is);
			axb(*upaux,*up,*uu, 3);
			aeb(*up,*upaux, 3);
		}    

		i4 = Nt-1;
		is = config->latmap(i1, i2, i3, i4);

		config->extract(*uu,3,is);
		poll += multtrace(*up,*uu, 3);
	}

	return poll/((double)3*Ns*Ns*Ns);
}

