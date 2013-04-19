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

#include "include/generic.hpp"

int Ns=0;
int Nt=0;

vector<vector<complex<double> > > pollev;

complex<double> localPoll(int is);

void gslhisto1d(string fhistofilename, int nbinsx, double *datax, int Ndatax);

int main(int argc, char *argv[]){
	if(argc<7){
		cout << "./cluster.x Ns Nt pollevin phaseout histoout nbinsx" << endl;
		return 1;
	}

	Ns=atoi(argv[1]);
	Nt=atoi(argv[2]);
	string fconfigname=argv[3];
	string foutname=argv[4];
	string fhistofilename=argv[5];
	int nbinsx=atoi(argv[6]);

	pollev.resize(Ns*Ns*Ns);
	for(unsigned int is=0; is<pollev.size();is++){
		pollev[is].resize(3);
	}

	if(readPollEvBinary(Ns, Ns, Ns, Nt, 3, pollev, fconfigname) != 0){
		cout << "ERROR: Problems with readPollEvBinary !" << endl;
		return 1;
	} 

/*	cout << "Performing a random Z_3 rotation... " << flush;
	config->z3rot();
	cout << "done!" << endl << endl; */
	
	cout << "Writing local Polyakov loop phases... " << flush;
	ofstream file;
	file.open(foutname.c_str());
	file << "# cnt Re(lPoll(cnt)) Im(lPoll(cnt)) arg(lPoll(cnt)) abs(lPoll(cnt))" << endl;
	complex<double> lpoll;

	double repoll[Ns*Ns*Ns], impoll[Ns*Ns*Ns], abspoll[Ns*Ns*Ns];

	for(int is=0;is<Ns*Ns*Ns;is++){
		lpoll=localPoll(is);
		file << is << " " << real(lpoll) << " " << imag(lpoll) << " " << arg(lpoll) << " " << abs(lpoll) << endl;

		repoll[is]=real(lpoll);
		impoll[is]=imag(lpoll);
		abspoll[is]=sqrt(real(lpoll)*real(lpoll)+imag(lpoll)*imag(lpoll));
	}
	file.close();
	cout << "done!" << endl << endl;

	cout << "Calculating histogram for AbsP(x)... " << endl;
	gslhisto1d(fhistofilename, nbinsx, abspoll, Ns*Ns*Ns);
	cout << "done!" << endl << endl;

	return 0;
}

complex<double> localPoll(int is){
	complex<double> poll;

	poll = pollev[is][0] + pollev[is][1] + pollev[is][2];
	
	return poll;
}

void gslhisto1d(string fhistofilename, int nbinsx, double *datax, int Ndatax){
	// Find x and y range
	double xmax=datax[0];
	double xmin=datax[1];
	for(int n=0; n<Ndatax; n++){
		if(datax[n]<xmin)
			xmin=datax[n];
		if(datax[n]>xmax)
			xmax=datax[n];
	}

	xmin=0;
	xmax=3;

	cout << "Ranges:" << endl;
	cout << "x_min = " << xmin << " x_max = " << xmax << endl;

	// Allocate gsl_histogram space
	gsl_histogram * h = gsl_histogram_alloc(nbinsx);
	gsl_histogram_set_ranges_uniform (h, xmin, xmax);

	int ret=0;
	for(int x=0;x<Ndatax;x++){
		ret = gsl_histogram_increment(h, datax[x]);
		/* if(ret == GSL_EDOM){
			cout << "Entry not added to histogram!" << endl;
			cout << "x = " << datax[x] << " y = " << datay[y] << endl;
		} */
	}

	cout << "Writing histogram data to file " << fhistofilename << " ..." << endl;
	ofstream fhisto;
	double lower, upper;
	fhisto.open(fhistofilename.c_str());
	for(int ibx=0;ibx<nbinsx; ibx++){
		gsl_histogram_get_range(h, ibx, &lower, &upper);
		fhisto << lower + (upper - lower)/(double)2 << " " << gsl_histogram_get(h, ibx) << endl;
	}

	fhisto.close();

	cout << "Sum over all bins = " << gsl_histogram_sum(h) << endl;

	// Deallocate gsl_histogram space
	gsl_histogram_free(h);
	h=0;
}
