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
#include "version.h"

#include "include/jackknife.h"

const int matrixdim=3;

int Ns=4;
int Nt=4;
int nbinsx=16;
int Nspace=Ns*Ns*Ns;

int nmeas=1;

vector<string> fevname;

vector<vector<complex<double> > > pollev;
vector<vector<double> > abshisto;
vector<vector<double> > arghisto;

#include "genpollhisto1d_jack_init.hpp"

complex<double> localPoll(int is);

void histoabs(vector<double> &histo, int nbinsx, double *datax, int Ndatax);
void writehistoabs(string fhistoname, vector<double> &histo, vector<double> &histoerr, int nbinsx);

int main(int argc, char *argv[]){

	// Handle command line information and initialize arrays and other stuff...
        if(init(argc, argv) != 0){
                cout << "Error in init()!" << endl;
                return 1;
        }  
        cout << "Starting histogram calculation... " << flush;
	for(int n=0;n<nmeas;n++){
		if(readPollEvBinary(Ns, Ns, Ns, Nt, 3, pollev, fevname[n]) != 0){
			cout << "ERROR: Problems with readPollEvBinary !" << endl;
			return 1;
		} 

	/*	cout << "Performing a random Z_3 rotation... " << flush;
		config->z3rot();
		cout << "done!" << endl << endl; */
	
		complex<double> lpoll;

		double repoll[Ns*Ns*Ns], impoll[Ns*Ns*Ns], abspoll[Ns*Ns*Ns];

		for(int is=0;is<Ns*Ns*Ns;is++){
			lpoll=localPoll(is);

			repoll[is]=real(lpoll);
			impoll[is]=imag(lpoll);
			abspoll[is]=sqrt(real(lpoll)*real(lpoll)+imag(lpoll)*imag(lpoll));
		}

		histoabs(abshisto[n], nbinsx, abspoll, Ns*Ns*Ns);

	/*	cout << "Calculating histogram for ArgP(x)... " << endl;
		gslhisto1d(arghisto[0], nbinsx, argpoll, Ns*Ns*Ns);
		cout << "done!" << endl << endl; */
	}
	cout << "done!" << endl;
	
	// Jackknifing...
	vector<double> mabshisto, mabshistoerr;
	mabshisto.resize(nbinsx);
	mabshistoerr.resize(nbinsx);
	
	double jackkarr[nmeas];
	
	cout << "Performing single elimination Jackknife... " << endl;
	for(int b=0;b<nbinsx;b++){
		for(int n=0;n<nmeas;n++){
			jackkarr[n]=0;
			for(int j=0;j<nmeas;j++){
				if(j!=n)
					jackkarr[n] += abshisto[j][b];
			}
			jackkarr[n] = jackkarr[n]/(double)(nmeas-1);
		}
		Jackknife(jackkarr, mabshisto[b], mabshistoerr[b], nmeas);
	}
	cout << "done!" << endl;
	
	cout << "Writing results to file... " << flush;
	stringstream fhistoname;
	fhistoname << "histo_" << Ns << "x" << Nt << ".res";
	writehistoabs(fhistoname.str(), mabshisto, mabshistoerr, nbinsx);
	cout << "done!" << endl;

	return 0;
}

complex<double> localPoll(int is){
	complex<double> poll;

	poll = pollev[is][0] + pollev[is][1] + pollev[is][2];
	
	return poll;
}

void histoabs(vector<double> &histo, int nbinsx, double *datax, int Ndatax){
	// Find x and y range

	double xmin=0;
	double xmax=3;

	// Allocate gsl_histogram space
	gsl_histogram * h = gsl_histogram_alloc(nbinsx);
	gsl_histogram_set_ranges_uniform (h, xmin, xmax);
	
	int ret=0;
	for(int x=0;x<Ndatax;x++){
		ret = gsl_histogram_increment(h, datax[x]);
	}
	
	double norm=0;
	double upper=0, lower=0;
	gsl_histogram_get_range(h, 1, &lower, &upper);
	double delta = upper-lower;
	for(int ibx=0;ibx<nbinsx;ibx++){
		norm += gsl_histogram_get(h, ibx)*delta;
	}

	for(int ibx=0;ibx<nbinsx; ibx++){
		histo[ibx] = gsl_histogram_get(h, ibx)/(double)norm;
	}

	// Deallocate gsl_histogram space
	gsl_histogram_free(h);
	h=0;
}

void writehistoabs(string fhistoname, vector<double> &histo, vector<double> &histoerr, int nbinsx){
	double xmin=0;
	double xmax=3;
	
	gsl_histogram * h = gsl_histogram_alloc(nbinsx);
	gsl_histogram_set_ranges_uniform (h, xmin, xmax);
	
	cout << "Writing histogram data to file " << fhistoname << " ..." << endl;
	ofstream fhisto;
	double lower, upper;
	fhisto.open(fhistoname.c_str());
	for(int ibx=0;ibx<nbinsx; ibx++){
		gsl_histogram_get_range(h, ibx, &lower, &upper);
		fhisto << lower + (upper - lower)/(double)2 << " " << histo[ibx] << " " << histoerr[ibx] << endl;
	}

	fhisto.close();
}
