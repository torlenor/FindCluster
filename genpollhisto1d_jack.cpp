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
complex<double> totalPoll();

void histoabs(vector<double> &histo, int nbinsx, double *datax, int Ndatax);
void writehistoabs(string fhistoname, vector<double> &histo, vector<double> &histoerr, int nbinsx);

void histoarg(vector<double> &histo, int nbinsx, double *datax, int Ndatax);
void writehistoarg(string fhistoname, vector<double> &histo, vector<double> &histoerr, int nbinsx);

int main(int argc, char *argv[]){

	// Handle command line information and initialize arrays and other stuff...
        if(init(argc, argv) != 0){
                cout << "Error in init()!" << endl;
                return 1;
        }  

	cout << endl << "Settings:" << endl;
	cout << "Lattice size = " << Ns << "x" << Nt << endl;
	cout << endl;

	complex<double> lpoll;
	double repoll[Ns*Ns*Ns], impoll[Ns*Ns*Ns], abspoll[Ns*Ns*Ns], argpoll[Ns*Ns*Ns];

        cout << "Starting histogram calculation... " << flush;
	for(int n=0;n<nmeas;n++){
		if(readPollEvBinary(Ns, Ns, Ns, Nt, 3, pollev, fevname[n]) != 0){
			cout << "ERROR: Problems with readPollEvBinary !" << endl;
			return 1;
		} 

		for(int is=0;is<Ns*Ns*Ns;is++){
			lpoll=localPoll(is);

			repoll[is]=real(lpoll);
			impoll[is]=imag(lpoll);
			abspoll[is]=sqrt(real(lpoll)*real(lpoll)+imag(lpoll)*imag(lpoll));
		}

		histoabs(abshisto[n], nbinsx, abspoll, Ns*Ns*Ns);

		// For arg Histogram we have to perform a phase rotation so that the
		// Polyakov loops are all in the same phase.
		//
		// To do that, we first calculate the total Polyakov loop, check the
		// phase and then rotate every local Polyakov loop by the correct phase.
		complex<double> tpoll;
		complex<double> phase(cos(2*M_PI/(double)3), sin(2*M_PI/(double)3));
		bool correctphase=false;
		int whatphase=0; // Is -1, 0 or 1
		while(! correctphase){
			tpoll=totalPoll();

			if( abs( arg(tpoll) + (double)whatphase*2*M_PI/(double)3) < M_PI/(double)3){
				correctphase=true;
			}else{
				for(int is=0;is<Ns*Ns*Ns;is++){
					pollev[is][0] = phase*pollev[is][0];
					pollev[is][1] = phase*pollev[is][1];
					pollev[is][2] = phase*pollev[is][2];
				}
			}
		}
		
		for(int is=0;is<Ns*Ns*Ns;is++){
			lpoll=localPoll(is);

			argpoll[is]=arg(lpoll);
		}
		histoarg(arghisto[n], nbinsx, argpoll, Ns*Ns*Ns);
	}
	cout << "done!" << endl;
	
	// Jackknifing...
	vector<double> mabshisto, mabshistoerr;
	mabshisto.resize(nbinsx);
	mabshistoerr.resize(nbinsx);
	
	double jackkarr[nmeas];
	
	cout << "Performing single elimination Jackknife (abs(P))... " << flush;
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

	vector<double> marghisto, marghistoerr;
	marghisto.resize(nbinsx);
	marghistoerr.resize(nbinsx);
	cout << "Performing single elimination Jackknife (arg(P))... " << flush;
	for(int b=0;b<nbinsx;b++){
		for(int n=0;n<nmeas;n++){
			jackkarr[n]=0;
			for(int j=0;j<nmeas;j++){
				if(j!=n)
					jackkarr[n] += arghisto[j][b];
			}
			jackkarr[n] = jackkarr[n]/(double)(nmeas-1);
		}
		Jackknife(jackkarr, marghisto[b], marghistoerr[b], nmeas);
	}
	cout << "done!" << endl;
	
	cout << "Writing results to file... " << flush;
	stringstream fhistoname;
	fhistoname << "histo_" << Ns << "x" << Nt << ".res";
	writehistoabs(fhistoname.str(), mabshisto, mabshistoerr, nbinsx);
	cout << "done!" << endl;
	
	cout << "Writing results to file (arg(P))... " << flush;
	stringstream farghistoname;
	farghistoname << "histo_arg_" << Ns << "x" << Nt << ".res";
	writehistoarg(farghistoname.str(), marghisto, marghistoerr, nbinsx);
	cout << "done!" << endl;

	return 0;
}

complex<double> localPoll(int is){
	complex<double> poll;

	poll = pollev[is][0] + pollev[is][1] + pollev[is][2];
	
	return poll;
}

complex<double> totalPoll(){
	complex<double> totalpoll=0;
	for(int is=0;is<Ns*Ns*Ns;is++){
		totalpoll += localPoll(is);
	}

	totalpoll = totalpoll/(double)(Ns*Ns*Ns);

	return totalpoll;

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
	
	ofstream fhisto;
	double lower, upper;
	fhisto.open(fhistoname.c_str());
	for(int ibx=0;ibx<nbinsx; ibx++){
		gsl_histogram_get_range(h, ibx, &lower, &upper);
		fhisto << lower + (upper - lower)/(double)2 << " " << histo[ibx] << " " << histoerr[ibx] << endl;
	}

	fhisto.close();
}

//
// arg(P) block
// 

void histoarg(vector<double> &histo, int nbinsx, double *datax, int Ndatax){
	// Find x and y range

	double xmin = -M_PI;
	double xmax =  M_PI;

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
	double delta = abs(upper-lower);
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

void writehistoarg(string fhistoname, vector<double> &histo, vector<double> &histoerr, int nbinsx){
	double xmin = -M_PI;
	double xmax =  M_PI;
	
	gsl_histogram * h = gsl_histogram_alloc(nbinsx);
	gsl_histogram_set_ranges_uniform (h, xmin, xmax);
	
	ofstream fhisto;
	double lower, upper;
	fhisto.open(fhistoname.c_str());
	for(int ibx=0;ibx<nbinsx; ibx++){
		gsl_histogram_get_range(h, ibx, &lower, &upper);
		fhisto << lower + abs(upper - lower)/(double)2 << " " << histo[ibx] << " " << histoerr[ibx] << endl;
	}

	fhisto.close();
}
