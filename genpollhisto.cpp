#include <iostream>
#include <fstream>
#include <complex>

#include <gsl/gsl_histogram.h>
#include <gsl/gsl_histogram2d.h>
#include <gsl/gsl_errno.h>

#include "include/ConfigData.hpp"

using namespace std;

#include "include/matrixfunctions.hpp"

int Ns=0;
int Nt=0;

ConfigData *config;

complex<double> localPoll(int i1, int i2, int i3);
complex<double> totalPoll();

void gslhisto(string fhistofilename, int nbinsx, int nbinsy, double *datax, double *datay, int Ndatax, int Ndatay);

int main(int argc, char *argv[]){
	if(argc<8){
		cout << "./cluster.x Ns Nt confin phaseout histoout nbinsx nbinsy" << endl;
		return 1;
	}

	Ns=atoi(argv[1]);
	Nt=atoi(argv[2]);
	string fconfigname=argv[3];
	string foutname=argv[4];
	string fhistofilename=argv[5];
	int nbinsx=atoi(argv[6]);
	int nbinsy=atoi(argv[7]);

	config = new ConfigData(Ns, Ns, Ns, Nt, 3);

	cout << "Reading the file..." << endl;
	if(config->readBinaryConfig2(fconfigname)){
		return 1;
	}

	complex<double> poll;

	poll = totalPoll();
	cout << "Total Polyakov loop P = " << poll << endl << endl;

	cout << "Performing a random Z_3 rotation... " << flush;
	config->z3rot();
	cout << "done!" << endl << endl;
	
	cout << "Writing local Polyakov loop phases... " << flush;
	ofstream file;
	file.open(foutname.c_str());
	file << "# cnt Re(lPoll(cnt)) Im(lPoll(cnt)) arg(lPoll(cnt)) abs(lPoll(cnt))" << endl;
	complex<double> lpoll;
	int cnt=0;

	double repoll[Ns*Ns*Ns], impoll[Ns*Ns*Ns];

	for(int i1=0;i1<Ns;i1++)
	for(int i2=0;i2<Ns;i2++)
	for(int i3=0;i3<Ns;i3++){
		cnt++;
		lpoll=localPoll(i1, i2, i3);
		file << cnt << " " << real(lpoll) << " " << imag(lpoll) << " " << arg(lpoll) << " " << abs(lpoll) << endl;
		// file << cnt << " " << arg(lpoll) << endl;
		repoll[cnt-1]=real(lpoll);
		impoll[cnt-1]=imag(lpoll);
	}
	file.close();
	cout << "done!" << endl << endl;

	cout << "Calculating histogram for ReP(x) and ImP(x)... " << endl;
	gslhisto(fhistofilename, nbinsx, nbinsy, repoll, impoll, Ns*Ns*Ns, Ns*Ns*Ns);
	cout << "done!" << endl << endl;

	return 0;
}

complex<double> localPoll(int i1, int i2, int i3){
	complex<double> poll;

	complex<double> up[3][3], uu[3][3], upaux[3][3];

	int t=0, i4=0, is=0;

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
	poll = multtrace(*up,*uu, 3);

	return poll/(double)3;
}

complex<double> totalPoll(){
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

void gslhisto(string fhistofilename, int nbinsx, int nbinsy, double *datax, double *datay, int Ndatax, int Ndatay){
	// Find x and y range
	double xmax=datax[0];
	double xmin=datax[1];
	for(int n=0; n<Ndatax; n++){
		if(datax[n]<xmin)
			xmin=datax[n];
		if(datax[n]>xmax)
			xmax=datax[n];
	}
	double ymax=datay[0];
	double ymin=datay[1];
	for(int n=0; n<Ndatay; n++){
		if(datay[n]<ymin)
			ymin=datay[n];
		if(datay[n]>ymax)
			ymax=datay[n];
	}

	cout << "Ranges:" << endl;
	cout << "x_min = " << xmin << " x_max = " << xmax << endl;
	cout << "y_min = " << ymin << " y_max = " << ymax << endl;

	// Allocate gsl_histogram space
	gsl_histogram2d * h = gsl_histogram2d_alloc(nbinsx, nbinsy);
	gsl_histogram2d_set_ranges_uniform (h, xmin, xmax, ymin, ymax);

	int ret=0;
	for(int x=0;x<Ndatax;x++)
		for(int y=0;y<Ndatay;y++){
			ret = gsl_histogram2d_increment(h, datax[x], datay[y]);
			/* if(ret == GSL_EDOM){
				cout << "Entry not added to histogram!" << endl;
				cout << "x = " << datax[x] << " y = " << datay[y] << endl;
			} */
		}

	cout << "Writing histogram data to file " << fhistofilename << " ..." << endl;
	ofstream fhisto;
	fhisto.open(fhistofilename.c_str());
	for(int ibx=0;ibx<nbinsx; ibx++){
		for(int iby=0;iby<nbinsy; iby++){
			fhisto << ibx << " " << iby << " " << gsl_histogram2d_get(h, ibx, iby) << endl;
		}
		fhisto << endl;
	}

	fhisto.close();

	cout << "Sum over all bins = " << gsl_histogram2d_sum(h) << endl;

	// Deallocate gsl_histogram space
	gsl_histogram2d_free(h);
	h=0;
}
