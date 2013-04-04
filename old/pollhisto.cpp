using namespace std;

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <limits>
#include <string>

long double **Create2D(int row, int col);
void Delete2D(long double **p, int row);

int main(int argv, char *argc[]){
	fstream fobs;
	fstream fmobs;
	
	string fobsname, fmobsname;
	long double *pl, *plhisto;
	
	int nmeas,nbeta,nbetacheck, dbeta;
        long double betafrom,betato,betadelta;
	
	string strtmp;
	
	if(argv<4)
		return 1;

	fobsname=string(argc[1]);
	
	fobs.open(fobsname.c_str());
	nmeas = 1000;
        
	fmobsname=string(argc[3]);
        
	int nbins=atoi(argc[2]);
	
        fmobs.open(fmobsname.c_str(),ios_base:: out);
	fmobs.precision(numeric_limits<long double>::digits10 + 1);
	cout.precision(numeric_limits<long double>::digits10 + 1);
	
	fmobs << "# Plaquette[i] | Hplaquette[i]" << endl;
	
	pl = new long double[nmeas];

	long double dtmp;

	int lineskip=9;

        for(int j=0;j<lineskip;j++){
                getline(fobs, strtmp);
        }

	// Plaquettes
	for(int j=0; j<nmeas; j++){
		// fobs >> dtmp >> pl[j] >> dtmp >> dtmp;
		fobs >> pl[j] >> dtmp >> dtmp >> dtmp >> dtmp >> dtmp;
	}
	
	/* ............................................................................... */
        plhisto = new long double[nbins];
        for(int i=0;i<nbins;i++){
        	plhisto[i]=0;
        }

        // search for smallest and largest value in plaquettes
        long double plmax=-10000,plmin=10000000;
        for(int i=0;i<nmeas;i++){
                if(pl[i]>plmax)
                        plmax=pl[i];
                if(pl[i]<plmin)
                        plmin=pl[i];
        }
        cout << "plmax=" << plmax << endl;
        cout << "plmin=" << plmin << endl;
	
	long double dpl;
	dpl=(plmax-plmin)/nbins;
	cout << "Delta Plaquette=" << dpl << endl;
	int ibin=0;
	for(int i=0;i<nmeas;i++){
		ibin=(pl[i]-plmin)/dpl;
		plhisto[ibin]++;
	}
	
	// normalization
	for(int i=0; i<nbins;i++){
		plhisto[i]=plhisto[i]/(nmeas*dpl);
	}
	
	// test normalization
	long double sumpl=0;
	for(int i=0; i<nbins;i++){
		sumpl=sumpl+plhisto[i]*dpl;
	}
	
	cout << "Sumpl=" << sumpl << endl;
	
	for(int i=0;i<nbins;i++){
		fmobs << plmin+i*dpl << " " << plhisto[i] << endl;
	}

        fmobs.close();
}

long double **Create2D(int row, int col)
{
   long double **p = new long double* [row];
   for (int j = 0; j < row; j ++)
      p[j] = new long double[col];
   return p;
}

void Delete2D(long double **p, int row)
{
   for (int j = 0; j < row; j ++)
      delete [] p[j];
   delete [] p;
}
