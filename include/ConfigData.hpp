// v1.1 - 2012-10-11 - 15:00
#ifndef CONFIGFUNCTIONS_H
#define CONFIGFUNCTIONS_H

#include <complex>
#include <vector>
#include <fstream>
#include <limits>
#include <cstdio>
#include <cstdlib>

#include <cstdlib>

#include "CArray.hpp"

using namespace std;

#define DIM 4 //Space-Time dimensions
#define MD 3 // SU(MD) Matrix dimension

class ConfigData{
	public:
		// Constructor
		ConfigData(int len1, int len2, int len3, int len4, int md, bool vb);
		~ConfigData();
		
		// Interfaces
		void extract(complex<double> *m, int mu, int nindex); // Returns the matrix @ (nindex,mu)
		void replace(complex<double> *m, int mu, int nindex); // Replaces the matrix @ (nindex,mu)

		void dumpConfig(); // Dumps the config to stdout
		void freeConfig(); // Creates a free config (diag(1,1,1) matrices on every link)

		int readConfig(string fconfigname); // Reads a config from file
		int readFConfig(string fconfigname); // Fortran config
		int readBinaryConfig(string fconfigname);
		int readBinaryConfig2(string fconfigname);
		int writeConfig(string fconfigname); // Reads a config from file
		int writeBinaryConfig(string fconfigname);
		int writeBinaryConfig2(string fconfigname);

		int latmap(int i1, int i2, int i3, int i4);
		int neib(int site, int mu);

		// MILC
		int MILCreadConfig(string fconfigname); // Reads a config from file
		int MILCreadBinaryConfig(string fconfigname);

		complex<double> U(int nindex, int mu, int c1, int c2);
		
		complex<double> calcPoll();
		complex<double> calcPlaq();

		// Function tests
		void runTests();

		// Switch to antiperiodic boundary conditions
		void antiperbc();

		// Perform random Z3 rotation
		void z3rot();

	private:
		// Variables which store basic dimensions
		int leng1, leng2, leng3, leng4, matrixdim, nsite;
		int **neibArray;

		bool verbose;
		
		// Stores the config
		// vector<vector<vector<vector< complex<double> > > > > A;
		CArray *A;

		// Used to create 2D arrays on heap (e.g., neib[][])
		int **Create2D(int row, int col);
		void Delete2D(int **p, int row);

		// Functions for checking the configuration
		void staplesum(complex<double> *S, int mu, int nindex);

		// Matrix functions
		void aeb(complex<double> *a, complex<double> *b);
		void apb(complex<double> *c, complex<double> *a);
		void amb(complex<double> *c,complex<double> *a, complex<double> *b);
		void capb(complex<double> *c, complex<double> *a, complex<double> *b);

		void axb(complex<double> *a,complex<double> *b, complex<double> *c);
		void axbdag(complex<double> *a,complex<double> *b, complex<double> *c);
		void adagxb(complex<double> *a,complex<double> *b, complex<double> *c);
		void adagxbdag(complex<double> *a,complex<double> *b, complex<double> *c);
		complex<double> multtrace(complex<double> *a, complex<double> *b);
		void za(complex<double> *c, double z, complex<double> *a);
		void adag(complex<double> *a);

		// Test functions
		void printMatrix(complex<double> *a, int md);
};

ConfigData::ConfigData(int len1, int len2, int len3, int len4, int md, bool vb=true){
	leng1=len1; leng2=len2; leng3=len3; leng4=len4; matrixdim=md;
	nsite=leng1*leng2*leng3*leng4;

	verbose=vb;

	A = new CArray(nsite, DIM, matrixdim, matrixdim);

	neibArray = Create2D(nsite, DIM*2);

	// Fills the neib array
        int i1p,i2p,i3p,i4p,i1m,i2m,i3m,i4m,is,isp1,isp2,isp3,isp4,ism1,ism2,ism3,ism4;
        for(int i1 = 0;i1<leng1;i1++){
                i1p = i1 + 1;
                i1m = i1 - 1;
                if (i1p == leng1) i1p = 0;
                if (i1m == -1) i1m = leng1-1;

                for(int i2 = 0;i2<leng2;i2++){
                        i2p = i2 + 1;
                        i2m = i2 - 1;
                        if (i2p == leng2) i2p = 0;
                        if (i2m == -1) i2m = leng2-1;

                        for(int i3 = 0;i3<leng3;i3++){
                                i3p = i3 + 1;
                                i3m = i3 - 1;
                                if (i3p == leng3) i3p = 0;
                                if (i3m == -1) i3m = leng3-1;

                                for(int i4 = 0;i4<leng4;i4++){
                                        i4p = i4 + 1;
                                        i4m = i4 - 1;
                                        if (i4p == leng4) i4p = 0;
                                        if (i4m == -1) i4m = leng4-1;

                                        // Compute the site address and the addresses of the sites shifted
                                        // by one unit in each direction

                                        is = i1 + i2*leng1 + i3*leng1*leng2 + i4*leng1*leng2*leng3;

                                        isp1 = i1p + i2*leng1 + i3*leng1*leng2 + i4*leng1*leng2*leng3;
                                        isp2 = i1 + i2p*leng1 + i3*leng1*leng2 + i4*leng1*leng2*leng3;
                                        isp3 = i1 + i2*leng1 + i3p*leng1*leng2 + i4*leng1*leng2*leng3;
                                        isp4 = i1 + i2*leng1 + i3*leng1*leng2 + i4p*leng1*leng2*leng3;

                                        ism1 = i1m + i2*leng1 + i3*leng1*leng2 + i4*leng1*leng2*leng3;
                                        ism2 = i1 + i2m*leng1 + i3*leng1*leng2 + i4*leng1*leng2*leng3;
                                        ism3 = i1 + i2*leng1 + i3m*leng1*leng2 + i4*leng1*leng2*leng3;
                                        ism4 = i1 + i2*leng1 + i3*leng1*leng2 + i4m*leng1*leng2*leng3;

                                        // Fill the neib array

                                        neibArray[is][0] = isp1;
                                        neibArray[is][1] = isp2;
                                        neibArray[is][2] = isp3;
                                        neibArray[is][3] = isp4;

                                        neibArray[is][4] = ism1;
                                        neibArray[is][5] = ism2;
                                        neibArray[is][6] = ism3;
                                        neibArray[is][7] = ism4;
                                }
                        }
                }
	}
}

ConfigData::~ConfigData(){
	Delete2D(neibArray, nsite);
}

complex<double> ConfigData::U(int nindex, int mu, int c1, int c2){
	#ifdef DEBUG
		if(nindex<0 || nindex>nsite-1 || c1<0 || c1>3 || c2<0 || c2>3 || mu<0 || mu>3){
			cout << "Error in ConfigData::U(int nindex, int mu, int c1, int c2)" << endl;
			return 0;
		}
	#endif
	return (*A)(nindex,mu,c1,c2);
}

void ConfigData::extract(complex<double> *m, int mu,int nindex){
        // Extracts the matrix at (x,mu) in A and writes it to m
        for(int i=0;i<matrixdim;i++)
		for(int j=0;j<matrixdim;j++){
			m[i*matrixdim + j]=(*A)(nindex,mu,i,j);
		}
	
	// memcpy( m, &(*A)(nindex,mu,0,0), entries*sizeof(complex<double>) );
}

void ConfigData::replace(complex<double> *m, int mu, int nindex){
        // Replaces the matrix at (x,mu) in A with m
        for(int i=0;i<matrixdim;i++)
		for(int j=0;j<matrixdim;j++){
			(*A)(nindex,mu,i,j)=m[i*matrixdim + j];
 	}
	
	// memcpy( &(*A)(nindex,mu,0,0),m, entries*sizeof(complex<double>) );
}

int ConfigData::latmap(int i1, int i2, int i3, int i4){
	int is = i1 + i2*leng1 + i3*leng1*leng2 + i4*leng1*leng2*leng3;

	return is;
}

int ConfigData::neib(int nindex, int mu){
	#ifdef DEBUG
		if(nindex>=nsite || nindex<0 || mu>=DIM*2 || mu<0){
			cout << "Error in ConfigData::neib(int nindex, int mu)" << endl;
			return -1;
		}else{ 
			return neibArray[nindex][mu];
		}
	#endif
	return neibArray[nindex][mu];
}

void ConfigData::antiperbc(){
	// Switches to antiperiodic boundary conditions
	if(verbose)
		cout << "Switching to anti-periodic boundary conditions... " << flush;
	int n=0;
	for(int x1=0;x1<leng1;x1++)
	for(int x2=0;x2<leng2;x2++)
	for(int x3=0;x3<leng3;x3++){
		n = x1 + x2*leng1 + x3*leng1*leng2
			+ (leng4-1)*leng1*leng2*leng3; // Not sure if this is correct!
		for(int i=0;i<matrixdim;i++)
		for(int j=0;j<matrixdim;j++){
			(*A)(n,3,i,j) = -(*A)(n,3,i,j);
		}
	}
	if(verbose)
		cout << "done!" << endl;
}

void ConfigData::dumpConfig(){
        std::complex<double> U[matrixdim][matrixdim];

        for(int n=0;n<nsite;n++){
                for(int mu=0;mu<DIM;mu++){
                        extract(*U,mu,n);
                        for(int i=0; i<matrixdim; i++){
                                for(int j=0; j<matrixdim; j++){
                                        cout << real(U[i][j]) << " " << imag(U[i][j]) << " ";
                                }
                        }
                        cout  << endl;
                }
        }
}

void ConfigData::freeConfig(){
	if(verbose)
		cout << "Using a free configuration... " << flush;
        for(int n=0;n<nsite;n++){
                for(int mu=0;mu<DIM;mu++){
                        //for(int i=0;i<matrixdim;i++)
                                // A[n][mu][i][i]=complex<double>(1,0);
			(*A)(n,mu,0,0)=complex<double>(1, 0);
			(*A)(n,mu,1,1)=complex<double>(1, 0);
			(*A)(n,mu,2,2)=complex<double>(1, 0);
                }
        }

	if(verbose)
		cout << "done!" << endl;
}

int ConfigData::readBinaryConfig(string fconfigname){
	int ERR = A->bread(fconfigname);

        complex<double> poll=calcPoll();
        complex<double> plaq=calcPlaq();

	if(verbose){
		cout << "Configfile " << fconfigname << " (binary) read." << endl << endl;
		cout << "Checks:" << endl;
		cout << "Polyakov loop: calculated=" << poll << endl;
		cout << "Plaquettes: calculated=" << plaq << endl << endl;
	}
	return ERR;
}

int ConfigData::writeBinaryConfig(string fconfigname){
	int ERR = A->bwrite(fconfigname);

	return ERR;
}

int ConfigData::writeBinaryConfig2(string fconfigname){
	int elems=0;
	complex<double> plaq, poll;
       	FILE* pFile;

	int nindex=nsite*matrixdim*matrixdim*DIM;

       	pFile = fopen(fconfigname.c_str(), "wb");
   	if (pFile == NULL) perror ("Error opening file");
	else{
		elems += fwrite(&leng1, sizeof(int), 1, pFile);
		elems += fwrite(&leng2, sizeof(int), 1, pFile);
		elems += fwrite(&leng3, sizeof(int), 1, pFile);
		elems += fwrite(&leng4, sizeof(int), 1, pFile);
		
		poll=calcPoll();
		plaq=calcPlaq();

		elems += fwrite(&poll, sizeof(std::complex<double>),1, pFile);
		elems += fwrite(&plaq, sizeof(std::complex<double>),1, pFile);

		elems += fwrite(&(*A)(0,0,0,0), sizeof(std::complex<double>), nindex, pFile);

		fclose(pFile);
	}
        return nindex + 6 - elems;
}

int ConfigData::readBinaryConfig2(string fconfigname){
	int elems=0;
	int nindex=nsite*matrixdim*matrixdim*DIM;
	
	int cleng1, cleng2, cleng3, cleng4;
	complex<double> pollref, plaqref, poll, plaq;
	
	double maxError=1E-10;
        
	FILE* pFile;
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

		elems += fread(&(*A)(0,0,0,0), sizeof(std::complex<double>), nindex, pFile);
		
		fclose(pFile);
		
		poll=calcPoll();
		plaq=calcPlaq();

		if(verbose){
		cout << "Configfile " << fconfigname << " read." << endl << endl;
		cout << "Checks:" << endl;
		cout << "Polyakov loop: file=" << pollref << " calculated=" << poll << " diff(abs)=" << abs(pollref-poll) << endl;
		cout << "Plaquettes: file=" << plaqref << " calculated=" << plaq << " diff(abs)=" << abs(plaqref-plaq) << endl << endl;
		}

		if(abs(pollref-poll)>maxError || abs(plaqref-plaq)>maxError)
			cout << "WARNING: Calculated values for Polyakov loop and/or Plaquettes differ!" << endl << endl;
	}

        return nindex + 6 - elems;
}


int ConfigData::readConfig(string fconfigname){
        // Reades the current lattice configuration into file with name
        // fconfigname
        complex<double> U[matrixdim][matrixdim];
	int cleng1, cleng2, cleng3, cleng4;
        double realnum, imagnum;
	complex<double> pollref, plaqref, poll, plaq;
        string strtmp;

	double maxError=1E-10;

        ifstream fconfig;
        fconfig.open(fconfigname.c_str());
	if(fconfig.is_open()!=true){
                cout  << "ERROR: File " << fconfigname <<  " to read configuration could not be opened!" << endl;
                return 1;
        } 

	fconfig >> cleng1 >> cleng2 >> cleng3 >> cleng4;

	// Check lattice dimensions
	if(leng1!=cleng1 || leng2!=cleng2 || leng3!=cleng3 || leng4!=cleng4){
		cout << "ERROR: lattice dimensions are different in config file and settings !" << endl;
		return 1;
	}

	double pollrefreal, pollrefimag, plaqrefreal, plaqrefimag;

	fconfig >> pollrefreal >> pollrefimag;
	fconfig >> plaqrefreal >> plaqrefimag;

	plaqref=complex<double>(plaqrefreal, plaqrefimag);
	pollref=complex<double>(pollrefreal, pollrefimag);

        for(int n=0;n<nsite;n++){
                for(int mu=0;mu<DIM;mu++){
                        for(int i=0; i<matrixdim; i++){
                                for(int j=0; j<matrixdim; j++){
                                        fconfig >> realnum >> imagnum;
					U[i][j]=complex<double>(realnum,imagnum);
                                }
                        }
			replace(*U,mu,n);
                }
        }
        
	/* for(int j=0;j<lineskip;j++){
                getline(fconfig, strtmp);
        } */

	// Calculation of Plaquette and Pollyakov loop
	poll=calcPoll();
	plaq=calcPlaq();

	if(verbose){
	cout << "Configfile " << fconfigname << " read." << endl << endl;
	cout << "Checks:" << endl;
	cout << "Polyakov loop: file=" << pollref << " calculated=" << poll << " diff(abs)=" << abs(pollref-poll) << endl;
	cout << "Plaquettes: file=" << plaqref << " calculated=" << plaq << " diff(abs)=" << abs(plaqref-plaq) << endl << endl;
	}

	if(abs(pollref-poll)>maxError || abs(plaqref-plaq)>maxError)
		cout << "WARNING: Calculated values for Polyakov loop and/or Plaquettes differ!" << endl << endl;

        fconfig.close();

        return 0;
}

int ConfigData::readFConfig(string fconfigname){
        // Reades the current lattice configuration into file with name
        // fconfigname
	// Fortran output
        complex<double> U[matrixdim][matrixdim];
        double realnum, imagnum;
	complex<double> pollref, plaqref, poll, plaq;
        string strtmp;

	double maxError=1E-10;

        ifstream fconfig;
        fconfig.open(fconfigname.c_str());
	if(fconfig.is_open()!=true){
                cout  << "ERROR: File " << fconfigname <<  " to read configuration could not be opened!" << endl;
                return 1;
        } 

	for(int i=0;i<28;i++){
		getline(fconfig, strtmp);
	}

        for(int n=0;n<nsite;n++){
                for(int mu=0;mu<DIM;mu++){
                        for(int i=0; i<matrixdim; i++){
                                for(int j=0; j<matrixdim; j++){
                                        fconfig >> realnum >> imagnum;
					U[i][j]=complex<double>(realnum,imagnum);
                                }
                        }
			replace(*U,mu,n);
                }
        }
        
	// Calculation of Plaquette and Pollyakov loop
	poll=calcPoll();
	plaq=calcPlaq();

	if(verbose){
	cout << "Configfile " << fconfigname << " read." << endl << endl;
	cout << "Checks:" << endl;
	cout << "Polyakov loop: file=" << pollref << " calculated=" << poll << " diff(abs)=" << abs(pollref-poll) << endl;
	cout << "Plaquettes: file=" << plaqref << " calculated=" << plaq << " diff(abs)=" << abs(plaqref-plaq) << endl << endl;
	}

	if(abs(pollref-poll)>maxError || abs(plaqref-plaq)>maxError)
		cout << "WARNING: Calculated values for Polyakov loop and/or Plaquettes differ!" << endl << endl;

        fconfig.close();

        return 0;
}

int ConfigData::writeConfig(string fconfigname){
        // Writes the current lattice configuration into file with name 
        // fconfigname
        complex<double> U[matrixdim][matrixdim], plaq, poll;

        ofstream fconfig;
        fconfig.open(fconfigname.c_str());
        if(fconfig.is_open()!=true){
                cout << "ERROR: File " << fconfigname <<  " to save configuration could not be   opened!" << endl;
                return 1;
        }   

        fconfig << leng1 << " " << leng2 << " " << leng3 << " " << leng4 << endl;

        fconfig.flags (std::ios::scientific);
        fconfig.precision(numeric_limits<double>::digits10 + 1); 

	poll=calcPoll();
	plaq=calcPlaq();

        fconfig << real(poll) << " " << imag(poll) << endl;
        fconfig << real(plaq) << " " << imag(plaq) << endl;
    
        for(int n=0;n<nsite;n++){
                for(int mu=0;mu<DIM;mu++){
                        extract(*U,mu,n);
                        for(int i=0; i<matrixdim; i++){
                                for(int j=0; j<matrixdim; j++){
                                        fconfig << real(U[i][j]) << " " << imag(U[i][j]) << " ";
                                }   
                        }   
                        fconfig << endl;
                }   
        }   
    
        fconfig << endl << "# SU(" << matrixdim  << ") Settings:" << endl;
        fconfig << "# Lattice size: " << leng1 << "x" << leng2 << "x" <<leng3 << "x" << leng4 <<    endl;
        fconfig << "# Plaquettes=" << plaq << endl;
        fconfig << "# Polyakov loop=" << poll << endl;
        fconfig << "# U[0][1] U[0][2] U[0][3] U[1][0] U[1][1] ... @ nsite=0" << endl;
        fconfig << "# U[0][1] U[0][2] U[0][3] U[1][0] U[1][1] ... @ nsite=1" << endl;

        fconfig.close();

        return 0;
}

// MILC import functions
int ConfigData::MILCreadConfig(string fconfigname){
        // Reades the current lattice configuration into file with name
        // fconfigname
	// Confifguration is generated using MILC
        complex<double> U[matrixdim][matrixdim];
	int cleng1, cleng2, cleng3, cleng4;
        double realnum, imagnum, pollrefreal, pollrefimag, dtmp, plaqrefreal;
	complex<double> pollref, plaqref, poll, plaq;
        string strtmp;
	
	double maxError=1E-7;

        ifstream fconfig;
        fconfig.open(fconfigname.c_str());
	if(fconfig.is_open()!=true){
                cout  << "ERROR: File " << fconfigname <<  " to read configuration could not be opened!" << endl;
                return 1;
        } 

	// Drop the first two lines
        getline(fconfig, strtmp);
        getline(fconfig, strtmp);

	// Lattice dimension
	fconfig >> cleng1 >> cleng2 >> cleng3 >> cleng4;

	// Check lattice dimensions
	if(leng1!=cleng1 || leng2!=cleng2 || leng3!=cleng3 || leng4!=cleng4){
		cout << "ERROR: lattice dimensions are different in config file and settings !" << endl;
		return 1;
	}

	int n;
        for(int t=0;t<leng4;t++)
        for(int z=0;z<leng3;z++)
        for(int y=0;y<leng2;y++)
        for(int x=0;x<leng1;x++){
                n = x + y*leng1 + z*leng1*leng2 + t*leng1*leng2*leng3;
                for(int mu=0;mu<DIM;mu++){
                        for(int i=0; i<matrixdim; i++){
                                for(int j=0; j<matrixdim; j++){
                                        fconfig >> realnum >> imagnum;
					U[i][j]=complex<double>(realnum,imagnum);
                                }
                        }
			replace(*U,mu,n);
                }
        }

	getline(fconfig, strtmp);
	getline(fconfig, strtmp);
	getline(fconfig, strtmp);
	fconfig >> pollrefreal >> pollrefimag >> dtmp >> dtmp >> plaqrefreal;

	pollref = complex<double>(pollrefreal/3.0, pollrefimag/3.0);
	plaqref = complex<double>(plaqrefreal/6.0, 0);
        
	// Calculation of Plaquette and Pollyakov loop
	poll=calcPoll();
	plaq=calcPlaq();

	double abspollerr=abs(pollref-poll);
	double absplaqerr=abs(real(plaqref)-real(plaq));

	if(verbose){
		cout << "Configfile " << fconfigname << " read." << endl << endl;
		cout << "Checks (note, MILC output only single precision):" << endl;
		cout << "Polyakov loop: file=" << pollref << " calculated=" << poll << " diff(abs)=" << abspollerr << endl;
		cout << "Plaquettes: file=" << real(plaqref) << " calculated=" << real(plaq) << " diff(abs)=" << absplaqerr << endl << endl;
	}

	if(abspollerr>maxError || absplaqerr>maxError)
		cout << "WARNING: Calculated values for Polyakov loop and/or Plaquettes differ!" << endl << endl;

        return 0;
}

int **ConfigData::Create2D(int row, int col)
{
        int **p = new int* [row];
        for (int j = 0; j < row; j ++)
                p[j] = new int[col];
        return p;
}

void ConfigData::Delete2D(int **p, int row)
{
        for (int j = 0; j < row; j ++)
                delete [] p[j];
        delete [] p;
}

/* ----------------------------------------------------------------------- */
/* -------- Observable calculations used for checking the config---------- */
/* ----------------------------------------------------------------------- */

void ConfigData::staplesum(complex<double> *S, int mu,int x){ 
        int xpmu,xpnu,xmnupmu,xmnu;
        complex<double> U1[matrixdim][matrixdim], U2[matrixdim][matrixdim], 
		U3[matrixdim][matrixdim], U12[matrixdim][matrixdim], U123[matrixdim][matrixdim];

        xpmu=neibArray[x][mu];

        for(int i=0;i<matrixdim;i++){
                for(int j=0;j<matrixdim;j++){
                        S[i*matrixdim + j]=complex<double>(0,0);
                }   
        }   
    
        for(int nu=0;nu<DIM;nu++){
                if(mu != nu){
                        xpnu=neibArray[x][nu];
                        extract(*U1,nu,xpmu);
                        extract(*U2,mu,xpnu);
                        extract(*U3,nu,x);

                        axbdag(*U12,*U1,*U2);
                        axbdag(*U123,*U12,*U3);
    
                        apb(S,*U123);

                        xmnu = neibArray[x][4+nu];
                        xmnupmu=neibArray[xmnu][mu];

                        extract(*U1,nu,xmnupmu);
                        extract(*U2,mu,xmnu);
                        extract(*U3,nu,xmnu);
    
                        adagxbdag(*U12,*U1,*U2);
                        axb(*U123,*U12,*U3);

                        apb(S,*U123);
                }   
        }   
}

complex<double> ConfigData::calcPoll(){
        // Calculates the Polyakov loop spatial average
        complex<double> poll, trace;
        poll=complex<double> (0,0);
        trace=complex<double> (0,0);

        complex<double> up[matrixdim][matrixdim], uu[matrixdim][matrixdim], upaux[matrixdim][matrixdim];
    
        int is0=0;
        int i4=0;

        for(int i1=0;i1<leng1;i1++){
                for(int i2=0;i2<leng2;i2++){
                        for(int i3=0;i3<leng3;i3++){
                                i4 = 0;
                                is0 = i1 + i2*leng1 + i3*leng1*leng2 + i4*leng1*leng2*leng3;

                                extract(*up,3,is0);

                                for(i4=1;i4<leng4-1;i4++){
                                        is0 = i1 + i2*leng1 + i3*leng1*leng2 + i4*leng1*leng2*leng3;
                                        extract(*uu,3,is0);
                                        axb(*upaux,*up,*uu);
                                        aeb(*up,*upaux);
                                }   

                                i4 = leng4-1;
                                is0 = i1 + i2*leng1 + i3*leng1*leng2 + i4*leng1*leng2*leng3;

                                extract(*uu,3,is0);
                                trace=multtrace(*up,*uu);

                                poll = poll + trace;
                        }   
                }   
        }   

        poll = poll/((double)matrixdim*leng1*leng2*leng3);

        return poll;
}

complex<double> ConfigData::calcPlaq(){
        // Calculates the plaquette spatial average over all 6N plaquettes
        int ispmu, ispnu;
        complex<double> sumplaqs = complex<double>(0,0), trace;
        complex<double> u1[matrixdim][matrixdim], u2[matrixdim][matrixdim],u3[matrixdim][matrixdim], u4[matrixdim][matrixdim], u23[matrixdim][matrixdim], u234[matrixdim][matrixdim];

        for(int is = 0;is<nsite;is++){
                for(int imu = 0;imu<DIM;imu++){
                        for(int inu = imu+1;inu<DIM;inu++){
                                ispmu = neibArray[is][imu];
                                ispnu = neibArray[is][inu];

                                extract(*u1,imu,is);
                                extract(*u2,inu,ispmu);
                                extract(*u3,imu,ispnu);
                                extract(*u4,inu,is);

                                axbdag(*u23,*u2,*u3);
                                axbdag(*u234,*u23,*u4);
                                trace=multtrace(*u1,*u234);

                                sumplaqs = sumplaqs + trace;
                        }
                }
        }

        sumplaqs=sumplaqs/((double)6*matrixdim*(double)nsite); // factor because sum over N lattice       points and md from trace

        return sumplaqs;
}

void ConfigData::z3rot(){
	// Init random number generator
	srand ( time(NULL) );
	
	// Performes a Z_3 center rotation (for SU(3) only)

	int i1, i2, i3, i4, nphase, is;
	double ran1=0, r=0;
	complex<double> uu[matrixdim][matrixdim], phase;

	// ................. 1st Volume ..............................	
	r=(double)rand()/((double)RAND_MAX+1);
	ran1=r*(double)leng1+0.5;
	i1 = floor(ran1);

	if(i1 == leng1)
		i1 = 0;

	r=(double)rand()/((double)RAND_MAX+1);
	nphase =floor(r*(double)3.0+0.5);
	phase = complex<double>(cos(2*M_PI*nphase/3.0),sin(2*M_PI*nphase/3.0));

	for(i2 = 0;i2<leng2;i2++){
		for(i3 = 0;i3<leng3;i3++){
			for(i4 = 0;i4<leng4;i4++){
				is = i1 + i2*leng1 + i3*leng1*leng2 + i4*leng1*leng2*leng3;
                                if(is>nsite){
                                        cout << "Error!!!" << endl;
                                        cout << "Loop 1st volume" << endl;
                                        cout << "i1=" << i1 << " i2=" << i2 << " i3=" << i3 << " i4=" << i4 << endl;
                                }
				extract(*uu,0,is);
				for(int i=0;i<matrixdim;i++)
 			       	for(int j=0;j<matrixdim;j++){
                			uu[i][j]=phase*uu[i][j];
        			}
				replace(*uu,0,is);
			}
		}
	}

	// ................. 2nd Volume ..............................	
	r=(double)rand()/((double)RAND_MAX+1);
	ran1=r*(double)leng2+0.5;
	i2 = floor(ran1);

	if(i2 == leng2)
		i2 = 0;

	r=(double)rand()/((double)RAND_MAX+1);
	nphase =floor(r*(double)3.0+0.5);
	phase = complex<double>(cos(2*M_PI*nphase/3.0),sin(2*M_PI*nphase/3.0));

	for(i1 = 0;i1<leng1;i1++){
		for(i3 = 0;i3<leng3;i3++){
			for(i4 = 0;i4<leng4;i4++){
				is = i1 + i2*leng1 + i3*leng1*leng2 + i4*leng1*leng2*leng3;
                                if(is>nsite){
                                        cout << "Error!!!" << endl;
                                        cout << "Loop 1st volume" << endl;
                                        cout << "i1=" << i1 << " i2=" << i2 << " i3=" << i3 << " i4=" << i4 << endl;
                                }
				extract(*uu,1,is);
				for(int i=0;i<matrixdim;i++)
 			       	for(int j=0;j<matrixdim;j++){
                			uu[i][j]=phase*uu[i][j];
        			}
				replace(*uu,1,is);
			}
		}
	}

	// ................. 3rd Volume ..............................	
	r=(double)rand()/((double)RAND_MAX+1);
	ran1=r*(double)leng3+0.5;
	i3 = floor(ran1);

	if(i3 == leng3)
		i3 = 0;

	r=(double)rand()/((double)RAND_MAX+1);
	nphase =floor(r*(double)3.0+0.5);
	phase = complex<double>(cos(2*M_PI*nphase/3.0),sin(2*M_PI*nphase/3.0));

	for(i1 = 0;i1<leng1;i1++){
		for(i2 = 0;i2<leng2;i2++){
			for(i4 = 0;i4<leng4;i4++){
				is = i1 + i2*leng1 + i3*leng1*leng2 + i4*leng1*leng2*leng3;
                                if(is>nsite){
                                        cout << "Error!!!" << endl;
                                        cout << "Loop 1st volume" << endl;
                                        cout << "i1=" << i1 << " i2=" << i2 << " i3=" << i3 << " i4=" << i4 << endl;
                                }
				extract(*uu,2,is);
				for(int i=0;i<matrixdim;i++)
 			       	for(int j=0;j<matrixdim;j++){
                			uu[i][j]=phase*uu[i][j];
        			}
				replace(*uu,2,is);
			}
		}
	}

	// ................. 4th Volume ..............................	
	r=(double)rand()/((double)RAND_MAX+1);
	ran1=r*(double)leng4+0.5;
	i4 = floor(ran1);

	if(i4 == leng4)
		i4 = 0;

	r=(double)rand()/((double)RAND_MAX+1);
	nphase =floor(r*(double)3.0+0.5);
	phase = complex<double>(cos(2*M_PI*nphase/3.0),sin(2*M_PI*nphase/3.0));

	for(i1 = 0;i1<leng1;i1++){
		for(i2 = 0;i2<leng2;i2++){
			for(i3 = 0;i3<leng3;i3++){
				is = i1 + i2*leng1 + i3*leng1*leng2 + i4*leng1*leng2*leng3;
                                if(is>nsite){
                                        cout << "Error!!!" << endl;
                                        cout << "Loop 1st volume" << endl;
                                        cout << "i1=" << i1 << " i2=" << i2 << " i3=" << i3 << " i4=" << i4 << endl;
                                }
				extract(*uu,3,is);
				for(int i=0;i<matrixdim;i++)
 			       	for(int j=0;j<matrixdim;j++){
                			uu[i][j]=phase*uu[i][j];
        			}
				replace(*uu,3,is);
			}
		}
	}
}

/* ----------------------------------------------------------------------- */
/* -------- Matrix functions --------------------------------------------- */
/* ----------------------------------------------------------------------- */

void ConfigData::za(complex<double> *c, double z, complex<double> *a){
	for(int i=0;i<matrixdim;i++)
	for(int j=0;j<matrixdim;j++){
		c[i*matrixdim + j]=complex<double>(z,0)*a[i*matrixdim + j];
	}
	
}

void ConfigData::aeb(complex<double> *a, complex<double> *b){
	for(int i=0;i<matrixdim;i++)
	for(int j=0;j<matrixdim;j++){
		a[i*matrixdim + j]=b[i*matrixdim + j];
	}
}

void ConfigData::amb(complex<double> *c,complex<double> *a, complex<double> *b){
	for(int i=0;i<matrixdim;i++)
	for(int j=0;j<matrixdim;j++){
		c[i*matrixdim + j]=a[i*matrixdim + j] - b[i*matrixdim + j];
	}
}

void ConfigData::apb(complex<double> *c, complex<double> *a){
	for(int i=0;i<matrixdim;i++)
	for(int j=0;j<matrixdim;j++){
		c[i*matrixdim + j]=a[i*matrixdim + j] + c[i*matrixdim + j];
	}
}

void ConfigData::capb(complex<double> *c, complex<double> *a, complex<double> *b){
	for(int i=0;i<matrixdim;i++)
	for(int j=0;j<matrixdim;j++){
		c[i*matrixdim + j]=a[i*matrixdim + j] + b[i*matrixdim + j];
	}
}

void ConfigData::axb(complex<double> *c, complex<double> *a, complex<double> *b){
	for(int i=0;i<matrixdim;i++)
	for(int j=0;j<matrixdim;j++){
		c[i*matrixdim + j]=0;
		for(int k=0;k<matrixdim;k++)
			c[i*matrixdim + j]=c[i*matrixdim + j] + a[i*matrixdim + k]*b[k*matrixdim + j];
	}
}

void ConfigData::axbdag(complex<double> *c, complex<double> *a, complex<double> *b){
	for(int i=0;i<matrixdim;i++)
	for(int j=0;j<matrixdim;j++){
		c[i*matrixdim + j]=0;
		for(int k=0;k<matrixdim;k++)
			c[i*matrixdim + j]=c[i*matrixdim + j] + a[i*matrixdim + k]*conj(b[j*matrixdim + k]);
	}
}

void ConfigData::adagxb(complex<double> *c, complex<double> *a, complex<double> *b){
	for(int i=0;i<matrixdim;i++)
	for(int j=0;j<matrixdim;j++){
		c[i*matrixdim + j]=0;
		for(int k=0;k<matrixdim;k++)
			c[i*matrixdim + j]=c[i*matrixdim + j]+conj(a[k*matrixdim + i])*b[k*matrixdim + j];
	}
}

void ConfigData::adagxbdag(complex<double> *c, complex<double> *a, complex<double> *b){
	for(int i=0;i<matrixdim;i++)
	for(int j=0;j<matrixdim;j++){
		c[i*matrixdim + j]=0;
		for(int k=0;k<matrixdim;k++)
			c[i*matrixdim + j]=c[i*matrixdim + j]+conj(a[k*matrixdim + i])*conj(b[j*matrixdim + k]);
	}
}

complex<double> ConfigData::multtrace(complex<double> *a, complex<double> *b){
	complex <double> tr;
	tr=complex<double>(0,0);

	for(int i=0;i<matrixdim;i++)
	for(int k=0;k<matrixdim;k++){
			tr=tr+a[i*matrixdim + k]*b[k*matrixdim + i];
	}
	return tr;
}

void ConfigData::adag(complex<double> *a){
       complex<double> tmp;
       for(int i=0;i<matrixdim;i++){
               a[i*matrixdim + i]=conj(a[i*matrixdim + i]);
       }
       
       tmp=a[0*matrixdim + 1];
       a[0*matrixdim + 1]=conj(a[1*matrixdim + 0]);
       a[1*matrixdim + 0]=conj(tmp);

       tmp=a[0*matrixdim + 2];
       a[0*matrixdim + 2]=conj(a[2*matrixdim + 0]);
       a[2*matrixdim + 0]=conj(tmp);
       
       tmp=a[1*matrixdim + 2];
       a[1*matrixdim + 2]=conj(a[2*matrixdim + 1]);
       a[2*matrixdim + 1]=conj(tmp);
}

/* ----------------------------------------------------------------------- */
/* -------- Test functions ----------------------------------------------- */
/* ----------------------------------------------------------------------- */

void ConfigData::runTests(){
	// Matrix operation tests
	complex<double> A[matrixdim][matrixdim], B[matrixdim][matrixdim], C[matrixdim][matrixdim];

	// Filling the test matrices
	A[0][0]=complex<double>( 1.2,-1.0); A[0][1]=complex<double>(-2.5, 5.1); A[0][2]=complex<double>(-6.1,-3.8);
	A[1][0]=complex<double>(-7.2, 6.1); A[1][1]=complex<double>(-8.2,-1.2); A[1][2]=complex<double>( 1.9, 1.6);
	A[2][0]=complex<double>( 3.2,-8.1); A[2][1]=complex<double>( 6.5, 2.1); A[2][2]=complex<double>( 3.5,-3.3);
	
	B[0][0]=complex<double>( 2.1, 2.4); B[0][1]=complex<double>( 1.7, 5.1); B[0][2]=complex<double>( 9.4,-6.1);
	B[1][0]=complex<double>(-6.2, 1.3); B[1][1]=complex<double>(-4.2,-3.2); B[1][2]=complex<double>(-3.9,-1.6);
	B[2][0]=complex<double>( 8.2, 1.6); B[2][1]=complex<double>(-4.5, 2.3); B[2][2]=complex<double>( 1.5, 3.3);

	// Starting tests
	cout << "---------------------------------------------" << endl;
	cout << "Config Data self test initiated..." << endl;
	cout << "---------------------------------------------" << endl;

	cout << "AxB:" << endl;
	cout << "Shall:" << endl;
	cout << "(-30.15,-75.01) (70.15,-5.93) (26.48,-58.44)\n(35.66,8.47) (-24.98,2.1) (-2.84,127.73)\n(17.11,-35.36) (18.01,-4.17) (-25.18,-107.65)" << endl;
	axb(*C, *A, *B);
	cout << "Is:" << endl;
	printMatrix(*C, matrixdim);
	cout << "---------------------------------------------" << endl;

	cout << endl << "AxBdag:" << endl;
	cout << "Shall:" << endl;
	cout << "(-12.28,-56.49) (15.31,-19.72) (9.53,-12.89)\n(-12.44,96.5) (80.88,-52.86) (-7.01,81.93)\n(62.07,-63.94) (-72.76,76.51) (-16.78,-112.44)" << endl;
	axbdag(*C, *A, *B);
	cout << "Is:" << endl;
	printMatrix(*C, matrixdim);
	cout << "---------------------------------------------" << endl;
	
	cout << endl << "AdagxB:" << endl;
	cout << "Shall:" << endl;
	cout << "(65.97,104.98) (-25.37,27.39) (13.77,60.1)\n(112.93,-41.63) (35.62,24.18) (-4.03,-5.95)\n(-8.21,38.39) (-66.19,-30.81) (-49.77,92.63)" << endl;
	adagxb(*C, *A, *B);
	cout << "Is:" << endl;
	printMatrix(*C, matrixdim);
	cout << "---------------------------------------------" << endl;
	
	cout << endl << "AdagxBdag:" << endl;
	cout << "Shall:" << endl;
	cout << "(-57.76,121.23) (18.18,-31.65) (61.34,51.88)\n(48.6,59.06) (17.48,22.18) (13.82,-48.96)\n(4.15,62.58) (20.97,-10.1) (-40.03,37.15)" << endl;
	adagxbdag(*C, *A, *B);
	cout << "Is:" << endl;
	printMatrix(*C, matrixdim);
	cout << "---------------------------------------------" << endl;

	cout << endl << "Tr(A*B):" << endl;
	cout << "Shall:" << endl;
	cout << "(-80.31,-180.56)" << endl;
	cout << "Is:" << endl;
	cout << multtrace(*A,*B) << endl;
	cout << "---------------------------------------------" << endl;

	cout << endl << "Adag:" << endl;
	cout << "Shall:" << endl;
	cout << "(1.2,1) (-7.2,-6.1) (3.2,8.1)\n(-2.5,-5.1) (-8.2,1.2) (6.5,-2.1)\n(-6.1,3.8) (1.9,-1.6) (3.5,3.3)" << endl;
	adag(*A);
	cout << "Is:" << endl;
	printMatrix(*A, matrixdim);

	cout << "---------------------------------------------" << endl;
	cout << "ConfigData self test finished!" << endl;
	cout << "---------------------------------------------" << endl;
}

void ConfigData::printMatrix(complex<double> *A, int md){
	for(int i=0;i<md;i++){
		for(int j=0;j<md;j++){
			cout << A[i*matrixdim + j] << " ";
		}
		cout << endl;
	} 
}

#endif /* CONFIGFUNCTIONS_H */

