#ifndef HELPER_HPP
#define HELPER_HPP

#include "optimatrixhandling.hpp"

// Helper functions to calculate various things, especially related to the
// computer methods elementary particle physics.
//
// Functions implemented in this file:
//
// long double fak(int k):
// 	Factorial calculation.
// T uincgamma0(T z, int sterms=150)
// 	Upper incomplete gamma function for any type of arguments with parameter a=0.
// complex<double> uincgammam1(complex<double> z)
// 	Upper incomplete gamma function for complex arguments with parameter a=-1.
// long double uincgamma(double a, double x, int sterms=200)
// 	Upper incomplete Gamma funtion for real arguments and a>0.
// T IntSimpson(T *data, double h, int nData)
// 	Calculates the integral over function values stored in array
// 	data[nData] with nData equidistant data points and a step length of
// 	h using an extended simpsons rule.
// int ix(int row, int col, int len)
// 	Gives back the linearized element number of the ij's entry of a
// 	len X len matrix.
// void FindPrivot(double *M, int len, int &i, int &j)
// 	Finds the largest off-diagonal element of a symmetric matrix given by
// 	M, dim(M)=len. Position will be given back in variables i and j.
// void jacobi(double *matrix, int len, double *eval, double *evec)
// 	Jacobi matrix eigenvalues and eigenvectors calculation
// 	matrix linearly addressed, where size is the dimension of the matrix
// 	results will be written in eval (sorted eigenvalues) and evec
//	(eigenvectors in rows). The original matrix will be overwritten by the
//	rotated matrix matrixp!
// T jbessel(int N, T z)
// 	Spherical Bessel functions 1st kind using recurrence relations.

#include <iostream>
#include <cmath>
#include <complex>

using namespace std;

#define eulerconst 0.57721566490153286060651209008240243104215933593992

// Factorial
long double fak(int k){
        long double faku=k;

        if(faku==0)
                return 1.0;

        for(unsigned long int i=k-1;i>0;i--)
                faku *= (long double)i;

        return faku;
}

// Upper incomplete gamma function for any type of arguments with parameter a=0
template <class T>
T uincgamma0(T z, int sterms=150){
        T sum = -eulerconst - log(z);

        for(int i=1;i<sterms;i++){
                sum -= pow(-1.0,(double)i)*pow(z,(double)i)/((double)fak(i)*i);
        }

        return sum;
}

// Upper incomplete gamma function for complex arguments with parameter a=-1
complex<double> uincgammam1(complex<double> z){
        complex<double> summ1 = -uincgamma0(z) + 1.0/z*exp(-z);

        return summ1;
}

// Gamma funtion for real arguments
long double uincgamma(double a, double x, int sterms=200){
        long double sum;
        sum=tgamma(a);

        for(int i=0;i<sterms;i++){
                sum -= pow(-1.0,i)*pow(x,(a+i))/((long double)fak(i)*(a+i));

                // cout << fak(i) << " " << sum << endl; 
        }

        return sum;
}

// Calculates the integral over function values stored in array data[nData]
// with nData equidistant data points and a step length of h using an
// extended simpsons rule.
template <class T>
T IntSimpson(T *data, double h, int nData){
        // extended Simpsons rule
        T intvalue=0;
        intvalue=data[0]/3.0;
        for(int k=1;k<nData-1;k++){
                if(k % 2 == 0){
                        intvalue=intvalue+2.0*data[k]/3.0;
                }else{
                        intvalue=intvalue+4.0*data[k]/3.0;
                }
        }
        intvalue=intvalue + data[nData-1]/3.0;
        return intvalue*h;
}

// Gives back the linearized element number of the ij's entry of a len X len matrix.
int ix(int row, int col, int len){
	return col + row*len;
}

// Finds the largest off-diagonal element of a symmetric matrix given by
// M, dim(M)=len. Position will be given back in variables i and j.
void FindPrivot(double *M, int len, int &i, int &j){
	double lnumber=0;

	for(int r=0;r<len;r++){
		for(int c=r+1; c<len; c++){
			if(abs( M[ix(r,c,len)] ) > abs(lnumber) ){
				lnumber=M[ix(r,c,len)];
				i=r;
				j=c;
			}
		}
	}

}

// Jacobi matrix eigenvalues and eigenvectors calculation
// matrix linearly addressed, where size is the dimension of the matrix
// results will be written in eval (sorted eigenvalues) and evec (eigenvectors).
// The original matrix will be overwritten by the rotated matrix matrixp!
void jacobi(double *matrix, int len, double *eval, double *evec){
	double c=0,s=0, theta=0;
	double matrixp[len*len];
	double evecp[len*len];
	int i=0,j=0;

	int iter = 0;
	int itermax = 100000;

	double eps=1E-16;
	double nullerr=1E5;

	// Set evec to unity
	for(int i=0;i<len;i++){
		for(int j=0;j<len;j++){
			if(i==j){
				evec[ ix(i,j,len) ] = 1;
			}else{
				evec[ ix(i,j,len) ] = 0;
			}
		}
	}

	while(nullerr > eps && iter < itermax){
		FindPrivot(matrix, len, i, j);
		nullerr = abs(matrix[ix(i,j,len)]);

		// Calculate rotation angle theta
		theta = atan( 2.0*matrix[ix(i,j,len)]/(matrix[ix(j,j,len)] - matrix[ix(i,i,len)]) )/2.0;

		// Build the new matrix matrixp and then write it into matrix
		c=cos(theta);
		s=sin(theta);	
		// A'_ii
		matrixp[ ix(i,i,len) ] = c*c*matrix[ ix(i,i,len) ] - 2.0*s*c*matrix[ ix(i,j,len) ] + 
						s*s*matrix[ ix(j,j,len) ];
		// A'_jj
		matrixp[ ix(j,j,len) ] = s*s*matrix[ ix(i,i,len) ] + 2.0*s*c*matrix[ ix(i,j,len) ] + 
						c*c*matrix[ ix(j,j,len) ];
		// A'_ij = A'_ji
		matrixp[ ix(i,j,len) ] = (c*c - s*s)*matrix[ ix(i,j,len) ] 
						+ s*c*(matrix[ ix(i,i,len)] - matrix[ ix(j,j,len) ]);
		matrixp[ ix(j,i,len) ] = (c*c - s*s)*matrix[ ix(i,j,len) ]
						+ s*c*(matrix[ ix(i,i,len)] - matrix[ ix(j,j,len) ]);
		// A'_ik = A'_ki
		for(int k=0;k<len;k++){
			if(k!=j && k!=i){
				matrixp[ ix(i,k,len) ] = c*matrix[ ix(i,k,len) ] - s*matrix[ ix(j,k,len) ];
				matrixp[ ix(k,i,len) ] = c*matrix[ ix(i,k,len) ] - s*matrix[ ix(j,k,len) ];
				
				matrixp[ ix(j,k,len) ] = s*matrix[ ix(i,k,len) ] + c*matrix[ ix(j,k,len) ];
				matrixp[ ix(k,j,len) ] = s*matrix[ ix(i,k,len) ] + c*matrix[ ix(j,k,len) ];
			}
			for(int l=0;l<len;l++){
				if( k!=j && k!=i && l!=j && l!=i ){
					matrixp[ ix(k,l,len) ] = matrix[ ix(k,l,len) ];
				}
			}
		}

		// Rotate eigenvector matrix
		// v'_rp
		for(int k=0;k<len;k++){
			evecp[ ix(k,i,len) ]  = c*evec[ ix(k,i,len) ]
						- s*evec[ ix(k,j,len) ];
		}
		// v'_rq
		for(int k=0;k<len;k++){
			evecp[ ix(k,j,len) ]  = s*evec[ ix(k,i,len) ]
						+ c*evec[ ix(k,j,len) ];
		}
		// v'_kl
		for(int k=0;k<len;k++){
			for(int l=0;l<len;l++){
				if( l!=j && l!=i ){
					evecp[ ix(k,l,len) ] = evec[ ix(k,l,len) ];
				}
			}
		}
		
		for(int n=0;n<len*len;n++){
			matrix[n]=matrixp[n];
			evec[n]=evecp[n];
		}

		iter++;
	}

	if(iter>=itermax)
		cout << "itermax in jacobi routine hit!!! iter=" << iter << " nullerr=" << nullerr << endl;

	// Write the eigenvalues into eval and sort it by absolute value (decending)
	for(int n=0;n<len;n++){
		eval[n]=matrixp[ ix(n,n,len) ];
	}

	cout << iter << " iterations needed ( maximum = " << itermax << " )" << endl;
}

// Spherical Bessel functions 1st kind using recurrence relations
// WARNING: Order N=0 needs special handling in upward recurrence!!!
template <class T>
T jbessel(int N, T z){
	double nulleps=1E-12;
	int M=2*(N + sqrt(40*N));
	if(abs(z) < nulleps){
		return 0.0;
	}
	if(abs(z) > N){ // Check if upward or downward recurrence relation
		T j[N+1];
		j[0]=sin(z)/z;
		j[1]=( j[0] - cos(z) )/z;

		for(int i=2;i<N+1;i++)
			j[i]=0;

		// Handle the case that n=0,1
		switch(N){
			case 0:	return j[0];
				break; // just to be sure...
			case 1: return j[1];
				break; // just to be sure...
		}

		// Handle the case that n>1 using recurrence relations
		for(int n=1;n<N;n++)
			j[n+1]=(2.0*n+1.0)*j[n]/z - j[n-1];
	
		if(j[N]!=j[N]){
			return 0;
		}else{
			return j[N];
		}
	}else{	// downward recurrence relation
		T norm;
		T j[M+1];
		T j0=sin(z)/z;
		T j1=( j0 - cos(z) )/z;

		for(int i=0;i<M+1;i++)
			j[i]=0;

		j[M]=0;
		j[M-1]=1;

		// Handle the case that n=0,1
		switch(N){
			case 0:	return j0;
				break; // just to be sure...
			case 1: return j1;
				break; // just to be sure...
		}

		// Handle the case that n>1 using recurrence relations
		for(int m=M-2;m>=1;m--)
			j[m-1]=(2.0*m+1.0)*j[m]/z - j[m+1];

		norm=j[0]/j0;

		for(int m=0;m<M+1;m++){
		 	j[m]=j[m]/norm;
		}
		
		if(j[N]!=j[N]){
			return 0;
		}else{
			return j[N];
		}
	}
}

// Calculates roots of the first N+1 bessel functions j_i(x) up to x=xmax
// (or up to the maxzeros root, whatever comes first).x are real numbers.
// The values are saved in the array zeros[order bessel][# of root]
template <class T>
void findBesselRoots(int N, int maxzeros, T xmax, T factor, T *zeros){
	double eps=1E-8; // desired precision for roots
	int iterlimit=100; // max number of bisection steps

	int ix=0;

	for(int l=0;l<N+1;l++){ // Loop for order of bessel function
		for(int zero=0;zero<maxzeros;zero++){
			// Set all root values to 0
			ix = zero + l*maxzeros;
			zeros[ix]=0;
		}
	}

	// We know the roots for order l=0, so put them in
	// in this case, we do not care for a xmax and fill all maxzeros roots
	// (CHECK IF GOOD IDEA)
	for(int zero=0;zero<maxzeros;zero++){
		ix = zero;
		zeros[ix]=M_PI*(zero+1)/factor;
//		if(zeros[ix]>xmax)
//			break; // If we reach xmax, stop the loop
	}

	// Now search all roots based on the fact, that we have to find one root
	// for l'=l+1 between two roots of l.
	int cnt=0;
	T xnew=0;
	T x1=0,x2=0;
	T j1=0,j2=0,jnew=0;
	int zero=0;

	for(int l=1;l<N+1;l++){
		for(int zero=0;zero<maxzeros-1;zero++){
			cnt=0;
			ix= zero + l*maxzeros;
			x1=zeros[zero + (l-1)*maxzeros];
			x2=zeros[zero+1 + (l-1)*maxzeros];
			if(abs(x2)<eps)
				break;
			while(abs(x1-x2)>eps && cnt<iterlimit){
				xnew=0.5*(x1+x2);

				j1=jbessel(l,x1*factor);
				j2=jbessel(l,x2*factor);
				jnew=jbessel(l,xnew*factor);

				if( (j1<0 && j2<0) || (j1>0 && j2>0)){
					cout << "Wrong start values, both are on the same side with respect to pole!" << endl;
					cout << "j1: " << j1 << endl;
					cout << "j2: " << j2 << endl;
					break;
				}   

				if(j1*jnew<0){
					x2=xnew;
				}else{
					x1=xnew;
				}   
				cnt++;
			}
			zeros[ix]=xnew;
			if(zeros[ix]>xmax)
				break;
		}
	}
}

#endif /* HELPER_HPP */
