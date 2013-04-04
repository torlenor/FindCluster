#ifndef OPTIMATRIXHANDLING_HPP
#define OPTIMATRIXHANDLING_HPP
#include <complex>

using namespace std;

// Functions which use templates, specify size with first parameter
template<class T>
void aeb(int size, T* a, T* b) {
        for(int i=0;i<size*size;i++){
                a[i]=b[i];
        }
}

template<class T>
void apb(int size, T* c, T* a) {
        for(int i=0;i<size*size;i++){
                c[i]+=a[i];
        }
}

template<class T,class Z>
void za(int size, T* c, Z z, T* a) {
        for(int i=0;i<size*size;i++){
                c[i]=z*a[i];
        }
}

template<class T>
void axb(int size, T* c, T* a, T* b) {
        for(int i=0;i<size;i++)
        for(int j=0;j<size;j++){
                c[i*size+j]=0;
                for(int k=0;k<size;k++)
                        c[i*size+j]+=a[i*size+k]*b[k*size+j];
                }
}

template<class T>
void capb(int size, T* c, T* a, T* b) {
        for(int i=0;i<size*size;i++){
		c[i]=a[i]+b[i];
        }
}

template<class T>
void amb(int size, T* c, T* a, T* b) {
        for(int i=0;i<size*size;i++){
		c[i]=a[i]-b[i];
        }
}

template<class T>
void axbdag(int size, T* c, T* a, T* b) {
        for(int i=0;i<size;i++)
        for(int j=0;j<size;j++){
                c[i*size+j]=0;
                for(int k=0;k<size;k++)
                        c[i*size+j]+=a[i*size+k]*conj(b[j*size+k]);
                }
}

template<class T>
void adagxb(int size, T* c, T* a, T* b) {
        for(int i=0;i<size;i++)
        for(int j=0;j<size;j++){
                c[i*size+j]=0;
                for(int k=0;k<size;k++)
                        c[i*size+j]+=conj(a[k*size+i])*b[k*size+j];
                }
}

template<class T>
void adagxbdag(int size, T* c, T* a, T* b) {
        for(int i=0;i<size;i++)
        for(int j=0;j<size;j++){
                c[i*size+j]=0;
                for(int k=0;k<size;k++)
                        c[i*size+j]+=conj(a[k*size+i]*b[j*size+k]);
                }
}

template<class T>
complex<double> multtrace(int size, T* a, T* b) {
	complex <double> tr;
	tr=complex<double>(0,0);

        for(int i=0;i<size;i++)
        for(int k=0;k<size;k++){
        	tr+=a[i*size+k]*b[k*size+i];
        }

	return tr;
}

// Functions for d3xd3 matrices

void adag(complex<double> a[3][3]){
	complex<double> tmp;
	for(int i=0;i<3;i++){
		a[i][i]=conj(a[i][i]);
	}
	
	tmp=a[0][1];
	a[0][1]=conj(a[1][0]);
	a[1][0]=conj(tmp);

	tmp=a[0][2];
	a[0][2]=conj(a[2][0]);
	a[2][0]=conj(tmp);
	
	tmp=a[1][2];
	a[1][2]=conj(a[2][1]);
	a[2][1]=conj(tmp);
}

#endif /*OPTIMATRIXHANDLING_HPP*/

