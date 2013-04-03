#ifndef MATRIXFUNCTIONS_H
#define MATRIXFUNCTIONS_H

#include <complex>

void za(complex<double> *c, double z, complex<double> *a, int matrixdim){
	for(int i=0;i<matrixdim;i++)
	for(int j=0;j<matrixdim;j++){
		c[i*matrixdim + j]=complex<double>(z,0)*a[i*matrixdim + j];
	}
}

void aeb(complex<double> *a, complex<double> *b, int matrixdim){
	for(int i=0;i<matrixdim;i++)
	for(int j=0;j<matrixdim;j++){
		a[i*matrixdim + j]=b[i*matrixdim + j];
	}
}

void amb(complex<double> *c,complex<double> *a, complex<double> *b, int matrixdim){
	for(int i=0;i<matrixdim;i++)
	for(int j=0;j<matrixdim;j++){
		c[i*matrixdim + j]=a[i*matrixdim + j] - b[i*matrixdim + j];
	}
}

void apb(complex<double> *c, complex<double> *a, int matrixdim){
	for(int i=0;i<matrixdim;i++)
	for(int j=0;j<matrixdim;j++){
		c[i*matrixdim + j]=a[i*matrixdim + j] + c[i*matrixdim + j];
	}
}

void capb(complex<double> *c, complex<double> *a, complex<double> *b, int matrixdim){
	for(int i=0;i<matrixdim;i++)
	for(int j=0;j<matrixdim;j++){
		c[i*matrixdim + j]=a[i*matrixdim + j] + b[i*matrixdim + j];
	}
}

void axb(complex<double> *c, complex<double> *a, complex<double> *b, int matrixdim){
	for(int i=0;i<matrixdim;i++)
	for(int j=0;j<matrixdim;j++){
		c[i*matrixdim + j]=0;
		for(int k=0;k<matrixdim;k++)
			c[i*matrixdim + j]=c[i*matrixdim + j] + a[i*matrixdim + k]*b[k*matrixdim + j];
	}
}

void axbdag(complex<double> *c, complex<double> *a, complex<double> *b, int matrixdim){
	for(int i=0;i<matrixdim;i++)
	for(int j=0;j<matrixdim;j++){
		c[i*matrixdim + j]=0;
		for(int k=0;k<matrixdim;k++)
			c[i*matrixdim + j]=c[i*matrixdim + j] + a[i*matrixdim + k]*conj(b[j*matrixdim + k]);
	}
}

void adagxb(complex<double> *c, complex<double> *a, complex<double> *b, int matrixdim){
	for(int i=0;i<matrixdim;i++)
	for(int j=0;j<matrixdim;j++){
		c[i*matrixdim + j]=0;
		for(int k=0;k<matrixdim;k++)
			c[i*matrixdim + j]=c[i*matrixdim + j]+conj(a[k*matrixdim + i])*b[k*matrixdim + j];
	}
}

void adagxbdag(complex<double> *c, complex<double> *a, complex<double> *b, int matrixdim){
	for(int i=0;i<matrixdim;i++)
	for(int j=0;j<matrixdim;j++){
		c[i*matrixdim + j]=0;
		for(int k=0;k<matrixdim;k++)
			c[i*matrixdim + j]=c[i*matrixdim + j]+conj(a[k*matrixdim + i])*conj(b[j*matrixdim + k]);
	}
}

complex<double> multtrace(complex<double> *a, complex<double> *b, int matrixdim){
	complex <double> tr;
	tr=complex<double>(0,0);

	for(int i=0;i<matrixdim;i++)
	for(int k=0;k<matrixdim;k++){
			tr=tr+a[i*matrixdim + k]*b[k*matrixdim + i];
	}
	return tr;
}

void adag(complex<double> *a, int matrixdim){
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

#endif //MATRIXFUNCTIONS_H
