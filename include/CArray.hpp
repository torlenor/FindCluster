// v0.3 - 2012-10-11 - 13:25
#ifndef CARRAY_H
#define CARRAY_H
#include <iostream>
#include <fstream>
#include <complex>
#include <limits>

class CArray{
	typedef unsigned int uint;

	public:
		CArray(uint len1);
		CArray(uint len1, uint len2);
		CArray(uint len1, uint len2, uint len3);
		CArray(uint len1, uint len2, uint len3, uint len4);
		CArray(uint len1, uint len2, uint len3, uint len4, uint len5);
		CArray(uint len1, uint len2, uint len3, uint len4, uint len5, uint len6);

		~CArray();
		
		void clear();

		int bwrite(const std::string &filename);
		int bread(const std::string &filename);
		int write(const std::string &filename);
		int read(const std::string &filename);

		std::complex<double> & operator()(uint n1);
		std::complex<double> & operator()(uint n1, uint n2);
		std::complex<double> & operator()(uint n1, uint n2, uint n3);
		std::complex<double> & operator()(uint n1, uint n2, uint n3, uint n4);
		std::complex<double> & operator()(uint n1, uint n2, uint n3, uint n4, uint n5);
		std::complex<double> & operator()(uint n1, uint n2, uint n3, uint n4, uint n5, uint n6);

		CArray& operator=(std::complex<double> C);
		CArray& operator+(std::complex<double> C);
		CArray& operator-(std::complex<double> C);

	private:
		uint DIM;
		std::complex<double> *A;
		uint *leng;
		uint *shift;
		uint nindex;
		void createStorage(uint N);

		void checkBounds(uint n1);
		void checkBounds(uint n1, uint n2);
		void checkBounds(uint n1, uint n2, uint n3);
		void checkBounds(uint n1, uint n2, uint n3, uint n4);
		void checkBounds(uint n1, uint n2, uint n3, uint n4, uint n5);
		void checkBounds(uint n1, uint n2, uint n3, uint n4, uint n5, uint n6);
};

CArray::~CArray(){
	delete [] A;
	A=0;
}

CArray& CArray::operator=(std::complex<double> C){
	for(uint n=0;n<nindex;n++){
		A[n] = C;
	}

	return *this;
}

CArray& CArray::operator+(std::complex<double> C){
	for(uint n=0;n<nindex;n++){
		A[n] += C;
	}

	return *this;
}

CArray& CArray::operator-(std::complex<double> C){
	for(uint n=0;n<nindex;n++){
		A[n] -= C;
	}

	return *this;
}

int CArray::bwrite(const std::string &filename){
       	FILE* pFile;
       	pFile = fopen(filename.c_str(), "wb");
       	int bytes = fwrite(A, sizeof(std::complex<double>),nindex, pFile);
       	fclose(pFile);

	return nindex - bytes;
}

int CArray::bread(const std::string &filename){
       	FILE* pFile;
       	pFile = fopen(filename.c_str(), "rb");
       	int bytes = fread(A, sizeof(std::complex<double>),nindex, pFile);
       	fclose(pFile);

	return nindex - bytes;
}

int CArray::write(const std::string &filename){
	std::ofstream file;
	file.open(filename.c_str());

	file << DIM << " ";
	
	for(uint d=0;d<DIM;d++){
		file << leng[d] << " ";
	}
	file << std::endl;

        file.flags (std::ios::scientific);
	file.precision(std::numeric_limits<double>::digits10 + 1);

	for(uint n=0;n<nindex;n++){
		file << real(A[n]) << " " << imag(A[n]) << " ";
	}

	file.close();

	return 0;
}

int CArray::read(const std::string &filename){
	uint dimfile;
	std::ifstream file;
	file.open(filename.c_str());
	
	file >> dimfile;

	if(dimfile != DIM){
		std::cout << "ERROR: File to read has wrong dimension!" << std::endl;
		return 1;
	}
	
	uint lengfile[dimfile];

	for(uint d=0;d<dimfile;d++){
		file >> lengfile[d];
		if(lengfile[d] != leng[d]){
			std::cout << "ERROR: File to read has wrong dimension size!" << std::endl;
			return 1;
		}
	}

	for(uint n=0;n<DIM;n++){
		file >> real(A[n]) >> imag(A[n]);
	}

	file.close();

	return 0;
}

CArray::CArray(uint len1){
	DIM=1;
	leng = new uint[1];
	leng[0]=len1;

	shift = new uint[1];
	shift[0]=1;

	createStorage(1);
}

CArray::CArray(uint len1, uint len2){
	DIM=2;
	leng = new uint[2];
	leng[0]=len1;
	leng[1]=len2;
	
	shift=new uint[DIM];
	shift[0]=leng[1];

	createStorage(2);
}

CArray::CArray(uint len1, uint len2, uint len3){
	DIM=3;
	leng = new uint[3];
	leng[0]=len1;
	leng[1]=len2;
	leng[2]=len3;
	
	shift=new uint[DIM];
	shift[0]=leng[2];
	shift[1]=leng[1]*leng[2];

	createStorage(3);
}

CArray::CArray(uint len1, uint len2, uint len3, uint len4){
	DIM=4;
	leng = new uint[4];
	leng[0]=len1;
	leng[1]=len2;
	leng[2]=len3;
	leng[3]=len4;
	
	shift=new uint[DIM];
	shift[0]=leng[3];
	shift[1]=leng[2]*leng[3];
	shift[2]=leng[1]*leng[2]*leng[3];

	createStorage(4);
}

CArray::CArray(uint len1, uint len2, uint len3, uint len4, uint len5){
	DIM=5;
	leng = new uint[5];
	leng[0]=len1;
	leng[1]=len2;
	leng[2]=len3;
	leng[3]=len4;
	leng[4]=len5;
	
	shift=new uint[DIM];
	shift[0]=leng[4];
	shift[1]=leng[3]*leng[4];
	shift[2]=leng[2]*leng[3]*leng[4];
	shift[3]=leng[1]*leng[2]*leng[3]*leng[4];

	createStorage(5);
}

CArray::CArray(uint len1, uint len2, uint len3, uint len4, uint len5, uint len6){
	DIM=6;
	leng = new uint[6];
	leng[0]=len1;
	leng[1]=len2;
	leng[2]=len3;
	leng[3]=len4;
	leng[4]=len5;
	leng[5]=len6;
	
	shift=new uint[DIM];
	shift[0]=leng[5];
	shift[1]=leng[4]*leng[5];
	shift[2]=leng[3]*leng[4]*leng[5];
	shift[3]=leng[2]*leng[3]*leng[4]*leng[5];
	shift[4]=leng[1]*leng[2]*leng[3]*leng[4]*leng[5];

	createStorage(6);
}
		
std::complex<double> & CArray::operator()(uint n1){
	checkBounds(n1);
	return A[n1];
}

std::complex<double> & CArray::operator()(uint n1, uint n2){
	checkBounds(n1, n2);
	return A[n2 + n1*shift[0]];
}

std::complex<double> & CArray::operator()(uint n1, uint n2, uint n3){
	checkBounds(n1, n2, n3);
	return A[n3 + n2*shift[0] + n1*shift[1]];
}

std::complex<double> & CArray::operator()(uint n1, uint n2, uint n3, uint n4){
	#ifdef DEBUG
		checkBounds(n1, n2, n3, n4);
	#endif
	return A[n4 + n3*shift[0] + n2*shift[1] + n1*shift[2]];
}

std::complex<double> & CArray::operator()(uint n1, uint n2, uint n3, uint n4, uint n5){
	#ifdef DEBUG
		checkBounds(n1, n2, n3, n4, n5);
	#endif
	return A[n5 + n4*shift[0] + n3*shift[1] + n2*shift[2] + n1*shift[3]];
}

std::complex<double> & CArray::operator()(uint n1, uint n2, uint n3, uint n4, uint n5, uint n6){
	#ifdef DEBUG
		checkBounds(n1, n2, n3, n4, n5, n6);
	#endif
	return A[n6 + n5*shift[0] + n4*shift[1] + n3*shift[2] + n2*shift[3] + n1*shift[4]];
}

void CArray::checkBounds(uint n1){
	if(n1 >= leng[0])
		std::cout << "ERROR: Out of array bounds" << std::endl;
}

void CArray::checkBounds(uint n1, uint n2){
	if(n1 >= leng[0] || n2 >= leng[1])
		std::cout << "ERROR: Out of array bounds" << std::endl;
}

void CArray::checkBounds(uint n1, uint n2, uint n3){
	if(n1 >= leng[0] || n2 >= leng[1] || n3 >= leng[2])
		std::cout << "ERROR: Out of array bounds" << std::endl;
}

void CArray::checkBounds(uint n1, uint n2, uint n3, uint n4){
	if(n1 >= leng[0] || n2 >= leng[1] || n3 >= leng[2] || n4 >= leng[3])
		std::cout << "ERROR: Out of array bounds" << std::endl;
}

void CArray::checkBounds(uint n1, uint n2, uint n3, uint n4, uint n5){
	if(n1 >= leng[0] || n2 >= leng[1] || n3 >= leng[2] || n4 >= leng[3] || n5 >= leng[4])
		std::cout << "ERROR: Out of array bounds" << std::endl;
}

void CArray::checkBounds(uint n1, uint n2, uint n3, uint n4, uint n5, uint n6){
	if(n1 >= leng[0] || n2 >= leng[1] || n3 >= leng[2] || n4 >= leng[3] || n5 >= leng[4] || n6 >= leng[5])
		std::cout << "ERROR: Out of array bounds" << std::endl;
}

void CArray::clear(){
	for(uint n=0;n<nindex;n++){
		A[n]=0;
	}
}

void CArray::createStorage(uint N){
	nindex=1;
	for(uint i=0;i<N;i++){
		nindex *= leng[i];
	}

	A = new std::complex<double>[nindex];

	clear();	
}

#endif // CARRAY_HPP
