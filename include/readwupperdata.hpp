#ifndef READWUPPERDATA_HPP
#define READWUPPERDATA_HPP

#include <complex>
#include <fstream>
#include <iostream>
#include <vector>
#include <string>

// Für jedes Gitter gibt es dann Nx*Ny*Nz*2*4+16 bytes
// im folgenden Format:
//
// struct file {
//   int Nx,Ny,Nz,Nt;
//         complex float L[Nz][Ny][Nx];
// };

int readWupperPollBinary(int leng1, int leng2, int leng3, int leng4, int matrixdim, std::vector<std::vector<std::complex<double> > > &pollev, std::string fevname){
	#ifdef DEBUG
  std::cout << "Reading 3d lattice with Polyakov loop evs... " << std::flush;
	#endif
	
	int elems=0;
  FILE* pFile;
  
  int cleng1, cleng2, cleng3, cleng4;
	unsigned int nindex=leng1*leng2*leng3;
	int is;

	if (nindex > pollev.size() * pollev[pollev.size()-1].size()) {
    std::cout << "ERROR: Something wrong in readPollEvBinary!" << std::endl;
    std::cout << "Want to read nindex = " << nindex << " entries, but has space for (only) " << pollev.size() * pollev[pollev.size()-1].size() << " entries!" << std::endl;
	}
	
 	pFile = fopen(fevname.c_str(), "rb");
  if (pFile == NULL) {
    perror ("Error opening file");
  } else {
		elems += fread(&cleng1, sizeof(int), 1, pFile);
		elems += fread(&cleng2, sizeof(int), 1, pFile);
		elems += fread(&cleng3, sizeof(int), 1, pFile);
		elems += fread(&cleng4, sizeof(int), 1, pFile);
		
		// Check lattice dimensions stored in file against settings
		if (leng1!=cleng1 || leng2!=cleng2 || leng3!=cleng3 || leng4!=cleng4) {
      std::cout << "ERROR: lattice dimensions are different in config file and settings !" << std::endl;
			return 1;
		}

    float re=0;
    float im=0;
		
		for (int i1=0; i1<leng1; i1++)
		for (int i2=0; i2<leng2; i2++)
		for (int i3=0; i3<leng3; i3++) {
			is = i1 + i2*leng1 + i3*leng1*leng2;
			elems += fread(&re, sizeof(float), 1, pFile);
			elems += fread(&im, sizeof(float), 1, pFile);
      pollev.at(is).at(0)=std::complex<double>(re,im);
      // We use only one entry as we get Polyakov loop values from Wuppertal instead of 3 eigenvalues
		}
		fclose(pFile);
	}

	#ifdef DEBUG
  std::cout << "done!" << std::endl;
	#endif

  return 2*nindex + 4 - elems; // factor 2 from complex and 4 are the 4 ints at the beginning
}

#endif // READWUPPERDATA_HPP
