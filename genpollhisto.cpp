#include <iostream>
#include <fstream>
#include <complex>

#include "include/ConfigData.hpp"

using namespace std;

#include "include/matrixfunctions.hpp"

int Ns=0;
int Nt=0;

ConfigData *config;

complex<double> localPoll(int i1, int i2, int i3);
complex<double> totalPoll();

int main(int argc, char *argv[]){
	if(argc<4){
		cout << "./cluster.x Ns Nt confin" << endl;
		return 1;
	}

	Ns=atoi(argv[1]);
	Nt=atoi(argv[2]);
	string fconfigname=argv[3];

	config = new ConfigData(Ns, Ns, Ns, Nt, 3);

	cout << "Reading the file..." << endl;
	config->readBinaryConfig2(fconfigname);

	complex<double> poll;

	poll = totalPoll();
	cout << "Total Polyakov loop P = " << poll << endl;

	cout << "Performing a random Z_3 rotation... " << flush;
	config->z3rot();
	cout << "done!" << endl;
	
	cout << "Writing local Polyakov loop phases... " << flush;
	ofstream file;
	file.open("phase.data");

	int cnt=0;
	for(int i1=0;i1<Ns;i1++)
	for(int i2=0;i2<Ns;i2++)
	for(int i3=0;i3<Ns;i3++){
		cnt++;
		file << cnt << " " << arg(localPoll(i1, i2, i3)) << endl;
	}

	file.close();

	cout << "done!" << endl;

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
