#ifndef FINDCLUSTER_INIT_HPP
#define FINDCLUSTER_INIT_HPP

#include <getopt.h>
#include <cstdio>
#include <string.h>

char texthelp[]="Usage: findcluster.x [OPTION] ... [POLLEVCONFIG/POLLEVCONFIGLIST]\n"
		"Finds cluster and performs calculations with it.\n" 
		"Uses Polyakov loop eigenvalues from file POLLEVCONFIG as input.\n"
		"If -n/--nmeas > 1, the input file is interpreted as a list of input files!\n"
		"\n"
		"Mandatory arguments to long options are mandatory for short options too.\n"
		"  -s, --Ns SSIZE             spatial lattice extent (default = 4)\n"
		"  -t, --Nt TSIZE             temporal lattice extent (default = 4)\n"
		"  -f, --fraction FRACTION    fraction (default = 1.0)\n"
		"  -n, --nmeas NMEAS          number of configurations (default = 1)\n"
		"  -d, --detail               write detailed output for every calculated configuration\n"
		"\n"  
		"  -h  --help                 display this help and exit\n"
		"  -v  --version              output version information and exit\n"
		"\n"
		"Exit status:\n"
		" 0  if OK,\n"
		" 1  if minor problems,\n"
		" 2  if serious trouble.\n"
		"\n"
		"Report bugs to hps@abyle.org\n";

int init(int &argc, char *argv[]){

	cout << endl;
	cout << "findcluster.x " << MAJOR_VERSION << "." << MINOR_VERSION << "." << REVISION_VERSION << endl
		<< "Finds cluster and performs calculations with it." << endl
		<< "Uses Polyakov loop eigenvalues as input." << endl << endl;


	cout << "Initializing... " << endl;

	if(argc<2){
		cout << endl << texthelp << endl;
		return 2;
	}

       	int c;
       	
       	while (1){
		static struct option long_options[] =
			{
			/* These options don't set a flag.
			We distinguish them by their indices. */
			{"Ns", required_argument, 0, 's'},
			{"Nt", required_argument, 0, 't'},
			{"fraction", required_argument, 0, 'f'},
			{"nmeas", required_argument, 0, 'n'},
			{"detail", no_argument, 0, 'd'},
			/* These options set a flag. */
			// {"free", no_argument, 0, 'f'},
			// {"u0", required_argument, 0, 0},
			{"help", no_argument, 0, 'h'},
			{"version", no_argument, 0, 'v'},
			{0, 0, 0, 0}
			};
			
		/* getopt_long stores the option index here. */
		int option_index = 0;

		c = getopt_long (argc, argv, "s:t:f:n:hvd",
		long_options, &option_index);

		/* Detect the end of the options. */
		if (c == -1)
			break;

		switch (c){
			case 0:
				/* If this option set a flag, do nothing else now. */
				/*if (long_options[option_index].flag != 0)
					break;
				printf ("option %s", long_options[option_index].name);
				if (optarg)
					printf (" with arg %s", optarg);
				printf ("\n");A */
		            //    if( strcmp( "u0", long_options[option_index].name ) == 0 ) {
			    //                u0 = atof(optarg);
			    //	}
				break;

			case 's':
				Ns = atoi(optarg);
				break;

			case 't':
				Nt = atoi(optarg);
				break;

			case 'f':
				fraction = atof(optarg);
				break;

			case 'n':
				nmeas = atoi(optarg);
				break;

			case 'd':
				detail = true;
				break;

			case 'v':
				cout << endl << "findcluster.x version " << MAJOR_VERSION << "." << MINOR_VERSION << "." << REVISION_VERSION << endl;
				abort();

			case 'h':
				cout << endl << texthelp << endl;
				abort();

			default:
				cout << endl << texthelp << endl;
				abort();
		}
	}

	string finname;
	/* Print any remaining command line arguments (not options). */
	if (optind < argc)
	{
		finname = argv[optind];
	}else{
		cout << "ERROR: No configuration file specified!" << endl;
	}
	
	delta0 = M_PI/3.0;
	delta = delta0*fraction;

	// Prepare for nmeas measurements
	clusterdata = new Clusterstruct[nmeas];
	obs = new Observablestruct[nmeas];
	fevname.resize(nmeas);

	if(nmeas>1){
		// Read finname file and fill fevname
		ifstream fin;
		fin.open(finname.c_str());
		if(fin.is_open()!=true){
                	cout  << "ERROR: File " << finname <<  " to read configuration filename list could not be opened!" << endl;
			throw 1;
        	}
        	
        	string strtmp;
        	int n=0;
        	while(n<nmeas && getline(fin, fevname.at(n)) ){
			n++;
		}
		if(n<nmeas){
			cout << "ERROR: Only found " << n << " names in " << finname << " !" << endl;
			return 1;
		}
	}else{
		fevname[0] = finname;
	}

	leng1=Ns; leng2=Ns; leng3=Ns; leng4=Nt;
	
	Nspace = Ns*Ns*Ns;

	neib.resize(Nspace);
	for(int is=0;is<Nspace;is++)
		neib[is].resize(6);

	fillNeib();

	#ifdef DEBUG
	cout << "Allocating 3d (Ns^3 * 3) lattice for Polyakov loop evs... " << flush;
	#endif
	pollev.resize(Nspace);
	for(int is=0;is<Nspace;is++)
		pollev[is].resize(matrixdim);
	#ifdef DEBUG
	cout << "done!" << endl;
	#endif
	
	#ifdef DEBUG
	cout << "Creating cluster data arrays..." << flush;
	#endif
	for(int n=0;n<nmeas;n++){
		(&clusterdata[n])->isinsector.resize(Nspace);
		// (&clusterdata[0])->clustersector will not be allocated here, but on the fly with push_back
		(&clusterdata[n])->isincluster.resize(Nspace);
		// (&clusterdata[0])->clustermembers will not be allocated here, but on the fly with push_back
	}
	#ifdef DEBUG
	cout << "done!" << endl;
	#endif
	
	return 0;
}

#endif // FINDCLUSTER_INIT_HPP
