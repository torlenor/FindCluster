#ifndef FINDCLUSTER_INIT_HPP
#define FINDCLUSTER_INIT_HPP

#include <cstdio>
#include <getopt.h>
#include <string.h>

#include "findcluster.h"

char texthelp[]="Usage: findcluster.x [OPTION] ... [POLLEVCONFIG/POLLEVCONFIGLIST]\n"
		"Finds cluster and performs calculations with it.\n" 
		"Uses Polyakov loop eigenvalues from file POLLEVCONFIG as input.\n"
		"If -n/--nmeas > 1, the input file is interpreted as a list of input files!\n"
		"\n"
		"Mandatory arguments to long options are mandatory for short options too.\n"
		"  -s, --Ns SSIZE             spatial lattice extent (default = 4)\n"
		"  -t, --Nt TSIZE             temporal lattice extent (default = 4)\n"
		"  -f, --fraction FRACTION    fraction (default = 0.0)\n"
		"  -n, --nmeas NMEAS          number of configurations (default = 1)\n"
		"  -d, --detail               write detailed output for every calculated configuration\n"
		"  -b, --boxes                performes the box counting calculation (expensive)\n"
		"  -a, --distance             performes the distance traveled calculation (expensive)\n"
		"  -r, --rsector RADIUS       use alternative sector classification with r = RADIUS\n"
		"  -o, --olddata              read old Polyakov loop ev data instead of wuppertal data (NOT ENTIRELY TESTED!!!)\n"
		"  -w, --writemeas            writes all measurements to file\n"
		"  -F, --fastmode             Fast Mode: Calculates only radius of largest cluster\n"
		"  --3d                       write 3dcluster data files\n"
		"\n"  
		"  -h  --help                 display this help and exit\n"
		"  -v  --version              output version information and exit\n"
		"\n"
    "Status output:\n"
    " r ... reading Polyakov loop data\n"
    " c ... identifying center cluster\n"
    " p ... identifying percolating cluster\n"
    " s ... sorting center cluster\n"
    " o ... calculating observables\n"
    " b ... calculating box counting observables\n"
		"\n"
		"Exit status:\n"
		" 0  if OK,\n"
		" 1  if minor problems,\n"
		" 2  if serious trouble.\n"
		"\n"
		"Report bugs to hps@abyle.org\n";

int init(int &argc, char *argv[]){

cout << endl;
cout << "findcluster.x " << MAJOR_VERSION << "." << MINOR_VERSION << "." << REVISION_VERSION << " ~ " << __DATE__ << " " << __TIME__ << endl 
  << endl 
  << "Finds cluster and performs calculations with it." << endl 
  << "Uses Polyakov loop eigenvalues as input." << endl 
  << endl;

	cout << "Initializing... " << endl;

	if (argc<2) {
		cout << endl << texthelp << endl;
		return 2;
	}

  // Set default options
  opt.Ns=4, opt.Nt=4;
  opt.matrixdim=1, opt.leng1=opt.Ns, opt.leng2=opt.Ns, opt.leng3=opt.Ns, opt.leng4=opt.Nt, opt.Nspace=opt.Ns*opt.Ns*opt.Ns;

  opt.nmeas=1;

  opt.wupperdata=true; // Controlls if we want to read wuppertal data, default true

  opt.usealternativesectors=false;
  opt.r=0;

  opt.do3d=false;

  opt.detail=false; // Controlls if we want detailed information for every configuration
  opt.doboxes=false; // Controlls if we want box counting calculations
  opt.dodistance=false; // Controlls if we want box counting calculations
  opt.doradius=true; // Controlls if we want radius calculation
  opt.domean=false; // Controlls if we want mean distance traveled calculations
  opt.writemeas=false; // Controlls if we want to writeout all measurements

  opt.fastmode=false; // Fast Mode: Only calc. the radius of largest cluster for f determination

  opt.fraction = 0.0;

	const double delta0 = M_PI/3.0;
	opt.delta = delta0*(1.0 - opt.fraction);

  int c;
    
  while (1) {
    static struct option long_options[] =
    {
      {"Ns", required_argument, 0, 's'},
      {"Nt", required_argument, 0, 't'},
      {"fraction", required_argument, 0, 'f'},
      {"nmeas", required_argument, 0, 'n'},
      {"rsector", required_argument, 0, 'r'},
      {"detail", no_argument, 0, 'd'},
      {"boxes", no_argument, 0, 'b'},
      {"distance", no_argument, 0, 'a'},
      {"olddata", no_argument, 0, 'o'},
      {"writemeas", no_argument, 0, 'w'},
      {"fastmode", no_argument, 0, 'F'},
      {"3d", no_argument, 0, 0},
      {"help", no_argument, 0, 'h'},
      {"version", no_argument, 0, 'v'},
      {0, 0, 0, 0}
    };
			
		/* getopt_long stores the option index here. */
		int option_index = 0;

		c = getopt_long(argc, argv, "s:t:f:n:r:hvdboawF",
		long_options, &option_index);

		/* Detect the end of the options. */
		if (c == -1)
			break;

		switch (c) {
			case 0:
        if (strcmp( "3d", long_options[option_index].name ) == 0)
          opt.do3d = true;
				break;

			case 's':
				opt.Ns = atoi(optarg);
				break;

			case 't':
				opt.Nt = atoi(optarg);
				break;

			case 'f':
				opt.fraction = atof(optarg);
				break;

			case 'n':
				opt.nmeas = atoi(optarg);
				break;

			case 'd':
				opt.detail = true;
				break;
			
			case 'b':
				opt.doboxes = true;
				break;
			
			case 'a':
				opt.dodistance = true;
				break;
			
      case 'o':
				opt.wupperdata = false;
        opt.matrixdim = 3;
				break;
      
      case 'w':
				opt.writemeas = true;
				break;
      
      case 'F':
				opt.fastmode = true;
				break;
			
			case 'r':
				opt.usealternativesectors = true;
				opt.r = atof(optarg);
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
	
	opt.delta = delta0*(1.0 - opt.fraction);

	// Prepare for nmeas measurements
	obs = new Observablestruct[opt.nmeas];

	opt.fevname.resize(opt.nmeas);

	if (opt.nmeas>1) {
		// Read finname file and fill fevname
		ifstream fin;
		fin.open(finname.c_str());
		if(fin.is_open()!=true) {
      cout  << "ERROR: File " << finname <<  " to read configuration filename list could not be opened!" << endl;
      throw 1;
    }

    string strtmp;
    int n=0;
    while (n<opt.nmeas && getline(fin, opt.fevname.at(n))) {
      n++;
		}

		if (n<opt.nmeas) {
			cout << "ERROR: Only found " << n << " names in " << finname << " !" << endl;
			return 1;
		}
	} else {
		opt.fevname[0] = finname;
	}

	opt.leng1=opt.Ns; opt.leng2=opt.Ns; opt.leng3=opt.Ns; opt.leng4=opt.Nt;
	
	opt.Nspace = opt.Ns*opt.Ns*opt.Ns;

	neib.resize(opt.Nspace);
	for (int is=0; is<opt.Nspace; is++)
		neib[is].resize(6);

	fillNeib();

	#ifdef DEBUG
	  cout << "Allocating 3d (Ns^3 * 3) lattice for Polyakov loop evs... " << flush;
	#endif
	pollev.resize(opt.Nspace);
	for (int is=0; is<opt.Nspace; is++)
		pollev[is].resize(opt.matrixdim);
	#ifdef DEBUG
	  cout << "done!" << endl;
	#endif

	if (opt.doboxes) {
		// Boxsize array
		for (int s=0; s<opt.Ns; s++) {
			if (opt.Ns % (s+1) == 0)
				boxsize.push_back(s+1);
		}
		
		for (unsigned int i=0; i<boxsize.size(); i++) {
			boxes.push_back(opt.Ns/boxsize[i]);
		}
	
		#ifdef DEBUG
		for (unsigned int i=0; i<boxsize.size(); i++)
			cout << "Size = " << boxsize[i] << " Nr. of boxes = " << boxes[i] << endl;
		#endif
	}
	#ifdef DEBUG
	cout << "done!" << endl;
	#endif
	
	return 0;
}

#endif // FINDCLUSTER_INIT_HPP
