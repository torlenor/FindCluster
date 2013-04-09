#ifndef FINDCLUSTER_INIT_HPP
#define FINDCLUSTER_INIT_HPP

#include <getopt.h>
#include <cstdio>
#include <string.h>

#define MAJOR_VERSION 0
#define MINOR_VERSION 0
#define REVISION_VERSION 0

char texthelp[]="Usage: 3dclusters.x [OPTION] ... [CLUSTERCONFIGLIST]\n"
		"Visualization of clusters.\n" 
		"Uses cluster configuration files from CLUSTERCONFIGLIST as input.\n"
		"\n"
		"Mandatory arguments to long options are mandatory for short options too.\n"
		"  -s, --Ns SSIZE             spatial lattice extent (default = 8)\n"
		"  -t, --Nt TSIZE             temporal lattice extent (default = 4)\n"
		"  -n, --nmeas NMEAS          number of configurations (default = 1)\n"
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

int parameterInit(int &argc, char *argv[]){
	cout << endl;
	cout << "3dclusters.x" << endl
		<< "Visualization of clusters." << endl

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
			{"nmeas", required_argument, 0, 'n'},
			/* These options set a flag. */
			// {"free", no_argument, 0, 'f'},
			// {"u0", required_argument, 0, 0},
			{"help", no_argument, 0, 'h'},
			{"version", no_argument, 0, 'v'},
			{0, 0, 0, 0}
			};
			
		/* getopt_long stores the option index here. */
		int option_index = 0;

		c = getopt_long (argc, argv, "s:t:n:hv",
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

			case 'n':
				nconfigs = atoi(optarg);
				break;

			case 'v':
				cout << endl << "3dcluster.x version " << MAJOR_VERSION << "." << MINOR_VERSION << "." << REVISION_VERSION << endl;
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

	leng1=Ns; leng2=Ns; leng3=Ns; leng4=Nt;

	cout << "done!" << endl;
	
	return 0;
}

#endif // FINDCLUSTER_INIT_HPP
