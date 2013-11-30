#ifndef FINDCLUSTER_H
#define FINDCLUSTER_H

// Vectors to store cluster information
struct Clusterstruct{
	vector<int> isinsector;  // lclusterdata.isinsector[is] stores the sector of the lattice point is
	vector<int> clustersector; // lclusterdata.clustersector[c] stores the sector of the cluster c
	vector<int> isincluster;  // lclusterdata.isincluster[is] stores the cluster of the lattice point is
	vector<vector<int> > clustermembers; // lclusterdata.clustermembers[c][i] stores the members/lattice points of the cluster c (0 < i < N_c)

	vector<int> percolatingclusters; // Percolating clusters
	vector<vector<int> > percolatingdirections;
	vector<int> clusterispercolating;
	
	vector<int> sortedcluster;
	vector<int> sortedrealcluster;
	vector<int> isinsortedcluster;

	vector<vector<int> > clusterisperiodic;

	int nsectm1, nsect0, nsectp1;

	vector<complex<double> > poll;
};

struct Observablestruct{
	int maxclustersize;
	int maxclusterid;
	int maxclustersector;

  double largestnonpercclustersize;
	double largestnonpercclusterid;
	double largetsnonpercclustersector;
	
	double avgclustersize;
	double avgclustersizeF;

	double rootmeansquaredistanceR;
	
	double avgclustersizenp;
	double avgclustersizeFnp;

	double cut;

	vector<vector<int> > numberofboxes;
	
	vector<vector<double> > centerofmass;
	vector<double> clusterradius;
	double largestclusterradius;
	double largestnpclusterradius;
	double avgclusterradius;
	double avgnpclusterradius;

	double percc;

	int largestclusterid;

	double area;
	double arealargestnonperccluster;
	double areaavgnonperccluster;

	double poll;
  
  vector<double> meanfreepath;
  double largestclustermeanfreepath;
  double largestnpclustermeanfreepath;
  double avgclustermeanfreepath;
  double avgnpclustermeanfreepath;
  
  vector<double> meanfreepathnew;
  double largestclustermeanfreepathnew;
  double largestnpclustermeanfreepathnew;
  
};

struct Resultstruct{
	double maxclustersize, maxclustersizeerr;
	double maxnonpercclustersize, maxnonpercclustersizeerr;
	double avgclustersize, avgclustersizeerr;
	double avgclusersizeFortunato, avgclusersizeFortunatoerr;
	double avgclustersizenp, avgclustersizenperr;
	double avgclustersizeFnp, avgclustersizeFnperr;
	double avgrootmeansquaredistance, avgrootmeansquaredistanceerr;
	double cut, cuterr;
	double totalperimeter, totalperimetererr;
	double largestnonpercperimeter, largestnonpercperimetererr;
	double avgnonpercperimeter, avgnonpercperimetererr;
	double largestclusterradius, largestclusterradiuserr;
	double largestnpclusterradius, largestnpclusterradiuserr;
	double avgclusterradius, avgclusterradiuserr;
	double avgnpclusterradius, avgnpclusterradiuserr;
	double polyakovloopaftercut, polyakovloopaftercuterr;
	double avgperccluster, avgpercclustererr;
	vector<double> largestclusterboxcount, largestclusterboxcounterr;
	vector<double> largestnonpercboxcount, largestnonpercboxcounterr;
  
  double largestclustermeanfreepath, largestclustermeanfreepatherr;
  double largestnpclustermeanfreepath, largestnpclustermeanfreepatherr;
  double avgclustermeanfreepath, avgclustermeanfreepatherr;
  double avgnpclustermeanfreepath, avgnpclustermeanfreepatherr;
  
  double largestclustermeanfreepathnew, largestclustermeanfreepathnewerr;
  double largestnpclustermeanfreepathnew, largestnpclustermeanfreepathnewerr;
};

// Stuff to find and categorize sectors/clusters for one configuration
void fillSectors(Clusterstruct &lclusterdata, double delta);
void fillSectorsAlt(Clusterstruct &lclusterdata, double r);
void findClusters(Clusterstruct &lclusterdata);
void checkClusters(Clusterstruct &lclusterdata);
void findPercolatingCluster(Clusterstruct &lclusterdata);

void writeConfigResultsstdout(Observablestruct &lobs, Clusterstruct &lclusterdata);
void writeClusterList(Clusterstruct &lclusterdata);

void calcExp();

void sortClusterSize(Clusterstruct &lclusterdata);
void cluster3doutput(Clusterstruct &clusterdata, string f3dname);

int latmap(int i1, int i2, int i3);
void getCoords(int is, int &i1, int &i2, int &i3);

void freeMem(Clusterstruct &lclusterdata);

void fillNeib();

#endif // FINDCLUSTER_H
