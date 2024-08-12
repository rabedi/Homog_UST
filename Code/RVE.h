#ifndef RVE__H
#define RVE__H
#include "SVE.h"


struct parameters
{
	int numSVE, idRVE, numLoad;
	string folderLocation;

};
class statStatSVE
{
public:
	STAT_VAR statName;
	vector<double> stat;
};
class statStatLoad
{
public:
	statStatLoad();
	LOAD_TYPE loadName;
	vector<statStatSVE> statStat_sve;
};
class statStatModMix
{
public:
	statStatModMix();
	vector<statStatSVE> statStat_sve;
};

class rve
{
public:
	rve();
	parameters params;
	int numSve_tot;
	vector<sve> SVEs;
	//vector<statLoad> statloads;
	vector<statStatLoad> statStat_load;
	statStatModMix statStat_ModMix;
	vector<int> nonExistingSVE;
	void set0step(int numSVE, double incrAngle, string FolderLoc);
	void readHMloads(int irve);
	void readHMloads2(int irve, int Nsve);
	void readHMloads3(int irve, int Nsve);
	void compStatStat_sve();
	vector<int> num2IJ(int sveID, int RVEpart_x);
	void writeStat_SVE();
	void writeSDGinput(int RVEpart_x, int RVEpart_y);
	void writeAllData(int RVEpart_x, int RVEpart_y);
	void writeAllSVEatTeta();
	void writeResults(int RVEpart_x, int RVEpart_y);
	void runRVE(int irve, double tanFrac, int Nsve, int incAngl, int RVEpart_x, int RVEpart_y);
	void writeSummAnal_rve();

};




#endif