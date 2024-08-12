#ifndef SVE__H
#define SVE__H
#include "common.h"
#include <fstream>

typedef enum{ tensile, shear, NUM_LT_TYPPE } LOAD_TYPE;

using namespace std;


class interfacePoint
{
public:
	vecDim posXY;
	vecDim normVec;
	double normAng_XY;
	double radius;
	void normalCalculation(const vecDim &center);

};

class interfGeom
{
public:
	interfGeom();
	int numIntPnts;
	vecDim center;
	vector<interfacePoint> intfPnts;
	void compGeomProp();
};
class geometry
{
public:
	int numIntf;
	vector<interfGeom> intfs;
	void initGeometry();
	void compGeomProp();
};

class appliedLoad
{
public:
	symTen stress;//has size
	//symTen strain;//has size
};

class intfSolution
{
public:
	int numPnts;
	vector<symTen> stressPnts;//should declare its size
	void initIntfSolution();

};

class FEsolution
{
public:
	vector<intfSolution> intfSol;//sould declare its size
	void initFEsolution(int numInt);
};
class solution
{
public:
	appliedLoad appLoad;
	FEsolution feSol;
	void readSol(string& fileLoc);
};

class stiffness
{
public:
	vector<vector<double>> realStiff;
	vector<vector<double>> symStiff;
	vector<vector<double>> antiSymStiff;
};

class criticalVal_intf
{
public:
	int pnt;
	double val;
};

class intrefacePostProc
{
public:
	int intfID;
	criticalVal_intf crStrng;
	vector<double> strengthPnts;
	void failCalc(const vector<solution> &sol, const vector<double> &LFsuperpose, const geometry &geom, double TanFac);
	double effStresCal(const vector<double> &traction, const vector<double> &normVec, double tangFrac);
};
class summIntf_intf
{
public:
	int pnt;
	int intrf;
	double strengh;
};
class tetaPostProc
{
public:
	double angle;
	vector<double> LFsuperpose;//superposition factors for auxLoad_XY
	vector<double> auxLoad_XY;//the target trial load in angle but in XY cord
	vector<intrefacePostProc> intfPP;//solution on each interface
	summIntf_intf criticalPntIntf;//worst case in all interfaces
	void compLFsuperpose(const vector<double> &auxload, const vector<solution> &sol);//doing superposition to find load factor of a tial load
	void compCriticalPntIntf(const vector<solution> &sol, const vector<double> &LFsuperpose, const geometry &geom, double tangFrac);
	void compIntfPP(const int &numIntf, const vector<double> &auxload, const vector<solution> &sol, const geometry &geom, double tanFrac);
};
class summPPIntf_teta
{
public:
	double teta;
	int pnt;
	int intf;
	double strength;
};
class summPP_teta
{
public:
	summPPIntf_teta intfaceSumm;
};

class angularProperty
{
public:
	angularProperty();
	double dteta;
	int numTeta;
	vector<tetaPostProc> postProc;//for each teta
	summPP_teta criticalAngle;//summary of all teta
	vector<double> statVar_teta;
	void summPostProc_teta();
	void compStat_teta();
	void compPostProc(const int &numIntf, const  vector<double> &auxloa, const vector<solution> &sol, const geometry &geom, double tanFrac);
	
};
class trialLoad
{
public:
	LOAD_TYPE loadname;
	void initTrialLoad();
	vector<double> auxLoad;
	angularProperty angProp;
	void set0step(double incAngl);
	void writeResults(const int iRVE, const int iSVE, double incAngl, string folderName);
};

class modeMixity
{
public:
	modeMixity();
	double IncrAngl;
	int Nteta;
	vector<double> modeMix_Intf;//at each teta
	vector<double> StatVar_Intf;
	void computModeMixity_intf(const vector<trialLoad> &trialLoads);
	void compStatModMix_Intf();
	void compModeMixity(const vector<trialLoad> &trialLoads);
	void writeResults(const int iRVE, const int iSVE, string folderName);
};

class sve
{
public:
	sve();
	
	int sveId;
	int numSol;
	int numTrialLoads;
	geometry geom;//geometry of SVE
	vector<solution> sol;//FE solutution from Outside
	stiffness stiff;//stiffness of SVE in X-Y cord
	vector<trialLoad> trialLoads;//testing loadesfor desired purpose (find strength at angles)
	modeMixity modemix;
	double bulkModu_sym;
	double bulkMod_real;
	bool stiffCalc(bool sympart,bool antiSympart);
	void set0step(double incAngle);
	void readSVEdata(int Irve, string folderLoc);
	void comptrialLoads(double tanFrac);
	void writeResults(const int irve, double incAngl, string folderName);
	bool compBulkMod();
	void runSVE(int irve, int isve, double tanFrac, double incAngl, string folderLoc);
	void deletData();
};

class criticalModeMixity_Intf
{
public:
	double teta;
	double value;
};























vector<double> readLine_d(ifstream& f1,int num);



#endif










