#ifndef CONFIH__H
#define CONFIH__H
#include <vector>

using namespace std;
class config
{
public:
	int numRVE;
	double* SVE_NumPart_x;
	double* SVE_NumPart_y;
	int* numSVE;
	double* incrAngle;
	double* contrastRatio;
	double* tangFac;
	double* intfStrength;
	int* microStruc_type;
	vector<string> RVEfolderName;
	void biuldBase();
	void read();
};
#endif CONFIH__H