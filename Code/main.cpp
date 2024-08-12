#include "params.h"
#include "RVE.h"
#include "CONFIG.h"
#include <iostream>
using std::cout;

fstream dbc;

int main()
{

	//----Start:creat a folder for output-----
	char tempC[MAX_SIZE_CHAR];
	string tempS = outPutMainRout;
	strcpy(tempC, tempS.c_str());
	CreateFolder(tempC);//creat the main root for output target 
	//----End:creat a folder for output-----


	config configProp;
	configProp.read();
	vector<rve> RVEs;
	int Nrve = configProp.numRVE;
	RVEs.resize(Nrve);
	for (int irve = 0; irve < Nrve; irve++)
	{
		int Nsve = *(configProp.numSVE + irve);
		double incAngl = *(configProp.incrAngle + irve);
		double tanFrac = *(configProp.tangFac + irve);
		int RVEpart_x = *(configProp.SVE_NumPart_x + irve);
		int RVEpart_y = *(configProp.SVE_NumPart_y + irve);

		
		RVEs[irve].set0step(Nsve, incAngl, configProp.RVEfolderName[irve]);
		
		RVEs[irve].runRVE(irve, tanFrac, Nsve, incAngl, RVEpart_x, RVEpart_y);

	}
	return 0;
}
