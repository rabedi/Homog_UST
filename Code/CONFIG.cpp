#include "CONFIG.h"
#include <fstream>
#include "params.h"
#include "CFEMTypes_Global.h"
#include "common.h"
using namespace std;

void config::read()
{
	ifstream file;
	file.open(configFile, ios::in);
	if (!file.is_open())
		THROW("Config File does not exist\n");
	
	string bufS;
	getline(file, bufS);
	file >> bufS;
	string fileRout = bufS;

	//----Start:creat a corresponding folder for output-----
	char tempC[MAX_SIZE_CHAR];
	string tempS = outPutMainRout + fileRout;
	//strcpy(tempC, tempS.c_str());
	//strcat(tempC, fileRout.c_str());
	strcpy(tempC, tempS.c_str());
	CreateFolder(tempC);//creat the main root for output target 
	//----End:creat a correspoding folder for output-----

	
	getline(file, bufS);
	getline(file, bufS);
	file >> numRVE;
	biuldBase();
	getline(file, bufS);
	getline(file, bufS);
	for (int iRVE = 0; iRVE < numRVE;iRVE++)
	{
		file >> bufS;
		RVEfolderName[iRVE] = (fileRout+bufS)+"//";

		//----Start:creat a corresponding folder for output-----
		char tempC2[MAX_SIZE_CHAR] ;
		string tempS2 = tempS+bufS;
		strcpy(tempC2, tempS2.c_str());
		CreateFolder(tempC2);//creat the main root for output target 
		//----End:creat a correspoding folder for output-----


		file >> *(microStruc_type+iRVE);
		file >> *(contrastRatio + iRVE);
		file >> *(SVE_NumPart_x + iRVE);
		file >> *(SVE_NumPart_y + iRVE);
		file >> *(incrAngle + iRVE);
		file >> *(tangFac + iRVE);
		file >> *(intfStrength + iRVE);
		*(numSVE + iRVE) = (*(SVE_NumPart_x + iRVE))* (*(SVE_NumPart_y + iRVE));
	}
	file.close();
}

void config::biuldBase()
{
	SVE_NumPart_x = new double[numRVE];
	SVE_NumPart_y = new double[numRVE];
	numSVE = new int[numRVE];
	incrAngle = new double[numRVE];
	contrastRatio = new double[numRVE];
	tangFac = new double[numRVE];
	intfStrength = new double[numRVE];
	microStruc_type = new int[numRVE];
	RVEfolderName.resize(numRVE);
}