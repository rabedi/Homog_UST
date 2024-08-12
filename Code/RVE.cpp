#include "RVE.h"
#include "params.h"
#include <iostream> //cout & ios
#include <sstream> //ostringstream
#include "CFEMTypes_Global.h"
rve::rve()
{
	/*
	params.numSVE = Nsve;
	params.numLoad = Nloading;
	SVEs.resize(params.numSVE);
	statStat_load.resize(NtriLoad);
	*/
	params.numSVE = -1;
	params.numLoad = -1;
	
}

void rve::set0step(int numSve, double incAngl, string folderLoc)
{
	params.numSVE = numSve;
	params.numLoad = Nloading;
	params.folderLocation = folderLoc;
	SVEs.resize(params.numSVE);
	for (int i = 0; i < numSve; i++)
	{
		SVEs[i].set0step(incAngl);
	}
	statStat_load.resize(NtriLoad);

}
void rve::readHMloads(int irve)
{
	ostringstream conv;
	conv << irve + 1;
	string fileLoc_HMloads = "RVE" + conv.str() + "//StressBoundaryOutput.txt";
	params.idRVE = irve + 1;
	//reading and assingning all loades to SVEs
	ifstream file;
	file.open(fileLoc_HMloads, ios::in);
	bool f1status = (file.is_open());


	//if (f1status == false)
	//THROW("Code can not find the boundary averaged stress file\n");
	int cont1, cont2, cont3;
	cont3 = 0;

	unsigned int a, b;
	for (int i = 0; i < params.numSVE*params.numLoad; i++)
	{
		{
			file >> a >> b;
			cont1 = (i + 1) / params.numLoad + 1 + cont3;
			if ((i + 1) % params.numLoad == 0)
				cont1--;
			if (a != cont1)
			{
				cont3++;
				cout << "There is a problem in file StressBoundaryOutput regard to SVE number" << cont1 << endl;
			}
			SVEs[b - 1].sveId = a;
			cont2 = (i + 1) % params.numLoad;
			if ((i + 1) % params.numLoad == 0)
				cont2 = params.numLoad;
			if (b != cont2)
				cout << "There is a problem in file StressBoundaryOutput regard to SVE number" << cont1 << "and boundary condition " << cont2 << endl;

			for (int j = 0; j < Ndim*(Ndim + 1) / 2; j++)
				file >> SVEs[a - 1].sol[b - 1].appLoad.stress.comp[j];// loads[b - 1].stress[j];
		}
	}
	file.close();
}

void rve::readHMloads2(int irve, int Nsve)
{
	ostringstream conv;
	conv << irve + 1;
	string fileLoc_HMloads = "RVE" + conv.str() + "//StressBoundaryOutput.txt";
	params.idRVE = irve + 1;
	//reading and assingning all loades to SVEs
	ifstream file;
	file.open(fileLoc_HMloads, ios::in);

	unsigned int a, b;
	string line;
	istringstream tempLine;

	int nSymTen = Ndim*(Ndim + 1) / 2;
	if (file.is_open())
	{
		int isve = 0, ill, tempID = 0, nTotNonExSVE = 0, nLocNonExSVE = 0;
		bool endForc = false;
		while (file.good() && isve < Nsve)
		{
			ill = 0;
			for (int iload = 0; iload < Nloading; iload++)
			{
				ill++;
				getline(file, line);
				tempLine = (istringstream)line;
				if (tempLine.rdbuf()->in_avail() == 0)
				{
					endForc = true;
					break;
				}

				tempLine >> a >> b;// > aa >> bb;
				if (b != ill)
					cout << "Error in reading HM load files\n";
				if (iload == 0)
					tempID = a;
				if (a != tempID)
					cout << "Error in reading HM files\n";
				for (int icomp = 0; icomp < nSymTen; icomp++)
					tempLine >> SVEs[isve].sol[iload].appLoad.stress.comp[icomp];
			}
			if (endForc == true)
				continue;
			SVEs[isve].sveId = a;
			if (isve == 0)
				nLocNonExSVE = a - 1;
			else
				nLocNonExSVE = a - SVEs[isve - 1].sveId - 1;
			if (nLocNonExSVE > 0)
			{
				nTotNonExSVE = nTotNonExSVE + nLocNonExSVE;
				nonExistingSVE.resize(nTotNonExSVE);
				for (int i = 0; i < nLocNonExSVE; i++)
					nonExistingSVE[nTotNonExSVE - i - 1] = a - i - 1;
			}


			isve++;

		}
		params.numSVE = isve;
		SVEs.resize(isve);
	}

	file.close();
}




void rve::readHMloads3(int irve, int Nsve)
{
	ostringstream conv,conv2;
	string fileLoc_HMloads, line,buf;
	conv << irve + 1;
	params.idRVE = irve + 1;
	//reading and assingning all loades to SVEs
	ifstream file;
	int a, b, nSymTen = Ndim*(Ndim + 1) / 2, ill, tempID = 0, SVECount = 0;
	istringstream tempLine;
	for (int isve = 0; isve < Nsve; isve++)
	{
		conv2.str("");
		conv2 << isve + 1;
		//fileLoc_HMloads = "RVE" + conv.str()+"\\SVE" + conv2.str() + ".whl";
		fileLoc_HMloads = params.folderLocation+ "SVE" + conv2.str() + ".whl";
		file.open(fileLoc_HMloads, ios::in);
		if (file.good())
		{
			
			for (int iload = 0; iload < Nloading; iload++)
			{
				getline(file, line);
				tempLine = (istringstream)line;
				for (int icomp = 0; icomp < nSymTen; icomp++)
					tempLine >> SVEs[SVECount].sol[iload].appLoad.stress.comp[icomp] >> buf;
			}
			SVEs[SVECount].sveId = isve+1;
			SVECount++;
		}
		else
		{
			nonExistingSVE.push_back(isve+1);
		}
		file.close();
	}
	params.numSVE = SVECount;
	SVEs.resize(SVECount);
}






void rve::compStatStat_sve()
{
	vector<double> tempVal;
	tempVal.resize(params.numSVE);
	for (int iload = 0; iload < NtriLoad; iload++)
	{
		statStat_load[iload].loadName = (LOAD_TYPE)iload;
		for (int is = 0; is < NUM_STAT_VAR; is++)
		{
			for (int i = 0; i < params.numSVE; i++)
				tempVal[i] = 0.0;
			for (int i = 0; i < params.numSVE; i++)
				tempVal[i] = SVEs[i].trialLoads[iload].angProp.statVar_teta[is];

			statStat_load[iload].statStat_sve[is].statName = (STAT_VAR)is;
			CompStatistics(tempVal, statStat_load[iload].statStat_sve[is].stat);
		}

	}




	for (int is = 0; is < NUM_STAT_VAR; is++)
	{
		for (int i = 0; i < params.numSVE; i++)
			tempVal[i] = 0.0;
		for (int i = 0; i < params.numSVE; i++)
			tempVal[i] = SVEs[i].modemix.StatVar_Intf[is];

		statStat_ModMix.statStat_sve[is].statName = (STAT_VAR)is;
		CompStatistics(tempVal, statStat_ModMix.statStat_sve[is].stat);
	}


}

void rve::runRVE(int irve, double tanFrac, int Nsve, int incAngl, int RVEpart_x, int RVEpart_y)
{
	readHMloads3(irve, Nsve);
	for (int isve = 0; isve < params.numSVE; isve++)
	{
		cout << "RVE\t" << irve << "\t" << "SVE\t" << isve << "\n";
		SVEs[isve].runSVE(irve, isve, tanFrac, incAngl, params.folderLocation);//computation over each SVE
	}
	compStatStat_sve();//computation over all SVEs
	writeResults( RVEpart_x, RVEpart_y);
}










void rve::writeResults(int RVEpart_x, int RVEpart_y)
{
	writeStat_SVE();
	writeSDGinput( RVEpart_x, RVEpart_y);
	writeAllData(RVEpart_x, RVEpart_y);
	writeAllSVEatTeta();
	writeSummAnal_rve();

}



statStatLoad::statStatLoad()
{
	statStat_sve.resize(NUM_STAT_VAR);

}
statStatModMix::statStatModMix()
{
	statStat_sve.resize(NUM_STAT_VAR);
}
vector<int> rve::num2IJ(int sveID, int RVEpart_x)
{
	//sveID shoul be in natural numbers [1-inf]
	//result will be in natural numbers [1-inf]
	vector<int> IJ(Ndim);
	IJ[1] = (int)sveID / RVEpart_x + 1;
	if (remainder(sveID, RVEpart_x) == 0.0) IJ[1] = IJ[1] - 1;
	IJ[0] = sveID - RVEpart_x * (IJ[1] - 1);
	return IJ;
}
void rve::writeSDGinput(int RVEpart_x, int RVEpart_y)
{
	ofstream file;
	ostringstream conv;
	string id;
	double scalFac = 1.0e9;//scale factor for Gpa
	conv << params.idRVE;
	//id = "OUTPUT//RVE" + conv.str() + "(SDGinput_SVE).OSDG";
	id = "OUTPUT//" + params.folderLocation + "(SDGinput_SVE).OSDG";
	file.open(id);
	file << "#NumGridPoint_X\t NumGridPoint_Y\n";
	file << RVEpart_x << "\t" << RVEpart_y << "\n";
	file << "#<field Type>\t <tensor Type>\t <tensor Order>\t <List of Tensor Dimensions for respective order>\t <number of independent components>\n";
	file << 1 << "\n";
	file << "Void\t tenFULL\t" << 1 << "\t" << 10 << "\t" << 10 << "\n";
	file << "#Xind\t Yind\t Nfstrength_mean\t Sfstrength_mean\t c00\t c10\t c11\t c20\t c21\t c22\t bulkMod\t betta_mean\n";
	//file << "#Xind\tYind\tNfstrength_mean\tSfstrength_mean\tc00\tc10\tc11\tc20\tc21\tc22\tbulkMod\tbetta_mean\n";
	vector<int> IJ;
	int NsymTen = Ndim*(Ndim + 1) / 2;
	for (int isve = 0; isve < params.numSVE; isve++)
	{
		IJ = num2IJ(SVEs[isve].sveId, RVEpart_x);
		for (int i = 0; i < Ndim; i++)
			file << IJ[i] << "\t";

		for (int i = 0; i < NtriLoad; i++)
			file << SVEs[isve].trialLoads[i].angProp.statVar_teta[mean] << "\t";

		for (int j = 0; j < NsymTen; j++)
			for (int i = 0; i <= j; i++)
				file << SVEs[isve].stiff.symStiff[i][j] / scalFac << "\t";

		file << SVEs[isve].bulkModu_sym / scalFac << "\t";
		file << SVEs[isve].modemix.StatVar_Intf[mean] << "\t";
		file << "\n";
	}


}

void rve::writeAllData(int RVEpart_x, int RVEpart_y)
{
	ofstream file;
	ostringstream conv;
	string id;
	double scalFac = 1.0e9;//scale factor for Gpa
	conv << params.idRVE;
	//id = "OUTPUT//RVE" + conv.str() + "(SDGinput_SVE).OSDG";
	id = "OUTPUT//" + params.folderLocation + "(AllFld_SVE).OAllFld";
	file.open(id);
	file << "#NumGridPoint_X\t NumGridPoint_Y\n";
	file << RVEpart_x << "\t" << RVEpart_y << "\n";
	file << "#<field Type>\t <tensor Type>\t <tensor Order>\t <List of Tensor Dimensions for respective order>\t <number of independent components>\n";
	file << 1 << "\n";
	file << "Void\t tenFULL\t" << 1 << "\t" << 13 << "\t" << 13 << "\n";
	file << "# Nfstrength normalized by mean value(<computed value>) GPa \n";
	file << 0<<endl;
	file << "# # Sfstrength normalized by mean value(<computed value>) GPa \n";
	file << 0<<endl;
	file << "#num of rows in data\n";
	file << params.numSVE << endl;
	file << "#Xind\t Yind\t Gind\t Nfstrength_mean\t Sfstrength_mean\t betta_mean\t c00\t c10\t c11\t c20\t c21\t c22\t bulkModulus\n";
	//file << "#Xind\tYind\tNfstrength_mean\tSfstrength_mean\tc00\tc10\tc11\tc20\tc21\tc22\tbulkMod\tbetta_mean\n";
	vector<int> IJ;
	int NsymTen = Ndim*(Ndim + 1) / 2;
	for (int isve = 0; isve < params.numSVE; isve++)
	{
		IJ = num2IJ(SVEs[isve].sveId, RVEpart_x);
		for (int i = 0; i < Ndim; i++)
			file << IJ[i] << "\t";

		file << SVEs[isve].sveId << "\t";

		for (int i = 0; i < NtriLoad; i++)
			file << SVEs[isve].trialLoads[i].angProp.statVar_teta[mean] << "\t";

		file << SVEs[isve].modemix.StatVar_Intf[mean] << "\t";

		for (int j = 0; j < NsymTen; j++)
			for (int i = 0; i <= j; i++)
				file << SVEs[isve].stiff.symStiff[i][j] / scalFac << "\t";

		file << SVEs[isve].bulkModu_sym / scalFac << "\t";
		file << "\n";
	}


}


void rve::writeStat_SVE()
{
	ofstream file;
	ostringstream conv;
	string id;
	conv << params.idRVE;

	for (int il = 0; il < NtriLoad; il++)
	{
		if (il == tensile)
		{
			//id = "OUTPUT//RVE" + conv.str() + "(Tstat_SVE).OTstat";
			id = "OUTPUT//" + params.folderLocation + "(Tstat_SVE).OTstat";
		}
			
		else if (il == shear)
		{
			//id = "OUTPUT//RVE" + conv.str() + "(Sstat_SVE).OSstat";
			id = "OUTPUT//" + params.folderLocation+ "(Sstat_SVE).OSstat";
		}
			
		file.open(id);
		file << "RVE\t Nsve\t NstatVar\t Load\n";
		file << params.idRVE << "\t" << params.numSVE << "\t" << NUM_STAT_VAR << "\t" << SVEs[0].trialLoads[il].loadname + 1 << "\n";
		file << "SVEid\t MIN\t MAX\t MEAN\t StdDiv\n";
		for (int isve = 0; isve < params.numSVE; isve++)
		{
			file << SVEs[isve].sveId << "\t";
			for (int is = 0; is < NUM_STAT_VAR; is++)
				file << SVEs[isve].trialLoads[il].angProp.statVar_teta[is] << "\t";
			file << "\n";
		}
		file << "Stat of Stat\n";
		for (int is = 0; is < NUM_STAT_VAR; is++)
		{
			for (int is2 = 0; is2 < NUM_STAT_VAR; is2++)
				file << statStat_load[il].statStat_sve[is].stat[is2] << "\t";
			file << "\n";
		}
		file.close();
	}



	//id = "OUTPUT//RVE" + conv.str() + "(MMstat_SVE).OMMstat";
	id = "OUTPUT//" + params.folderLocation + "(MMstat_SVE).OMMstat";
	file.open(id);
	file << "RVE\t Nsve\t NstatVar ModeMixity\n";
	file << params.idRVE << "\t" << params.numSVE << "\t" << NUM_STAT_VAR << "\t" << "\n";
	file << "SVEid\t MIN\t MAX\t MEAN\t StdDiv\n";
	for (int isve = 0; isve < params.numSVE; isve++)
	{
		file << SVEs[isve].sveId << "\t";
		for (int is = 0; is < NUM_STAT_VAR; is++)
			file << SVEs[isve].modemix.StatVar_Intf[is] << "\t";
		file << "\n";
	}
	file << "Stat of Stat\n";
	for (int is = 0; is < NUM_STAT_VAR; is++)
	{
		for (int is2 = 0; is2 < NUM_STAT_VAR; is2++)
			file << statStat_ModMix.statStat_sve[is].stat[is2] << "\t";
		file << "\n";
	}
	file.close();


}

void rve::writeAllSVEatTeta()
{
	ofstream file;
	ostringstream conv, conv2, conv3;
	string id;
	conv << params.idRVE;
	
	int Nteta = SVEs[0].trialLoads[0].angProp.numTeta;
	for (int il = 0; il < NtriLoad; il++)
	{
		conv2.str("");
		conv2 << il + 1;
		for (int it = 0; it < Nteta; it++)
		{
			conv3.str("");
			conv3 << SVEs[0].trialLoads[il].angProp.postProc[it].angle;
			if (il == 0)
			{
				//id = "OUTPUT//RVE" + conv.str() + "Teta" + conv3.str() + "(Tstrn_SVE).OSveSTeta";
				id = "OUTPUT//" + params.folderLocation + "Teta" + conv3.str() + "(Tstrn_SVE).OSveSTeta";
			}
				
			else
			{
				//id = "OUTPUT//RVE" + conv.str() + "Teta" + conv3.str() + "(Sstrn_SVE).OSveSTeta";
				id = "OUTPUT//" + params.folderLocation + "Teta" + conv3.str() + "(Sstrn_SVE).OSveSTeta";
			}
				

			
			file.open(id);
			file << "data type\n";
			if (il == 0)
				file << "Tensile\n";
			else
				file << "Shear\n";
			
			file << "RVE\t Nsve\n";
			file << params.idRVE << "\t" <<  params.numSVE << "\n";
				file << "Teta\n";
				file << SVEs[0].trialLoads[il].angProp.postProc[it].angle<<"\n";
				file << "SVE\t Strength\n";
			for (int is = 0; is < params.numSVE; is++)
			{
				file << SVEs[is].sveId << "\t" << SVEs[is].trialLoads[il].angProp.postProc[it].criticalPntIntf.strengh<<"\n";
			}
			file.close();
		}

	}




	for (int it = 0; it < Nteta; it++)
	{
		conv3.str("");
		conv3 << SVEs[0].trialLoads[0].angProp.postProc[it].angle;
		//id = "OUTPUT//RVE" + conv.str() + "Teta" + conv3.str() + "(ModeMixity_SVE).OSveSTeta";
		id = "OUTPUT//" + params.folderLocation + "Teta" + conv3.str() + "(ModeMixity_SVE).OSveSTeta";
		file.open(id);
		file << "data type\n";
		file << "ModeMixity\n";
		file << "RVE\t Nsve\n";
		file << params.idRVE << "\t" << params.numSVE << "\n";
		file << "Teta\n";
		file << SVEs[0].trialLoads[0].angProp.postProc[it].angle << "\n";
		file << "SVE\t Strength\n";
		for (int is = 0; is < params.numSVE; is++)
		{
			file << SVEs[is].sveId << "\t" << SVEs[is].modemix.modeMix_Intf[it] << "\n";
		}
		file.close();
	}

}

void rve::writeSummAnal_rve()
{
	ofstream file;
	ostringstream conv, conv2, conv3;
	string id;
	conv << params.idRVE;
	//id = "OUTPUT//RVE" + conv.str() + "SummaryAnalysis.OSumAnlys";
	id = "OUTPUT//" + params.folderLocation + "SummaryAnalysis.OSumAnlys";
	file.open(id);
	file << "number of active SVEs\t" << params.numSVE<<"\n";
	int tNumIncl = 0;
	for (int i = 0; i < params.numSVE; i++)
	{
		if (SVEs[i].geom.numIntf <= 0)
			THROW("Error in finding correct active SVEs");
		tNumIncl += SVEs[i].geom.numIntf;
	}
	file << "Total number of inclusions in the RVE\t" << tNumIncl << "\n";
	file << "Average number of inclusions in active SVES\t" << tNumIncl / params.numSVE << "\n";

		

}