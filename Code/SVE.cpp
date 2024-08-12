#include "SVE.h"
#include "params.h"
#include <sstream>
#include <iostream>
#include "CFEMTypes_Global.h"
#include "common.h"

// Bahador's code is designed to compute strengths(theta) (stress loading along theta direction)
// Turn on the flag below (0 -> 1) to compute material failure initiation stretches (strain-based resistance)
#define COMPUTE_STRETCH	1

sve::sve()
{
	/*
	numSol = Nloading;
	numTrialLoads = NtriLoad;
	int Nsym = Ndim*(Ndim + 1) / 2;
	sol.resize(numSol);
	sveId = -1;
	trialLoads.resize(numTrialLoads);
	for (int i = 0; i < numSol; i++)
	{
		trialLoads[i].set0step(incAngl);
	}
	*/
}

void geometry::initGeometry()
{
	intfs.resize(numIntf);
}
interfGeom::interfGeom()
{
	numIntPnts = -1;
	//intfPnts.resize(numIntPnts);
}
void FEsolution::initFEsolution(int numInt)
{
	intfSol.resize(numInt);
}

void trialLoad::initTrialLoad()
{
	int Nsym = Ndim*(Ndim + 1) / 2;
	auxLoad.resize(Nsym);
	for (int i = 0; i < Nsym; i++)
		auxLoad[i] = 0.0;
	if (loadname == tensile)
	{
		auxLoad[1] = 1.0;
		//auxLoad[0] = 1.0;
	}
		
	else if (loadname == shear)
		auxLoad[2] = 1.0;
	else
		cout << "The Load Type in trial load is not implemented\n";
}

angularProperty::angularProperty()
{
	//dteta = IncrAngl;
	//numTeta = 180 / dteta + 1;
	//postProc.resize(numTeta);
	dteta = -1;
	numTeta =-1;
	
}
void tetaPostProc::compLFsuperpose(const vector<double> &auxload, const vector<solution> &sol)
{
	int Nsym = Ndim*(Ndim + 1) / 2;
	LFsuperpose.resize(Nsym);
	for (int i = 0; i < Nsym; i++)
		LFsuperpose[i] = 0.0;
	vector<double> StresVec_XY;
	CompStressRotation(auxload, -angle*PI/180, auxLoad_XY);//angle should be in radian!!!

#if DB_CALCULATIONS
	dbc << "auxLoad_XY\t" << auxLoad_XY[0] << '\t' << auxLoad_XY[1] << '\t' << auxLoad_XY[2] << '\n';
#endif

	LFsuperpose.resize(Nsym);
	for (int i = 0; i < Nsym; i++)
		LFsuperpose[i] = 0.0;

	if (Ndim == 2)
	{
		double a11, a12, a13, a21, a22, a23, a31, a32, a33;
		double b11, b12, b13, b21, b22, b23, b31, b32, b33;//inverse components
		double det;
		a11 = sol[0].appLoad.stress.comp[0]; a21 = sol[0].appLoad.stress.comp[1]; a31 = sol[0].appLoad.stress.comp[2];
		a12 = sol[1].appLoad.stress.comp[0]; a22 = sol[1].appLoad.stress.comp[1]; a32 = sol[1].appLoad.stress.comp[2];
		a13 = sol[2].appLoad.stress.comp[0]; a23 = sol[2].appLoad.stress.comp[1]; a33 = sol[2].appLoad.stress.comp[2];

#if DB_CALCULATIONS
		dbc << "LC0Ave-S11\t" << a11 << "\tS22\t" << a21 << "\tS12\t" << a31 << '\n';
		dbc << "LC1Ave-S11\t" << a12 << "\tS22\t" << a22 << "\tS12\t" << a32 << '\n';
		dbc << "LC2Ave-S11\t" << a13 << "\tS22\t" << a23 << "\tS12\t" << a33 << '\n';
#endif

		//2019/09/26 -> computing stretches
#if COMPUTE_STRETCH
		double e = 1e-4;
		a11 = e;	a21 =  e;	a31 = 0.0;
		a12 = e;	a22 = -e;	a32 = 0.0;
		a13 = 0.0;	a23 = 0.0;	a33 = e;
#endif

#if DB_CALCULATIONS
#if COMPUTE_STRETCH
		dbc << "stretch\n";
#else
		dbc << "strength\n";
#endif
		dbc << "KMat-Column1-comp1\t" << a11 << "\tcomp2\t" << a21 << "\tcomp3\t" << a31 << '\n';
		dbc << "KMat-Column2-comp1\t" << a12 << "\tcomp2\t" << a22 << "\tcomp3\t" << a32 << '\n';
		dbc << "KMat-Column3-comp1\t" << a13 << "\tcomp2\t" << a23 << "\tcomp3\t" << a33 << '\n';
#endif

		det = a11*a22*a33 + a21*a32*a13 + a31*a12*a23 - a11*a32*a23 - a31*a22*a13 - a21*a12*a33;
		if (det == 0.0)
			cout << "WARNING :: SINGULAR INVERSION";
		b11 = (1.0 / det)*(a22*a33 - a23*a32); b12 = (1.0 / det)*(a13*a32 - a12*a33); b13 = (1.0 / det)*(a12*a23 - a13*a22);
		b21 = (1.0 / det)*(a23*a31 - a21*a33); b22 = (1.0 / det)*(a11*a33 - a13*a31); b23 = (1.0 / det)*(a13*a21 - a11*a23);
		b31 = (1.0 / det)*(a21*a32 - a22*a31); b32 = (1.0 / det)*(a12*a31 - a11*a32); b33 = (1.0 / det)*(a11*a22 - a12*a21);
		LFsuperpose[0] = b11*auxLoad_XY[0] + b12*auxLoad_XY[1] + b13*auxLoad_XY[2];
		LFsuperpose[1] = b21*auxLoad_XY[0] + b22*auxLoad_XY[1] + b23*auxLoad_XY[2];
		LFsuperpose[2] = b31*auxLoad_XY[0] + b32*auxLoad_XY[1] + b33*auxLoad_XY[2];

#if DB_CALCULATIONS
		dbc << "LoadFactors-comp1\t" << LFsuperpose[0] << "\tcomp2\t" << LFsuperpose[1] << "\tcomp3\t" << LFsuperpose[2] << '\n';
#endif
	}

}
void sve::readSVEdata(int Irve,string folderLoc)
{
	static double close2ZeroVal = 1e-40;
	ostringstream convert, convert1, convert2;
	ifstream f1;
	bool f1status;
	convert << (Irve + 1); convert1 << sveId ;
	string Id, buf, Id0;
	int iSVE, Iload, Nincl, Iincl, aa;
	double bufd;
	//Id0 = "RVE" + convert.str() + "\\SVE" + convert1.str();
	Id0 = folderLoc+ "SVE" + convert1.str();
	int NsymTen = Ndim*(Ndim + 1) / 2;
	bool oldFileFormat = true;
	for (int iload = 0; iload < Nloading; iload++)//loop over Type of loading
	{
		convert2.str("");
		convert2 << iload + 1;
		Id = Id0 + "BC" + convert2.str() + ".PTH";
		f1.open(Id, ios::in);
		if (f1.is_open())
		{
			iSVE = -1; Iload = -1; Nincl = -1;
			f1 >> buf;
			if (buf == "X")
			{
				oldFileFormat = false;
				if (iload == 0)
				{
					Nincl = 0;
					fstream f2;
					f2.open(Id.c_str(), ios::in);
					string buf2 = "";
					while (!f2.eof())
					{
						f2 >> buf2;
						if (buf2 == "Coord1")
							++Nincl;
					}
					f2.close();
				}
			}
			if (oldFileFormat)
				f1 >> iSVE >> buf >> Iload >> buf >> Nincl;//Isve,Iload,Nincl
			else
			{
				Iload = iload;
				if (iload != 0)
					Nincl = geom.numIntf;
			}

			if (iload == 0)
			{
				geom.numIntf = Nincl;
				geom.initGeometry();
			}

			sol[iload].feSol.initFEsolution(Nincl);
			if (oldFileFormat)
				f1 >> buf;//begining of the line (Incl)

			for (int iintf = 0; iintf < Nincl; iintf++)//loop over interfaces
			{
				sol[iload].feSol.intfSol[iintf].initIntfSolution();

				if ((iload == 0) || (!oldFileFormat))
				{
					//read everything
					if (oldFileFormat)
					{
						f1 >> Iincl >> buf;//Iincl
						//1012
						if (Iincl != (iintf + 1))
						{
							cout << "Iincl\t" << Iincl << '\n';
							cout << "iintf\t" << iintf << '\n';
							THROW("(Iincl != (iintf + 1)), this may not be a bug, is an inclusion missing?!\n");
						}
					}
					else
					{
						Iincl = iintf + 1;
						f1 >> buf;//Y
					}
					for (int idim = 0; idim < Ndim; idim++)//read center point
					{
						f1 >> geom.intfs[iintf].center.comp[idim];
						if (oldFileFormat)
							if (idim != Ndim - 1)
								f1 >> buf;
					}
					if (!oldFileFormat)
						for (int k = 0; k < 6; ++k) f1 >> buf;

					f1 >> buf;
					int ip = -1;
					//1012 Bug - buf[0] can be - in this case it won't work!
					while (isDigit(buf) && buf != "\0" && !f1.eof())
					{
						if (!oldFileFormat)
							f1 >> buf;
						bool valid = true;
						ip++;
						geom.intfs[iintf].intfPnts.resize(ip + 1);
						sol[iload].feSol.intfSol[iintf].stressPnts.resize(ip + 1);
						for (int idim = 0; idim < Ndim; idim++)
						{
							if (buf != "NoValue")
								geom.intfs[iintf].intfPnts[ip].posXY.comp[idim] = stod(buf);
							else
							{
								geom.intfs[iintf].intfPnts[ip].posXY.comp[idim] = close2ZeroVal;
								valid = false;
							}
							if (oldFileFormat)
								f1 >> buf;//,
							f1 >> buf;//num
						}
						for (int icomp = 0; icomp < NsymTen; icomp++)
						{
							if (buf != "NoValue")
								sol[iload].feSol.intfSol[iintf].stressPnts[ip].comp[icomp]=stod(buf);
							else
							{
								sol[iload].feSol.intfSol[iintf].stressPnts[ip].comp[icomp] = close2ZeroVal;
								valid = false;
							}
							if (icomp != NsymTen - 1)
							{
								if (oldFileFormat)
									f1 >> buf;//,
								f1 >> buf;//num
							}
						}
						if (!valid)
						{
							ip--;
							geom.intfs[iintf].intfPnts.resize(ip + 1);
							sol[iload].feSol.intfSol[iintf].stressPnts.resize(ip + 1);
						}
						f1 >> buf;//begining of a new lin
					}
					geom.intfs[iintf].numIntPnts = ip+1;
					//1012
//					if (geom.intfs[iintf].numIntPnts != 51)
//						THROW("geom.intfs[iintf].numIntPnts != 51, this may not be a bug, but in current data there are 51 points. Just make sure if this error happens the actual number of points in the file is not 0 or another off number and then deactivate this message\n");
					sol[iload].feSol.intfSol[iintf].numPnts = ip + 1;
					aa=2;
				}
				else
				{
					f1 >> Iincl >> buf;//Iincl
					for (int idim = 0; idim < Ndim; idim++)
					{
						f1 >> bufd;
						if (idim != Ndim - 1)
							f1 >> buf;
					}

					f1 >> buf;
					int ip = -1;
					//1012
//					while (isdigit(buf[0]) && buf != "\0" && !f1.eof())
					while (isDigit(buf) && buf != "\0" && !f1.eof())
					{
						ip++;
						sol[iload].feSol.intfSol[iintf].stressPnts.resize(ip + 1);
						for (int idim = 0; idim < Ndim; idim++)
						{
							//geom.intfs[iintf].intfPnts[ip].posXY.comp[idim] = stod(buf);
							f1 >> buf;//,
							f1 >> buf;//num
						}
						for (int icomp = 0; icomp < NsymTen; icomp++)
						{
							sol[iload].feSol.intfSol[iintf].stressPnts[ip].comp[icomp] = stod(buf);
							if (icomp != NsymTen - 1)
							{
								f1 >> buf;//,
								f1 >> buf;//num
							}
						}
						f1 >> buf;//begining of a new lin

					}
					sol[iload].feSol.intfSol[iintf].numPnts = ip + 1;
					//1012
/*
					if (sol[iload].feSol.intfSol[iintf].numPnts != 51)
					{
						cout << "sol[" << iload << "].feSol.intfSol[" << iintf << "].numPnts\t" << sol[iload].feSol.intfSol[iintf].numPnts << '\n';
						THROW("sol[iload].feSol.intfSol[iintf].numPnts != 51, this may not be a bug, but in current data there are 51 points. Just make sure if this error happens the actual number of points in the file is not 0 or another off number and then deactivate this message\n");
					}
*/
				}

			}


		}
		else
			cout << "Error: the following file does not esixt->" << Id << "\n";
		f1.close();
	}

}

bool sve::stiffCalc(bool sympart, bool antiSympart)
{
	int n = Ndim*(Ndim + 1) / 2;
	stiff.realStiff.resize(n);
	for (int i = 0; i < n; i++)
		stiff.realStiff[i].resize(n);
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			stiff.realStiff[i][j] = 0.0;
	vector<vector<double>> stressMat, stainMat, invStrn;
	stressMat = stiff.realStiff;
	stainMat = stiff.realStiff;
	invStrn = stiff.realStiff;
	for (int il = 0; il < Nloading; il++)
		for (int ic = 0; ic < n; ic++)
			stressMat[ic][il] = sol[il].appLoad.stress.comp[ic];

	double eps0 = 0.0001;
	stainMat[0][0] = eps0; stainMat[0][1] = eps0; stainMat[0][2] = 0;
	stainMat[1][0] = eps0; stainMat[1][1] = -eps0; stainMat[1][2] = 0;
	stainMat[2][0] = 0; stainMat[2][1] = 0; stainMat[2][2] = 2.0 * eps0;



	bool invR = inversMat(stainMat, invStrn);
	if (invR == false)
	{
		cout << "ERROR IN COMPUTING INVERS_IN STIFFNESS";
		return false;
	}
	bool I = matrixMultiplication(stressMat, invStrn, stiff.realStiff);
	if (invR == false)
	{
		cout << "ERROR IN COMPUTING INVERS_MAT MUKTi";
		return false;
	}
	if (sympart == true)
	{
		stiff.symStiff.resize(n);
		for (int i = 0; i < n; i++)
			stiff.symStiff[i].resize(n);
		for (int i = 0; i < n; i++)
			for (int j = 0; j < n; j++)
				stiff.symStiff[i][j] = 0.0;

		for (int i = 0; i < n; i++)
			for (int j = 0; j < n; j++)
				stiff.symStiff[i][j] = 0.5*(stiff.realStiff[i][j] + stiff.realStiff[j][i]);
	}
	if (antiSympart == true)
	{
		stiff.antiSymStiff.resize(n);
		for (int i = 0; i < n; i++)
			stiff.antiSymStiff[i].resize(n);
		for (int i = 0; i < n; i++)
			for (int j = 0; j < n; j++)
				stiff.antiSymStiff[i][j] = 0.0;
		for (int i = 0; i < n; i++)
			for (int j = 0; j < n; j++)
				stiff.antiSymStiff[i][j] = stiff.realStiff[i][j] - 0.5*(stiff.realStiff[i][j] + stiff.realStiff[j][i]);
	}
}

void intfSolution::initIntfSolution()
{
	//numPnts = NpintIncl;
	numPnts = -1;
	//stressPnts.resize(numPnts);
}
vector<double> readLine_d(ifstream& f1, int num)
{
	vector<double> data;
	data.resize(num);
	for (int i = 0; i < num; i++)
		f1 >> data[i];
	return data;
}
void interfacePoint::normalCalculation(const vecDim &center)
{
	double dy = posXY.comp[1] - center.comp[1];
	double dx = posXY.comp[0] - center.comp[0];
	radius = sqrt(dy*dy + dx*dx);
	normVec.comp[0] = dx / radius;
	normVec.comp[1] = dy / radius;
	normAng_XY = atan2(dy, dx);
	normAng_XY = normAng_XY*180.0 / PI;
}



void intrefacePostProc::failCalc(const vector<solution> &sol, const vector<double> &LFsuperpose, const geometry &geom, double tangFrac)
{
	vector<double> stressSupPos;
	int NsymTen = (Ndim + 1)*Ndim / 2;
	stressSupPos.resize(NsymTen);
	double crtNormStrength = 1e200;
	int crtPnt = -1;
	int numPnt = sol[0].feSol.intfSol[intfID].numPnts;
	strengthPnts.resize(numPnt);
	for (int iPnt = 0; iPnt < numPnt; iPnt++)//loop over points of the interface
	{
		strengthPnts[iPnt] = 0.0;
		//doing superposition on selected point
		for (int i = 0; i < NsymTen; i++)
			stressSupPos[i] = 0.0;
		for (int icomp = 0; icomp < NsymTen; icomp++)
		{
			for (int isol = 0; isol < sol.size(); isol++)
			{
				//cout << sol[isol].feSol.intfSol[intfID].stressPnts[iPnt].comp[icomp]<<"\n";
				stressSupPos[icomp] = stressSupPos[icomp] + sol[isol].feSol.intfSol[intfID].stressPnts[iPnt].comp[icomp]
					* LFsuperpose[isol];
			}
		}
		//find the traction on the interface at the selected point
		//1.transform from vector to tensor
		vector<vector<double>> stressTen;
		stressTen.resize(Ndim);
		for (int i = 0; i < Ndim; i++)
			stressTen[i].resize(Ndim);
		stressTen[0][0] = stressSupPos[0];
		stressTen[0][1] = stressSupPos[2];
		stressTen[1][0] = stressSupPos[2];
		stressTen[1][1] = stressSupPos[1];
		vector<double> nVec;
		nVec.resize(Ndim);
		for (int i = 0; i < Ndim; i++)
			nVec[i] = geom.intfs[intfID].intfPnts[iPnt].normVec.comp[i];
		vector<double> traction;
		//2.find the traction
		bool buf = MultiplicationMAT_VEC(stressTen, nVec, traction);
		//computing effective stress
		double effStrs = effStresCal(traction, nVec,tangFrac);
		//double normStrength = effStrs / IntfStrength;
		double normStrength = IntfStrength / effStrs;
		if (crtNormStrength>normStrength)//higher value of strength means strongrt and so we should find the min
		{
			crtNormStrength = normStrength;
			crtPnt = iPnt;
		}
		strengthPnts[iPnt] = normStrength;
	}
	crStrng.pnt = crtPnt;//worse case over aal points of the interface
	crStrng.val = crtNormStrength;//worse case over all points of the interface

}

double intrefacePostProc::effStresCal(const vector<double> &traction, const vector<double> &normVec, double TanFac)
{
	//StresVecNorm should be normal to a plane of consideration
	double SnP = 0.0, St = 0.0, Sn = 0.0;
	int n = traction.size();
	for (int i = 0; i < n; i++)
		Sn = Sn + traction[i] * normVec[i];
	vector<double> tangVec;
	tangVec.resize(n);
	for (int i = 0; i < n; i++)
	{
		tangVec[i] = traction[i] - Sn*normVec[i];
	}
	for (int i = 0; i < n; i++)
	{
		St = St + tangVec[i] * tangVec[i];
	}
	St = sqrt(St);

	SnP = Sn;
	if (SnP < 0)
		SnP = 0;
	double EffStress = sqrt(SnP*SnP + St*St*TanFac*TanFac);
	return EffStress;
}
void tetaPostProc::compCriticalPntIntf(const vector<solution> &sol, const vector<double> &LFsuperpose, const geometry &geom, double tangFrac)
{
	double val = 1e100;
	int pnt = -1, intf = -1;
	for (int i = 0; i < intfPP.size(); i++)//loop over interfaces
	{
		if (val > intfPP[i].crStrng.val)//worse case is the weakeset case
		{
			val = intfPP[i].crStrng.val;
			pnt = intfPP[i].crStrng.pnt;
			intf = i;
		}
	}
	criticalPntIntf.intrf = intf;
	criticalPntIntf.pnt = pnt;
	criticalPntIntf.strengh = val;
#if DB_CALCULATIONS
	dbc << "_Interface_\t" << intf << "\t_pnt_\t" << pnt << "\t_strngth_\t" << val << '\n';
	dbc << "//////////\n";
	dbc << "Crd\t" << geom.intfs[intf].intfPnts[pnt].posXY.comp[0] << "\t" << geom.intfs[intf].intfPnts[pnt].posXY.comp[1] << '\n';
	dbc << "center\t" << geom.intfs[intf].center.comp[0] << "\t" << geom.intfs[intf].center.comp[1] << "\n";
	vector<double> del(2); double tmp, R = 0;
	for (int i = 0; i < 2; ++i)
	{
		tmp = geom.intfs[intf].intfPnts[pnt].posXY.comp[i] - geom.intfs[intf].center.comp[i];
		del[i] = tmp;
		R += tmp * tmp;
	}
	R = sqrt(R);
	double compAngle = atan2(del[1], del[0]);
	vector<double> n(2);
	n[0] = geom.intfs[intf].intfPnts[pnt].normVec.comp[0];
	n[1] = geom.intfs[intf].intfPnts[pnt].normVec.comp[1];
	dbc << "computedR\t" << R << "\tcompRelAngleR\t" << compAngle << "\tcompRelAngleD\t" << compAngle * 180 / PI << "\tcompNormal\t" << cos(compAngle) << "\t" << sin(compAngle) << '\n';
	dbc << "angle\t" << geom.intfs[intf].intfPnts[pnt].normAng_XY << "\tNormal\t" << n[0] << "\t" << n[1] << '\n';
	double test = n[0] * cos(compAngle) + n[1] * sin(compAngle);
	dbc << "Stresses\n";
	vector<double> strsVoigt(3);
	strsVoigt[0] = 0.0;
	strsVoigt[1] = 0.0;
	strsVoigt[2] = 0.0;
	for (int isol = 0; isol < sol.size(); ++isol)
	{
		dbc << "Sln" << isol << "(\t";
		for (int j = 0; j < sol[isol].feSol.intfSol[intf].stressPnts[pnt].comp.size(); ++j)
		{
			tmp = sol[isol].feSol.intfSol[intf].stressPnts[pnt].comp[j];
			strsVoigt[j] += tmp * LFsuperpose[isol];
			dbc << tmp << "\t";
		}
		dbc << ")\n";
	}
	double sxx = strsVoigt[0];
	double syy = strsVoigt[1];
	double sxy = strsVoigt[2];
	double tracx = sxx * n[0] + sxy * n[1];
	double tracy = sxy * n[0] + syy * n[1];
	double sn = tracx * n[0] + tracy * n[1];
	double tracShearx = tracx - sn * n[0];
	double tracSheary = tracy - sn * n[1];
	double st = sqrt(tracShearx * tracShearx + tracSheary * tracSheary);
	double snp = 0;
	if (sn > 0)
		snp = sn;
	double effectiveStress = sqrt(snp * snp + tangFrac * tangFrac * st * st);

	dbc << "LFs\t" << LFsuperpose[0] << '\t' << LFsuperpose[1] << '\t' << LFsuperpose[2] << '\n';
	dbc << "sxx\t" << sxx << "\tsyy\t" << syy << "\tsxy\t" << sxy << '\n';
	dbc << "tracx\t" << tracx << "\ttracy\t" << tracy << '\n';
	dbc << "tracShearx\t" << tracShearx << "\ttracSheary\t" << tracSheary << '\n';
	dbc << "sn\t" << sn << "\tst\t" << st << '\n';
	dbc << "snp\t" << snp << "\ttangFrac\t" << tangFrac << '\n';
	dbc << "sEff\t" << effectiveStress << '\n';
	double compStrength = 1.0 / effectiveStress;

	double maxStrength = Max(compStrength, criticalPntIntf.strengh);
	if (fabs(maxStrength) < 1e-40)
		maxStrength = 1.0;
	double relError = fabs(compStrength - criticalPntIntf.strengh) / maxStrength;
	dbc << "compStrength\t" << compStrength << "\trelError\t" << relError << "\n";

	dbc.flush();
	if (fabs(test - 1.0) > 1e-3)
	{
		cout << "n\t" << n[0] << '\t' << n[1] << '\n';
		cout << "nCompted\t" << cos(compAngle) << '\t' << sin(compAngle) << '\n';
		cout << "n.nComputed\t" << test << '\n';
		THROW("Normal not computed correctly\n");
	}
	if (relError > 1e-3)
	{
		cout << "compStrength\t" << compStrength << "\trelError\t" << relError << "\n";
		THROW("Strength not computed correctly\n");
	}
#endif
}

void angularProperty::summPostProc_teta()
{
	double crtStrngIntf = 1e200, crtTeta = -1;
	int crtPnt = -1, crtIntf = -1;
	for (int i = 0; i < numTeta; i++)//loop over all angles
	{
		if (crtStrngIntf > postProc[i].criticalPntIntf.strengh)
		{
			crtStrngIntf = postProc[i].criticalPntIntf.strengh;
			crtIntf = postProc[i].criticalPntIntf.intrf;
			crtPnt = postProc[i].criticalPntIntf.pnt;
			crtTeta = postProc[i].angle;
		}
	}
	criticalAngle.intfaceSumm.intf = crtIntf;
	criticalAngle.intfaceSumm.pnt = crtPnt;
	criticalAngle.intfaceSumm.teta = crtTeta;
	criticalAngle.intfaceSumm.strength = crtStrngIntf;
}

void modeMixity::computModeMixity_intf(const vector<trialLoad> &trialLoads)
{
	modeMix_Intf.resize(Nteta);
	//computing mode mixity at each point
	for (int iteta = 0; iteta < Nteta; iteta++)
	{
		modeMix_Intf[iteta] = trialLoads[0].angProp.postProc[iteta].criticalPntIntf.strengh
			/ trialLoads[1].angProp.postProc[iteta].criticalPntIntf.strengh;

	}

}
void modeMixity::compStatModMix_Intf()
{

	CompStatistics(modeMix_Intf, StatVar_Intf);
}
void angularProperty::compStat_teta()
{
	vector<double> tempVal;
	tempVal.resize(numTeta);
	for (int i = 0; i < numTeta; i++)//loop over all angles
	{
		tempVal[i] = postProc[i].criticalPntIntf.strengh;
	}
	CompStatistics(tempVal, statVar_teta);
}
void sve::set0step(double incAngle)
{
	numSol = Nloading;
	numTrialLoads = NtriLoad;
	int Nsym = Ndim*(Ndim + 1) / 2;
	sol.resize(numSol);
	sveId = -1;
	trialLoads.resize(numTrialLoads);
	for (int i = 0; i < numTrialLoads; i++)
	{
		trialLoads[i].set0step(incAngle);
	}
	modemix.IncrAngl = incAngle;
	modemix.Nteta = 180 / incAngle + 1;;
}
void sve::runSVE(int irve, int isve, double tanFrac, double incAngl, string folderLoc)
{
#if	DB_CALCULATIONS
	string fileName = "dbCalcRVE";
	string buf;
	ostringstream convert, convert2;
	convert << irve;
	buf = convert.str();
	fileName += buf;
	fileName += "_SVE_";
	convert2 << isve;
	buf = convert2.str();
	fileName += buf;
	fileName += ".txt";
	dbc.open(fileName.c_str(), ios::out);
#endif

	bool sympart, antiSympart, res;
	readSVEdata(irve, folderLoc);
	sympart = true;
	antiSympart = true;
	res = stiffCalc(sympart, antiSympart);//compute stiffness
	res = compBulkMod();//compute bulk modulus
	geom.compGeomProp();//compute geometry properties
	comptrialLoads(tanFrac);
	modemix.compModeMixity(trialLoads);
	writeResults(irve + 1, incAngl,  folderLoc);
	deletData();

#if	DB_CALCULATIONS
	dbc.close();
#endif
}
void interfGeom::compGeomProp()
{
	for (int ip = 0; ip < numIntPnts; ip++)
	{
		intfPnts[ip].normalCalculation(center);
	}
}

void geometry::compGeomProp()
{
	for (int ii = 0; ii < numIntf; ii++)
	{
		intfs[ii].compGeomProp();
	}
}
void sve::comptrialLoads(double tanFrac)
{
	
	trialLoads[0].loadname = tensile;
	trialLoads[1].loadname = shear;
	for (int il = 0; il < numTrialLoads; il++)
	{
#if DB_CALCULATIONS
		if (il == 0)
			dbc << "TENSILE\n++++++++++++++++++++++++++++++++++++++\n";
		else
			dbc << "SHEAR\n++++++++++++++++++++++++++++++++++++++\n";
		dbc << "tanFrac\t" << tanFrac << '\n';
#endif
		trialLoads[il].initTrialLoad();
		trialLoads[il].angProp.compPostProc(geom.numIntf, trialLoads[il].auxLoad, sol, geom, tanFrac);
	}
}
void tetaPostProc::compIntfPP(const int &numIntf, const vector<double> &auxload, const vector<solution> &sol, const geometry &geom, double tanFrac)
{
	intfPP.resize(numIntf);
	compLFsuperpose(auxload, sol);
	for (int ii = 0; ii < numIntf; ii++)//loop over interfaces
	{
		intfPP[ii].intfID = ii;
		intfPP[ii].failCalc(sol, LFsuperpose, geom,tanFrac);
	}
	compCriticalPntIntf(sol, LFsuperpose, geom, tanFrac);
}
void angularProperty::compPostProc(const int &numIntf, const  vector<double> &auxload, const vector<solution> &sol, const geometry &geom, double tanFrac)
{
	for (int i = 0; i < numTeta; i++)
	{
		postProc[i].angle = dteta*i;
#if DB_CALCULATIONS
		dbc << "\n=========\nanglei" << i << "\tDegree\t" << postProc[i].angle << '\n';
#endif
		postProc[i].compIntfPP(numIntf, auxload, sol, geom, tanFrac);
	}
	compStat_teta();
	summPostProc_teta();
}

modeMixity::modeMixity()
{
	//Nteta = 180 / IncrAngl + 1;
}

void modeMixity::compModeMixity(const vector<trialLoad> &trialLoads)
{
	 computModeMixity_intf(trialLoads);
	 compStatModMix_Intf();
}

void modeMixity::writeResults(const int iRVE, const int iSVE, string folderName)
{
	//iRVE and iSVE should be in natural numbers 
	ofstream file;
	ostringstream conv, conv2;
	string id;
	conv << iRVE; conv2 << iSVE;
	//id = "OUTPUT//RVE" + conv.str() + "SVE" + conv2.str() + "(ModeMixity_Angle).OMM";
	id = "OUTPUT//"+folderName+ "SVE" + conv2.str() + "(ModeMixity_Angle).OMM";
	file.open(id);
	file << "RVE\t SVE\t dteta\t Nteta\n";
	file << iRVE << "\t" <<iSVE << "\t" << IncrAngl << "\t" << Nteta << "\n";
	file << "Angle\t ModeMixity\n";
	for (int i = 0; i < Nteta; i++)
	{
		file << i*IncrAngl << "\t" << modeMix_Intf[i] << "\n";
	}
	file << "number of Stat Var\n" << NUM_STAT_VAR << "\n";
	file << "Min\t Max\t Mean\t StdDiv\n";
	for (int i = 0; i < NUM_STAT_VAR; i++)
		file << StatVar_Intf[i] << "\t";
	file << "\n";
	file.close();
}

void trialLoad::set0step(double incAngl)
{
	angProp.dteta = incAngl;
	angProp.numTeta = 180 / angProp.dteta + 1;
	angProp.postProc.resize(angProp.numTeta);
}
void trialLoad::writeResults(const int iRVE, const int iSVE, double incAngl,string folderName)
{
	//iRVE, itrLoad, and iSVE should be in natural numbers 
	fstream file;
	ostringstream conv, conv2, conv3;
	string id;

	conv << iRVE; conv2 << iSVE; conv3 << loadname;
	if (loadname == tensile)
	{
		//id = "OUTPUT//RVE" + conv.str() + "SVE" + conv2.str() + "(TStrn_Angle).OTS";
		id = "OUTPUT//" +folderName+ "SVE" + conv2.str() + "(TStrn_Angle).OTS";

	}
		
	else if (loadname == shear)
	{
		//id = "OUTPUT//RVE" + conv.str() + "SVE" + conv2.str() + "(SStrn_Angle).OSS";
		id = "OUTPUT//" + folderName+ "SVE" + conv2.str() + "(SStrn_Angle).OSS";
	}
	else
		cout << "Errror in output of angular properties\n";
	
	file.open(id.c_str(), ios::out);
	if (!file.is_open())
	{
		cout << id << '\n';
		THROW("cannot open file\n");
	}
	file << "RVE\t SVE\t dteta\t Nteta\n";
	file << iRVE << "\t" << iSVE << "\t" << incAngl << "\t" << angProp.numTeta << "\n";
	if (loadname == tensile)
		file << "Angle\t TensileStrength\n";
	else if (loadname == shear)
		file << "Angle\t ShearStrength\n";
	else
		cout << "Errror in output of angular properties\n";
	
	for (int i = 0; i < angProp.numTeta; i++)
	{
		file << i*incAngl << "\t" << angProp.postProc[i].criticalPntIntf.strengh << "\n";
	}
	file << "number of Stat Var\n" << NUM_STAT_VAR << "\n";
	file << "Min\t Max\t Mean\t StdDiv\n";
	for (int i = 0; i < NUM_STAT_VAR; i++)
		file << angProp.statVar_teta[i] << "\t";
	file << "\n";
	file.close();
}
void sve::writeResults(const int irve, double incAngl, string folderName)
{
	//irve should be a natural number
	//sveId should be in a natural number
	for (int il = 0; il < numTrialLoads; il++)
	{
		trialLoads[il].writeResults(irve, sveId, incAngl, folderName);
	}
	modemix.writeResults(irve, sveId,folderName);
}
void sve::deletData()
{
	geom.intfs.clear();
	//stiff.antiSymStiff.clear();
	//stiff.symStiff.clear();
	//stiff.realStiff.clear();
	sol.clear();
}

bool sve::compBulkMod()
{
	int Nsymten = Ndim*(Ndim + 1) / 2;
	vector<double> Hstrs(Nsymten);//hydrostatic stress
	vector<double> Hstrn(Nsymten);
	Hstrs[0] = 1.0;//unit load
	Hstrs[1] = 1.0;
	vector<vector<double>> invStiff;
	invStiff.resize(Nsymten);
	for (int i = 0; i < Nsymten; i++)
	{ 
		invStiff[i].resize(Nsymten);
		for (int j = 0; j < Nsymten; j++)
			invStiff[i][j] = 0;
	}
		
	bool invR = inversMat(stiff.realStiff, invStiff);
	if (invR == false)
	{
		cout << "ERROR in computin bulk modulus";
		return false;
	}
	bool I1 = MultiplicationMAT_VEC(invStiff, Hstrs, Hstrn);
	bulkMod_real = 1.0 / (Hstrn[0] + Hstrn[1]);

	for (int i = 0; i < Nsymten; i++)
		for (int j = 0; j < Nsymten; j++)
			invStiff[i][j] = 0;
	invR = inversMat(stiff.symStiff, invStiff);
	if (invR == false)
	{
		cout << "ERROR in computin bulk modulus";
		return false;
	}

	bool I2 = MultiplicationMAT_VEC(invStiff, Hstrs, Hstrn);
	bulkModu_sym = 1.0 / (Hstrn[0] + Hstrn[1]);
	return I1&&I2;
}