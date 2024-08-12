#include "common.h"
#include "params.h"
#include <Windows.h>
vecDim::vecDim()
{
	comp.resize(Ndim);
	for (int i = 0; i < Ndim; i++)
		comp[i] = 0.0;
}
symTen::symTen()
{
	int n = Ndim*(Ndim + 1) / 2;
	comp.resize(n);
	for (int i = 0; i < n; i++)
		comp[i] = 0.0;
}

bool inversMat(vector<vector<double>> mat, vector<vector<double>>& invMat)
{
	int n = mat.size();
	if (n != 3)
		return false;
	switch (n)
	{
	case 3:
	{
		double a11, a12, a13, a21, a22, a23, a31, a32, a33;
		a11 = mat[0][0]; a21 = mat[1][0]; a31 = mat[2][0];
		a12 = mat[0][1]; a22 = mat[1][1]; a32 = mat[2][1];
		a13 = mat[0][2]; a23 = mat[1][2]; a33 = mat[2][2];
		double det = a11*a22*a33 + a21*a32*a13 + a31*a12*a23 - a11*a32*a23 - a31*a22*a13 - a21*a12*a33;
		if (det == 0)
			return false;
		invMat[0][0] = (1.0 / det)*(a22*a33 - a23*a32); invMat[0][1] = (1.0 / det)*(a13*a32 - a12*a33);
		invMat[0][2] = (1.0 / det)*(a12*a23 - a13*a22); invMat[1][0] = (1.0 / det)*(a23*a31 - a21*a33);
		invMat[1][1] = (1.0 / det)*(a11*a33 - a13*a31); invMat[1][2] = (1.0 / det)*(a13*a21 - a11*a23);
		invMat[2][0] = (1.0 / det)*(a21*a32 - a22*a31); invMat[2][1] = (1.0 / det)*(a12*a31 - a11*a32);
		invMat[2][2] = (1.0 / det)*(a11*a22 - a12*a21);
		return true;
	}
	break;
	}
}

bool matrixMultiplication(vector<vector<double>> mat1, vector<vector<double>> mat2, vector<vector<double>>& mat)
{
	int r1, r2, c1, c2;
	r1 = mat1.size(); c1 = mat1[0].size();
	r2 = mat2.size(); c2 = mat2[0].size();
	if (c1 != r2)
		return false;
	for (int i = 0; i < r1; i++)
	{
		for (int j = 0; j < c2; j++)
		{
			double ij = 0;
			for (int k = 0; k < c1; k++)
				ij = ij + mat1[i][k] * mat2[k][j];
			mat[i][j] = ij;
		}
	}
}

void CompStressRotation(const vector<double>& StresVec0, const double Teta, vector<double>& StresVec1)
{
	//Alredy it works for 2D!
	//for 2D case::Stress vector=[Snn, Stt, Snt]
	//Teta is positive counter clockwisely!
	//Teta should be define in Radian!
	unsigned int N;
	N = StresVec0.size();
	StresVec1.resize(N);
	for (int i = 0; i < N; i++)
		StresVec1[i] = 0.0;
	if (Ndim == 2)
	{
		StresVec1[0] = (StresVec0[0] + StresVec0[1]) / 2.0 + (StresVec0[0] - StresVec0[1]) / 2.0*cos(2.0*Teta) + StresVec0[2] * sin(2.0*Teta);
		StresVec1[1] = StresVec0[0] + StresVec0[1] - StresVec1[0];
		StresVec1[2] = -(StresVec0[0] - StresVec0[1]) / 2.0*sin(2.0*Teta) + StresVec0[2] * cos(2.0*Teta);
	}

}

bool MultiplicationMAT_VEC(const vector<vector<double>>& mat1, const vector<double>& vec1, vector<double>& vec)
{
	int r1, r2, c1, c2;
	r1 = mat1.size(); c1 = mat1[0].size();
	r2 = vec1.size();
	if (c1 != r2)
		return false;
	vec.clear();
	vec.resize(r1);
	for (int i = 0; i < r1; i++)
	{
		for (int j = 0; j < c1; j++)
		{
			vec[i] = vec[i] + mat1[i][j] * vec1[j];
		}
	}
}



void CompStatistics(const vector<double>& vec,vector<double> &StatVar)
{
	StatVar.resize(NUM_STAT_VAR);
	unsigned int n;
	n = vec.size();
	double l, s, x2;
	double mean0 = 0.0; 
	l = 1.0e-160; 
	s = 1.0 / l; x2 = 0.0;
	for (unsigned int i = 0; i < n; i++)
	{
		if (l < vec[i]) l = vec[i];
		if (s > vec[i]) s = vec[i];
		mean0 = mean0 + vec[i];
		x2 = x2 + vec[i] * vec[i];
	}
	StatVar[min] = s;
	StatVar[max] = l;
	StatVar[mean] = mean0 / n;
	StatVar[stdDiv] = sqrt(x2 / n - StatVar[mean] * StatVar[mean]);
}

void CreateFolder(const char* path)
{

	if (!CreateDirectory(path, NULL))
	{
		return;
	}
}


bool isDigit(string& buf)
{
	return (isdigit(buf[0]) || ((buf[0] == '-') && (isdigit(buf[1]))));
}
