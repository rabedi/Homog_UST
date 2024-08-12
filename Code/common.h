#ifndef COMMON__H
#define COMMON__H
#include <vector>
#include <cstring>
#include <string>
#include <sstream>
#include <cstdlib>


#define MAX(a,b) (a<b)?b:a
#define MIN(a,b) (a>b)?b:a
#define PI 3.14159265359
typedef enum{ min, max, mean, stdDiv, NUM_STAT_VAR } STAT_VAR;

using namespace std;


class vecDim
{
public:
	vecDim();
	vector<double> comp;
};
class symTen
{
public:
	symTen();
	vector<double> comp;
};


bool inversMat(vector<vector<double>> mat, vector<vector<double>>& invMat);
bool matrixMultiplication(vector<vector<double>> mat1, vector<vector<double>> mat2, vector<vector<double>>& mat);
void CompStressRotation(const vector<double>& StresVec0, const double Teta, vector<double>& StresVec1);
bool MultiplicationMAT_VEC(const vector<vector<double>>& mat1, const vector<double>& vec1, vector<double>& vec);
void CompStatistics(const vector<double>& vec, vector<double> &StatVar);

void CreateFolder(const char* path);

template<class T>
bool fromString(string& name, T& dat)
{
	istringstream ss(name);
	ss >> dat;
	if (ss.fail())
		return false;
	return true;

}

template<class T>
void toString(T& dat, string& name)
{
	ostringstream convert;
	convert << dat;
	name = convert.str();
}

bool isDigit(string& buf);

#endif