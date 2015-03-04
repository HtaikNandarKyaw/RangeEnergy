#pragma once
#include <math.h>

class CRangeEnergy
{
public:
	CRangeEnergy(void);
	~CRangeEnergy(void);

	double function0(double Mass,double Range,int Z,double D,double r);
	double function1(double Mass,double KE,int Z,double D,double r);
	double Rs_function1(double LR);
	double Rs_function2(double LR);
	double Rs_function3(double LR);
	double FunctionCz(int Z, double B);
	double FunctionRs(double KEM, double MKEM);
	double FunctionRsRwRatio(double Rs);


};

