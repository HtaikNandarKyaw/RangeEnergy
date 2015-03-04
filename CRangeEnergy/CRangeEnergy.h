#pragma once
#include <math.h>

class CRangeEnergy
{
public:
	CRangeEnergy(void);
	~CRangeEnergy(void);

	double function0(double Mass,double Range,int Z,double D,double r);
	double function1(double Mass,double KE,int Z,double D,double r);
	double CRangeEnergy::Rs_function1(double LR);
	double CRangeEnergy::Rs_function2(double LR);
	double CRangeEnergy::Rs_function3(double LR);

	double rate_grobal,Rp_grobal,Rs_grobal;


};

