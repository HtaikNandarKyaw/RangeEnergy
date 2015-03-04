#include <stdio.h>
#include <math.h>

#include "CRangeEnergy.h"


static const double Mp = 938.272;//proton mass
static const double LMp = log(Mp);//log(proton mass)
static const double D0 = 3.815;//density of standard emulsion
static const double r = 0.884;//



CRangeEnergy::CRangeEnergy(void)
{
}


CRangeEnergy::~CRangeEnergy(void)
{
}


double CRangeEnergy::GetKineticEnergyFromRange(double Mass,double Range,int Z,double densityEM)
{
	if( Range==0.0 ) return 0.0;

	double R = 0.00;
	double R0 = 0.00;
	double KE = 0.00;
	double KE0 = 0.00;
	double dKE = 10.0;

	while(dKE>0.00005){

		(Range > R0)?  KE = KE + dKE : 	KE = KE - dKE;

		if(KE<=0) KE = 0.0;

		double return_R = GetRangeFromKineticEnergy(Mass,KE,Z,densityEM);

		if((return_R>=Range && Range>=R0) || (return_R<=Range && Range<=R0))
		{
			dKE = dKE/10.0;
			KE = KE0+(KE-KE0)*(Range-R0)/(return_R-R0); 
			if(KE<=0) KE = 0.0;
			R0 = GetRangeFromKineticEnergy(Mass,KE,Z,densityEM);
		}
		else
		{
			R0 = return_R;
		}

		KE0 = KE;
	}//while

	//double E = Mass+KE;
	//double P = sqrt(E*E-Mass*Mass);

	return KE;
}



double CRangeEnergy::Rs_function1(double LR)
{

	double LK =
		-2.288460778
		+1.382747508*LR
		-0.439300692*LR*LR
		+0.162697682*LR*LR*LR
		-0.037735480*LR*LR*LR*LR
		+0.005152047*LR*LR*LR*LR*LR
		-0.000373872*LR*LR*LR*LR*LR*LR
		+0.000010917*LR*LR*LR*LR*LR*LR*LR;//fitted by Mishina

	return LK;
}


double CRangeEnergy::Rs_function2(double LR)
{
	
	double LK =
		12.499454326
		-12.637449190*LR
		+5.296813187*LR*LR
		-1.163641812*LR*LR*LR
		+0.151898030*LR*LR*LR*LR
		-0.011803694*LR*LR*LR*LR*LR
		+0.000505820*LR*LR*LR*LR*LR*LR
		-0.000009219*LR*LR*LR*LR*LR*LR*LR;//fitted by Mishina

	return LK;

}


double CRangeEnergy::Rs_function3(double LR)
{

	double LK =
		-0.52629642
		+0.31555326*LR
		+0.021856192*LR*LR
		+0.0012217823*LR*LR*LR
		-0.00026892371*LR*LR*LR*LR
		+0.00001057489*LR*LR*LR*LR*LR;//fitted by Mishina

	return LK;
}



double CRangeEnergy::FunctionRs(double KE, double Mass){

	double KEM = KE/Mass;// KE/mass
	double MKEM = log(KE*Mp/Mass);//log10(K.E. * protonmass / ThisMass)

	double dd = 0.00001;//step for italation
	double d0;
	double y0;
	double Rs;


	if(KEM < 0.0001)
	{

		Rs = 479.210*pow(KEM,0.675899);

	}
	else if(MKEM < 1.930606146327)
	{

		d0 = 3.0000; 
		y0 = Rs_function1(d0);
		while(abs(MKEM-y0)>0.00001)
		{
			(MKEM > y0)?  d0 = d0 + dd	 : 	d0 = d0 - dd;
			y0 = Rs_function1(d0);
		}
		Rs = exp(d0);

	}
	else if(MKEM < 37.156634656805)
	{

		d0 = 6.0000; 
		y0 = Rs_function2(d0);
		while(abs(MKEM-y0)>0.00001)
		{
			(MKEM > y0)?  d0 = d0 + dd	 : 	d0 = d0 - dd;
			y0 = Rs_function2(d0);
		}
		Rs = exp(d0);
	}
	else
	{

		d0 = 10.0000;
		y0 = Rs_function3(d0);
		while(abs(MKEM-y0)>0.00001)
		{
			(MKEM > y0)?  d0 = d0 + dd	 : 	d0 = d0 - dd;
			y0 = Rs_function3(d0);
		}
		Rs = exp(d0);
	}

	return Rs;

}



double CRangeEnergy::FunctionRsRwRatio(double Rs)
{
	const double LRs = log(Rs);

	double rate =
		-0.107714711
		-0.332543998*LRs
		+0.141029694*LRs*LRs
		-0.044679440*LRs*LRs*LRs
		+0.008162611*LRs*LRs*LRs*LRs
		-0.000830409*LRs*LRs*LRs*LRs*LRs
		+0.000044038*LRs*LRs*LRs*LRs*LRs*LRs
		-0.000000951*LRs*LRs*LRs*LRs*LRs*LRs*LRs;//fitted by D.Tovee and Gajewski

	return  exp(rate);
}



double CRangeEnergy::FunctionCz(int Z, double beta)
{
	if(Z==1) return 0.0;

	const double FX = 137.0*beta/Z;

	if(FX<=0.5)//regionI: a*FX^b
	{
		return  0.168550736771407*pow(FX,1.90707106569386);
	}
	else if(FX<=2.51)//regionII: polinominal7
	{
		return
			0.002624371
			-0.081622520*FX
			+0.643381535*FX*FX
			-0.903648583*FX*FX*FX
			+0.697505012*FX*FX*FX*FX
			-0.302935572*FX*FX*FX*FX*FX
			+0.067662990*FX*FX*FX*FX*FX*FX
			-0.006004180*FX*FX*FX*FX*FX*FX*FX;//by Mishina
	}
	else//regionIII: constant
	{
		return 0.217598079611354;
	}

}



double CRangeEnergy::GetRangeFromKineticEnergy(double Mass, double KE, int Z, double densityEM)
{
	if(KE <= 0.0) return 0.0;

	const double CPS = 1;//correction factor for 1st term
	const double CPM = 1;//correction factor for 2nd term
	const double CF = 1;//correction factor


	double Rs = FunctionRs(KE, Mass);//proton range in standard emulsion
	double ratio = FunctionRsRwRatio(Rs);// Rs/Rw ratio
	double F = densityEM/D0 + ((r*(D0-densityEM))/(r*D0-1.0))*ratio;//factor for proton-range
	double Rp = Rs/F;//proton range in this emulsion


	double E = Mass + KE;//total energy
	double P = sqrt(E*E-Mass*Mass);//momentum norm
	double beta = P/E;//beta of particle
	double Cz = FunctionCz(Z, beta);


	double R1 = CPS * (Mass/Mp) / (Z*Z) * Rp;
	double R2 = CPM * (Mass/Mp) * pow(Z,2.0/3.0) * Cz;
	double R = (R1+R2)/CF;

	return R;
}

