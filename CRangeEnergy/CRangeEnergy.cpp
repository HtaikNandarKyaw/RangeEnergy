#include <stdio.h>
#include <math.h>

#include "CRangeEnergy.h"

#pragma warning(disable : 4996)


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
		+0.000010917*LR*LR*LR*LR*LR*LR*LR;

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
		-0.000009219*LR*LR*LR*LR*LR*LR*LR;

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
		+0.00001057489*LR*LR*LR*LR*LR;

	return LK;
}



double CRangeEnergy::FunctionRs(double KE, double Mass){

	double KEM = KE/Mass;// KE/mass
	//double LKEM = log10(KEM);//log10(K.E. / ThisMass)
	double MKEM = log(KE*Mp/Mass);//log10(K.E. * protonmass / ThisMass)

	double dd = 0.00001;//step for italation
	double d0;//
	double y0;//
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
		-0.000000951*LRs*LRs*LRs*LRs*LRs*LRs*LRs;

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
			-0.006004180*FX*FX*FX*FX*FX*FX*FX;
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



/*
int main()
{
	int Z,A,S;
	double Mass;
	double Range;

loop:

	double M[4][19][8] ={{{938.272,0.,0.,0.,0.,0.,0.,0.},
	{1875.613,0.,0.,0.,0.,0.,0.,0.},
	{2808.922,2808.392,0.,0.,0.,0.,0.,0.},
	{0.,3727.380,0.,0.,0.,0.,0.,0.},
	{0.,4667.845,4667.624,0.,0.,0.,0.,0.},
	{0.,5605.541,5601.520,5605.298,0.,0.,0.,0.},
	{0.,6545.550,6533.836,6534.186,0.,0.,0.,0.},
	{0.,7482.542,7471.369,7454.852,7472.321,0.,0.,0.},
	{0.,0.,8406.871,8392.753,8393.310,8409.295,0.,0.},
	{0.,0.,0.,9325.507,9324.440,9327.580,0.,0.},
	{0.,0.,10285.846,10264.570,10252.550,10254.022,0.,0.},
	{0.,0.,0.,11200.922,11188.746,11174.866,11191.693,0.},
	{0.,0.,0.,12142.542,12123.434,12109.485,12111.195,12128.443},
	{0.,0.,0.,0.,13062.023,13040.874,13040.207,13044.841},
	{0.,0.,0.,0.,0.,13979.222,13968.939,13971.182},
	{0.,0.,0.,0.,0.,14914.536,14906.014,14895.084},
	{0.,0.,0.,0.,0.,0.,0.,15830.506},
	{0.,0.,0.,0.,0.,0.,0.,0.},
	{0.,0.,0.,0.,0.,0.,0.,0.}},
	{{0.,0.,0.,0.,0.,0.,0.,0.},
	{0.,0.,0.,0.,0.,0.,0.,0.},
	{2991.166,0.,0.,0.,0.,0.,0.,0.},
	{3922.565,3921.685,0.,0.,0.,0.,0.,0.},
	{0.,4839.943,0.,0.,0.,0.,0.,0.},
	{0.,5779.348,5778.807,0.,0.,0.,0.,0.},
	{0.,0.,6711.623,6715.821,0.,0.,0.,0.},
	{0.,7654.073,7642.719,7643.029,0.,0.,0.,0.},
	{0.,0.,8578.552,8563.825,8579.714,0.,0.,0.},
	{0.,0.,0.,9499.326,9500.103,0.,0.,0.},

	{0.,0.,0.,0.,10429.883,0.,0.,0.},

	{0.,0.,0.,0.,11356.863,11358.905,0.,0.},

	{0.,0.,0.,0.,12293.059,12278.859,0.,0.},

	{0.,0.,0.,0.,0.,13212.998,13214.708,0.},

	{0.,0.,0.,0.,0.,0.,14142.300,0.},

	{0.,0.,0.,0.,0.,0.,0.,15074.365},

	{0.,0.,0.,0.,0.,0.,0.,0.},

	{0.,0.,0.,0.,0.,0.,0.,16931.689},

	{0.,0.,0.,0.,0.,0.,0.,0.}},

	{{0.,0.,0.,0.,0.,0.,0.,0.},

	{0.,0.,0.,0.,0.,0.,0.,0.},

	{0.,0.,0.,0.,0.,0.,0.,0.},

	{4106.719,0.,0.,0.,0.,0.,0.,0.},

	{5036.208,5034.978,0.,0.,0.,0.,0.,0.},

	{0.,5952.506,0.,0.,0.,0.,0.,0.},

	{0.,6890.851,6889.990,0.,0.,0.,0.,0.},

	{0.,0.,7821.726,7826.344,0.,0.,0.,0.},

	{0.,8762.596,8751.602,8751.872,0.,0.,0.,0.},

	{0.,0.,9685.735,9672.798,9687.107,0.,0.,0.},

	{0.,0.,0.,10605.899,10606.896,0.,0.,0.},

	{0.,0.,0.,11538.653,11535.326,0.,0.,0.},

	{0.,0.,0.,0.,12461.176,12463.788,0.,0.},

	{0.,0.,0.,0.,13397.372,13382.852,0.,0.},

	{0.,0.,0.,0.,0.,14316.511,14318.221,0.},

	{0.,0.,0.,0.,0.,15247.900,15244.393,0.},

	{0.,0.,0.,0.,0.,0.,0.,16177.548},

	{0.,0.,0.,0.,0.,0.,0.,0.},

	{0.,0.,0.,0.,0.,0.,0.,18032.872}},

	{{139.570,0.,0.,0.,0.,0.,0.,0.},

	{134.977,0.,0.,0.,0.,0.,0.,0.},

	{493.677,0.,0.,0.,0.,0.,0.,0.},

	{497.648,0.,0.,0.,0.,0.,0.,0.},

	{1869.4,0.,0.,0.,0.,0.,0.,0.},

	{1115.683,0.,0.,0.,0.,0.,0.,0.},

	{1189.37,0.,0.,0.,0.,0.,0.,0.},

	{1192.449,0.,0.,0.,0.,0.,0.,0.},

	{1314.83,0.,0.,0.,0.,0.,0.,0.},

	{1321.31,0.,0.,0.,0.,0.,0.,0.},

	{1672.45,0.,0.,0.,0.,0.,0.,0.},

	{0.,0.,0.,0.,0.,0.,0.,0.},

	{0.,0.,0.,0.,0.,0.,0.,0.},

	{0.,0.,0.,0.,0.,0.,0.,0.},

	{0.,0.,0.,0.,0.,0.,0.,0.},

	{0.,0.,0.,0.,0.,0.,0.,0.},

	{0.,0.,0.,0.,0.,0.,0.,0.},

	{0.,0.,0.,0.,0.,0.,0.,0.},

	{0.,0.,0.,0.,0.,0.,0.,0.}}};








	//	Z = 1;
	double D;
	double r=0.884;
	//	Mass = 938.272;



	printf("\n");
	printf(" Z = ");
	scanf_s("%d",&Z);
	printf(" A = ");
	scanf_s("%d",&A);
	printf(" S = ");
	scanf_s("%d",&S);


	if(S!=2)
	{
		Mass = M[S][A-1][Z-1];
	}
	else
	{
		printf(" Mass = ");
		scanf_s("%lf",&Mass);
	}


	printf(" Density [g/cm^3] = ");

	scanf_s("%lf",&D);

	printf(" Range = ");

	scanf_s("%lf",&Range);




	// Range = 37.3;
	// double Range_up = 69.4;//2495.9;
	// double Range_down = 68.2;//2488.1;


	// Z=1;
	// A=2;
	// S=0;

	//for(Z = 1 ; Z < 4 ; Z += 1)
	//	for(A = 1 ; A < 8 ; A += 1)
	//for(Range = 1.0 ;Range <  300.0; Range += 1.0)

	//	{
	//	{
	//	Mass =M[S][A-1][Z-1];
	//	
	//	if(Mass == 0)
	//	{
	//	continue;
	//	}
	//	double Range_up = Range + 68.8;
	//	if(Range_up > 6370.2)
	//	{
	//	return 0;
	//	}


	double answer = function0(Mass,Range,Z,D,r);
	// double answer1 = function0(Mass,Range_up,Z,D,r);
	//   double answer2 = function0(Mass,Range_down,Z,D,r);
	double total_Energy = Mass+answer;
	// double total_Energy1 = Mass+answer1;
	//   double total_Energy2 = Mass+answer2;
	double momentum = sqrt(total_Energy*total_Energy-Mass*Mass);
	// double momentum1 = sqrt(total_Energy1*total_Energy1-Mass*Mass);
	//   double momentum2 = sqrt(total_Energy2*total_Energy2-Mass*Mass);
	//   double err_KE = abs(answer1-answer2)/2;
	//   double err_momentum = abs(momentum1-momentum2)/2;
	//  double delta_momentum = abs(momentum1-momentum);



	printf(" Mass = %.3f MeV/c^2   KE = %.3f MeV   P = %.3f MeV/c \n",Mass,answer,momentum);
	// printf("%d  %d  %.3f  %.3f  %.3f  %.3f  %.3f\n",Z,A,Mass,answer,err_KE,momentum,err_momentum);
	goto loop;

	// getchar();
	//  FILE *file;
	//	file = fopen("result.txt","a");
	//	fprintf(file,"%.1f  %.3f\n",Range,answer);
	//	fprintf(file,"%d  %d  %.3f  %.3f  %.3f  %.3f  %.3f\n",Z,A,Mass,answer,err_KE,momentum,err_momentum);
	//	fclose(file);
	//	}
	//	}


	return 0;
} 
*/