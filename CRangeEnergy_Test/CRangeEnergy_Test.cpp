#include <assert.h>

#include "..\CRangeEnergy\Range_Energy_Relation.cpp"
#include "..\CRangeEnergy\CRangeEnergy.h"


#include <cmath>
#include <limits>

//http://stackoverflow.com/questions/17333/most-effective-way-for-float-and-double-comparison
bool AreSame(double a, double b) {
	return std::fabs(a - b) < 1E-10;//std::numeric_limits<double>::epsilon();
}



int main(void){


	double Mass[5] = {938.272, 3727.380, 4667.845, 12109.485, 15074.365};//tekitou na mass
	int Z[5] = {1,2,3,6,7};
	double densityEM = 3.6;
	
	for(int i=0; i<5; i++){
		for(double Range = 0.1;  Range<10000; Range+=5.0){
			CRangeEnergy cre;
			double a = cre.GetKineticEnergyFromRange(Mass[i], Range*1.0, Z[i], densityEM);
			double b = function0(Mass[i], Range, Z[i], densityEM, 0.884);//r=0.884;
			assert(AreSame( a, b));
		}
		printf("i=%d OK\n",i);
	}

	return 0;
}