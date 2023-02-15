#include <con2020.h>
#include <stdio.h>
#define _USE_MATH_DEFINES
#include <math.h>

int main() {
	printf("Testing C++\n");
	double x = 10.0;
	double y = 20.0;
	double z = 5.0;
	double Bx, By, Bz;

	Con2020Field(x,y,z,&Bx,&By,&Bz);

	printf("B = [%5.1f,%5.1f,%5.1f] at [%4.1f,%4.1f,%4.1f]\n",Bx,By,Bz,x,y,z);

	double thetamm, dthetamm, thetaoc, dthetaoc;
	thetamm = con2020.GetThetaMM();
	dthetamm = con2020.GetdThetaMM();
	thetaoc = con2020.GetThetaOC();
	dthetaoc = con2020.GetdThetaOC();
	printf("thetamm = %f\n",thetamm);
	printf("dthetamm = %f\n",dthetamm);
	printf("thetaoc = %f\n",thetaoc);
	printf("dthetaoc = %f\n",dthetaoc);

	thetamm = con2020.GetThetaMM();
	dthetamm = con2020.GetdThetaMM();
	thetaoc = con2020.GetThetaOC();
	dthetaoc = con2020.GetdThetaOC();
	printf("thetamm = %f\n",thetamm);
	printf("dthetamm = %f\n",dthetamm);
	printf("thetaoc = %f\n",thetaoc);
	printf("dthetaoc = %f\n",dthetaoc);

}
