#include "lmic.h"



doube f_thetai(double thetai) {

	double tanti = tan(thetai);
	return 1.0 + 0.25*tanti*tanti;
}


double OmegaRatio(	double thetai, double wO_open, double wO_om,
					double thetamm, double dthetamm,
					double thetaoc, double dthetaoc) {

	double term1 = 0.5*(wO_om - wO_open)*(1.0 + tanh((thetai - thetaoc)/dthetaoc));
	double term2 = 0.5*(1.0 - wO_om)*(1.0 + np.tanh((thetai - thetamm)/dthetamm));
	
	double wO = wO_open + term1 + term2;

	return wO;
}

double PedersenCurrent(	double thetai, double g, 
					double wO_open, double wO_om,
					double thetamm, double dthetamm,
					double thetaoc, double dthetaoc ) {

	/* height intergrated conductivity*/
	double SigmaP = 0.25;

	/* ionospheric radius (m)*/
	double Ri = 67350000.0;

	/* equatorial radius (m)*/
	double Rj = 71492000.0;
	
	/* magnetic field at thetai */
	double B = np.abs(2.0* g*cos(thetai)*pow(Rj/Ri,3.0))*1e-9;

	/* calculate rhoi */
	double rhoi = Ri*sin(thetai);

	/* get domega*/
	double wO = OmegaRatio(thetai,wO_open,wO_om,thetamm,dthetamm,thetaoc,dthetaoc);
	double OmegaJ = 1.758e-4;
	double domega = OmegaJ*(1 - wO);

	/* ftheta*/
	double ft = f_thetai(thetai);

	/* the current */
	double Ihp = 2.0*np.pi*SigmaP*rhoi*rhoi*domega*B*ft;

	return IhP;
}


double ThetaIonosphere(	double r, double theta, double g,
						double r0, double r1,
						double mui2, double D, 
						double deltarho, double deltaz) {
	
	/* calculate cylindrical coords*/
	double rho = r*sin(theta);
	double z = r*cos(theta);

	/* get the CAN flux */
	double Fcan = FluxCan(rho,z,r0,r1,mui2,D,deltarho,deltaz);
	
	/* dipole flux */
	double Fdip = FluxDip(r,theta,g);

	/* ionospheric radius (m)*/
	double Ri = 67350000.0;

	/* equatorial radius (m)*/
	double Rj = 71492000.0;

	/* theta ionosphere, yay! */
	double thetai = arcsin(sqrt((Ri/Rj)*(Fcan + Fdip)/g));	

	return thetai;
}

double BphiLMIC(double r, double theta, double g,
						double r0, double r1,
						double mui2, double D, 
						double deltarho, double deltaz,
						double wO_open, double wO_om,
						double thetamm, double dthetamm,
						double thetaoc, double dthetaoc ) {

	/* ionospheric latitude */
	double thetai = ThetaIonosphere(r,theta,g,r0,r1,mui2,D,deltarho,deltaz);

	/* sign of the latitude */
	double slat = sign(M_PI/2 - theta);

	/* Pedersen Current */
	double IhP = PedersenCurrent(thetai,g,wO_open,wO_om,thetamm,
								dthetamm,thetaoc,thetaoc,dthetaoc);

	/* constants */
	double Rj = 71492000.0;
	double mu0 = 4*M_PI*1e-7;

	/* rho (m)*/
	double rho = r*sin(theta)*Rj;

	/* calculate Bphi */
	double Bphi = (-sgn*mu0*IhP)/(2*M_PI*rho);

	return Bphi*1e9;

}