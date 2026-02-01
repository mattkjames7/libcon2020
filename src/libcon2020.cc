#include "libcon2020.h"



namespace con2020 {

Con2020 con2020inst; /* default model object */

void Con2020FieldArray(int n, double *p0, double *p1, double *p2,
			double *B0, double *B1, double *B2) {

	/* could create a separate function for default model */
	con2020inst.Field(n,p0,p1,p2,B0,B1,B2);
						
}

void Con2020Field(double p0, double p1, double p2,
			double *B0, double *B1, double *B2) {

	/* could create a separate function for default model */
	con2020inst.Field(p0,p1,p2,B0,B1,B2);
						
}

void GetCon2020Params(double *mui, double *irho, double *r0, double *r1,
				double *d, double *xt, double *xp, char *eqtype,
				bool *Edwards, bool *ErrChk, bool *CartIn, bool *CartOut, 
				bool  *smooth, double *DeltaRho, double *DeltaZ,
				double *g, char *azfunc, double *wO_open, double *wO_om,
				double *thetamm, double *dthetamm, double *thetaoc, double *dthetaoc) {

	mui[0] = con2020inst.GetAzCurrentParameter();
	irho[0] = con2020inst.GetRadCurrentParameter();
	r0[0] = con2020inst.GetR0();
	r1[0] = con2020inst.GetR1();
	d[0] = con2020inst.GetCSHalfThickness();
	xt[0] = con2020inst.GetCSTilt();
	xp[0] = con2020inst.GetCSTiltAzimuth();
	Edwards[0] = con2020inst.GetEdwardsEqs();
	ErrChk[0] = con2020inst.GetErrCheck();
	CartIn[0] = con2020inst.GetCartIn();
	CartOut[0] = con2020inst.GetCartOut();
	con2020inst.GetEqType(eqtype);	
	smooth[0] = con2020inst.GetSmooth();
	DeltaRho[0] = con2020inst.GetDeltaRho();
	DeltaZ[0] = con2020inst.GetDeltaZ();
	/* new LMIC parameters */
	g[0] = con2020inst.GetG();
	con2020inst.GetAzimuthalFunc(azfunc);
	wO_open[0] = con2020inst.GetOmegaOpen();
	wO_om[0] = con2020inst.GetOmegaOM();
	thetamm[0] = con2020inst.GetThetaMM();
	dthetamm[0] = con2020inst.GetdThetaMM();
	thetaoc[0] = con2020inst.GetThetaOC();
	dthetaoc[0] = con2020inst.GetdThetaOC();

	
}
void SetCon2020Params(double mui, double irho, double r0, double r1,
				double d, double xt, double xp, const char *eqtype,
				bool Edwards, bool ErrChk, bool CartIn, bool CartOut, 
				bool smooth, double DeltaRho, double DeltaZ,
				double g, const char *azfunc, double wO_open, double wO_om,
				double thetamm, double dthetamm, double thetaoc, double dthetaoc) {

	con2020inst.SetAzCurrentParameter(mui);
	con2020inst.SetRadCurrentParameter(irho);
	con2020inst.SetR0(r0);
	con2020inst.SetR1(r1);
	con2020inst.SetCSHalfThickness(d);
	con2020inst.SetCSTilt(xt);
	con2020inst.SetCSTiltAzimuth(xp);
	con2020inst.SetEdwardsEqs(Edwards);
	con2020inst.SetErrCheck(ErrChk);
	con2020inst.SetCartIn(CartIn);
	con2020inst.SetCartOut(CartOut);
	con2020inst.SetEqType(eqtype);
	con2020inst.SetSmooth(smooth);
	con2020inst.SetDeltaRho(DeltaRho);
	con2020inst.SetDeltaZ(DeltaZ);

	/*set LMIC parameters */
	con2020inst.SetG(g);
	con2020inst.SetAzimuthalFunc(azfunc);
	con2020inst.SetOmegaOpen(wO_open);
	con2020inst.SetOmegaOM(wO_om);
	con2020inst.SetThetaMM(thetamm);
	con2020inst.SetdThetaMM(dthetamm);
	con2020inst.SetThetaOC(thetaoc);
	con2020inst.SetdThetaOC(dthetaoc);
}

void Con2020AnalyticField(	int n, double a, 
							double *rho, double *z, 
							double *Brho, double *Bz) {

	/* define a few required variables */
	int i;

	for (i=0;i<n;i++) {
		con2020inst.AnalyticField(a,rho[i],z[i],&Brho[i],&Bz[i]);
	}

}


void Con2020AnalyticFieldSmooth(	int n, double a, 
							double *rho, double *z, 
							double *Brho, double *Bz) {

	/* define a few required variables */
	int i;

	for (i=0;i<n;i++) {
		con2020inst.AnalyticFieldSmooth(a,rho[i],z[i],&Brho[i],&Bz[i]);
	}

}
} //namespace con2020

extern "C" {
	
void Con2020FieldArray(
	int n, double *p0, double *p1, double *p2,
	double *B0, double *B1, double *B2) {
	con2020::Con2020FieldArray(n, p0, p1, p2, B0, B1, B2);
}

void Con2020Field(
	double p0, double p1, double p2,
	double *B0, double *B1, double *B2) {
	con2020::Con2020Field( p0, p1, p2, B0, B1, B2);
}

void GetCon2020Params(
	double *mui, double *irho, double *r0, double *r1,
	double *d, double *xt, double *xp, char *eqtype,
	bool *Edwards, bool *ErrChk, bool *CartIn, bool *CartOut, 
	bool  *smooth, double *DeltaRho, double *DeltaZ,
	double *g, char *azfunc, double *wO_open, double *wO_om,
	double *thetamm, double *dthetamm, double *thetaoc, double *dthetaoc) {
	con2020::GetCon2020Params( mui, irho, r0, r1,
		 d, xt, xp, eqtype,
		 Edwards, ErrChk, CartIn, CartOut, 
		 smooth, DeltaRho, DeltaZ,
		 g, azfunc, wO_open, wO_om,
		 thetamm, dthetamm, thetaoc, dthetaoc);
}

void SetCon2020Params(
	double mui, double irho, double r0, double r1,
	double d, double xt, double xp, const char *eqtype,
	bool Edwards, bool ErrChk, bool CartIn, bool CartOut, 
	bool smooth, double DeltaRho, double DeltaZ,
	double g, const char *azfunc, double wO_open, double wO_om,
	double thetamm, double dthetamm, double thetaoc, double dthetaoc) {
	con2020::SetCon2020Params( mui, irho, r0, r1,
		 d, xt, xp, eqtype,
		 Edwards, ErrChk, CartIn, CartOut, 
		 smooth, DeltaRho, DeltaZ,
		 g, azfunc, wO_open, wO_om,
		 thetamm, dthetamm, thetaoc, dthetaoc);
}

void Con2020AnalyticField(
	int n, double a, 
	double *rho, double *z, 
	double *Brho, double *Bz) {
	con2020::Con2020AnalyticField(
		n, a, rho, z, Brho, Bz);
}

void Con2020AnalyticFieldSmooth(
	int n, double a, 
	double *rho, double *z, 
	double *Brho, double *Bz) {
	con2020::Con2020AnalyticFieldSmooth(
		n, a, rho, z, Brho, Bz);
}

} // extern "C"