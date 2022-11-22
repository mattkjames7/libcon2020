#ifndef __LIBCON2020_H__
#define __LIBCON2020_H__
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

#define LIBCON2020_VERSION_MAJOR 0
#define LIBCON2020_VERSION_MINOR 0
#define LIBCON2020_VERSION_PATCH 1

#define _USE_MATH_DEFINES
#define deg2rad M_PI/180.0;

/* these wrappers can be used to get the magnetic field vectors */
void Con2020FieldArray(int n, double *p0, double *p1, double *p2,
					double *B0, double *B1, double *B2);
	
void Con2020Field(double p0, double p1, double p2,
			double *B0, double *B1, double *B2);


void GetCon2020Params(double *mui, double *irho, double *r0, double *r1,
					double *d, double *xt, double *xp, char *eqtype,
					bool *Edwards, bool *ErrChk, bool *CartIn, bool *CartOut, 
					bool *smooth, double *DeltaRho, double *DeltaZ,
					double *g, char *azfunc, double *wO_open, double *wO_oc,
					double *thetamm, double *dthetamm, double *thetaoc, double *dthetaoc);
						
	
void SetCon2020Params(double mui, double irho, double r0, double r1,
					double d, double xt, double xp, const char *eqtype,
					bool Edwards, bool ErrChk, bool CartIn, bool CartOut, 
					bool smooth, double DeltaRho, double DeltaZ,
					double g, const char *azfunc, double wO_open, double wO_oc,
					double thetamm, double dthetamm, double thetaoc, double dthetaoc);

void Con2020AnalyticField(	int n, double a, 
							double *rho, double *z, 
							double *Brho, double *Bz);

void Con2020AnalyticFieldSmooth(	int n, double a, 
							double *rho, double *z, 
							double *Brho, double *Bz);

#endif