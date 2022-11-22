#ifndef __LIBCON2020_H__
#define __LIBCON2020_H__
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

#define LIBCON2020_VERSION_MAJOR 0
#define LIBCON2020_VERSION_MINOR 1
#define LIBCON2020_VERSION_PATCH 0

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


/***************************************************************
*
*   NAME : ScalarPotentialSmallRho(rho,z,a,mui2,D)
*
*   DESCRIPTION : Calcualte the small rho approximation
*       of the scalar potential accoring to Edwards et al.,
*       2001 (equation 8).
*
*   INPUTS : 
*       double  rho     Cylindrical rho coordinate (in disc 
*                       coordinate system, Rj)
*       double  z       z-coordinate, Rj
*       double  a       inner edge of semi-infinite current 
*                       sheet, Rj
*       double mui2     mu_0 I_0 /2 parameter (default 139.6 nT)
*       double D        Current sheet half-thickness, Rj
*
***************************************************************/
	double ScalarPotentialSmallRho( double rho, double z, double a,
									double mui2, double D);


/***************************************************************
*
*   NAME : ScalarPotentialLargeRho(rho,z,a,mui2,D,deltaz)
*
*   DESCRIPTION : Calcualte the large rho approximation
*       of the scalar potential accoring to Edwards et al.,
*       2001 (equation 12).
*
*   INPUTS : 
*       double  rho     Cylindrical rho coordinate (in disc 
*                       coordinate system, Rj)
*       double  z       z-coordinate, Rj
*       double  a       inner edge of semi-infinite current 
*                       sheet, Rj
*       double mui2     mu_0 I_0 /2 parameter (default 139.6 nT)
*       double D        Current sheet half-thickness, Rj
*       double deltaz   Scale length over which to smooth 4th
*                        term of the equation
*
***************************************************************/
	double ScalarPotentialLargeRho( double rho, double z, double a,
									double mui2, double D, double deltaz);


/***************************************************************
*
*   NAME : ScalarPotential(rho,z,a,mui2,D,deltarho,deltaz)
*
*   DESCRIPTION : Calculate the small/large rho approximation
*       of the scalar potential accoring to Edwards et al.,
*       2001 (equations 8 & 12).
*
*   INPUTS : 
*       double  rho     Cylindrical rho coordinate (in disc 
*                       coordinate system, Rj)
*       double  z       z-coordinate, Rj
*       double  a       inner edge of semi-infinite current 
*                       sheet, Rj
*       double mui2     mu_0 I_0 /2 parameter (default 139.6 nT)
*       double D        Current sheet half-thickness, Rj
*       double deltarho Scale length to smoothly transition from
*                       small to large rho approx
*       double deltaz   Scale length over which to smooth 4th
*                        term of the equation
*
***************************************************************/
	double ScalarPotential( double rho, double z, double a,
							double mui2, double D, 
							double deltarho, double deltaz);

	/*************************************************************
	*
	*	NAME: f_theta(thetai)
	*
	*	DESCRIPTION: Equation 5 of Cowley et al., 2008
	*
	*	INPUTS:
	*		double thetai	colatitude of the ionospheric footprint
	*						in radians!
	*
	*	RETURNS:
	*		double f_theta	1 + 0.25*tan^2 thetai
	*
	*************************************************************/
	double f_thetai(double thetai);

	/*************************************************************
	*
	*	NAME: OmegaRatio(thetai,wO_open,wO_om,thetamm,dthetamm,
	*						thetaoc,dthetaoc)
	*
	*	DESCRIPTION: Ratio of the angular velocity mapped to
	*		thetai to the planetary rotation. Equation 15 of 
	*		Cowley et al., 2008.
	*
	*	INPUTS:
	*		double thetai	colatitude of the ionospheric footprint
	*						in radians!
	*		double wO_open	angular velocity ratio of open flux to
	*						planetary spin
	*		double wO_om	angular velocity ratio of outer magnetosphere
	*						to planetary spin
	*		double thetamm	ionospheric footprint latitude of the 
	*						middle magnetosphere (where plasma 
	*						goes from rigid corotation to subcorotation)
	*						in radians.
	*		double dthetamm	width of the middle magnetosphere in radians.
	*		double thetaoc	ionospheric latitude of the open-closed field
	*						line boundary, in radians.
	*		double dthetaoc	width of the open-closed field line boundary,
	*						in radians.
	*
	*	RETURNS:
	*		double wO		Ratio of plasma angular veloctiy to Jupiter
	*						spin.
	*
	*************************************************************/
	double OmegaRatio(	double thetai, double wO_open, double wO_om,
						double thetamm, double dthetamm,
						double thetaoc, double dthetaoc);

	/*************************************************************
	*
	*	NAME: PedersenCurrent(thetai,g,wO_open,wO_om,thetamm,dthetamm,
	*						thetaoc,dthetsoc)
	*
	*	DESCRIPTION: Calculate the Pedersen current which maps to a
	*		given ionospheric latitude using equation 6 of Cowley et
	*		al., 2008.
	*
	*	INPUTS:
	*		double thetai	colatitude of the ionospheric footprint
	*						in radians!
	*		double g		dipole coefficient, nT.
	*		double wO_open	angular velocity ratio of open flux to
	*						planetary spin
	*		double wO_om	angular velocity ratio of outer magnetosphere
	*						to planetary spin
	*		double thetamm	ionospheric footprint latitude of the 
	*						middle magnetosphere (where plasma 
	*						goes from rigid corotation to subcorotation)
	*						in radians.
	*		double dthetamm	width of the middle magnetosphere in radians.
	*		double thetaoc	ionospheric latitude of the open-closed field
	*						line boundary, in radians.
	*		double dthetaoc	width of the open-closed field line boundary,
	*						in radians.
	*	RETURNS:
	*		double Ihp		Ionospheric Pedersen current.
	*
	*************************************************************/
	double PedersenCurrent(	double thetai, double g, 
						double wO_open, double wO_om,
						double thetamm, double dthetamm,
						double thetaoc, double dthetaoc );				

	/*************************************************************
	*
	*	NAME: ThetaIonosphere(r,theta,g,r0,r1,mui2,D,deltarho,deltaz)
	*
	*	DESCRIPTION: Use the flux functions of the CAN model and a 
	*		dipole field to map the current position to a position
	*		on the ionosphere.
	*
	*	INPUTS:
	*		double r		radial coordinate, Rj.
	*		double theta	colatitude, radians.
	*		double g		dipole coefficient, nT.
	*		double r0		Inner edge of the current sheet, Rj.
	*		double r1		Outer edge of the current sheet, Rj.
	*		double mui2		current parameter, nT.
	*		double D		half-thickness of the current sheet, Rj.
	*		double deltarho	scale distance of the smoothing between
	*						inner and outer approximations, Rj.
	*		double deltaz	scale distance to smooth across the
	*						+/-D boundary, Rj.
	*
	*	RETURNS:
	*		double thetai	Ionospheric latitude in radians.
	*
	*
	*************************************************************/
	double ThetaIonosphere(	double r, double theta, double g,
							double r0, double r1,
							double mui2, double D, 
							double deltarho, double deltaz);

	/*************************************************************
	*
	*	NAME: BphiLMIC(r,theta,g,r0,r1,mui2,D,deltarho,deltaz,
	*					wO_open,wO_om,thetamm,dthetamm,
	*					thetaom,dthetaom)
	*
	*	DESCRIPTION: Calculate the azimuthal field using the LMIC 
	*		model.
	*
	*	INPUTS:
	*		double r		radial coordinate, Rj.
	*		double theta	colatitude, radians.
	*		double g		dipole coefficient, nT.
	*		double r0		Inner edge of the current sheet, Rj.
	*		double r1		Outer edge of the current sheet, Rj.
	*		double mui2		current parameter, nT.
	*		double D		half-thickness of the current sheet, Rj.
	*		double deltarho	scale distance of the smoothing between
	*						inner and outer approximations, Rj.
	*		double deltaz	scale distance to smooth across the
	*						+/-D boundary, Rj.
	*		double wO_open	angular velocity ratio of open flux to
	*						planetary spin
	*		double wO_om	angular velocity ratio of outer magnetosphere
	*						to planetary spin
	*		double thetamm	ionospheric footprint latitude of the 
	*						middle magnetosphere (where plasma 
	*						goes from rigid corotation to subcorotation)
	*						in radians.
	*		double dthetamm	width of the middle magnetosphere in radians.
	*		double thetaoc	ionospheric latitude of the open-closed field
	*						line boundary, in radians.
	*		double dthetaoc	width of the open-closed field line boundary,
	*						in radians.
	*
	*	RETURNS:
	*		double Bphi		Azimuthal field, nT.
	*
	*************************************************************/
	double BphiLMIC(double r, double theta, double g,
							double r0, double r1,
							double mui2, double D, 
							double deltarho, double deltaz,
							double wO_open, double wO_om,
							double thetamm, double dthetamm,
							double thetaoc, double dthetaoc );

	/*************************************************************
	*
	*	NAME: BphiIonosphere(thetai,g,wO_open,wO_om,thetamm,dthetamm,
	*					thetaom,dthetaom)
	*
	*	DESCRIPTION: Calculate the ionospheric azimuthal field using the LMIC 
	*		model.
	*
	*	INPUTS:
	*		double thetai	ionospheric colatitude, radians.
	*		double g		dipole coefficient, nT.
	*		double wO_open	angular velocity ratio of open flux to
	*						planetary spin
	*		double wO_om	angular velocity ratio of outer magnetosphere
	*						to planetary spin
	*		double thetamm	ionospheric footprint latitude of the 
	*						middle magnetosphere boundary (where plasma 
	*						goes from rigid corotation to subcorotation)
	*						in radians.
	*		double dthetamm	width of the middle magnetosphere boundary
	*						in radians.
	*		double thetaoc	ionospheric latitude of the open-closed field
	*						line boundary, in radians.
	*		double dthetaoc	width of the open-closed field line boundary,
	*						in radians.
	*
	*	RETURNS:
	*		double Bphi		Azimuthal field, nT.
	*
	*************************************************************/
	double BphiIonosphere( 	double thetai, double g,
							double wO_open, double wO_om,
							double thetamm, double dthetamm,
							double thetaoc, double dthetaoc );


#endif