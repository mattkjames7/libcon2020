#include "bessel.h"

namespace con2020 {


/* J0 polynomials as defined in the Cephes library */
static double PP0[7] = {
  7.96936729297347051624E-4,
  8.28352392107440799803E-2,
  1.23953371646414299388E0,
  5.44725003058768775090E0,
  8.74716500199817011941E0,
  5.30324038235394892183E0,
  9.99999999999999997821E-1,
};
static double PQ0[7] = {
  9.24408810558863637013E-4,
  8.56288474354474431428E-2,
  1.25352743901058953537E0,
  5.47097740330417105182E0,
  8.76190883237069594232E0,
  5.30605288235394617618E0,
  1.00000000000000000218E0,
};

static double QP0[8] = {
-1.13663838898469149931E-2,
-1.28252718670509318512E0,
-1.95539544257735972385E1,
-9.32060152123768231369E1,
-1.77681167980488050595E2,
-1.47077505154951170175E2,
-5.14105326766599330220E1,
-6.05014350600728481186E0,
};
static double QQ0[7] = {
  6.43178256118178023184E1,
  8.56430025976980587198E2,
  3.88240183605401609683E3,
  7.24046774195652478189E3,
  5.93072701187316984827E3,
  2.06209331660327847417E3,
  2.42005740240291393179E2,
};



static double DR1 = 5.78318596294678452118E0;
static double DR2 = 3.04712623436620863991E1;

static double RP0[4] = {
-4.79443220978201773821E9,
 1.95617491946556577543E12,
-2.49248344360967716204E14,
 9.70862251047306323952E15,
};
static double RQ0[8] = {
 4.99563147152651017219E2,
 1.73785401676374683123E5,
 4.84409658339962045305E7,
 1.11855537045356834862E10,
 2.11277520115489217587E12,
 3.10518229857422583814E14,
 3.18121955943204943306E16,
 1.71086294081043136091E18,
};

static double RP1[4] = {
-8.99971225705559398224E8,
 4.52228297998194034323E11,
-7.27494245221818276015E13,
 3.68295732863852883286E15,
};
static double RQ1[8] = {
/* 1.00000000000000000000E0,*/
 6.20836478118054335476E2,
 2.56987256757748830383E5,
 8.35146791431949253037E7,
 2.21511595479792499675E10,
 4.74914122079991414898E12,
 7.84369607876235854894E14,
 8.95222336184627338078E16,
 5.32278620332680085395E18,
};

static double PP1[7] = {
 7.62125616208173112003E-4,
 7.31397056940917570436E-2,
 1.12719608129684925192E0,
 5.11207951146807644818E0,
 8.42404590141772420927E0,
 5.21451598682361504063E0,
 1.00000000000000000254E0,
};
static double PQ1[7] = {
 5.71323128072548699714E-4,
 6.88455908754495404082E-2,
 1.10514232634061696926E0,
 5.07386386128601488557E0,
 8.39985554327604159757E0,
 5.20982848682361821619E0,
 9.99999999999999997461E-1,
};

static double QP1[8] = {
 5.10862594750176621635E-2,
 4.98213872951233449420E0,
 7.58238284132545283818E1,
 3.66779609360150777800E2,
 7.10856304998926107277E2,
 5.97489612400613639965E2,
 2.11688757100572135698E2,
 2.52070205858023719784E1,
};
static double QQ1[7] = {
/* 1.00000000000000000000E0,*/
 7.42373277035675149943E1,
 1.05644886038262816351E3,
 4.98641058337653607651E3,
 9.56231892404756170795E3,
 7.99704160447350683650E3,
 2.82619278517639096600E3,
 3.36093607810698293419E2,
};



static double Z1 = 1.46819706421238932572E1;
static double Z2 = 4.92184563216946036703E1;

/* pi related constants */
double PIO4   =  7.85398163397448309616E-1;    /* pi/4 */
double SQ2OPI =  7.9788456080286535587989E-1;  /* sqrt( 2/pi ) */
double TWOOPI =  6.36619772367581343075535E-1; /* 2/pi */
double THPIO4 =  2.35619449019234492885;

/***********************************************************************
 * NAME : j0(x)
 * 
 * DESCRIPTION : Fucntion to calculate an estimate of the Bessel 
 * function j0 using code based on the Cephes C library
 * (Cephes Mathematical Functions Library, http://www.netlib.org/cephes/).
 * 
 * INPUTS :
 * 		double 	x	position to calculate J0 at.
 * 
 * RETURNS :
 * 		double j	j0 function evaluated at x.
 * 
 * ********************************************************************/
double j0(double x) {
	double w, z, p, q, xn;

	if( x < 0 ) {
		x = -x;
	}
	
	if( x <= 5.0 ) {
		z = x * x;
		if( x < 1.0e-5 ) {
			return (1.0 - z/4.0);
		}
		p = (z - DR1) * (z - DR2);
		p = p * polyeval( z, RP0, 3)/pol1eval( z, RQ0, 8 );
		return p;
	}

	w = 5.0/x;
	q = 25.0/(x*x);
	p = polyeval( q, PP0, 6)/polyeval( q, PQ0, 6 );
	q = polyeval( q, QP0, 7)/pol1eval( q, QQ0, 7 );
	xn = x - PIO4;
	p = p * cos(xn) - w * q * sin(xn);
	return (p*SQ2OPI/sqrt(x));
}


/***********************************************************************
 * NAME : j0(n,x,j)
 * 
 * DESCRIPTION : Fucntion to calculate an estimate of the Bessel 
 * function j0 using code based on the Cephes C library
 * (Cephes Mathematical Functions Library, http://www.netlib.org/cephes/).
 * 
 * INPUTS :
 * 		int 	n	Number of elements in x
 * 		double 	*x	position to calculate J0 at.
 * 
 * OUTPUTS :
 * 		double *j	j0 function evaluated at x.
 * 
 * ********************************************************************/
void j0(int n, double *x, double *j) {
	
	int i;
	for (i=0;i<n;i++) {
		j[i] = j0(x[i]);
	}	
}

/***********************************************************************
 * NAME : j0(n,x,multx,j)
 * 
 * DESCRIPTION : Fucntion to calculate an estimate of the Bessel 
 * function j0 using code based on the Cephes C library
 * (Cephes Mathematical Functions Library, http://www.netlib.org/cephes/).
 * 
 * INPUTS :
 * 		int 	n	Number of elements in x
 * 		double 	*x	position to calculate J0(x*multx) at.
 * 		double multx	Constant to multiply x by
 * 
 * OUTPUTS :
 * 		double *j	j0 function evaluated at x*multx.
 * 
 * ********************************************************************/
void j0(int n, double *x, double multx, double *j) {
	int i;
	for (i=0;i<n;i++) {
		j[i] = j0(multx*x[i]);
	}
}

/***********************************************************************
 * NAME : j1(x)
 * 
 * DESCRIPTION : Fucntion to calculate an estimate of the Bessel 
 * function j1 using code based on the Cephes C library
 * (Cephes Mathematical Functions Library, http://www.netlib.org/cephes/).
 * 
 * INPUTS :
 * 		double 	x	position to calculate J1 at.
 * 
 * RETURNS :
 * 		double j	j1 function evaluated at x.
 * 
 * ********************************************************************/
double j1(double x) {
	double w, z, p, q, xn;

	w = x;
	if ( x < 0 ) {
		w = -x;
	}
	
	if( w <= 5.0 ) {
		z = x * x;	
		w = polyeval( z, RP1, 3 ) / pol1eval( z, RQ1, 8 );
		w = w * x * (z - Z1) * (z - Z2);
		return w ;
		}

	w = 5.0/x;
	z = w * w;
	p = polyeval( z, PP1, 6)/polyeval( z, PQ1, 6 );
	q = polyeval( z, QP1, 7)/pol1eval( z, QQ1, 7 );
	xn = x - THPIO4;
	p = p * cos(xn) - w * q * sin(xn);
	return (p * SQ2OPI / sqrt(x));
}

/***********************************************************************
 * NAME : j1(n,x,j)
 * 
 * DESCRIPTION : Fucntion to calculate an estimate of the Bessel 
 * function j1 using code based on the Cephes C library
 * (Cephes Mathematical Functions Library, http://www.netlib.org/cephes/).
 * 
 * INPUTS :
 * 		int 	n	Number of elements in x
 * 		double 	*x	position to calculate J1 at.
 * 
 * OUTPUTS :
 * 		double *j	j1 function evaluated at x.
 * 
 * ********************************************************************/
void j1(int n, double *x, double *j) {
	int i;
	for (i=0;i<n;i++) {
		j[i] = j1(x[i]);
	}
}

/***********************************************************************
 * NAME : j1(n,x,multx,j)
 * 
 * DESCRIPTION : Fucntion to calculate an estimate of the Bessel 
 * function j1 using code based on the Cephes C library
 * (Cephes Mathematical Functions Library, http://www.netlib.org/cephes/).
 * 
 * INPUTS :
 * 		int 	n	Number of elements in x
 * 		double 	*x	position to calculate J1(x*multx) at.
 * 		double multx	Constant to multiply x by
 * 
 * OUTPUTS :
 * 		double *j	j1 function evaluated at x*multx.
 * 
 * ********************************************************************/
void j1(int n, double *x, double multx, double *j) {
	
	int i;
	for (i=0;i<n;i++) {
		j[i] = j1(multx*x[i]);
	}
}


/***********************************************************************
 * Bessel functions below are from 
 * 
 * https://doi.org/10.1007/s40314-020-01238-z
 * 
 * and 
 * 
 * https://doi.org/10.1007/s40314-020-01238-z
 * 
***********************************************************************/



void j0m(int n, double *x, double *j) {
	
	/* define the parameters for this approx. */
	double lambda = 0.865;
	double q = 0.7172491568;
	double p0 = 0.6312725339;
	double pt0 = 0.4308049446;
	double p1 = 0.3500347951;
	double pt1 = 0.4678202347;
	double p2 = -0.06207747907;
	double pt2 = 0.04253832927;
	double l4 = lambda*lambda;
	l4 = l4*l4;
	
	/* some temporary variables */
	double x2, A, B, C, l4x2p1, rtl4x2p1, rtrtl4x2p1;
	
	/* loop through each point */
	int i;
	for (i=0;i<n;i++) {
		/* these temporary variables are used in multiple places so
		 * defining them once should hopefully speed things up a bit */
		x2 = x[i]*x[i];
		l4x2p1 = 1.0 + l4*x2;	
		rtl4x2p1 = sqrt(l4x2p1);
		rtrtl4x2p1 = sqrt(rtl4x2p1);
		
		A = 1.0/(rtrtl4x2p1*(1.0 + q*x2));
		B = p0 + p1*x2 + p2*rtl4x2p1;
		C = (pt0 + pt1*x2)*rtl4x2p1 + pt2*x2;
		
		j[i] = A*(B*cos(x[i]) + C*sin(x[i])/x[i]);
	}
		
}

void j0m(int n, double *x, double multx, double *j) {
	
	/* define the parameters for this approx. */
	double lambda = 0.865;
	double q = 0.7172491568;
	double p0 = 0.6312725339;
	double pt0 = 0.4308049446;
	double p1 = 0.3500347951;
	double pt1 = 0.4678202347;
	double p2 = -0.06207747907;
	double pt2 = 0.04253832927;
	double l4 = lambda*lambda;
	l4 = l4*l4;
	
	/* some temporary variables */
	double mx, x2, A, B, C, l4x2p1, rtl4x2p1, rtrtl4x2p1;
	
	/* loop through each point */
	int i;
	for (i=0;i<n;i++) {
		/* these temporary variables are used in multiple places so
		 * defining them once should hopefully speed things up a bit */
		mx = x[i]*multx;
		x2 = mx*mx;
		l4x2p1 = 1.0 + l4*x2;	
		rtl4x2p1 = sqrt(l4x2p1);
		rtrtl4x2p1 = sqrt(rtl4x2p1);
		
		A = 1.0/(rtrtl4x2p1*(1.0 + q*x2));
		B = p0 + p1*x2 + p2*rtl4x2p1;
		C = (pt0 + pt1*x2)*rtl4x2p1 + pt2*x2;
		
		j[i] = A*(B*cos(mx) + C*sin(mx)/mx);
	}
		
}

void j1m(int n, double *x, double *j) {
	
	/* parameters from the paper */
	double lambda = 0.1;
	double q1 = 0.4120981204;
	double q2 = 0.006571619275;
	double P0 = -0.7763224930;
	double p0 = 1.776322448;
	double P1 = -0.03147133771;
	double p1 = 0.2250803518;
	double P2 = -(2*pow(lambda,1.5))*q2/sqrt(M_PI);
	double p2 = (2*sqrt(lambda))*q2/sqrt(M_PI);
	double l2 = lambda*lambda;
	
	/* temporary variables */
	double x2, x4, A, B, C, ac, rtac, bc;
	
	/* loop through each value */
	int i;
	for (i=0;i<n;i++) {
		x2 = x[i]*x[i];
		x4 = x2*x2; 
		ac = sqrt(1 + l2*x2);
		rtac = sqrt(ac);
		A = 1.0/(2.0*rtac);
		bc = (1 + q1*x2 + q2*x4);
		B = (p0 + p1*x2 + p2*x4)/bc;
		C = (x[i]/ac)*(P0 + P1*x2 + P2*x4)/bc;
		j[i] = A*(B*sin(x[i]) + C*cos(x[i]));
	}

}

void j1m(int n, double *x, double multx, double *j) {
	
	/* parameters from the paper */
	double lambda = 0.1;
	double q1 = 0.4120981204;
	double q2 = 0.006571619275;
	double P0 = -0.7763224930;
	double p0 = 1.776322448;
	double P1 = -0.03147133771;
	double p1 = 0.2250803518;
	double P2 = -(2*pow(lambda,1.5))*q2/sqrt(M_PI);
	double p2 = (2*sqrt(lambda))*q2/sqrt(M_PI);
	double l2 = lambda*lambda;
	
	/* temporary variables */
	double mx, x2, x4, A, B, C, ac, rtac, bc;
	
	/* loop through each value */
	int i;
	for (i=0;i<n;i++) {
		mx = multx*x[i];
		x2 = mx*mx;
		x4 = x2*x2; 
		ac = sqrt(1 + l2*x2);
		rtac = sqrt(ac);
		A = 1.0/(2.0*rtac);
		bc = (1 + q1*x2 + q2*x4);
		B = (p0 + p1*x2 + p2*x4)/bc;
		C = (mx/ac)*(P0 + P1*x2 + P2*x4)/bc;
		j[i] = A*(B*sin(mx) + C*cos(mx));
	}

}


} // end namespace con2020	
