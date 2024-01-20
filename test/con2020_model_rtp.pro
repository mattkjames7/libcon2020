FUNCTION con2020_model_rtp, eq_type, r_rj, colat_rads, elong_rads, use_these_params
  ;% ======
  ;% CON2020_MODEL_RTP (Spherical)
  ;% ======
  ;% Code to calculate the perturbation magnetic field produced by the Connerney et al. 1981 (CAN) current sheet, which is
  ;% represented by a finite disk of current.
  ;%  This disk has variable parameters including (among others) the current density, and current sheet inner edge, outer
  ;%   edge and thickness.
  ;%  The disk is centered on the magnetic equator (shifted in longitude and tilted as specified by model parameters xp__cs_rhs_azimuthal_angle_of_tilt_degs and xt__cs_tilt_degs)
  ;%  This 2020 version includes a radial current per Connerney et al. (2020),
  ;%   https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2020JA028138
  ;%  For more details about the model and the development of this code please see the PDF at 
  ;%   https://github.com/marissav06/con2020_idl/blob/main/con2020_final_code_documentation_sept13_2021.pdf
  ;%
  ;% Use in one of the following ways:
  ;%  Use default current sheet model parameter structure:  B = con2020_model_rtp(eq_type, r_rj, colat_rads, elong_rads)
  ;%  Use your own current sheet model parameter structure: B = con2020_model_rtp(eq_type, r_rj, colat_rads, elong_rads, use_these_params)
  ;%  Obtain the default model parameters:             params = con2020_model_rtp('default_values')
  ;%  Then you can edit the structure and use, e.g. params.r1__outer_rj = 51.400001,
  ;%    then B = con2020_model_rtp(eq_type, r_rj, colat_rads, elong_rads, params)
  ;%
  ;% Required inputs (System III Spherical, right handed):
  ;%  eq_type - equation type: 'integral', 'analytic' or 'hybrid',
  ;%   or set to 'default_values' to return a structure of all default values.
  ;%  r_rj       - radial distance, in Rj.                    Value(s) must be 0 <      r_rj   <  200.
  ;%  colat_rads - colatitude, in radians.                    Value(s) must be 0 <= colat_rads <=  pi.
  ;%  elong_rads - East longitude, right handed, in radians.  Value(s) must be 0 <= elong_rads <= 2pi.
  ;% r_rj, colat_rads and elong_rads can be scalars or 1D arrays (nx1), but only one eq_type.
  ;%
  ;% Unless an option structure is provided it will default to parameters from Connerney et al., 2020.
  ;% Optional input of a structure: use_these_params
  ;% with the structure fields:
  ;%  use_these_params.mu_i_div2__current_parameter_nT         - mu0i0/2 term (current sheet field parameter), in nT
  ;%  use_these_params.i_rho__radial_current_MA                - radial current term from Connerney et al., 2020 (set this to zero to turn radial currents off as in Connerney et al. 1981)
  ;%  use_these_params.r0__inner_rj                            - inner edge of current disk in Rj
  ;%  use_these_params.r1__outer_rj                            - outer edge of current disk in Rj
  ;%  use_these_params.d__cs_half_thickness_rj                 - current sheet half thickness in Rj
  ;%  use_these_params.xt__cs_tilt_degs                        - current sheet tilt in degrees
  ;%  use_these_params.xp__cs_rhs_azimuthal_angle_of_tilt_degs - current sheet longitude (right handed) in degrees
  ;%  use_these_params.error_check                             - 1 to check that inputs are valid (Default),
  ;%                                                             or set to 0 to skip input checks (faster).
  ;%
  ;% To retrieve the default value parameters, run
  ;% params = con2020_model_rtp('default_values')
  ;% then you may edit the values in the structure, then use as B = con2020_model_rtp(eq_type, r_rj, colat_rads, elong_rads, params)
  ;%
  ;% Outputs:
  ;%  B - Spherical Magnetic field vector from current sheet model, [Br, Btheta, Bphi], units of nT.
  ;%
  ;% This code can take a hybrid approach to calculating the current sheet field, using the integral equations in some regions
  ;% and the analytic equations in others.
  ;% Following Connerney et al. 1981, figure A1, and Edwards et al. (2001), figure 2, the choice of integral vs. analytic
  ;% equations is most important near rho = r0 and z = 0.
  ;% By default, this hybrid method uses the analytic equations everywhere except |Z| < D*1.5 and |Rho-R0| < 2
  ;%    Analytic equations:
  ;%        For the analytic equations, we use the equations provided by Edwards et al. 2001:
  ;%         https://doi.org/10.1016/S0032-0633(00)00164-1
  ;%        Other analytic approximations to the CAN sheet equations are provided in Connerney et al., 1981
  ;%         https://doi.org/10.1029/JA086iA10p08370
  ;%    Integral equations:
  ;%        For the integral equations we use the Bessel functions from Connerney et al. 1981, eqs. 14, 15, 17, 18
  ;%        We do not integrate lambda from zero to infinity, but vary the integration limit depending on the value of the
  ;%        Bessel functions.
  ;%        For computational speed, the inner disk edge calculations use the integral equations,
  ;%        but the outer disk edge calculations use the anlytical equations.
  ;%
  ;% Updates:
  ;% by Marissa Vogt (mvogt@bu.edu), March 2021,
  ;% RJ Wilson did some speedups and re-formatting of lines, also March 2021
  ;%
  ;% Converted to MATLAB by Marty Brennan, June 2021
  ;% RJ Wilson did some reformatting, June 2021, and added
  ;% int_tabulated_rjw2_sub as a subfunction, rather than separate file int_tabulated_rjw2.m
  ;% which was then replaced by some in-line code instead of calling the subfunction
  ;% RJ Wilson split initial Matlab and IDL code in to Cartesian and a Spherical wrapper code and updated this help text,
  ;% in August 2021, to make con2020_model_xyz and con2020_model_rtp.
  ;% RJ Wilson renamed i_rho__radial_current_density_nT to i_rho__radial_current_intensity_MA in June 2022.
  ;% RJ Wilson renamed i_rho__radial_current_intensity_MA to i_rho__radial_current_MA and mu_i_div2__current_density_nT to mu_i_div2__current_parameter_nT in November 2022.
  ;% RJ Wilson put in a fix to prevent Infinities/NaNs in Bphi1 when rho1 = 0, and updated the Integral equation comment above to note that the outer edge always uses analytical equations, even if integral is chosen.

  ON_ERROR, 2 ; % Exit code if an error in main, don't stop in code - no Matlab equivalent, just delete line in Matlab
  FORWARD_FUNCTION con2020_model_xyz ; telling IDL this is a function, in case con2020_model_xyz has not been compiled yet.

  IF STRCMP(eq_type,'default_values',FOLD_CASE = 1) THEN RETURN, con2020_model_xyz(eq_type) ; case insensitive. % Not yet checked if eq_type is a character string, IDL doesn't care

  ;% Covert to Doubles, and rename input variables so as not to over-write them (only an IDL issue)
  r_in  = DOUBLE(r_rj      ); % units Rj
  theta = DOUBLE(colat_rads); % SYSIII Colat (radians)
  phi   = DOUBLE(elong_rads); % SYSIII ELong or Azimuth (radians)

  is_use_these_params = KEYWORD_SET(use_these_params)

  ;% Very basis error check that theta and phi are probably in radians
  WHILE 1 DO BEGIN
    IF (is_use_these_params EQ 1) THEN BEGIN
      IF use_these_params.error_check EQ 1 THEN BREAK ;% skip error check that values are probably radians.
    ENDIF
    ;% Not checking value of r_rj here.  Cartesian Code will check x, y and z are all within +/-200.
    IF N_ELEMENTS(theta) EQ 1 THEN BEGIN
      IF   (theta LT 0d) OR (   theta GT    !DPI) THEN MESSAGE,'Input colat_rads must be in radians and from 0 to  pi only!'
      IF   (phi   LT 0d) OR (   phi   GT 2d*!DPI) THEN MESSAGE,'Input elong_rads must be in radians and from 0 to 2pi only!'
      BREAK
    ENDIF
    ;% assume vector, not scalar
    min_theta = MIN(theta, MAX=max_theta)
    min_phi   = MIN(phi  , MAX=max_phi  )
    IF (min_theta LT 0) OR (max_theta GT    !DPI) THEN MESSAGE,'Input colat_rads must be in radians and from 0 to  pi only!'
    IF (min_phi   LT 0) OR (max_phi   GT 2d*!DPI) THEN MESSAGE,'Input elong_rads must be in radians and from 0 to 2pi only!'
    BREAK ;% escape WHILE loop!
  ENDWHILE

  ;% Convert to cartesian coordinates and rotate into magnetic longitude
  ;% (x,y,z) are the shifted (phi) coordinates
  sin_theta = sin(theta)
  cos_theta = cos(theta)
  sin_phi   = sin(phi)
  cos_phi   = cos(phi)

  x = r_in * sin_theta * cos_phi
  y = r_in * sin_theta * sin_phi
  z = r_in * cos_theta

  ;% Do calculation in cartesian
  IF (is_use_these_params EQ 0) THEN $
    bxyz = con2020_model_xyz(eq_type, x, y, z) $
  ELSE $
    bxyz = con2020_model_xyz(eq_type, x, y, z, use_these_params) ; don't want to pass use_these_params if undefined

  ;% Convert to spherical coordinates and return
  RETURN, [$
    [ bxyz[*,0]*cos_phi*sin_theta + bxyz[*,1]*sin_phi*sin_theta + bxyz[*,2]*cos_theta], $; br
    [ bxyz[*,0]*cos_phi*cos_theta + bxyz[*,1]*sin_phi*cos_theta - bxyz[*,2]*sin_theta], $; bt
    [-bxyz[*,0]*sin_phi           + bxyz[*,1]*cos_phi                                ]  $; bp
    ]; % size n x 3
END
