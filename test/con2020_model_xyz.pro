  ;% ======
  ;% CON2020_MODEL_XYZ (Cartesian)
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
  ;%  Use default current sheet model parameter structure:  B = con2020_model_xyz(eq_type, x_rj, y_rj, z_rj)
  ;%  Use your own current sheet model parameter structure: B = con2020_model_xyz(eq_type, x_rj, y_rj, z_rj, use_these_params)
  ;%  Obtain the default model parameters:             params = con2020_model_xyz('default_values')
  ;%  Then you can edit the structure and use, e.g. params.r1__outer_rj = 51.400001,
  ;%    then B = con2020_model_xyz(eq_type, x_rj, y_rj, z_rj, params)
  ;%
  ;% Required inputs (System III Cartesian, right handed; see con2020_model_rtp.pro for spherical version):
  ;%  eq_type - equation type: 'integral', 'analytic' or 'hybrid',
  ;%   or set to 'default_values' to return a structure of all default values.
  ;%  x_rj      - SYSIII x position, in Rj, Values must be -200 < x_rj < 200.
  ;%  y_rj      - SYSIII y position, in Rj, Values must be -200 < x_rj < 200.
  ;%  z_rj      - SYSIII z position, in Rj, Values must be -200 < x_rj < 200.
  ;% x_rj, y_rj and z_rj can be scalars or 1D arrays (nx1), but only one eq_type.
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
  ;% params = con2020_model_xyz('default_values')
  ;% then you may edit the values in the structure, then use as B = con2020_model_xyz(eq_type, x_rj, y_rj, z_rj, params)
  ;%
  ;% Outputs:
  ;%  B - Cartesian Magnetic field vector from current sheet model, [Bx, By, Bz], units of nT.
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
  ;% by Marissa Vogt, March 2021,
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
  ;% RJ Wilson added the _con2020_model_xyz_params_to_string function to the end of the pro file.  May be useful for Field Line Tracing codes at a later date, but otherwise not called.

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; Function con2020_model_xyz begins around line 250. Sub-functions are listed first (in order to compile appropriately in IDL).
  ; Function _con2020_model_xyz_params_to_string begins around line 190.
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

FUNCTION _con2020_model_xyz_analytic, rho1, z1, rho1_sq, d__cs_half_thickness_rj, r, mu_i_div2__current_parameter_nT, scalar_input
  COMPILE_OPT HIDDEN
  ON_ERROR, 2 ;Return to caller if an error occurs.

  ;% Analytic equations
  ;% Connerney et al. 1981's equations for the field produced by a semi-infinite disk of thickness D, inner edge R0, outer edge R1 -
  ;%  see their equations A1 through A9
  ;% the analytic equations for Brho and Bz vary depending on the region with respect to the current disk

  ;% Doing these 3 equations on the whole array to save getting confused by indices,  Will swap to just required indices later
  z1md = z1-d__cs_half_thickness_rj; % Matt's zmd
  z1pd = z1+d__cs_half_thickness_rj; % Matt's zpd

  r_sq = r*r; % Matt's a2

  IF scalar_input THEN BEGIN
    IF rho1 LT r THEN BEGIN
      ind_LT   = 0
      n_ind_LT = 1
      ;% ind_GE = []
      n_ind_GE = 0
    ENDIF ELSE BEGIN
      ;% ind_LT = []
      n_ind_LT = 0
      ind_GE   = 0
      n_ind_GE = 1
    ENDELSE
    brho1 = !VALUES.D_NaN
    bz1   = !VALUES.D_NaN
  ENDIF ELSE BEGIN
    ind_LT = WHERE(rho1 LT r, n_ind_LT, COMPLEMENT = ind_GE, NCOMPLEMENT = n_ind_GE, NULL = 1); ind_LT is Matt's small approx., ind_GE is Matt's large approx
    brho1  = DBLARR(N_ELEMENTS(rho1), NOZERO=1) * !VALUES.D_NaN;% NaNs
    bz1    = brho1;% NaNs
  ENDELSE

  IF (n_ind_LT NE 0) THEN BEGIN ;% Matt's small approx
    zmd2 = z1md[ind_LT]*z1md[ind_LT]
    zpd2 = z1pd[ind_LT]*z1pd[ind_LT]
    f1_sq = zmd2 + r_sq
    f2_sq = zpd2 + r_sq
    f1 = SQRT(f1_sq)
    f2 = SQRT(f2_sq)
    f1_cubed = f1_sq*f1
    f2_cubed = f2_sq*f2

    ;% calculate the terms which make equations 9a and 9b
    rhoov2   = rho1[ind_LT]/2d
    rho2ov4  = rhoov2  * rhoov2
    rho3ov16 = rho2ov4 * rhoov2/2d

    ;% these bits are used to form 9a
    f3 = (r_sq - 2d*zmd2)/(f1_sq*f1_sq*f1)
    f4 = (r_sq - 2d*zpd2)/(f2_sq*f2_sq*f2)

    terma0 = rhoov2  *(1d/f1 - 1d/f2)
    terma1 = rho3ov16*(f3    - f4   )
    brho1[ind_LT] = mu_i_div2__current_parameter_nT*(terma0 + terma1)

    ;% now equation 9b
    termb0 = ALOG((z1pd[ind_LT] + f2)/(z1md[ind_LT] + f1))
    termb1 = rho2ov4*(z1pd[ind_LT]/f2_cubed - z1md[ind_LT]/f1_cubed)
    bz1[  ind_LT] = mu_i_div2__current_parameter_nT*(termb0 + termb1)
  ENDIF
  IF (n_ind_GE NE 0) THEN BEGIN
    zmd2 = z1md[ind_GE]*z1md[ind_GE]
    zpd2 = z1pd[ind_GE]*z1pd[ind_GE]
    f1_sq = zmd2 + rho1_sq[ind_GE]
    f2_sq = zpd2 + rho1_sq[ind_GE]
    f1 = SQRT(f1_sq)
    f2 = SQRT(f2_sq)
    f1_cubed = f1_sq*f1
    f2_cubed = f2_sq*f2

    ;%equation 13a
    terma0 = (1d/rho1[ind_GE])*(f1 - f2)
    terma1 = (rho1[ind_GE]*r_sq/4d)*(1d/f2_cubed - 1d/f1_cubed)
    ;% from Python z1_clip = z1.clip(max=D,min=-D)
    IF scalar_input THEN BEGIN
      IF   z1 GT  d__cs_half_thickness_rj THEN z1_clip =  d__cs_half_thickness_rj ELSE $
        IF z1 LT -d__cs_half_thickness_rj THEN z1_clip = -d__cs_half_thickness_rj ELSE z1_clip = z1
    ENDIF ELSE BEGIN
      z1_clip = z1[ind_GE]
      z1_clip(WHERE(z1_clip GT  d__cs_half_thickness_rj, NULL =1)) =  d__cs_half_thickness_rj
      z1_clip(WHERE(z1_clip LT -d__cs_half_thickness_rj, NULL =1)) = -d__cs_half_thickness_rj
    ENDELSE
    terma2 = (2d/rho1[ind_GE])*z1_clip
    brho1[ind_GE] = mu_i_div2__current_parameter_nT*(terma0 + terma1 + terma2)

    ;%equation 13b - same as before
    bz1[ind_GE] = mu_i_div2__current_parameter_nT*(ALOG((z1pd[ind_GE]+f2)/(z1md[ind_GE]+f1)) + (r_sq/4d)*((z1pd[ind_GE]/f2_cubed) - (z1md[ind_GE]/f1_cubed)))
  ENDIF

  RETURN, [[brho1], [bz1]]
END
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;FUNCTION _con2020_model_int_tabulated_rjw2_sub, X, F
;  ; Numerical Integration based on Trapezoidal rule only https://www.math24.net/trapezoidal-rule
;  ; Assumes equally spaced x-data
;  COMPILE_OPT HIDDEN
;  ON_ERROR, 2 ;Return to caller if an error occurs.
;  Xsegments = N_ELEMENTS(X) - 1L
;  h = (X[Xsegments] - X[0]) / DOUBLE(Xsegments)
;  RETURN, h * (TOTAL(F) - (F[0] + F[Xsegments])*0.5d )
;END
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
FUNCTION _con2020_model_xyz_params_to_string, params
  ; Writes out the Con2020 code parameters structure to string, that could be executed later.
  ;
  ; To change the default parameters of the Con2020 model, first get the parameters with:
  ;   params = con2020_model_xyz('default_values')
  ; Then change any existing values as you wish, e.g.
  ;   params.I_RHO__RADIAL_CURRENT_MA = 20.2d ; don't forget the d for double!
  ; to write out the new structure to a string that could be executed later:
  ;   struct_string = _con2020_model_xyz_params_to_string(params)
  ; which would give the string:
  ;   {MU_I_DIV2__CURRENT_PARAMETER_NT:139.6d,R0__INNER_RJ:7.8d,R1__OUTER_RJ:51.4d,D__CS_HALF_THICKNESS_RJ:3.6d,XT__CS_TILT_DEGS:9.3d,XP__CS_RHS_AZIMUTHAL_ANGLE_OF_TILT_DEGS:155.8d,I_RHO__RADIAL_CURRENT_MA:20.2d,ERROR_CHECK:1}
  ;
  ; The Con2020 code itself does not use this function. However, this may be useful later
  ; for Field Line Tracing codes, etc.

  COMPILE_OPT HIDDEN
  ON_ERROR, 2

  IF (SIZE(params, /TYPE) NE 8) THEN MESSAGE,'ERROR: Expected an input structure'

  T     = TAG_NAMES(params)
  nT_m1 = N_ELEMENTS(T) - 1L

  cmd  = "{"
  FOR z = 0L, nT_m1 DO BEGIN
    IF ISA(params.(z), /NUMBER, /SCALAR) EQ 0 THEN MESSAGE,'ERROR: Fields are not all scalar numbers.'
    IF ISA(params.(z), /COMPLEX        ) EQ 1 THEN MESSAGE,'ERROR: Fields are not all scalar real numbers: some are complex.'

    cmd += T[z] + ':'

    IF (ISA(params.(z),'DOUBLE') EQ 1) THEN BEGIN

      cmd += STRTRIM(     params.(z)  ,2)
      ; Remove trailing zeroes if stuff after decimal place
      WHILE (STRCMP(STRMID(cmd, 0, 1, /REVERSE_OFFSET),'0') EQ 1) DO BEGIN
        cmd = STRMID(cmd, 0L, STRLEN(cmd)-1L )
        IF (STRLEN(cmd) LT 2) THEN MESSAGE,'ERROR: Had issues removing trailing zeros.'
      ENDWHILE
      ; If number ends with a decimal point, remove it.  i.e, 20.00000 becomes 20, nto 20.
      IF (STRCMP(STRMID(cmd, 0, 1, /REVERSE_OFFSET),'.') EQ 1) THEN cmd = STRMID(cmd, 0L, STRLEN(cmd)-1L )
      ; Add IDL symbol for double
      cmd += 'd' ; d for double

    ENDIF ELSE BEGIN ; must be a byte, for ERROR_CHECK

      IF (ISA(params.(z),'BYTE') EQ 0) THEN MESSAGE,'ERROR: If not a Double, was expecting a Byte for ERROR_CHECK'
      cmd += STRTRIM( FIX(params.(z)) ,2)

    ENDELSE

    IF (z NE nT_m1) THEN cmd += ','

  ENDFOR
  cmd += '}'

  RETURN, cmd
END
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  Con2020 function begins here!
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
FUNCTION con2020_model_xyz, eq_type, x_rj, y_rj, z_rj, use_these_params
  ;% ======
  ;% CAN_SHEET_VARIABLE_2020_HYBRID (Cartesian)
  ;% ======
  ;% Code to calculate the perturbation magnetic field produced by the Connerney et al. 1981 (CAN) current sheet, which is
  ;% represented by a finite disk of current.
  ;%
  ;% More details are provided in comments at top of this file.

  ON_ERROR, 2 ; % Exit code if an error in main, don't stop in code - no Matlab equivalent, just delete line in Matlab
  FORWARD_FUNCTION _con2020_model_xyz_analytic ; telling IDL this is a function, in case _con2020_model_xyz_analytic has not been compiled yet.
  ; FORWARD_FUNCTION _con2020_model_int_tabulated_rjw2_sub ; no longer used, but function kept in code above for reference.

  ;% RJW  I keep sin and cos type commands lower case (IDL doesn't care, Matlab will), less to edit when copied later
  ;% RJW '%' is the comment symbol in Matlab, again thinking ahead.

  IF KEYWORD_SET(use_these_params) THEN BEGIN
    IF (ISA(use_these_params,'STRUCT') EQ 0) THEN MESSAGE,'Must be a structure of terms to use in code'
    IF N_ELEMENTS(TAG_NAMES(use_these_params)) NE 8 THEN MESSAGE,'ERROR: Expecting 8 fields in the structure, not '+STRTRIM(N_ELEMENTS(TAG_NAMES(use_these_params)),2)
    mu_i_div2__current_parameter_nT         = DOUBLE(use_these_params.mu_i_div2__current_parameter_nT           )
    r0__inner_rj                            = DOUBLE(use_these_params.r0__inner_rj                              )
    r1__outer_rj                            = DOUBLE(use_these_params.r1__outer_rj                              )
    d__cs_half_thickness_rj                 = DOUBLE(use_these_params.d__cs_half_thickness_rj                   )
    xt__cs_tilt_degs                        = DOUBLE(use_these_params.xt__cs_tilt_degs                          )
    xp__cs_rhs_azimuthal_angle_of_tilt_degs = DOUBLE(use_these_params.xp__cs_rhs_azimuthal_angle_of_tilt_degs   )
    i_rho__radial_current_MA                = DOUBLE(use_these_params.i_rho__radial_current_MA                  )
    error_check = DOUBLE(use_these_params.error_check) ;

  ENDIF ELSE BEGIN
    ; Make sure all of these numbers are doubles!
    Default_values = {$
      mu_i_div2__current_parameter_nT         : 139.6d  , $ ;% current sheet field parameter (nT)
      r0__inner_rj                            :   7.8d  , $ ;% inner radius (Rj)
      r1__outer_rj                            :  51.4d  , $ ;% outer radius (Rj)
      d__cs_half_thickness_rj                 :   3.6d  , $ ;% half-height  (Rj)
      xt__cs_tilt_degs                        :   9.3d  , $ ;% current sheet tilt (Deg.)
      xp__cs_rhs_azimuthal_angle_of_tilt_degs : 155.8d  , $ ;% current sheet longitude (right handed) (Deg.), Table 1 xp = 204.2 but that value is in left handed SIII
      i_rho__radial_current_MA                :  16.7d  , $ ;% radial current term from Connerney et al., 2020, in mega-Amps.
      ;% NOTE: The default value (16.7 MA) is the average value from Connerney et al 2020. This value was shown to vary from one
      ;% pass to the next, where Table 2 (units of MA) provides radial current intensity values for 23 of the first 24 perijoves
      ;% (units mistakenly given as nT in the article text).
      error_check : 1b } ;% input error check: 1 = yes, 0 = no

    IF STRCMP(eq_type,'default_values',FOLD_CASE=1) THEN BEGIN
      PRINT,'Returning structure of Default terms used in code.'
      RETURN,Default_values
    ENDIF

    mu_i_div2__current_parameter_nT         = Default_values.mu_i_div2__current_parameter_nT
    r0__inner_rj                            = Default_values.r0__inner_rj
    r1__outer_rj                            = Default_values.r1__outer_rj
    d__cs_half_thickness_rj                 = Default_values.d__cs_half_thickness_rj
    xt__cs_tilt_degs                        = Default_values.xt__cs_tilt_degs
    xp__cs_rhs_azimuthal_angle_of_tilt_degs = Default_values.xp__cs_rhs_azimuthal_angle_of_tilt_degs
    i_rho__radial_current_MA                = Default_values.i_rho__radial_current_MA
    error_check = Default_values.error_check
  ENDELSE

  Deg2Rad = !DPI/180d;

  N_input = N_ELEMENTS(x_rj)
  scalar_input = (N_input EQ 1) ;% scalar or not

  eq_type = STRLOWCASE(eq_type) ;% make lower case

  ;% RJW Still need to check inputs (if set) are scalar values, and set to doubles.  (Exception of equation type)
  IF error_check[0] THEN BEGIN
    IF (ISA(error_check, NUMBER=1, SCALAR=1) EQ 0) THEN MESSAGE,'ERROR: Field  error_check must be a scalar number'
    IF (ISA(mu_i_div2__current_parameter_nT        , NUMBER=1, SCALAR=1) EQ 0) THEN MESSAGE,'ERROR: Field     mu_i_div2__current_parameter_nT     must be a scalar number'
    IF (ISA(r0__inner_rj                           , NUMBER=1, SCALAR=1) EQ 0) THEN MESSAGE,'ERROR: Field               r0__inner_rj              must be a scalar number'
    IF (ISA(r1__outer_rj                           , NUMBER=1, SCALAR=1) EQ 0) THEN MESSAGE,'ERROR: Field               r1__outer_rj              must be a scalar number'
    IF (ISA(d__cs_half_thickness_rj                , NUMBER=1, SCALAR=1) EQ 0) THEN MESSAGE,'ERROR: Field          d__cs_half_thickness_rj        must be a scalar number'
    IF (ISA(xt__cs_tilt_degs                       , NUMBER=1, SCALAR=1) EQ 0) THEN MESSAGE,'ERROR: Field             xt__cs_tilt_degs            must be a scalar number'
    IF (ISA(xp__cs_rhs_azimuthal_angle_of_tilt_degs, NUMBER=1, SCALAR=1) EQ 0) THEN MESSAGE,'ERROR: Field xp__cs_rhs_azimuthal_angle_of_tilt_degs must be a scalar number'
    IF (ISA(i_rho__radial_current_MA               , NUMBER=1, SCALAR=1) EQ 0) THEN MESSAGE,'ERROR: Field         i_rho__radial_current_MA        must be a scalar number'

    ;IF (error_check NE 0) AND (error_check NE 1)   THEN MESSAGE,'ERROR: Field  error_check must be 0 or 1'
    ;% if here error_check must be scalar and not 0
    IF error_check NE 1 THEN MESSAGE,'ERROR: Field  error_check must be 0 or 1'
    error_check = 1b
    IF mu_i_div2__current_parameter_nT LE 0d           THEN MESSAGE,'mu_i_div2__current_parameter_nT must be GT 0'
    IF r0__inner_rj                    LE 0d           THEN MESSAGE,                   'r0__inner_rj must be GT 0'
    IF r1__outer_rj                    LE r0__inner_rj THEN MESSAGE,                   'r1__outer_rj must be GT r0__inner_rj'
    IF d__cs_half_thickness_rj         LE 0d           THEN MESSAGE,        'd__cs_half_thickness_rj must be GT 0'


    IF (ISA(eq_Type, STRING=1, SCALAR=1) EQ 0) THEN MESSAGE,'ERROR: First argument equation_type must be a scalar string and can only be ''hybrid'' (default), ''analytic'' or ''integral'' (or ''default_values'')'
    IF (TOTAL(STRCMP(eq_type,['integral','analytic','hybrid'])) EQ 0) THEN MESSAGE,'ERROR: First argument equation_type  can only be ''hybrid'' (default), ''analytic'' or ''integral'' (or ''default_values'')' ;% must be lower case here

    ;% Check inputs x_rj, y_rj and z_rj are all numbers, and same size (also scalar or 1D only)
    ; N_DIMENSION=1 is 1 if 1D, and 0 if scalar, hence the GT 1 below
    IF (ISA(x_rj      , NUMBER=1) EQ 0) OR (SIZE(x_rj      , N_DIMENSIONS=1) GT 1) THEN MESSAGE,'ERROR: Second argument x_rj must be a scalar number or 1D array of numbers';
    IF (ISA(y_rj      , NUMBER=1) EQ 0) OR (SIZE(y_rj      , N_DIMENSIONS=1) GT 1) THEN MESSAGE,'ERROR: Third  argument y_rj must be a scalar number or 1D array of numbers';
    IF (ISA(z_rj      , NUMBER=1) EQ 0) OR (SIZE(z_rj      , N_DIMENSIONS=1) GT 1) THEN MESSAGE,'ERROR: Fourth argument z_rj must be a scalar number or 1D array of numbers';
    IF (N_input NE N_ELEMENTS(y_rj)) THEN MESSAGE,'ERROR: Second argument x_rj must be the same size as 3rd argument y_rj';
    IF (N_input NE N_ELEMENTS(z_rj)) THEN MESSAGE,'ERROR: Second argument x_rj must be the same size as 4th argument z_rj';
    IF scalar_input THEN BEGIN
      IF (x_rj LE -200d) OR (x_rj GE 200d) OR (y_rj LE -200d) OR (y_rj GE 200d) OR (z_rj LE -200d) OR (z_rj GE 200d) THEN $
        MESSAGE,'ERROR: Positions x_rj, y_rj and z_rj must all be in units of Rj, and each >-200 and <200 only, and not outside that range (did you use km instead?)'
    ENDIF ELSE BEGIN
      min_x = MIN(x_rj, MAX=max_x)
      min_y = MIN(y_rj, MAX=max_y)
      min_z = MIN(z_rj, MAX=max_z)
      IF (min_x LE -200d) OR (max_x GE 200d) THEN MESSAGE,'ERROR: Second argument, Position x_rj, must be in units of Rj and >-200 and <200 only, and not outside that range (did you use km instead?)'
      IF (min_x LE -200d) OR (max_y GE 200d) THEN MESSAGE,'ERROR: Third  argument, Position y_rj, must be in units of Rj and >-200 and <200 only, and not outside that range (did you use km instead?)'
      IF (min_z LE -200d) OR (max_z GE 200d) THEN MESSAGE,'ERROR: Fourth argument, Position z_rj, must be in units of Rj and >-200 and <200 only, and not outside that range (did you use km instead?)'
    ENDELSE
  ENDIF

  xp__cs_rhs_azimuthal_angle_of_tilt_degs -= 180.0d ;shift needed because of the way we define the rotation angle
  dipole_shift = xp__cs_rhs_azimuthal_angle_of_tilt_degs * Deg2Rad; % xp__cs_rhs_azimuthal_angle_of_tilt_degs is longitude of the current sheet tilt (roughly the dipole longitude). dipole_shift used here and at end of code
  theta_cs     = xt__cs_tilt_degs * Deg2Rad ; % dipole tilt is xt__cs_tilt_degs
  cos_dipole_shift = cos(dipole_shift)
  sin_dipole_shift = sin(dipole_shift)
  cos_theta_cs     = cos(theta_cs)
  sin_theta_cs     = sin(theta_cs)

  ;% Covert to Doubles, and rename input variables so as not to over-write them (only an IDL issue)
  xSYSIII = DOUBLE(x_rj)
  ySYSIII = DOUBLE(y_rj)
  zSYSIII = DOUBLE(z_rj)

  ;% Two rotations are needed to align the IAU_JUPITER frame with the MAG frame
  ;% first by (if using defaults) +155.2 degrees about Z, second by -9.3 degrees about Y.
  ;% so first rotate around Z
  x0 =  xSYSIII*cos_dipole_shift + ySYSIII*sin_dipole_shift
  y0 = -xSYSIII*sin_dipole_shift + ySYSIII*cos_dipole_shift
  z0 =  zSYSIII; % no chance as we rotate about Z
  ;% Second rotate around Y
  x1 =  x0*cos_theta_cs + z0*sin_theta_cs
  y1 =  y0; % no change as we rotate about Y %RJW - NOT NEEDED??? - But used in ATAN equivalent code near end of function (could just use y)
  z1 = -x0*sin_theta_cs + z0*cos_theta_cs
  ;% Calculate rho
  rho1_sq = x1*x1 + y1*y1
  rho1 = sqrt(rho1_sq) ; %cylindrical radial distance


  abs_z1 = abs(z1)

  ;% Decide whether to use integral equations or analytic equations
  do_integral = BYTARR(N_input); preallocate array of right size, each element 0b
  CASE 1 OF ; Only check first letter for speed, since a, h and i are unique for analytic, hybrid and integral
    STRCMP(eq_type,'analytic') : do_integral[0] = 0b ;% Really do nothing as it's already all zeros, so just setting element to 0 (from the already 0), just to have a command here for the analytic case
    STRCMP(eq_type,'hybrid'  ) : do_integral[WHERE((abs_z1 LE d__cs_half_thickness_rj*1.5d) AND (abs(rho1-r0__inner_rj) LE 2d),NULL=1)] = 1b ;% scalar form is IF ((abs_z1 LE d__cs_half_thickness_rj*1.5d) and (abs(rho1-r0__inner_rj) LE 2.0d)) THEN do_integral = 1 ELSE do_integral = 0;
    STRCMP(eq_type,'integral') : do_integral[*] = 1b
    ELSE : MESSAGE,'ERROR: case statement has unrecognized string - was your equation_type lower case?' ;% if no_error_check, this will still catch bad inputs for equation_type
  ENDCASE


  IF scalar_input THEN BEGIN
    brho1 = 0d
    bz1   = 0d
    IF do_integral[0] EQ 0 THEN BEGIN
      n_ind_analytic = 1
      n_ind_integral = 0
      ind_analytic   = [0]
      ind_integral   = !NULL
    ENDIF ELSE BEGIN
      n_ind_analytic = 0
      n_ind_integral = 1
      ind_analytic   = !NULL
      ind_integral   = [0]
    ENDELSE
  ENDIF ELSE BEGIN
    ind_analytic = WHERE(do_integral EQ 0b,n_ind_analytic,COMPLEMENT=ind_integral, NCOMPLEMENT=n_ind_integral,NULL=1)
    ;% preallocate arrays with NaNs (For safety)
    brho1    = DBLARR(N_input,NOZERO=1); make array of right size, don't care what the values are (hence NOZERo)
    brho1[*] = !Values.D_NAN ;% set all to NaN
    bz1      = brho1         ;% copy of Brho of array of NaNs
  ENDELSE



  IF (n_ind_integral NE 0) THEN BEGIN
    ;% Integral equations - Brho and Bz eqs. vary depending on region with respect to the disk
    dlambda_brho    = 1d-4  ;% default step size for Brho function
    dlambda_bz      = 5d-5  ;% default step size for Bz function

    check1 = abs(abs_z1[ind_integral] -  d__cs_half_thickness_rj      )
    check2 =     abs_z1[ind_integral] LE d__cs_half_thickness_rj*1.1d  ;% Byte 1 or 0 (true or false)
    FOR zcase = 1L,6L DO BEGIN
      CASE zcase OF
        1 : BEGIN
          check3 = (check2 EQ 1b) AND (check1 GE 0.7d);
          lambda_max_brho =   4d ;% default integration limit for Brho function;
          lambda_max_bz   = 100d
        END
        2 : BEGIN
          check3 = (check2 EQ 0b) AND (check1 GE 0.7d);
          lambda_max_brho =   4d ;% default integration limit for Brho function;
          lambda_max_bz   =  20d ;% Small Z or default in else
        END
        3 : BEGIN
          check3 = (check2 EQ 1b) AND (check1 GE 0.1d) AND (check1 lt 0.7d);
          lambda_max_brho =  40d ;% Z > D             ;
          lambda_max_bz   = 100d
        END
        4 : BEGIN
          check3 = (check2 EQ 0b) AND (check1 GE 0.1d) AND (check1 lt 0.7d);
          lambda_max_brho =  40d ;% Z > D             ;
          lambda_max_bz   =  20d ;% Small Z or default in else
        END
        5 : BEGIN
          check3 = (check2 EQ 1b) AND (check1 LT 0.1d);
          lambda_max_brho = 100d ;% Z very close to D
          lambda_max_bz   = 100d
        END
        6 : BEGIN
          check3 = (check2 EQ 0b) AND (check1 LT 0.1d);
          lambda_max_brho = 100d ;% Z very close to D;
          lambda_max_bz   =  20d ;% Small Z or default in else
        END
        ELSE: MESSAGE,'Should never get to this else case statement'
      ENDCASE

      IF scalar_input THEN BEGIN
        IF (check3[0] EQ 0b) THEN CONTINUE
        ind_case   = 0L
        ;n_ind_case = 1L; % do not need this line in scalar case
        n_ind_case_minus_1 = 0
      ENDIF ELSE BEGIN
        ;ind_case = WHERE(check3 EQ 1b, n_ind_case,NULL = 1)
        ind_case = WHERE(check3      , n_ind_case,NULL = 1)
        if n_ind_case EQ 0 THEN CONTINUE ; Can not do if NOT n_ind_case THEN CONTINUE, since IDL NOT 5 = -6, etc.
        n_ind_case_minus_1 = n_ind_case -1L
      ENDELSE

      lambda_int_brho = DINDGEN(lambda_max_brho/dlambda_brho - 1d, INCREMENT=dlambda_brho, START=dlambda_brho)
      lambda_int_bz   = DINDGEN(lambda_max_bz  /dlambda_bz   - 1d, INCREMENT=dlambda_bz  , START=dlambda_bz  )

      beselj_rho_r0_0   = beselj(lambda_int_brho*r0__inner_rj,0); % Only 6 sets of values
      beselj_z_r0_0     = beselj(lambda_int_bz  *r0__inner_rj,0); % Only 6 sets of values

      FOR zi = 0L, n_ind_case_minus_1 DO BEGIN
        ind_for_integral = ind_integral[ind_case[zi]];% sub-indices of sub-indices!

        beselj_rho_rho1_1 = beselj(lambda_int_brho*rho1[ind_for_integral]    ,1)
        beselj_z_rho1_0   = beselj(lambda_int_bz  *rho1[ind_for_integral]    ,0)

        IF (abs_z1[ind_for_integral] gt d__cs_half_thickness_rj) THEN BEGIN ;% Connerney et al. 1981 eqs. 14 and 15
          brho_int_funct = beselj_rho_rho1_1*beselj_rho_r0_0*sinh(d__cs_half_thickness_rj*lambda_int_brho)*exp(-abs_z1[ind_for_integral]*lambda_int_brho)/lambda_int_brho
          bz_int_funct   = beselj_z_rho1_0  *beselj_z_r0_0  *sinh(d__cs_half_thickness_rj*lambda_int_bz  )*exp(-abs_z1[ind_for_integral]*lambda_int_bz  )/lambda_int_bz

          ;brho1[ind_for_integral] = mu_i_div2__current_parameter_nT*2d * _con2020_model_int_tabulated_rjw2_sub(lambda_int_brho,brho_int_funct)
          ;% Can use Trapezoidal rule approx
          brho1[ind_for_integral] = mu_i_div2__current_parameter_nT*2d * dlambda_brho * (TOTAL(brho_int_funct) - (brho_int_funct[0] + brho_int_funct[-1])*0.5d )

          IF z1[ind_for_integral] LT 0 THEN brho1[ind_for_integral] = -brho1[ind_for_integral]

        ENDIF ELSE BEGIN ;% Connerney et al. 1981 eqs. 17 and 18
          brho_int_funct = beselj_rho_rho1_1*beselj_rho_r0_0*(     sinh(z1[ind_for_integral]*lambda_int_brho)*exp(-d__cs_half_thickness_rj*lambda_int_brho))/lambda_int_brho
          bz_int_funct   = beselj_z_rho1_0  *beselj_z_r0_0  *(1d - cosh(z1[ind_for_integral]*lambda_int_bz  )*exp(-d__cs_half_thickness_rj*lambda_int_bz  ))/lambda_int_bz

          ;brho1[ind_for_integral] = mu_i_div2__current_parameter_nT*2d * _con2020_model_int_tabulated_rjw2_sub(lambda_int_brho,brho_int_funct)
          ;% Can use Trapezoidal rule approx
          brho1[ind_for_integral] = mu_i_div2__current_parameter_nT*2d * dlambda_brho * (TOTAL(brho_int_funct) - (brho_int_funct[0] + brho_int_funct[-1])*0.5d )
        ENDELSE

        ;bz1[ind_for_integral]   = mu_i_div2__current_parameter_nT*2d * _con2020_model_int_tabulated_rjw2_sub(lambda_int_bz  ,bz_int_funct  )
        ;% Can use Trapezoidal rule approx
        bz1[ind_for_integral]   = mu_i_div2__current_parameter_nT*2d * dlambda_bz * (TOTAL(bz_int_funct) - (bz_int_funct[0] + bz_int_funct[-1])*0.5d )

      ENDFOR

      IF scalar_input THEN BREAK ;% If only scalar, stop doing for loops once we found the one we want
    ENDFOR

  ENDIF

  ;% Work out the finite sheet here - re-use some bits for the final analy
  brho_finite_AND_bz_finite = _con2020_model_xyz_analytic( rho1, z1, rho1_sq, d__cs_half_thickness_rj, r1__outer_rj, mu_i_div2__current_parameter_nT, scalar_input)
  ;brho_finite = brho_finite_AND_bz_finite[*,0]; we'll just use brho_finite_AND_bz_finite later
  ;bz_finite   = brho_finite_AND_bz_finite[*,1]

  ;% Do now near Analytic bit
  IF (n_ind_analytic NE 0) THEN BEGIN
    brho1_AND_bz1         = _con2020_model_xyz_analytic( rho1[ind_analytic], z1[ind_analytic], rho1_sq[ind_analytic], d__cs_half_thickness_rj, r0__inner_rj, mu_i_div2__current_parameter_nT, scalar_input)
    brho1[ind_analytic] = brho1_AND_bz1[*,0]
    bz1[  ind_analytic] = brho1_AND_bz1[*,1]
  ENDIF

  ;% Check that brho1 and bz1 do not contain NaNs or Infs
  IF error_check THEN BEGIN
    IF TOTAL(FINITE(brho1)) NE N_input THEN MESSAGE,'ERROR: some brho1 values are NaN or Infs'
    IF TOTAL(FINITE( bz1 )) NE N_input THEN MESSAGE,'ERROR: some  bz1  values are NaN or Infs'
  ENDIF


  ;% New to CAN2020 (not included in CAN1981): radial current produces an azimuthal field, so Bphi is nonzero
  ;% bphi1 = 2.7975d*i_rho__radial_current_MA/rho1
  ;% Above could give Infinities if rho1 = 0, or NaN is rho1 = 0 and i_rho__radial_current_MA=0 and rho1 = 0
  ;% deal with this below, separately for scalar_input or vector
  IF scalar_input THEN BEGIN
    IF rho1 NE 0d THEN BEGIN
      bphi1 = 2.7975d*i_rho__radial_current_MA/rho1
      
      IF abs_z1 LT d__cs_half_thickness_rj  THEN bphi1 =  bphi1 * abs_z1 / d__cs_half_thickness_rj

      IF     z1 GT     0d                   THEN bphi1 = -bphi1
    ENDIF ELSE BEGIN ;% if rho1 == 0, a very very rare occurence
      bphi1 = 0d ;% Remove Infinity or NaN
    ENDELSE
  ENDIF ELSE BEGIN
    ind_rho1_ne_0 = WHERE(rho1 NE 0d, NULL = 1) ;% we use this again later

    bphi1 = DBLARR(N_input) ;% set bphi1 to zeros, the rho1 == 0 case
    bphi1[ind_rho1_ne_0] = 2.7975d*i_rho__radial_current_MA/rho1[ind_rho1_ne_0]

    ind = WHERE(abs_z1 LT d__cs_half_thickness_rj, NULL=1)
    IF N_ELEMENTS(ind) NE 0 THEN bphi1[ind] =  bphi1[ind] * abs_z1[ind] / d__cs_half_thickness_rj

    ind = WHERE(    z1 GT 0d     , NULL=1)
    IF N_ELEMENTS(ind) NE 0 THEN bphi1[ind] = -bphi1[ind]
  ENDELSE


  ;% Account for finite nature of current sheet by subtracting the field values
  ;% calculated brho_finite and bz_finite earlier
  ;brho1 = brho1 - brho_finite;
  ;bz1   = bz1   - bz_finite  ;
  brho1 -= brho_finite_AND_bz_finite[*,0]
  bz1   -= brho_finite_AND_bz_finite[*,1]


  ;% brho1, bphi1, and bz1 here are the ultimately calculated brho and bz values from the CAN model
  ;% the remaining calculations just rotate the field back into SIII

  ;% Calculate 'magnetic longitude' and convert the field into cartesian coordinates
  ;% rho1 is alwas positive, from above, but could be 0, and don't want a divide by zero
  IF scalar_input THEN BEGIN
    IF rho1 NE 0d THEN BEGIN
      cos_phi1 = x1/rho1
      sin_phi1 = y1/rho1
      bx1 = brho1*cos_phi1 - bphi1*sin_phi1
      by1 = brho1*sin_phi1 + bphi1*cos_phi1
    ENDIF ELSE BEGIN ;% else rho1==0
      bx1 = 0d
      by1 = 0d
    ENDELSE
  ENDIF ELSE BEGIN
    ;% ind_rho1_ne_0 was defined in the radial current section
    IF N_ELEMENTS(ind_rho1_ne_0) EQ N_input THEN BEGIN 
      cos_phi1 = x1/rho1
      sin_phi1 = y1/rho1
      bx1 = brho1*cos_phi1 - bphi1*sin_phi1
      by1 = brho1*sin_phi1 + bphi1*cos_phi1
    ENDIF ELSE BEGIN
      bx1 = DBLARR(N_input) ;% set bx1 and by1 to zeros, the rho1 == 0 case
      by1 = bx1
      cos_phi1 = x1[ind_rho1_ne_0]/rho1[ind_rho1_ne_0]
      sin_phi1 = y1[ind_rho1_ne_0]/rho1[ind_rho1_ne_0]
      bx1[ind_rho1_ne_0] = brho1[ind_rho1_ne_0]*cos_phi1 - bphi1[ind_rho1_ne_0]*sin_phi1
      by1[ind_rho1_ne_0] = brho1[ind_rho1_ne_0]*sin_phi1 + bphi1[ind_rho1_ne_0]*cos_phi1
    ENDELSE
  ENDELSE

  ;% Now convert back to SYSIII

  ;% Rotate back by current sheet tilt amount, into coordinate system that is aligned with Jupiter's spin axis
  bx = bx1*cos_theta_cs - bz1*sin_theta_cs
  ;by = by1; % just using by1 below
  bz = bx1*sin_theta_cs + bz1*cos_theta_cs

  ;% Finally, shift back to SIII longitude
  cos_xp = cos(dipole_shift)
  sin_xp = sin(dipole_shift)
  bx2 = bx *cos_xp - by1*sin_xp; % used sin(-a) = asin(a) & cos(-a) = cos(a)
  by2 = by1*cos_xp + bx *sin_xp
  ;% bz2 = bz so just using bz below

  RETURN, [[bx2],[by2],[bz]]
END
