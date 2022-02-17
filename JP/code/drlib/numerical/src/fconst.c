#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif

typedef struct {
    Mchar	*name;
    Mfloat	value;
    Mchar	*units;
} Const_table;

#ifdef ANSI
static void     l_cunit(Mfloat x, Mchar *xunits, Mfloat *y, Mchar *yunits);
static Mint	l_comp_const(Mvoid *t1, Mvoid *t2);
static Mint	l_iicsr(Mchar *s1, Mchar *s2);
static void 	l_c3nit (Mchar *unit, Mint *ipow, Mfloat *cvt);
static Mint     l_comp_conv(Mvoid *t1, Mvoid *t2);
#else
static void     l_cunit();
static Mint	l_comp_const();
static Mint	l_iicsr();
static void 	l_c3nit ();
static Mint     l_comp_conv();
#endif

#define	INCH	            0.0254e0
#define	ATM	            1.01325e+5
#define AU	            1.496e+11
#define POUND	            0.45359237
#define CLIGHT              2.997924580e+8
#define STANDARD_GRAVITY    9.80665

/*
 -----------------------------------------------------------------------
   IMSL Name:  CONST/DCONST (Single/Double precision version)
 
   Computer:   $COMPUTER/$PRECISION
 
   Revised:    January 19, 1990
 
   Purpose:    Return the value of various mathematical and physical
               constants.
 
   Usage:      CONST(NAME)
 
   Arguments:
      NAME   - Character string containing the name of the desired
               constant.  (Input)
               See remark 3 for a list of valid constants.
      CONST  - Value of the constant.  (Output)
 
   Remarks:
   1. The case of the character string in NAME does not matter.  The
      names 'PI', 'Pi', 'pI' and 'pi' are equivalent.
 
   2. The units of the physical constants are in SI units, (meter-
      kilogram-second).
 
   3. The names allowed are as follows:
        Name             Description            Value          Reference
      AMU              Atomic mass unit       1.6605402E-27 kg       [1]
      ATM              Standard atm pressure  1.01325E+5 N/m^2      E[2]
      AU               Astronomical unit      1.496E+11 m            [ ]
      Avogadro         Avogadro's number      6.0221367E+23 1/mole   [1]
      Boltzman         Boltzman's constant    1.380658E-23 J/K       [1]
      C                Speed of light         2.997924580E+8 m/s    E[1]
      Catalan          Catalan's constant     0.915965...           E[3]
      E                Base of natural logs   2.718...              E[3]
      ElectronCharge   Electron change        1.60217733E-19 C       [1]
      ElectronMass     Electron mass          9.1093897E-31 kg       [1]
      ElectronVolt     Electron volt          1.60217733E-19 J       [1]
      Euler            Euler's constant gamma 0.577...              E[3]
      Faraday          Faraday constant       9.6485309E+4 C/mole    [1]
      FineStructure    Fine structure         7.29735308E-3          [1]
      Gamma            Euler's constant       0.577...              E[3]
      Gas              Gas constant           8.314510 J/mole/K      [1]
      Gravity          Gravitational constant 6.67259E-11N*m^2/kg^2  [1]
      Hbar             Planck constant / 2 pi 1.05457266E-34 J*s     [1]
      PerfectGasVolume Std vol ideal gas      2.241383E-2 m^3/mole   [*]
      Pi               Pi                     3.141...              E[3]
      Planck           Planck's constant h    6.6260755E-34 J*s      [1]
      ProtonMass       Proton mass            1.6726231E-27 kg       [1]
      Rydberg          Rydberg's constant     1.0973731534E+7 /m     [1]
      SpeedLight       Speed of light         2.997924580E+8 m/s    E[1]
      StandardGravity  Standard g             9.80665 m/s^2         E[2]
      StandardPressure Standard atm pressure  1.01325E+5 N/m^2      E[2]
      StefanBoltzmann  Stefan-Boltzman        5.67051E-8 W/K^4/m^2   [1]
      WaterTriple      Triple point of water  2.7316E+2 K           E[2]
 
   Keywords:   Utilities; Mathematical constants; Physical constants;
               Avogadro's number; Atomic mass unit; Astronomical unit;
               Boltzman's constant; Speed of light; Catalan's constant;
               Electron charge; Electron mass; Electron volt; Gamma;
               Euler's constant; Faraday's constant; Fine structure
               constant; Gas constant; Gravitational constant;
               Planck's constant; Protron mass; Stefan-Boltzman's
               constant; Pi; Triple point of water
 
   GAMS:       C19
 
   Chapter:    MATH/LIBRARY Utilities
 
   Copyright:  1986 by IMSL, Inc.  All Rights Reserved.
 
   Warranty:   IMSL warrants only that IMSL testing has been applied
               to this code.  No other warranty, expressed or implied,
               is applicable.
 
 -----------------------------------------------------------------------
*/

static Const_table	lv_const_table[] = {
    {"AMU",		1.6605402e-27,	    "kg"},
    {"ATM",		ATM,		    "N/m^2"},
    {"AU",		AU,		    "m"},
    {"AVOGADRO",	6.0221367e+23,	    "1"},
    {"BOLTZMAN",	1.380658e-23,	    "J/K"},
    {"C",		CLIGHT,		    "m/s"},
    {"CATALAN",		0.9159655941772190150546035149323841107741, "1"},
    {"E",		2.7182818284590452353602874713526624977572, "1"},
    {"ELECTRONCHARGE",	1.60217733e-19,	    "C"},
    {"ELECTRONMASS",	9.1093897e-31,	    "kg"},
    {"ELECTRONVOLT",	1.60217733e-19,	    "J"},
    {"EULER",		0.5772156649015328606065120900824024310422, "1"},
    {"FARADAY",		9.6485309e+4,	    "C"},
    {"FINESTRUCTURE",	7.29735308e-3,	    "1"},
    {"GAMMA",		0.5772156649015328606065120900824024310422, "1"},
    {"GAS",		8.314510,	    "J/mole/K"},
    {"GRAVITY",		6.67259e-11,	    "N*m^2/kg^2"},
    {"HBAR",		1.05457266e-34,	    "J*s"},
    {"PERFECTGASVOLUME",2.241383e-2,	    "m^3/mole"},
    {"PI",		3.1415926535897932384626433832795028841972, "1"},
    {"PLANCK",		6.6260755e-34,	    "J*s"},
    {"PROTONMASS",	1.6726231e-27,	    "kg"},
    {"RYDBERG",		1.0973731534e+7,    "1/m"},
    {"SPEEDLIGHT",	CLIGHT,		    "m/s"},
    {"STANDARDGRAVITY",	STANDARD_GRAVITY,   "m/s^2"},
    {"STANDARDPRESSURE",ATM,		    "N/m^2"},
    {"STEFANBOLTZMANN",	5.67051e-8,	    "W/K^4/m^2"},
    {"WATERTRIPLE",	2.7316e+2,	    "K"}
};

#ifdef ANSI
Mfloat imsl_f_constant(Mchar *name, Mchar *units)
#else
Mfloat imsl_f_constant(name, units)
    Mchar   *name;
    Mchar   *units;
#endif 
{
    Const_table	    key;
    Mchar	    capnam[50];
    Mchar	    *from;
    Mchar	    *to;
    Const_table	    *node;
    Mfloat	    value;
#ifdef DOUBLE
	imsl_e1psh("imsl_d_constant");
#else
	imsl_e1psh("imsl_f_constant");
#endif
                                /* Convert NAME to upper case CAPNAM */
				/* deleting spaces and underscores */
    for (from=name, to=capnam; *from != '\0'; from++)
	if (islower(*from))
	    *to++ = imsl_toupper(*from);
	else if (*from != ' '  &&  *from != '_')
	    *to++ = *from;
    *to = '\0';
				/* Binary search of CTABLE for CAPNAM */
    key.name = capnam;
#if defined(COMPUTER_HP93C) || defined(COMPUTER_PMXUX) || defined(COMPUTER_ALFAC_IEEE)
    node = (Const_table*) bsearch((void*)&key, (void*)lv_const_table,
		sizeof(lv_const_table)/sizeof(Const_table),
		sizeof(Const_table), (int (*)())l_comp_const);
#else
#ifdef ANSI
    node = (Const_table*) bsearch((void*)&key, (void*)lv_const_table,
		sizeof(lv_const_table)/sizeof(Const_table),
		sizeof(Const_table), (int (*) (const void*, const void*))l_comp_const);
#else
    node = (Const_table*) bsearch((void*)&key, (void*)lv_const_table,
		sizeof(lv_const_table)/sizeof(Const_table),
		sizeof(Const_table), l_comp_const);
#endif
#endif

    if (node == NULL) {
	    /* The argument NAME = %(L1) is illegal. */
        imsl_e1stl (1, name);
        imsl_ermes (IMSL_TERMINAL, IMSL_BAD_CONST_NAME);
				/* Return not-a-number on error */
        value = imsl_amach(6);
        goto RETURN_CONSTANT;
    }
    value = node->value;

    if (units != NULL)
	l_cunit (value, node->units, &value, units);
RETURN_CONSTANT:
#ifdef DOUBLE
	imsl_e1pop("imsl_d_constant");
#else
	imsl_e1pop("imsl_f_constant");
#endif
    return value;
}


#ifdef ANSI
static Mint l_comp_const(Mvoid *t1, Mvoid *t2)
#else
static Mint l_comp_const(t1, t2)
    Mvoid    *t1;
    Mvoid    *t2;
#endif
{
    Const_table   *s1 = (Const_table*) t1;
    Const_table   *s2 = (Const_table*) t2;

    return strcmp(s1->name, s2->name);
}


/* Structured by FOR_STRUCT, v0.2, on 07/19/90 at 16:17:50
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  CUNIT/DCUNIT (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    March 6, 1990

    Purpose:    Convert X in units XUNITS to Y in units YUNITS.

    Usage:      CALL CUNIT (X, XUNITS, Y, YUNITS)

    Arguments:
       X      - Value to be converted.  (Input)
       XUNITS - Character string containing the name of the units for X.
                (Input)
                See remarks for a description of units allowed.
       Y      - Value in YUNITS corresponding to X in XUNITS.  (Output)
       YUNITS - Character string containing the name of the units for Y.
                (Input)
                See remarks for a description of units allowed.

    Remarks:
    1. Strings XUNITS and YUNITS have the form U1*U2*...*Um/V1/.../Vn,
       where Ui and Vi are the names of basic units or are the names
       of basic units raised to a power.  Examples are,
       'METER*KILOGRAM/SECOND', 'M*KG/S', 'METER', or 'M/KG^2'.

    2. The case of the character string in XUNITS and YUNITS does not
       matter.  The names 'METER', 'Meter' and 'meter' are equivalent.

    3. If XUNITS is 'SI', then X is assumed to be in the standard
       international units corresponding to YUNITS.  Similarly,
       if YUNITS is 'SI', then Y is assumed to be in the standard
       international units corresponding to XUNITS.

    4. The basic unit names allowed are as follows.
       Units of time
          day, hour = hr, min = minute, s = sec = second, year
       Units of frequency
          Hertz = Hz
       Units of mass
          AMU, g = gram, lb = pound, ounce = oz, slug
       Units of distance
          Angstrom, AU, feet = foot, in = inch, m = meter = metre,
          micron, mile, mill, parsec, yard
       Units of area
          acre
       Units of volume
          l = liter = litre
       Units of force
          dyne, N = Newton, poundal
       Units of energy
          BTU(thermochemical), Erg, J = Joule
       Units of work
          W = watt
       Units of pressure
          ATM = atomosphere, bar, Pascal
       Units of temperature
          degC = Celsius, degF = Fahrenheit, degK = Kelvin
       Units of viscosity
          poise, stoke
       Units of charge
          Abcoulomb, C = Coulomb, statcoulomb
       Units of current
          A = ampere, abampere, statampere,
       Units of voltage
          Abvolt, V = volt
       Units of magnetic induction
          T = Tesla, Wb = Weber
       Other units
          1, farad, mole, Gauss, Henry, Maxwell, Ohm
       The following metric prefixes may be used with the above units.
       Note that the one or two letter prefixes may only be used with
       one letter unit abbreviations.
          A  = atto  = 1.E-18
          F  = femto = 1.E-15
          P  = pico  = 1.E-12
          N  = nano  = 1.E-9
          U  = micro = 1.E-6
          M  = milli = 1.E-3
          C  = centi = 1.E-2
          D  = deci  = 1.E-1
          DK = deca  = 1.E+2
          K  = kilo  = 1.E+3
               myria = 1.E+4 (no single letter prefix; M means milli)
               mega  = 1.E+6 (no single letter prefix; M means milli)
          G  = giga  = 1.E+9
          T  = tera  = 1.E+12

    5. Informational error
       Type Code
         3   8  A conversion of units of mass to units of
                force was required for consistency.

    Keyword:    Utilities

    GAMS:       C19

    Chapter:    MATH/LIBRARY Utilities

    Copyright:  1986 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#define	NCLASS	5

#ifdef ANSI
static void l_cunit(Mfloat x, Mchar *xunits, Mfloat *y, Mchar *yunits)
#else
static void l_cunit(x, xunits, y, yunits)
	Mfloat          x;
	Mchar           *xunits;
	Mfloat          *y;
	Mchar           *yunits;
#endif
{
    Mint	    xsi, ysi;
    Mint            i, icntr, ipwnew[NCLASS], ipwold[NCLASS], j;
    Mfloat          cvtnew, cvtold;

    imsl_e1psh("l_cunit");
        /* Initialize y to NaN */
    *y = imsl_amach(6);

    /* Check for YUNITS = 'SI' */
    ysi = (l_iicsr(yunits, "SI") == 0);
    if (ysi) {
        cvtnew = F_ONE;
    } else {
        /* Convert YUNITS into SI system */
        l_c3nit(yunits, ipwnew, &cvtnew);
        if (imsl_n1rty(1) != 0)
            goto L_9000;
    }
    /* Check for XUNITS = 'SI' */
    xsi = (l_iicsr(xunits, "SI") == 0);
    if (xsi) {
        cvtold = F_ONE;
    } else {
        /* Convert XUNITS into SI system */
        l_c3nit(xunits, ipwold, &cvtold);
        if (imsl_n1rty(1) != 0)
            goto L_9000;
    }
    if (!(xsi || ysi)) {
        /* Initialize substitution counter */
        icntr = 0;
        /* Check compatibility */
        for (i = 1; i <= NCLASS; i++) {
            if (ipwold[i - 1] != ipwnew[i - 1]) {
                /*
                 * Substitute kg*m/s^2 for kg in XUNITS
                 */
                ipwold[0] += 1;
                ipwold[2] -= 2;
                /* Increment substitution counter */
                icntr += 1;
                /* Check compatibility */
        L_10:
                for (j = 1; j <= NCLASS; j++) {
                    if (ipwold[j - 1] != ipwnew[j - 1]) {
                        if (icntr == 2) {
			/* The units XUNITS = %(L1) and YUNITS  = %(L2) */
			/* are not compatible. */
                            imsl_e1stl(1, xunits);
                            imsl_e1stl(2, yunits);
                            imsl_ermes(5, IMSL_INCOMPATIBLE_UNITS);
                            goto L_9000;
                        }
                        /*
                         * Reset XUNITS as before the
                         * substitution
                         */
                        ipwold[0] -= 1;
                        ipwold[2] += 2;
                        /*
                         * Substitute kg*m/s^2 for
                         * kg in YUNITS
                         */
                        ipwnew[0] += 1;
                        ipwnew[2] -= 2;
                        /*
                         * Increment substitution
                         * counter
                         */
                        icntr += 1;
                        goto L_10;
                    }
                }
                /* Print warning message */
		    /* A conversion of units of mass to units of */
		    /* force was required for consistency.	 */
                imsl_ermes(IMSL_WARNING, IMSL_MASS_TO_FORCE);
                if (icntr == 1) {
                    *y = (STANDARD_GRAVITY * x * cvtold) / cvtnew;
                    goto L_9000;
                } else {
                    *y = (x * cvtold) / (STANDARD_GRAVITY * cvtnew);
                    goto L_9000;
                }
            }
        }
    }
    *y = (x * cvtnew) / cvtold;

L_9000:
    imsl_e1pop("l_cunit");
    return;
}                /* end of function */


        /* compare s1 and s2 with */
        /* regard for case      */
        /* return 0 if equal      */
#ifdef ANSI
static Mint l_iicsr(Mchar *s1, Mchar *s2)
#else
static Mint l_iicsr(s1, s2)
    Mchar   *s1;
    Mchar   *s2;
#endif
{
    Mchar   *p;
    Mchar   *q;

    for (p=s1, q=s2;  *p != '\0';  p++, q++) {
	if (*p == *q) continue;
	if ((abs(*p-*q) == 32) && isalpha(*p) && isalpha(*q)) continue;
	return 1;
    }
    return 0;
}


typedef struct {
    Mchar	*name;
    Mint	length;
    Mint	mass;
    Mint	time;
    Mint	current;
    Mint	temperature;
    Mfloat	convert;
} Conv_table;


/* Structured by FOR_STRUCT, v0.2, on 07/19/90 at 17:30:46
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  C3NIT/DC3NIT (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    March 6, 1990

    Purpose:    Return the value of conversion factor from UNIT to
                SI units.

    Usage:      CALL C3NIT (UNIT, IPOW, CVT)

    Arguments:
       UNIT   - Character string containing the name of the unit.
                (Input)
       IPOW   - Vector of length NCLASS containing the powers of the
                basic units (classes) in UNIT as follows.  (Output)
                     IPOW(1) = Length  (meters)
                     IPOW(2) = Mass    (kilograms)
                     IPOW(3) = Time    (seconds)
                     IPOW(4) = Current (amperes)
                     IPOW(5) = Temperature (Celcius)
       CVT    - Value of the convergence factor between a quantity in
                UNIT units to a quantity in SI units.  (Output)

    Chapter:    MATH/LIBRARY Utilities

    Copyright:  1986 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */

static Conv_table    lv_conv_table[] = {
    /* name	    l   m   t   i   T   factor */
    {"1",	    0,  0,  0,  0,  0,  1.0},
    {"A",	    0,  0,  0,  1,  0,  1.0},
    {"ABAMPERE",    0,  0,  0,  1,  0,  10.0},
    {"ABCOULOMB",   0,  0,  1,  1,  0,  10.0},
    {"ABVOLT",	    2,  1, -3, -1,  0,  1.0e-8},
    {"ACRE",	    2,  0,  0,  0,  0,  (5280.0*12.0*INCH)*(5280.0*12.0*INCH)
								    /640.0},
    {"AMPERE",	    0,  0,  0,  1,  0,  1.0},
    {"AMU",	    0,  1,  0,  0,  0,  1.0},
    {"ANGSTROM",    1,  0,  0,  0,  0,  1.0e-10},
    {"ATM",	   -1,  1, -2,  0,  0,  ATM},
    {"ATMOSPHERE", -1,  1, -2,  0,  0,  ATM},
    {"AU",	    1,  0,  0,  0,  0,  AU},
    {"BAR",	   -1,  1, -2,  0,  0,  1.0e+5},
    {"BTU",	    2,  1, -2,  0,  0,  1054.35},
    {"C",	    0,  0, -1,  1,  0,  1.0},
    {"CELSIUS",	    0,  0,  0,  0,  1,  1.0},
    {"COULOMB",	    0,  0,  1,  1,  0,  1.0},
    {"DAY",	    0,  0,  1,  0,  0,  3600.0*24.0},
    {"DEGC",	    0,  0,  0,  0,  1,  1.0},
    {"DEGF",	    0,  0,  0,  0,  1,  5.0/9.0},
    {"DEGK",	    0,  0,  0,  0,  1,  1.0},
    {"DYNE",	    1,  1, -2,  0,  0,  1.0e-5},
    {"ERG",	    2,  1, -2,  0,  0,  1.0e-7},
    {"FAHRENHEIT",  0,  0,  0,  0,  1,  5.0/9.0},
    {"FARAD",	   -2, -1,  4,  2,  0,  1.0},
    {"FEET",	    1,  0,  0,  0,  0,  12.0*INCH},
    {"FOOT",	    1,  0,  0,  0,  0,  12.0*INCH},
    {"G",	    0,  1,  0,  0,  0,  1.0e-3},
    {"GAUSS",	    1, -2, -1,  0,  0,  1.0e-4},
    {"GRAM",	    0,  1,  0,  0,  0,  1.0e-3},
    {"HENRY",	    2,  1, -2, -2,  0,  1.0},
    {"HERTZ",	    0,  0, -1,  0,  0,  1.0},
    {"HOUR",	    0,  0,  1,  0,  0,  3600.0},
    {"HR",	    0,  0,  1,  0,  0,  3600.0},
    {"HZ",	    0,  0, -1,  0,  0,  1.0},
    {"IN",	    1,  0,  0,  0,  0,  INCH},
    {"J",	    2,  1, -2,  0,  0,  1.0},
    {"JOULE",	    2,  1, -2,  0,  0,  1.0},
    {"KELVIN",	    0,  0,  0,  0,  1,  1.0},
    {"L",	    0,  0,  3,  0,  0,  1.0e-3},
    {"LB",	    0,  1,  0,  0,  0,  POUND},
    {"LITER",	    3,  0,  0,  0,  0,  1.0e-3},
    {"LITRE",	    3,  0,  0,  0,  0,  1.0e-3},
    {"M",	    1,  0,  0,  0,  0,  1.0},
    {"MAXWELL",	    2,  1, -2, -1,  0,  1.0e-8},
    {"METER",	    1,  0,  0,  0,  0,  1.0},
    {"METRE",	    1,  0,  0,  0,  0,  1.0},
    {"MICRON",	    1,  0,  0,  0,  0,  1.0e-6},
    {"MILE",	    1,  0,  0,  0,  0,  12.0*5280.0*INCH},
    {"MILL",	    1,  0,  0,  0,  0,  1.0e-3*INCH},
    {"MIN",	    0,  0,  1,  0,  0,  60.0},
    {"MINUTE",	    0,  0,  1,  0,  0,  60.0},
    {"MOLE",	    0,  0,  0,  0,  0,  1.0},
    {"N",	    1,  1, -2,  0,  0,  1.0},
    {"NEWTON",	    1,  1, -2,  0,  0,  1.0},
    {"OHM",	    2,  1, -3, -2,  0,  1.0},
    {"OUNCE",	    0,  1,  0,  0,  0,  POUND/16.0},
    {"OZ",	    0,  1,  0,  0,  0,  POUND/16.0},
				/* tan(2.0*PI/360.0/3600.0) = 0.48e-5 */
    {"PARSEC",	    1,  0,  0,  0,  0,  AU/0.48481368111333441675396429e-5}, 
    {"PASCAL",	   -1,  1, -2,  0,  0,  1.0},
    {"POINT",	    1,  0,  0,  0,  0,  0.013837*INCH},
    {"POISE",	   -1,  1, -1,  0,  0,  0.1},
    {"POUND",	    0,  1,  0,  0,  0,  POUND},
    {"POUNDAL",	    1,  1, -2,  0,  0,  0.138255},
    {"S",	    0,  0,  1,  0,  0,  1.0},
    {"S",	    0,  0,  1,  0,  0,  1.0},
    {"SEC",	    0,  0,  0,  0,  0,  1.0},
    {"SECOND",	    0,  0,  1,  0,  0,  1.0},
    {"SLUG",	    0,  1,  0,  0,  0,  1.45939E+1},
    {"STATAMPERE",  0,  0,  0,  1,  0,  0.1/CLIGHT},
    {"STATCOULOMB", 0,  0,  1,  1,  0,  0.1/CLIGHT},
    {"STOKE",	    2,  0, -1,  0,  0,  1.0e-4},
    {"T",	    0,  1, -2, -1,  0,  1.0},
    {"TESLA",	    0,  1, -2, -1,  0,  1.0},
    {"V",	    2,  1, -3, -1,  0,  1.0},
    {"VOLT",	    2,  1, -3, -1,  0,  1.0},
    {"W",	    2,  1, -3,  0,  0,  1.0},
    {"WATT",	    2,  1, -3,  0,  0,  1.0},
    {"WB",	    2,  1, -2, -1,  0,  1.0},
    {"WEBER",	    2,  1, -2, -1,  0,  1.0},
    {"YARD",	    1,  0,  0,  0,  0,  36.0*INCH},
    {"YEAR",	    0,  0,  1,  0,  0,  3.1536e+7},
};

typedef struct {
    Mchar   *name;
    char    length;
    Mfloat  value;
} Prefix_table;

#define ENTRY(NAME)	NAME, sizeof(NAME)-1
static Prefix_table  lv_prefix_table[] = {
    { ENTRY("ATTO"),	1.0e-18},
    { ENTRY("FEMTO"),	1.0e-15},
    { ENTRY("PICO"),	1.0e-12},
    { ENTRY("NANO") ,	1.0e-9},
    { ENTRY("MICRO"),	1.0e-6},
    { ENTRY("MILLI"),	1.0e-3},
    { ENTRY("CENTI"),	1.0e-2},
    { ENTRY("DECI"),	1.0e-1},
    { ENTRY("DECA"),	1.0e+2},
    { ENTRY("KILO"),	1.0e+3},
    { ENTRY("MYRIA"),   1.0e+4},
    { ENTRY("MEGA"),    1.0e+6},
    { ENTRY("GIGA"),	1.0e+9},
    { ENTRY("TERA"),	1.0e+12},
    { NULL, 0,          0.0},
};
#undef ENTRY


#ifdef ANSI
static void l_c3nit (Mchar *unit, Mint *ipow, Mfloat *cvt)
#else
static void l_c3nit (unit, ipow, cvt)
    Mchar	*unit;
    Mint	ipow[5];
    Mfloat	*cvt;
#endif
{
#define MAX_TERMS  10
    Mchar	    input[256];
    Mchar	    *token;
    Conv_table	    key;
    Conv_table	    *node;
    Mchar	    *terminator;
    Mchar	    *t;
    Mchar	    *s;
    Prefix_table    *prefix;
    Mint	    k;
    Mint	    n_terms;
    Mchar	    old_terminator;
    struct {
	Mchar	*name;
	Mint	power;
	Mfloat  factor;
    } term[MAX_TERMS];

    imsl_e1psh("l_c3nit");

			/* copy unit into input and convert to uppercase */
    for (t=input, s=unit;  *s != '\0';  t++, s++)
        *t = (isalpha(*s) && islower(*s)) ? (*s+'A'-'a') : *s;
    *t = '\0';

    ipow[0] = ipow[1] = ipow[2] = ipow[3] = ipow[4] = 0;
    *cvt   = F_ONE;
				/* parse input */
    for (k = 0;  k < MAX_TERMS;  k++) {
	term[k].power = 1;
	term[k].factor = F_ONE;
    }

    terminator = unit;
    old_terminator = '*';
    for (k = 0;  ; k++) {
	term[k].name = strtok((k==0) ? input : NULL, "*/^");
	if (term[k].name == NULL) break;
	terminator += strlen(term[k].name);
	if (old_terminator == '^') {
	    term[k-1].power *= atoi(term[k].name);
	    k--;
	} else if (old_terminator == '/') {
	    term[k-1].power *= -1;
	}
        old_terminator = terminator[0];
	terminator++;
    }
    n_terms = k;

    for (k = 0;  k < n_terms;  k++) {
	token = term[k].name;
        if (strlen(token) == 2) {
            Mchar           *p;
            Mint            exp;
            static Mchar    *list = "AFPNUMCK GTD";
            p = strchr(list, *token);
            if (p != NULL) {
                exp = -18 + 3*(p-list);
                if (exp == 0)  exp = -2;    /* C = centi */
                if (exp == 15) exp = -1;    /* D = deci */
		term[k].name++;
                term[k].factor *= imsl_fi_power(F_TEN, exp);
            }
        } else if (strlen(token)==3 && token[0]=='D' && token[1]=='K') {
	    term[k].name += 2;
	    term[k].factor *= 100.;	/* DK = deca = 100 */
	} else {
	    for (prefix=lv_prefix_table;  prefix->length != 0; prefix++) {
		if (strncmp(token, prefix->name, prefix->length) == 0) {
		    term[k].name   += prefix->length;
		    term[k].factor *= prefix->value;
		    break;
		}
	    }
        }
    }

    for (k = 0;  k < n_terms;  k++) {
	token = term[k].name;
				/* Binary search */
	key.name = token;
#if defined(COMPUTER_HP93C) || defined(COMPUTER_PMXUX) || defined(COMPUTER_ALFAC_IEEE)
	node = (Conv_table*) bsearch((void*)&key, (void*)lv_conv_table,
				     sizeof(lv_conv_table)/sizeof(Conv_table),
				     sizeof(Conv_table), 
				     (int (*)())l_comp_conv);
#else
#ifdef ANSI
	node = (Conv_table*) bsearch((void*)&key, (void*)lv_conv_table,
				     sizeof(lv_conv_table)/sizeof(Conv_table),
				     sizeof(Conv_table), 
			       	     (int (*) (const void*, const void*))l_comp_conv);
#else
	node = (Conv_table*) bsearch((void*)&key, (void*)lv_conv_table,
				     sizeof(lv_conv_table)/sizeof(Conv_table),
				     sizeof(Conv_table), l_comp_conv);
#endif
#endif
	if (node == NULL) {
		    /* The unit %(L1) is illegal. */
	    imsl_e1stl(1, token);
	    imsl_ermes(IMSL_TERMINAL, IMSL_ILLEGAL_UNIT);
	    goto RETURN;
	}
	term[k].factor *= node->convert;
	ipow[0] += term[k].power*node->length;
	ipow[1] += term[k].power*node->mass;
	ipow[2] += term[k].power*node->time;
	ipow[3] += term[k].power*node->current;
	ipow[4] += term[k].power*node->temperature;
	if (term[k].power != 1  &&  term[k].factor != F_ONE)
	    *cvt *= imsl_fi_power(term[k].factor,term[k].power);
	else
	    *cvt *= term[k].factor;
    }
RETURN:
    imsl_e1pop("l_c3nit");
}


#ifdef ANSI
static Mint l_comp_conv(Mvoid *t1, Mvoid *t2)
#else
static Mint l_comp_conv(t1, t2)
    Mvoid    *t1;
    Mvoid    *t2;
#endif
{
    Conv_table   *s1 = (Conv_table *)t1;
    Conv_table   *s2 = (Conv_table *)t2;
    return strcmp(s1->name, s2->name);
}
