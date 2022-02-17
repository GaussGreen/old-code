#ifndef NAG_G13MESG
#define NAG_G13MESG  

/* Output messages for Chapter g13 */

#define NM_NO_MESG 0
#define NM_ITC_S_D 1
#define NM_PHI_SER 2
#define NM_THETA_SER 3
#define NM_SPHI_SER 4
#define NM_STHETA_SER 5
#define NM_OMEGA_SI_SER 6
#define NM_OMEGA_SER 7
#define NM_DELTA_SER 8
#define NM_CONST_SER 9

#define NM_PHI 10
#define NM_THETA 11
#define NM_SPHI 12
#define NM_STHETA 13
#define NM_OMEGA_SI 14
#define NM_OMEGA 15
#define NM_DELTA_G13 16
#define NM_CONST_G13 17

#define NM_REQ_WA 18
#define NM_REQ_MWA 19


/* Option field list output */
#define NM_FUNCT_TITLE 20
#define NM_NUM_NSERIES 21
#define NM_CFIXED_G13BEC_VAL 22
#define NM_CFIXED_G13BJC_VAL 23
#define NM_LIKELIHOOD 24
#define NM_G13BEC_ALPHA 25
#define NM_G13BEC_BETA 26
#define NM_G13BEC_DELTA 27
#define NM_G13BEC_GAMMA 28
#define NM_G13BEC_PRINT_LEVEL 29
#define NM_G13BEC_OUTFILE 30

/* Final Solution */
#define NM_G13BEC_ITER 31
#define NM_G13BEC_SOLN_TITLE 32
#define NM_G13BEC_SOLN 33
#define NM_G13BEC_RESD 34
#define NM_G13BEC_FUNCT_VAL 35
#define NM_G13BEC_DF 36

#ifdef NAG_MESG
char *nag_g13mesg[] =
{
": Dummy message for Chapter g13\n",

/* Output when monitoring the course of the optimization */
"NM_ITC_S_D: Iter =%4ld     Residual =%15.6e     Objf =%15.6e\n",
"NM_PHI_SER: \nphi        series %3ld%15.6e\n",
"NM_THETA_SER: theta      series %3ld%15.6e\n",
"NM_SPHI_SER: sphi       series %3ld%15.6e\n",
"NM_STHETA_SER: stheta     series %3ld%15.6e\n",
"NM_OMEGA_SI_SER: omega/si   series %3ld%15.6e\n",
"NM_OMEGA_SER: omega      series %3ld%15.6e\n",
"NM_DELTA_SER: delta      series %3ld%15.6e\n",
"NM_CONST_SER: constant   series %3ld%15.6e\n\n",

"NM_PHI: \nphi                  %15.6e\n",
"NM_THETA: theta                %15.6e\n",
"NM_SPHI: sphi                 %15.6e\n",
"NM_STHETA: stheta               %15.6e\n",
"NM_OMEGA_SI: omega/si             %15.6e\n",
"NM_OMEGA: omega                %15.6e\n",
"NM_DELTA_G13: delta                %15.6e\n",
"NM_CONST_G13: constant             %15.6e\n\n",

"NM_REQ_WA: \n\nThe required minimum dimension of wa is %6ld\n",
"NM_REQ_MWA: \n\nThe required minimum dimension of mwa is %6ld\n",

/* Option field list output */
"NM_FUNCT_TITLE: \nParameters to %s\n____________________\n\n",
"NM_NUM_NSERIES: nseries...................... %3ld\n\n",
"NM_CFIXED_G13BEC_VAL:     cfixed.....................%s\n",
"NM_CFIXED_G13BJC_VAL: cfixed.....................%s\n",
"NM_LIKELIHOOD: criteria.....%s",
"NM_G13BEC_ALPHA: alpha...................%9.2e",
"NM_G13BEC_BETA:     beta....................%9.2e\n",
"NM_G13BEC_DELTA: delta...................%9.2e",
"NM_G13BEC_GAMMA:     gamma...................%9.2e\n",
"NM_G13BEC_PRINT_LEVEL: print_level..%s\n",
"NM_G13BEC_OUTFILE: outfile................ %9s\n\n",

/* Final Solution */
"NM_G13BEC_ITER: \n\nThe number of iterations carried out is %4ld\n\n",
"NM_G13BEC_SOLN_TITLE: The final values of the parameters and their standard\
 deviations are\n\n   i            para[i]                 sd\n",
"NM_G13BEC_SOLN: %4ld%20.6f%20.6f\n",
"NM_G13BEC_RESD: \n\nThe residual sum of squares =  %15.6e\n\n",
"NM_G13BEC_FUNCT_VAL: The objective function =  %15.6e\n\n",
"NM_G13BEC_DF: The degrees of freedom = %8.2f\n",
""
};
#else
extern char *nag_g13mesg[];
#endif

#endif  /* not NAG_G13MESG */
