#ifndef NAG_G02MESG
#define NAG_G02MESG  

#define NM_G02_MON_BANNER 0
#define NM_G02_IT_DELTA 1
#define NM_G02_J_THETA 2
#define NM_G02_RVAL 3
#define NM_G02_MAT_BANNER 4
#define NM_G02_THETA_BANNER 5
#define NM_G02_BLANK 6

/* defines for g02gbz */
#define NM_G02_IT_DEV 7
#define NM_G02_WT_EQU_SING 8
#define NM_G02_BETA_ARR 9

/*define for g02haw */
#define NM_G02_IT_MON_WT 10
#define NM_G02_IT_S 11
#define NM_G02_LAB_A 12
#define NM_G02_ROW_LAB 13
#define NM_G02_J 14
#define NM_G02_A 15

/*defines for g02hax */
#define NM_G02_IT_MON_THETA 16
#define NM_G02_IT_MON_THETA_LAB1 17
#define NM_G02_MON_NSJTR 18
#define NM_G02_MON_JTR 19

#ifdef NAG_MESG
char *nag_g02mesg[] =
{
"NM_G02_MON_BANNER:                     ** Iteration Monitoring **",
"NM_G02_IT_DELTA: Iteration %16ld  Max Delta = %13.5e",
"NM_G02_J_THETA:  %16ld    %13.5e",
"NM_G02_RVAL:    %13.5e",
"NM_G02_MAT_BANNER: Matrix %s",
"NM_G02_THETA_BANNER:                 I       theta(I)",
"NM_G02_BLANK: \n",

/* Messages for g02gbz */
"NM_G02_IT_DEV: Iteration %ld   Deviance = %13.5e\n",
"NM_G02_WT_EQU_SING: Weighted Least-square equations are singular\n",
"NM_G02_BETA_ARR: beta[%3ld] = %13.5e\n",

/* Messages for g02haw */
"NM_G02_IT_MON_WT: ** Iteration monitoring for weights **\n\n",
"NM_G02_IT_S: Iteration  %5ld  max(abs(s(i,j))) = %13.5e\n",
"NM_G02_LAB_A:        A\n",
"NM_G02_ROW_LAB: Row\n",
"NM_G02_J:  %2ld",
"NM_G02_A: %10.2e%c",

/* Messages for g02hax */
"NM_G02_IT_MON_THETA:\n                 ** Iteration monitoring for theta **\n\n",
"NM_G02_IT_MON_THETA_LAB1:  iteration       sigma       j       theta        \
rs\n\n",
"NM_G02_MON_NSJTR:    %5ld   %13.5e   %3ld  %13.5e  %13.5e\n",
"NM_G02_MON_JTR:                            %3ld  %13.5e  %13.5e\n", 
""
};
#else
extern char *nag_g02mesg[];
#endif

#endif  /* not NAG_G02MESG */
