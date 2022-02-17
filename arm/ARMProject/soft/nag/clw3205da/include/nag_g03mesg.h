/* <nag_g03mesg.h>
 *
 * Copyright 1997 Numerical Algorithms Group.
 *
 * Mark 5, 1997.
 *
 */

#ifndef NAG_G03MESG
#define NAG_G03MESG  

/* Output messages for Chapter g03 */

#define NM_G03_NO_MESG 0
#define NM_G03CAZ_IT_FUN 1
#define NM_G03CAZ_CRIT 2
#define NM_G03CAZ_HEAD 3
#define NM_G03CAZ_RES_1 4
#define NM_G03CAZ_RES_2 5
#define NM_G03CAZ_RES_3 6
#define NM_G03CAZ_RES_4 7
/* END OF DEFINES */


#ifdef NAG_MESG
char *nag_g03mesg[] =
{
": Dummy message for Chapter g03\n",

/* Result output */
"NM_G03CAZ_IT_FUN: \nIterations performed = %ld,   function evaluations = %ld\n",
"NM_G03CAZ_CRIT: Criterion = %16.6e\n",
"NM_G03CAZ_HEAD: \n%13sVariable    Standardized\n%25sCommunalities\n",
"NM_G03CAZ_RES_1:  %16ld%10s%7.4f\n",
"NM_G03CAZ_RES_2:  %16ld%10s%7.4f * at lower bound\n",
"NM_G03CAZ_RES_3:  %16ld%10s%7.4f * at upper bound\n",
"NM_G03CAZ_RES_4:  %16ld%10s%7.4f * held constant\n",
/* END OF MESSAGE STRINGS */

""
};
#else
extern char *nag_g03mesg[];
#endif

#endif  /* not NAG_E04MESG */








