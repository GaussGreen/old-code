#ifndef SRT_H_LOGTYPES_H
#define SRT_H_LOGTYPES_H

/*<%%STA-----------------------------------------------------------------
  TYPE            :SrtBasicTmInf 
  AUTHOR          :Ka Lok Chau
  DESCRIPTION     :information needed for the mdl at a particular time
  DEFINITION      :

<%%END---------------------------------------------------------------------*/
/*<%%STA*/
typedef struct{
  double df;                    /* B(0,t) */
}SrtBasicTmInf;


/*<%%STA-----------------------------------------------------------------
  TYPE            :SrtLogTmInf 
  AUTHOR          :Ka Lok Chau
  DESCRIPTION     :information needed for the mdl at a particular time
  DEFINITION      :

<%%END---------------------------------------------------------------------*/
/*<%%STA*/
typedef struct{
  double df;                    /* B(0,t) */
  double int_sig2_dt;           /* cum vol from node to next node */
  double sqrt_int_sig2_dt;      /* sqrt of cum vol from node to next node */
  double int_sig_dt ;		/* integral of sig*dt from node to next node */
  double inv_exp_int;

}SrtLogTmInf;


/*<%%STA-----------------------------------------------------------------
  TYPE            :SrtLogTreNdInf 
  AUTHOR          :O. Van Eyseren
  DESCRIPTION     :information about a particular node in a trin tree
  DEFINITION      :

<%%END---------------------------------------------------------------------*/

/*<%%STA*/
typedef struct {
	SrtSample 	cur_sam;
	SrtSample 	drift_sam;
   	double 	p[3];
	double 	df;  		/* df to next time point */  
	long	sons_index;
        double 	var_at_sam;
      } SrtLogTreNdInf;
/*<%%END*/


/*<%%STA-----------------------------------------------------------------
  TYPE            :SrtLogStpInf 
  AUTHOR          :E.Auld
  DESCRIPTION     :information about a particular time step in a log tree
  DEFINITION      :

<%%END---------------------------------------------------------------------*/

/*<%%STA*/
typedef struct{
	long max_x_index;
	long min_x_index;
	double logxmin; /*used in lognormal models */
	double xmin;    /*used in normal models */
	double u;	/*multiplicative or additive*/
	}SrtLogTreInf;
/*<%%END*/


/*<%%STA
Centering convention --- node closest in x or closest in log x
*/
#define CLOSESTINX 0
#define CLOSESTINLNX 1
/*<%%END*/

#endif
