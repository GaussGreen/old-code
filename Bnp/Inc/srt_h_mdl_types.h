#ifndef SRT_H_MDLTYPES_H
#define SRT_H_MDLTYPES_H

/********************************************************************
	srt_h_mdl_TYPES.H
	STRUCTURES FOR REPRESENTING
	BOTH FINANCIAL mdlS (IN PARTICULAR, INTEREST RATE mdlS),
	AND THE WAYS IN WHICH THEY ARE DISCRETIZED)
 ********************************************************************/

#define IS_TREE_CHEYETTE_2D(M) ( (int) (M) == CHEY )

#define IS_TREE_CHEYETTE_1D(M) ( (int) (M) == CHEY1D )

#define IS_TREE_SLICED(M) ( (int)(M)==CHEY )


/*******************************************************************/




/*<%%STA******************************************************************
  TYPE            :Sample 
  DESCRIPTION     :contains information representing a certain state
	at some point in the future for some model
  DEFINITION      :

<%%END
******************************************************************************/
/*<%%STA*/
#define SRTMAXSTOCK 2
typedef struct Sample{
      double short_rate;   /* the instantaneous short rate */
      double phi;          /* Cheyette's v^2=the cumulative volatility fn */
      double pv_money_mkt; /* 1/rolling money-mkt a/c in (A.31) */
      double s[SRTMAXSTOCK];/* value(s) of underlying(s) in 1(2) fact BS mdl */
                     } Sample;
/*<%%END*/



                                       
/********************************************************************/



/*** for communicating with a model **/

typedef enum SrtStrEleType{SRTDVAL,SRTIVAL,SRTSVAL,SRTEVAL}SrtStrEleType;
#define SRTSTRELEBUFSZ 100
typedef struct SrtStrEle
{
	SrtStrEleType type;
	char name[SRTSTRELEBUFSZ];
	char sval[SRTSTRELEBUFSZ];
	double dval;
	long ival;
	long offset;
	double lub;
	double upb;
}SrtStrEle;


/*****************************************************************/



#endif	
