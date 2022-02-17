/* ========================================================================== */
 
#ifndef		BGMTYPES_H
#define		BGMTYPES_H

/* Type of Vol */
typedef enum
{
	STRIKE,		 
	STDEV	
}STRIKE_TYPE;


/* Type of Vol */
typedef enum
{
	SLIDING,		
	FLAT	
}SMILE_MOVE_TYPE;

/* Type of Display for Fwd Vol Matrix */
typedef enum
{
	D_TS,		/* Std and default -> term structure of volatility at time t */
	D_CUM		/* Cumulative equivalent volatility for t up to the expery date */
}DISPLAY_TYPE;

#endif