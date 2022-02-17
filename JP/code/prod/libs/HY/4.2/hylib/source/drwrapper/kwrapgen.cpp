/*------------------------------------------------------------------
C FILE:         kwrapgen.c

CREATED BY:     Leda Braga - April 1996
MODIFIED :      Antonio Paras - Feb 1999 took away spot vol calibration
                                         added function to read FXVolatility File.


PURPOSE:        General data conversion functions to be used 
                in the Kapital DR wrapper in connection with
                ALIB based executable.
---------------------------------------------------------------------- */
#include <ctype.h>
#include <string.h>
#include <math.h>

#include "bastypes.h"
#include "cgeneral.h" 
#include "cfileio.h"
#include "cerror.h"  
#include "cmemory.h"
#include "macros.h"
#include "convert.h"
#include "dateconv.h" 
#include "strutil.h"
#include "ldate.h"
#include "bastypes.h"
#include "gtomat.h"
#include "interp.h"
#include "tcurve.h"
#include "ratelist.h"
#include "datelist.h"           /* GtoNewDateListFromDates     */   
#include "dbllist.h"            /* GtoNewDoubleList            */
#include "streamcf.h"           /* TStreamFixed, Float         */ 
#include "cashflow.h"           /* GtoNewEmptyCFL              */
#include "date_sup.h"           /* TDateInterval   */

#if ! defined(_WIN32)
#include "Legacy.h"
#endif
#include "eqsettle.h"           /* TEquitySettlement           */
#include "eqdiv.h"              /* TDividendList               */
#include "gtobf.h"              /* GtoDIVIDENDRATE_AMOUNT etc. */ 

#include "kwrapgen.h"           /* Prototype consistency       */
  
#define  IS_COMMENT_LINE(a) ((a)[0] IS '#')



 
    
    

/*----------------------------------------------------------------------------
FUNCTION:    DrKapWrapReadEquitySta

CREATED BY:  Neil Yang -  Febuary 2000

DESCRIPTION: Kapital wrapper function to read the equaity1.sta file to get divident info  
             
             Restrictions are applied concerning the format of the data
             file.
----------------------------------------------------------------------------*/
GTO_EXPORT (int) DrKapWrapReadEquitySta(char   *filename,       /*(I) Usually "zero.dat" */
										long   *numPoints,      /* (O) */
									    TDate  **dividentDates, /* (O) */
									    double **dividentRates) /* (O) */
{
    static char routine[]="DrKapWrapReadEquitySta";
    int        status = FAILURE;    
    
    FILE     *fp   = NULL;   /* File pointer               */
    
	TDate    *dDates = NULL;
	double   *dRates = NULL;
    double    fNumItems;     /* Number of items in zc      */  
	long       iNumItems;

	char      type[10];
    long       dateDR;        /* Date in YYYYMMDD format    */
    TDate      dateLIB;       /* Date in TDate format       */
    
    double     zRate;         /* Zero rate                  */
    
    int idx;
	char buffer[BUFFER_SIZE];
    
    /* Open data file   */
    fp = fopen(filename, "r");     
    if (fp IS NULL)
    {
        GtoErrMsg("%s: File %s cannot be opened.\n",  
                  routine,
                  filename);
        goto done;
    }

    /* find the first "#" indicating the start of the divident info */

	do 
	{
		if (fgets(buffer, BUFFER_SIZE , fp) == NULL)
		{
			GtoErrMsg("%s: unable to find the start of the divident info.\n",routine);
		}

	}
    while( !(IS_COMMENT_LINE(buffer)) );
    
    /* read in num of points */
    if (fscanf(fp,
                 "%lf \n", 
                 &fNumItems) != 1)
	{  
       GtoErrMsg("%s: Unable to read # of points in %s.\n", routine,
                 filename);
	}

	iNumItems = int(fNumItems);
	dDates = NEW_ARRAY(TDate,iNumItems);
	dRates = NEW_ARRAY(double,iNumItems);

	 /* Main loop */
    for (idx=0; idx<iNumItems; idx++)
    {
        if (fscanf(fp,
                      "%ld %lf %s\n", 
                      &dateDR,
                      &zRate,
					  type) != 3)
        {   
            GtoErrMsg("%s: Not enough divident points in file.\n",
                      routine);
            goto done;
        }
        
        if (DrDateToTDate(dateDR,
                          &dateLIB) IS FAILURE)
        {
            goto done;
        }  
        
        dDates[idx] = dateLIB;
        dRates[idx] = zRate/100.0;
    } 

	*numPoints = iNumItems;
	*dividentDates = dDates;
	*dividentRates = dRates;

	status = SUCCESS;

 done:
  
    if(fp ISNT NULL) fclose(fp);
    
    if (status IS FAILURE)
    {   
        free(dDates);
		free(dRates);
        
        return (status);
    }  
    else
    {
        return(status);
    }
	
}
/*----------------------------------------------------------------------------
FUNCTION:    DrKapWrapReadEquityDyn

CREATED BY:  Neil Yang -  Febuary 2000

DESCRIPTION: Kapital wrapper function to read the equaity1.dyn file to get value date,
			 price, and vol curve (i.e. TCurve format)  
             
             vol curve returned is checked for date order and for non-negative
             zero rates. 
             
             Restrictions are applied concerning the format of the data
             file.
----------------------------------------------------------------------------*/
GTO_EXPORT (int) DrKapWrapReadEquityDyn(char   *filename,   /*(I) Usually "zero.dat" */
									  TDate  valueDate,  /* (I) */
									  double *stockPrice, /* (O) */
									  double *correlation, /* (O) */
									  TCurve **volCurve)    /* (O) */
{
    static char routine[]="DrKapWrapReadEquityDyn";
    int        status = FAILURE;    
    
    FILE     *fp   = NULL;   /* File pointer               */
    
    TCurve    *vc   = NULL;   /* Zcurve to be returned      */
    int        fNumItems;     /* Number of items in zc      */  
   
    
    long       dateDR;        /* Date in YYYYMMDD format    */
    TDate      dateLIB;       /* Date in TDate format       */
    long       valueDateDR;   /* V date in YYYYMMDD format  */
    TDate      valueDateLIB;  /* V date in TDate format     */
    
    double     zRate;         /* Zero rate                  */
    
    int idx;

    
    /* Open data file   */
    fp = fopen(filename, "r");     
    if (fp IS NULL)
    {
        GtoErrMsg("%s: File %s cannot be opened.\n",  
                  routine,
                  filename);
        goto done;
    }

    /* File should contain value date, MM basis, an- */
    /* nual or semi-annual, swap basis, number of zc */
    /* points, and  zc dates+rates. One comment line */
    /* is  expected in  between  input sections, but */
    /* not partway through an array of data.         */
    
    /* 1 - Skip comment line, get + convert value date */
    if (DrFindAndSkipComLine(fp,
                              filename) IS FAILURE)
    {
        goto done;
    }
    if (fscanf(fp,
                 "%ld \n", 
                 &valueDateDR) != 1)
    {  
       GtoErrMsg("%s: Unable to read value date in %s.\n",
                 routine,
                 filename);
       goto done;
    }
    if (DrDateToTDate(valueDateDR,
					&valueDateLIB) IS FAILURE)
    {
        goto done;
    } 
 
//	*valueDate=valueDateLIB;
	/* 2 - Skip comment line and get equity price */
    if (DrFindAndSkipComLine(fp,
                              filename) IS FAILURE)
    {
        goto done;
    }
    if (fscanf(fp,
                 "%lf \n", 
                 stockPrice) != 1)
    {
       goto done;
    }

	/* 3 - Skip comment line and get correlation */
    if (DrFindAndSkipComLine(fp,
                              filename) IS FAILURE)
    {
        goto done;
    }
    if (fscanf(fp,
                 "%lf \n", 
                 correlation) != 1)
    {
       goto done;
    }


	/* 4 - Skip comment line and get number of zc items */
    if (DrFindAndSkipComLine(fp,
                              filename) IS FAILURE)
    {
        goto done;
    }
    if (fscanf(fp,
                 "%d \n", 
                 &fNumItems) != 1)
    {
       goto done;
    }
    
    /* 5 - Skip comment, make and populate zero curve */
    vc  = GtoNewTCurve(valueDate, 
                       fNumItems,
                       1.0,         
                       GTO_ACT_365F); 
    if (vc IS (TCurve *)NULL)  
    {
        goto done;
    }                      

    if (DrFindAndSkipComLine(fp,
                              filename) IS FAILURE)
    {
        goto done;
    }
    /* Main loop */
    for (idx=0; idx<fNumItems; idx++)
    {
        if (fscanf(fp,
                      "%ld %lf \n", 
                      &dateDR,
                      &zRate) != 2)
        {   
            GtoErrMsg("%s: Not enough zero curve points in file.\n",
                      routine);
            goto done;
        }
        
        if (DrDateToTDate(dateDR,
                          &dateLIB) IS FAILURE)
        {
            goto done;
        }  
        
        vc->fArray[idx].fDate = dateLIB;
        vc->fArray[idx].fRate = zRate/100.0;
    } 
    
    if (DrKapWrapCheckTCurve(vc,
                              "Zero curve") IS FAILURE)
    {
        goto done;
    }


	*volCurve = vc;   
    status = SUCCESS;  
    
  done:
  
    if(fp ISNT NULL) fclose(fp);
    
    if (status IS FAILURE)
    {   
        GtoFreeTCurve(vc);
        GtoErrMsg("%s: Failed.\n",routine); 
        return (status);
    }  
    else
    {
        return(status);
    }
}
  



/*----------------------------------------------------------------------------
FUNCTION:    DrKapWrapGetZC

CREATED BY:  Leda Braga -  April 1996

DESCRIPTION: Kapital wrapper function to read the zero.dat file (in DR
             format) into a Library zero curve (i.e. TCurve format)  
             
             ZC returned is checked for date order and for non-negative
             zero rates. 
             
             Restrictions are applied concerning the format of the data
             file.
----------------------------------------------------------------------------*/
GTO_EXPORT(TCurve *)DrKapWrapGetZC(char  *filename) /*(I) Usually "zero.dat" */
{
    static char routine[]="DrKapWrapGetZC";
    int        status = FAILURE;    
    
    FILE     *fp   = NULL;   /* File pointer               */
    
    TCurve    *zc   = NULL;   /* Zcurve to be returned      */
    int        fNumItems;     /* Number of items in zc      */  
    int        mmBasis;
    char       swapBasis[BUFFER_SIZE];
    char       annualSemi;
    
    long       dateDR;        /* Date in YYYYMMDD format    */
    TDate      dateLIB;       /* Date in TDate format       */
    long       valueDateDR;   /* V date in YYYYMMDD format  */
    TDate      valueDateLIB;  /* V date in TDate format     */
    
    double     zRate;         /* Zero rate                  */
    
    int idx;

    
    /* Open data file   */
    fp = fopen(filename, "r");     
    if (fp IS NULL)
    {
        GtoErrMsg("%s: File %s cannot be opened.\n",  
                  routine,
                  filename);
        goto done;
    }

    /* File should contain value date, MM basis, an- */
    /* nual or semi-annual, swap basis, number of zc */
    /* points, and  zc dates+rates. One comment line */
    /* is  expected in  between  input sections, but */
    /* not partway through an array of data.         */
    
    /* 1 - Skip comment line, get + convert value date */
    if (DrFindAndSkipComLine(fp,
                              filename) IS FAILURE)
    {
        goto done;
    }
    if (fscanf(fp,
                 "%ld \n", 
                 &valueDateDR) != 1)
    {  
       GtoErrMsg("%s: Unable to read value date in %s.\n",
                 routine,
                 filename);
       goto done;
    }
    if (DrDateToTDate(valueDateDR,
                      &valueDateLIB) IS FAILURE)
    {
        goto done;
    } 
    
    /* 2 - Skip comment line and get money market basis */
    if (DrFindAndSkipComLine(fp,
                             filename) IS FAILURE)
    {
        goto done;
    }
    if (fscanf(fp,
                 "%d \n", 
                 &mmBasis) != 1)
    {
       goto done;
    }  
    if ((mmBasis ISNT 360) && (mmBasis ISNT 365)) 
    {
        GtoErrMsg("%s: Unrecognised MM basis.\n", routine);
        goto done;
    }
    
    /* 3 - Skip comment line and get annual/semiannual */
    if (DrFindAndSkipComLine(fp,
                             filename) IS FAILURE)
    {
        goto done;
    }
    if (fscanf(fp,
                 "%c \n", 
                 &annualSemi) != 1)
    {
       goto done;
    } 
    annualSemi = (char)toupper(annualSemi);
    if((annualSemi ISNT 'A') && (annualSemi ISNT 'S')) 
    {
        GtoErrMsg("%s: Unrecognised char for annual/semi.\n",routine);
        goto done;
    }
    
    /* 4 - Skip comment line and get swap basis */
    if (DrFindAndSkipComLine(fp,
                              filename) IS FAILURE)
    {
        goto done;
    }
    if (fscanf(fp,
                 "%s \n", 
                 swapBasis) != 1)
    {
       goto done;
    }
    if ( (strcmp(swapBasis,"ACT")) &&
         (strcmp(swapBasis,"360")) &&
         (strcmp(swapBasis,"365")) )
    {
        GtoErrMsg("%s: Unrecognised swap basis.\n", routine);
        goto done;
    }
    
    /* 5 - Skip comment line and get number of zc items */
    if (DrFindAndSkipComLine(fp,
                              filename) IS FAILURE)
    {
        goto done;
    }
    if (fscanf(fp,
                 "%d \n", 
                 &fNumItems) != 1)
    {
       goto done;
    }
    
    /* 6 - Skip comment, make and populate zero curve */
    zc  = GtoNewTCurve(valueDateLIB, 
                       fNumItems,
                       1.0,         
                       GTO_ACT_365F); 
    if (zc IS (TCurve *)NULL)  
    {
        goto done;
    }                      

    if (DrFindAndSkipComLine(fp,
                              filename) IS FAILURE)
    {
        goto done;
    }
    /* Main loop */
    for (idx=0; idx<fNumItems; idx++)
    {
        if (fscanf(fp,
                      "%ld %lf \n", 
                      &dateDR,
                      &zRate) != 2)
        {   
            GtoErrMsg("%s: Not enough zero curve points in file.\n",
                      routine);
            goto done;
        }
        
        if (DrDateToTDate(dateDR,
                          &dateLIB) IS FAILURE)
        {
            goto done;
        }  
        
        zc->fArray[idx].fDate = dateLIB;
        zc->fArray[idx].fRate = zRate/100.0;
    } 
    
    if (DrKapWrapCheckTCurve(zc,
                              "Zero curve") IS FAILURE)
    {
        goto done;
    }
   
    status = SUCCESS;  
    
  done:
  
    if(fp ISNT NULL) fclose(fp);
    
    if (status IS FAILURE)
    {   
        GtoFreeTCurve(zc);
        GtoErrMsg("%s: Failed.\n",routine); 
        return (NULL);
    }  
    else
    {
        return(zc);
    }
}
  
 


/*----------------------------------------------------------------------------
FUNCTION:    DrKapWrapGetFXVol

CREATED BY:  Antonio Paras -- Feb 1999 based on DrKapWrapGetVC

DESCRIPTION: Kapital wrapper function to read the fxvolatility.dat file (in 
             DR format) into a Library vol curve (i.e. TCurve format) and return
             the fxspot value.
             
             Restrictions are applied concerning the format of the data
             file.
----------------------------------------------------------------------------*/
GTO_EXPORT(TCurve *) DrKapWrapGetFXVol
                     (char      *filename,   /*(I) Usually "FXVolatility.dat"     */
                      double    *FXSpotRate) /*(O) FX Spot Rate  */
{
    static char routine[]="DrKapWrapGetFXVol";
    int        status = FAILURE;    
        
    FILE      *fp   = NULL;   /* File pointer               */
    
    TDate      vcValueDate;   /* Value Date */
    TCurve    *vc   = NULL;   /* Vcurve to be returned      */
    int        fNumItems;     /* Number of items in zc      */ 
    char       freqChar;      /* Freq: A, S, Q, M           */  
    double     freqDouble;    /* As above, but in double    */
    
    long       dateDR;        /* Date in YYYYMMDD format    */
    TDate      dateLIB;       /* Date in TDate format       */
    
    double     bvol;          /* Base volatility            */
    
    int idx;
    
    /* Open data file   */
    fp = fopen(filename, "r");     
    if (fp IS NULL)
    {
        GtoErrMsg("%s: File %s cannot be opened.\n",  
                  routine,
                  filename);
        goto done;
    }

    /* File  should  contain the value date, FX Spot */
    /* rate, base vol frequency,                     */
    /* the number of base vols and bvols dates+rates.*/
    /* One comment line is expected in between input */
    /* sections, but not partway through an array.   */
    
    /* 1 - Skip comment line and get value Date*/
    if (DrFindAndSkipComLine(fp,
                              filename) IS FAILURE)
    {
        goto done;
    }

    if (fscanf(fp,
                  "%ld\n", 
                  &dateDR) != 1)
    {   
        GtoErrMsg("%s: Value Date not present.\n",
                  routine);
        goto done;
    }
    
    if (DrDateToTDate(dateDR,
                      &vcValueDate) IS FAILURE)
    {
        goto done;
    }  

    /* 2 - Skip comment line and get FX Spot Rate*/
    if (DrFindAndSkipComLine(fp,
                              filename) IS FAILURE)
    {
        goto done;
    }

    if (fscanf(fp,
                  "%lf\n", 
                  FXSpotRate) != 1)
    {   
        GtoErrMsg("%s: Couldn't read FX Spot Rate.\n",
                  routine);
        goto done;
    }
    
    /* 3 - Skip comment line and get+convert bvol freq*/
    if (DrFindAndSkipComLine(fp,
                              filename) IS FAILURE)
    {
        goto done;
    }
    if (fscanf(fp,
                 "%c \n", 
                 &freqChar) != 1)
    {
       goto done;
    }
    
    switch ((char)toupper(freqChar))
    {
        case 'A':
            freqDouble = 1.0;
            break;
        case 'S':
            freqDouble = 2.0;
            break;
        case 'Q':
            freqDouble = 4.0;
            break;
        case 'M':
            freqDouble = 12.0;
            break; 
        default:
            GtoErrMsg("%s: Unrecognised vol frequency (%c).\n",
                      routine,
                      freqChar);
            goto done;  
    }
        
    /* 4 - Skip comment line and get number of vc items */
    if (DrFindAndSkipComLine(fp,
                              filename) IS FAILURE)
    {
        goto done;
    }
    if (fscanf(fp,
                 "%d \n", 
                 &fNumItems) != 1)
    {
       goto done;
    }
    
    /* 5 - Skip comment, make and populate vol curve */
    vc  = GtoNewTCurve(vcValueDate, 
                       fNumItems,
                       freqDouble,      
                       GTO_ACT_365F); 
    if (vc IS (TCurve *)NULL)  
    {
        goto done;
    }                      

    if (DrFindAndSkipComLine(fp,
                              filename) IS FAILURE)
    {
        goto done;
    }
    /* Main loop */
    for (idx=0; idx<fNumItems; idx++)
    {
        if (fscanf(fp,
                      "%ld %lf \n", 
                      &dateDR,
                      &bvol) != 2)
        {   
            GtoErrMsg("%s: Not enough basevol points in file.\n",
                      routine);
            goto done;
        }
        
        if (DrDateToTDate(dateDR,
                          &dateLIB) IS FAILURE)
        {
            goto done;
        }  
        
        vc->fArray[idx].fDate = dateLIB;
        vc->fArray[idx].fRate = bvol/100.0;
    }
    
    if (DrKapWrapCheckTCurve(vc,
                              "Basevol curve") IS FAILURE)
    {
        goto done;
    }
    
    status = SUCCESS;  
    
  done:
  
    if(fp ISNT NULL) fclose(fp);
    
    if (status IS FAILURE)
    {   
        GtoFreeTCurve(vc);
        GtoErrMsg("%s: Failed.\n",routine);
        return(NULL); 
    }
    else
    {                          
        return (vc);
    }
}
    
    
    
    
    
/*----------------------------------------------------------------------------
FUNCTION:    DrKapWrapGetVC

CREATED BY:  Leda Braga -  May 1996

DESCRIPTION: Kapital wrapper function to read the basevol.dat file (in 
             DR format) into a Library vol curve (i.e. TCurve format)
             
             Since a value date is not given in the DR format basevol.
             dat, this must be provided externally.   
             
             Restrictions are applied concerning the format of the data
             file.
----------------------------------------------------------------------------*/
GTO_EXPORT(TCurve *) DrKapWrapGetVC
                     (char  *filename,      /*(I) Usually "basevol.dat"     */
                      TDate  vcValueDate)   /*(I) Not given in the file!    */

{
    static char routine[]="DrKapWrapGetVC";
    int        status = FAILURE;    
    
    FILE     *fp   = NULL;   /* File pointer               */
    
    TCurve    *vc   = NULL;   /* Vcurve to be returned      */
    int        fNumItems;     /* Number of items in zc      */ 
    char       freqChar;      /* Freq: A, S, Q, M           */  
    double     freqDouble;    /* As above, but in double    */
    
    long       dateDR;        /* Date in YYYYMMDD format    */
    TDate      dateLIB;       /* Date in TDate format       */
    
    double     bvol;          /* Base volatility            */
    
    int idx;
    
    /* Open data file   */
    fp = fopen(filename, "r");     
    if (fp IS NULL)
    {
        GtoErrMsg("%s: File %s cannot be opened.\n",  
                  routine,
                  filename);
        goto done;
    }

    /* File  should  contain the base vol frequency, */
    /* the number of base vols and bvols dates+rates.*/
    /* One comment line is expected in between input */
    /* sections, but not partway through an array.   */
    
    /* 1 - Skip comment line and get+convert bvol freq*/
    if (DrFindAndSkipComLine(fp,
                              filename) IS FAILURE)
    {
        goto done;
    }
    if (fscanf(fp,
                 "%c \n", 
                 &freqChar) != 1)
    {
       goto done;
    }
    
    switch ((char)toupper(freqChar))
    {
        case 'A':
            freqDouble = 1.0;
            break;
        case 'S':
            freqDouble = 2.0;
            break;
        case 'Q':
            freqDouble = 4.0;
            break;
        case 'M':
            freqDouble = 12.0;
            break; 
        default:
            GtoErrMsg("%s: Unrecognised vol frequency (%c).\n",
                      routine,
                      freqChar);
            goto done;  
    }
        
    /* 2 - Skip comment line and get number of vc items */
    if (DrFindAndSkipComLine(fp,
                              filename) IS FAILURE)
    {
        goto done;
    }
    if (fscanf(fp,
                 "%d \n", 
                 &fNumItems) != 1)
    {
       goto done;
    }
    
    /* 3 - Skip comment, make and populate vol curve */
    vc  = GtoNewTCurve(vcValueDate, 
                       fNumItems,
                       freqDouble,      
                       GTO_ACT_365F); 
    if (vc IS (TCurve *)NULL)  
    {
        goto done;
    }                      

    if (DrFindAndSkipComLine(fp,
                              filename) IS FAILURE)
    {
        goto done;
    }
    /* Main loop */
    for (idx=0; idx<fNumItems; idx++)
    {
        if (fscanf(fp,
                      "%ld %lf \n", 
                      &dateDR,
                      &bvol) != 2)
        {   
            GtoErrMsg("%s: Not enough basevol points in file.\n",
                      routine);
            goto done;
        }
        
        if (DrDateToTDate(dateDR,
                          &dateLIB) IS FAILURE)
        {
            goto done;
        }  
        
        vc->fArray[idx].fDate = dateLIB;
        vc->fArray[idx].fRate = bvol/100.0;
    }
    
    if (DrKapWrapCheckTCurve(vc,
                              "Basevol curve") IS FAILURE)
    {
        goto done;
    }
    
    status = SUCCESS;  
    
  done:
  
    if(fp ISNT NULL) fclose(fp);
    
    if (status IS FAILURE)
    {   
        GtoFreeTCurve(vc);
        GtoErrMsg("%s: Failed.\n",routine);
        return(NULL); 
    }
    else
    {                          
        return (vc);
    }
}
    
    
    
    
    
/*----------------------------------------------------------------------------
FUNCTION:    DrKapWrapGetSwpvolGrid

CREATED BY:  Leda Braga -  May 1996

DESCRIPTION: Kapital wrapper function to read the swapvol.dat file(in 
             DR format) into an  array of doubles.  The final grid as 
             exported by this function has vols decimal (rather  than
             percent) format, i.e. 0.10 for 10% vol.   

             Restrictions are applied concerning the format of the data
             file.
----------------------------------------------------------------------------*/

GTO_EXPORT(int ) DrKapWrapGetSwpvolGrid
               (char           *filename,       /*(I) Usually "swapvol.dat" */
                int            *numSwpvolExp,   /*(O) Number of expiries    */
                int            *numSwpvolMat,   /*(O) Number of maturities  */
                double        **swpvolExpPtr,   /*(O) Array of exp in mths  */
                double        **swpvolMatPtr,   /*(O) Array of mats in yrs  */
                double        **swpvolGridPtr)  /*(O) Swapvol grid(not in %)*/
{       
    static char routine[] = "DrKapWrapGetSwpvolGrid";
    int status = FAILURE;
    
    
    FILE     *fp   = NULL;       /* File pointer               */
    
    double    *swpvolExpLocal  = NULL;  /* Local vars to ensure */
    double    *swpvolMatLocal  = NULL;  /* protection  of  the  */
    double    *swpvolGridLocal = NULL;  /* pointers  passed in  */
    
    int numExp;    /* Nb of expiries in grid   */
    int numMat;    /* Nb of maturities in grid */
    
    int idx , jdx;
    
    /* Open data file   */
    fp = fopen(filename, "r");     
    if (fp IS NULL)
    {
        GtoErrMsg("%s: File %s cannot be opened.\n",  
                  routine,
                  filename);
        goto done;
    }

    
    /* 1 - Skip comment line and number of expiries */
    if (DrFindAndSkipComLine(fp,
                              filename) IS FAILURE)
    {
        goto done;
    }
    if (fscanf(fp,
                 "%d \n", 
                 &numExp) != 1)
    {
        goto done;   
    }
    /* 2 - Skip comment line and number of maturities */
    if (DrFindAndSkipComLine(fp,
                              filename) IS FAILURE)
    {
        goto done;
    }
    if (fscanf(fp,
                 "%d \n", 
                 &numMat) != 1)
    {
        goto done;   
    }   
    
    /* 3 - Size is known, allocate memory */
    swpvolExpLocal  = NEW_ARRAY(double,(numExp));
    if (swpvolExpLocal IS NULL)
    {
        goto done;
    }
    swpvolMatLocal  = NEW_ARRAY(double,(numMat));
    if (swpvolMatLocal IS NULL)
    {
        goto done;
    }
    swpvolGridLocal = NEW_ARRAY(double,numExp*numMat);
    if (swpvolGridLocal IS NULL)
    {
        goto done;
    }
    
    /* 4 - Skip comment line and read arrays */
    if (DrFindAndSkipComLine(fp,
                              filename) IS FAILURE)
    {
        goto done;
    }  
    for (idx=0; idx<numMat; idx++)
    {
        if (fscanf(fp,
                      "%lf", 
                      &(swpvolMatLocal[idx])) != 1)
        {
            goto done;   
        }
    }    
    fscanf(fp,"\n"); 
    
    for (jdx=0; jdx<numExp; jdx++)
    {
        if (fscanf(fp,
                      "%lf", 
                      &(swpvolExpLocal[jdx])) != 1)
        {
            goto done;   
        } 
    
        for (idx=0; idx<numMat; idx++)
        {
            if (fscanf(fp,
                         "%lf", 
                          &(swpvolGridLocal[jdx*numMat+idx])) != 1)
            {
                goto done;   
            }
            swpvolGridLocal[jdx*numMat+idx] /= 100.0;
        } 
        fscanf(fp,"\n"); /* Move to next line */
    }
    
    /* All set so, before we can be done, pass the pointers back to caller */
    *swpvolExpPtr  = swpvolExpLocal;
    *swpvolMatPtr  = swpvolMatLocal;  
    *swpvolGridPtr = swpvolGridLocal;
    *numSwpvolExp  = numExp;
    *numSwpvolMat  = numMat;
    
    status = SUCCESS;
    
  done:  
  
    if(fp ISNT NULL) fclose(fp);
  
    if (status IS FAILURE)
    {
        FREE_ARRAY(swpvolExpLocal);
        FREE_ARRAY(swpvolMatLocal);
        FREE_ARRAY(swpvolGridLocal);
        GtoErrMsg("%s: Failed.\n",routine);
    }
    return(status);


}  /* End of DrKapWrapGetSwpvolGrid() */


/*----------------------------------------------------------------------------
FUNCTION:    DrKapWrapGetTwoFacPars

CREATED BY:  Leda Braga -  October 1996

DESCRIPTION: Kapital wrapper function to read the interest rate two fac-
             tor model parameters.

             Restrictions are applied concerning the format of the data
             file.

             Returns: SUCCESS or FAILURE
----------------------------------------------------------------------------*/
GTO_EXPORT(int ) DrKapWrapGetTwoFacPars
                (char     *filename,    /* (I) "2factormeanreversion.dat"   */
                 double   *meanRev1,    /* (O) Mean rev of 1st factor       */
                 double   *meanRev2,    /* (O) Mean rev of 2nd factor       */
                 double   *correl,      /* (O) Correlation between factors  */
                 double   *weight)      /* (O) Relative weight of factors   */


{
    static char routine[]="DrKapWrapGetTwoFacPars";
    int        status = FAILURE;    
    
    FILE     *fp   = NULL;   /* File pointer               */
    
   
    
    /* Open data file   */
    fp = fopen(filename, "r");     
    if (fp IS NULL)
    {
        GtoErrMsg("%s: File %s cannot be opened.\n",  
                  routine,
                  filename);
        goto done;
    }

    
    /*  Mean reversion 1 */
    if (DrFindAndSkipComLine(fp,
                              filename) IS FAILURE)
    {
        goto done;
    }
    if (fscanf(fp,
                 "%lf \n", 
                 meanRev1) != 1)
    {  
       GtoErrMsg("%s: Unable to read first mean rev coefficient in %s.\n",
                 routine,
                 filename);
       goto done;
    }



    /* Mean reversion 2 */
    if (DrFindAndSkipComLine(fp,
                              filename) IS FAILURE)
    {
        goto done;
    }
    if (fscanf(fp,
                 "%lf \n", 
                 meanRev2) != 1)
    {  
       GtoErrMsg("%s: Unable to read second mean rev coefficient in %s.\n",
                 routine,
                 filename);
       goto done;
    }


    /* Correlation */
    if (DrFindAndSkipComLine(fp,
                              filename) IS FAILURE)
    {
        goto done;
    }
    if (fscanf(fp,
                 "%lf \n", 
                 correl) != 1)
    {  
       GtoErrMsg("%s: Unable to read correlation in %s.\n",
                 routine,
                 filename);
       goto done;
    }


    /* Relative weight */
    if (DrFindAndSkipComLine(fp,
                              filename) IS FAILURE)
    {
        goto done;
    }
    if (fscanf(fp,
                 "%lf \n", 
                 weight) != 1)
    {  
       GtoErrMsg("%s: Unable to read relative weight in %s.\n",
                 routine,
                 filename);
       goto done;
    }


      
    status = SUCCESS;  
    
  done:
  
    if(fp ISNT NULL) fclose(fp);
    
    if (status IS FAILURE)
    {   
        GtoErrMsg("%s: Failed.\n",routine); 
    }  
    return(status);
 
}
  


/*----------------------------------------------------------------------------
FUNCTION:    DrKapWrapGetOneFacPars

CREATED BY:  Leda Braga -  October 1996

DESCRIPTION: Kapital wrapper function to read the interest rate one fac-
             tor model parameters.

             Restrictions are applied concerning the format of the data
             file.

             Returns: SUCCESS or FAILURE
----------------------------------------------------------------------------*/
GTO_EXPORT(int ) DrKapWrapGetOneFacPars
                (char     *filename,    /* (I) Usually "meanreversion.dat"  */
                 double   *meanRev1)    /* (O) Mean rev of 1st factor       */
                 
{
    static char routine[]="DrKapWrapGetOneFacPars";
    int        status = FAILURE;    
    
    FILE     *fp   = NULL;   /* File pointer               */
    
   
    
    /* Open data file   */
    fp = fopen(filename, "r");     
    if (fp IS NULL)
    {
        GtoErrMsg("%s: File %s cannot be opened.\n",  
                  routine,
                  filename);
        goto done;
    }

    
    /*  Mean reversion coefficient */
    if (DrFindAndSkipComLine(fp,
                              filename) IS FAILURE)
    {
        goto done;
    }
    if (fscanf(fp,
                 "%lf \n", 
                 meanRev1) != 1)
    {  
       GtoErrMsg("%s: Unable to read first mean rev coefficient in %s.\n",
                 routine,
                 filename);
       goto done;
    }

      
    status = SUCCESS;  
    
  done:
  
    if(fp ISNT NULL) fclose(fp);
    
    if (status IS FAILURE)
    {   
        GtoErrMsg("%s: Failed.\n",routine); 
    }  
    return(status);
 
}   /* End of DrKapWrapGetOneFacPars() */
  
 
   
                                                                              
/*----------------------------------------------------------------------------
FUNCTION:    DrDateToTDate

CREATED BY:  Leda Braga -  May 1996

DESCRIPTION: Kapital wrapper convenience function to convert a DR format
             date to a TDate.
----------------------------------------------------------------------------*/
GTO_EXPORT(int ) DrDateToTDate(long   DRDate,
                             TDate *libDate)

{
    static char routine[] = "DrDateToTDate";
    
    TMonthDayYear mdy;
    
    mdy.year  =  DRDate/10000; /* Truncation forced */
    mdy.month = (DRDate - 10000*mdy.year)/100;
    mdy.day   = (DRDate - 10000*mdy.year - 100*mdy.month);
    
    if (GtoMDYToDate(&mdy,
                     libDate) IS FAILURE)
    {
        GtoErrMsg("%s: Failed.\n",routine);
        return(FAILURE);
    }
    
    return(SUCCESS);
    
}       


/*----------------------------------------------------------------------------
FUNCTION:    DrFindAndSkipComLine

CREATED BY:  Leda Braga -  May 1996

DESCRIPTION: Kapital wrapper convenience function which finds and skips 
             a comment line  and flags an error if line is not found or 
             is not a comment.
 
                            FILE MUST BE OPEN!
----------------------------------------------------------------------------*/
GTO_EXPORT(int ) DrFindAndSkipComLine(FILE *fp,
                                  char  *filename)  
{
    static char routine[] = "DrFindAndSkipComLine";
    
    char buffer[BUFFER_SIZE];
    
    if (fgets(buffer, BUFFER_SIZE , fp) != NULL)
    {
        if (!(IS_COMMENT_LINE(buffer)))
        {
            buffer[10] = 0;  /* print first 10 chars */
            GtoErrMsg("%s: Comment line expected in %s.'%s ...' \n",
                      routine,
                      filename,
                      buffer);
            return(FAILURE);
        }
    }
    else
    {
        GtoErrMsg("%s: Line expected and not found in %s.\n",
                  routine,
                  filename);
        return(FAILURE);
    }    
    
    return(SUCCESS);
    
}     




/*----------------------------------------------------------------------------
FUNCTION:    DRLabelToTDateInterval

CREATED BY:  Leda Braga -  May 1996

DESCRIPTION: Kapital  wrapper convenience to convert a DR date interval to
             a TDateInterval. It allows for "2Y" and "2A" as valid strings

----------------------------------------------------------------------------*/
GTO_EXPORT(int ) DRLabelToTDateInterval
                (char             *matLabel,
                 TDateInterval    *matInterval)
                 
{
    static char routine[] = "DRLabelToTDateInterval";
    int status = FAILURE;
    
    int length;    /* Length of string */
    
    /* Convert from Y to A since this is the only exception */
    /* format which is not listed in the Lib standard       */
    
    length = strlen(matLabel);
    if(toupper(matLabel[length-1]) IS 'Y')
    { 
        matLabel[length-1] = 'A';
    }
    
    if (GtoStringToDateInterval(matLabel,
                                "Maturity of eq vc",
                                matInterval) IS FAILURE)
    {
        goto done;
    }
     
    status = SUCCESS;
    
  done:
    if (status IS FAILURE)
    {
        GtoErrMsg("%s: Failed.\n", routine);
    }
    return(status); 
    
    
}




/*----------------------------------------------------------------------------
FUNCTION:    DrKapWrapCheckTCurve

CREATED BY:  Leda Braga -  May 1996

DESCRIPTION: Kapital  wrapper convenience check a TCurve for dates in
             ascending order and for non-negative rates (either zeros
             or vols). 

----------------------------------------------------------------------------*/
GTO_EXPORT(int ) DrKapWrapCheckTCurve(TCurve     *curve,
                                        char       *label)
                 
{
    static char routine[] = "DrKapWrapCheckTCurve";
    int status = FAILURE;
    
    int idx;
    
    /* For each point in the curve, check ascending order */
    /* of dates and that values are non-negative          */
    for (idx=0; idx<curve->fNumItems; idx++)
    { 
        if(curve->fArray[idx].fRate < 0)
        {
            GtoErrMsg("%s: Negative value in %s.\n",
                      routine,
                      label);
            goto done;
        }
        if (idx>0)
        {
            if (curve->fArray[idx].fDate < curve->fArray[idx-1].fDate)
            {
                GtoErrMsg("%s: Dates not in ascending order in %s.\n",
                          routine,
                          label);
                goto done;
            }
        }
    } 
    
     
    status = SUCCESS;
    
  done:
    if (status IS FAILURE)
    {
        GtoErrMsg("%s: Failed.\n", routine);
    }
    return(status); 
    
    
} /* End of DrKapWrapCheckTCurve() */




/*----------------------------------------------------------------------------
FUNCTION:    DrKapWrapDCCtoString

CREATED BY:  David Fung -  September 1997

DESCRIPTION: Kapital wrapper convenience routine to convert a clib day count
             convention to a string

----------------------------------------------------------------------------*/
GTO_EXPORT(int ) DrKapWrapDCCtoString(
    long clibDCC,       /* (I) CLib Daycount conv identifier */
    char *DCCstring)    /* (O) string output mem alloc'd on entry */
{
    int status = FAILURE;
    static char routine[] = "DrKapWrapDCCtoString";

    if (DCCstring IS NULL) goto done;

    switch (clibDCC)
    {
    case GTO_ACT_365F:
        strcpy(DCCstring,"ACT/365F");
        break;
    case GTO_ACT_365:
        strcpy(DCCstring,"ACT/365");
        break;
    case GTO_ACT_360:
        strcpy(DCCstring,"ACT/360");
        break;
    case GTO_B30_360:
        strcpy(DCCstring,"30/360");
        break;
    default:
        GtoErrMsg("%s: Unrecognised clib DCC. (%ld).\n",routine,clibDCC);
        goto done;
    }

    status = SUCCESS;

done:
    if (status IS FAILURE)
    {
        GtoErrMsg("%s: Failed.\n", routine);
    }
    return(status); 
}   

/*----------------------------------------------------------------------------
FUNCTION:    DrKapWrapGetRateConv

CREATED BY:  David Fung -  September 1997

DESCRIPTION: Kapital wrapper convenience routine to get standard rate
             conventions from a zero curve env file

----------------------------------------------------------------------------*/
GTO_EXPORT(int ) DrKapWrapGetRateConv(
    char *fileName,     /* (I)  Kapital zero curve env file */
    long *stdMMDcc,     /* (O)  standard MM Rate Convention */
    long *stdSWDcc,     /* (O)  standard Swap Rate Convention */
    long *stdSWFreq)    /* (O)  standard Swap Rate Freq */    
{
    int         status = FAILURE;
    static char routine[] = "DrKapWrapGetRateConv";

    FILE      *fp;
    char        auxChar;
    long        auxLong;
    char        buf[BUFFER_SIZE];

    /* Open data file   */
    fp = fopen(fileName, "r");     
    if (fp IS NULL)
    {
        GtoErrMsg("%s: File %s cannot be opened.\n", routine, fileName);
        goto done;
    }

    /* skip the value date entry */
    if (DrFindAndSkipComLine(fp, fileName) IS FAILURE) goto done;
    if (fgets(buf,BUFFER_SIZE,fp) == NULL) goto done;

    /* read MM DCC */
    if (DrFindAndSkipComLine(fp, fileName) IS FAILURE) goto done;
    if (fscanf(fp, "%ld \n", &auxLong) != 1) goto done;

    switch (auxLong)
    {
    case 360L:
        *stdMMDcc = GTO_ACT_360;
        break;
    case 365L:
        *stdMMDcc = GTO_ACT_365F;
        break;
    default:
        GtoErrMsg("%s: Unrecognised MM Dcc (%ld) in file (%s).\n",
                  routine,
                  auxLong,
                  fileName);
        goto done;
    }

    /* read Swap Freq */
    if (DrFindAndSkipComLine(fp, fileName) IS FAILURE) goto done;
    if (fscanf(fp, "%c \n", &auxChar) != 1) goto done;

    switch (auxChar)
    {
    case 'S':
        *stdSWFreq = 2L;
        break;
    case 'A':
        *stdSWFreq = 1L;
        break;
    default:
        GtoErrMsg("%s: Unrecognised Swap Freq (%c) in file (%s).\n",
                  routine,
                  auxChar,
                  fileName);
        goto done;
    }

    /* read Swap DCC */
    if (DrFindAndSkipComLine(fp, fileName) IS FAILURE) goto done;
    if (fscanf(fp, "%s \n", buf) != 1) goto done;

    if (strstr(buf,"ACT") ISNT NULL)
    {
        *stdSWDcc = GTO_B30_360;
    }
    else if (strstr(buf,"360") ISNT NULL)
    {
        *stdSWDcc = GTO_ACT_360;
    }
    else if (strstr(buf,"365") ISNT NULL)
    {
        *stdSWDcc = GTO_ACT_365;
    }
    else
    {
        GtoErrMsg("%s: Unrecognised Swap DCC (%s) in file (%s).\n",
                  routine,
                  buf,
                  fileName);
        goto done;
    }

    status = SUCCESS;

done:

    if(fp ISNT NULL) fclose(fp);

    if (status IS FAILURE)
    {
        GtoErrMsg("%s: Failed.\n", routine);
    }
    return(status); 
} /* DrKapWrapGetRateConv */  


/*----------------------------------------------------------------------------
FUNCTION:    DrKapWrapConvDCC

CREATED BY:  David Fung -  September 1997

DESCRIPTION: Kapital wrapper convenience routine to convert DR wrapper
             day count convention to clib convention

----------------------------------------------------------------------------*/
GTO_EXPORT(int ) DrKapWrapConvDCC(
    char  drDCC,                /* (I) Daycount conv in DR wrapper format */
    long *clibDCC)              /* (I) CLib Daycount conv identifier      */
{
    int status = FAILURE;
    static char routine[] = "DrKapWrapConvDCC";

    switch (drDCC)
    {
    case '5':
        *clibDCC = GTO_ACT_365F;
        break;
    case '0':
        *clibDCC = GTO_ACT_360;
        break;
    case '3':
        *clibDCC = GTO_B30_360;
        break;
    default:
        GtoErrMsg("%s: Invalid DR DCC (%c).\n",routine,drDCC);
        goto done;
    }

    status = SUCCESS;

done:
    if (status IS FAILURE)
    {
        GtoErrMsg("%s: Failed.\n", routine);
    }
    return(status); 
}   


/*----------------------------------------------------------------------------
FUNCTION:    DrKapWrapDatesInList

CREATED BY:  David Fung -  September 1997

DESCRIPTION: Kapital wrapper convenience routine to check if a list(2) of 
             dates is a subset of list (1)
             Both lists must be sorted in ascending order on entry 
             with no repeats
             list (1) = fullList
             list (2) = subList
----------------------------------------------------------------------------*/
GTO_EXPORT(int ) DrKapWrapDatesInList(
    long   nbFullList,          /* (I) nb dates in full list */
    TDate *fullList,
    long   nbSubList,           /* (I) nb dates in sub list */
    TDate *subList)
{
    long  fIdx;      /* index for fullList */
    long  sIdx;      /* index for subList */
    TDate prevDate;

    if (nbSubList>nbFullList)
    {
        return(FALSE);
    }

    /* for each subList[sIdx], find fIdx s.t.  */
    /* fullList[fIdx] has the same date        */

    fIdx = 0L;
    prevDate = -999L;

    for (sIdx=0L; sIdx<nbSubList; sIdx++)
    {   
        /* check for repeats */
        if (subList[sIdx] IS prevDate) return (FALSE);

        /* find fIdx as above */
        for (; fullList[fIdx] ISNT subList[sIdx]; fIdx++)
        {
            if (fIdx >= nbFullList-1) return (FALSE);

            if (fullList[fIdx] > subList[sIdx]) return (FALSE);
        }

        prevDate = subList[sIdx];
    }

    return (TRUE);

}   


/*----------------------------------------------------------------------------
FUNCTION:    DrKapWrapCalcAmortPrns

CREATED BY:  David Fung -  September 1997

DESCRIPTION: Kapital wrapper convenience routine to calculate the current
             notional at each date in a payment date list given the info
             from a Dr Wrapper Amortisation schedule

----------------------------------------------------------------------------*/
GTO_EXPORT(int ) DrKapWrapCalcAmortPrns(
    double    origNotional,    /* (I) original notional amount */
    double    origNotionalPct, /* (I) original notional %      */
    long      nbAmortDates,    /* (I) nb amort dates in the schedule */
    TDate    *amortDates,      /* (I) amort dates in the schedule */
    double   *amortPcts,       /* (I) amort % in the schedule */
    long      nbTargetDates,   /* (I) nb of target dates */
    TDate    *targetDates,     /* (I) target list of dates  */
    double  **currNotionals)   /* (O) notional for each date in target list */
{
    int status = FAILURE;
    static char routine[] = "DrKapWrapCalcAmortPrns";

    double *currNotionalsLocal = NULL;

    long amortIdx;
    long targetIdx;
    double sumAmortPct;
    double currAmt;     /* = curr notional for a target date */


    /* check that sum of amrotPcts = origNotionalPct */

    sumAmortPct = 0.0;
    for (amortIdx=0;amortIdx<nbAmortDates;amortIdx++) 
    {
        sumAmortPct += amortPcts[amortIdx];
    }
    if (!IS_ALMOST_ZERO((sumAmortPct-origNotionalPct)/100.0))
    {
        GtoErrMsg("%s: sum of amort%% (%lf)  <> original amort%% (%lf)\n", 
                  routine,
                  sumAmortPct,
                  origNotionalPct);
        goto done;
    }


    /* allocate memory for the current notional list */

    if ((currNotionalsLocal = NEW_ARRAY(double, nbTargetDates)) IS NULL)
    {
        GtoErrMsg("%s: Cannot allocate memory for notional list.\n",routine);
        goto done;
    }


    /* for each targetIdx, amortIdx = min(i) s.t. */
    /* amortDates[i] >= targetDates[targetIdx] */

    currAmt = origNotional * origNotionalPct / 100.0;
    amortIdx = 0L;
    for (targetIdx=0; targetIdx<nbTargetDates; targetIdx++)
    {
        /* update amortIdx = min(i) as above */

        for(;;amortIdx++)
        {
            if (amortIdx>=nbAmortDates)
            {
                GtoErrMsg("%s ...\n",routine);
                goto done;
            }

            if (amortDates[amortIdx]>=targetDates[targetIdx]) break;
            currAmt -= origNotional * amortPcts[amortIdx] / 100.0;
        }

        /* store the notional */

        currNotionalsLocal[targetIdx] = currAmt;
    }

    *currNotionals = currNotionalsLocal;

    status = SUCCESS;
    
  done:

    if (status IS FAILURE)
    {
        FREE_ARRAY(currNotionalsLocal);
        GtoErrMsg("%s: Failed.\n", routine);
    }
    return(status); 
    
}   


/*----------------------------------------------------------------------------
FUNCTION:    DrKapWrapExpandStepUpSchedule

CREATED BY:  David Fung -  September 1997

DESCRIPTION: Kapital wrapper convenience routine to calculate the current 
             rate level at each date in the target date list given the compact
             date/rate step-up schedule in a DR Wrapper format.


----------------------------------------------------------------------------*/
GTO_EXPORT(int ) DrKapWrapExpandStepUpSchedule(
    long     nbSchDates,          /* (I) nb of dates in the schedule */
    TDate   *schDates,            /* (I) schedule dates */
    double  *schRates,            /* (I) schedule rate levels */
    TBoolean isOnOrAfter,         /* (I) T => apply rate on/after sch date 
                                         F => apply only after sch date    */
    long     nbTargetDates,       /* (I) nb of target dates */
    TDate   *targetDates,         /* (I) expanded list of dates */
    double **targetRates)         /* (O) expanded list of rates */
{
    int status = FAILURE;
    static char routine[] = "DrKapWrapExpandStepUpSchedule";

    double  *targetRatesLocal = NULL;

    long schIdx = 0;
    long targetIdx = 0;
    long i;

    /* allocate memory for the targetRates */
    if ((targetRatesLocal = NEW_ARRAY(double, nbTargetDates)) IS NULL)
    {
        GtoErrMsg("%s: Cannot allocate memory for rate list.\n",routine);
        goto done;
    }

    /* for each targetIdx, schIdx = max(i) s.t. */
    /* if isOnOrAfter, then schDates[i] <= targetDates[targetIdx] */
    /*           otherwise, schDates[i] <  targetDates[targetIdx] */

    for (targetIdx=0; targetIdx<nbTargetDates; targetIdx++)
    {
        if (schDates[schIdx] > targetDates[targetIdx])
        {
            GtoErrMsg("%s: Incorrect step-up schedule. \
                       schIdx(%ld) > targetIdx(%ld).\n",
                       schIdx, 
                       targetIdx);
            goto done;
        }

        /* update schIdx = max(i) as above */        
        for (i=schIdx+1; ; i++)
        {
            if (i>=nbSchDates)
            {
                schIdx = nbSchDates-1;
                break;
            }

            if (((isOnOrAfter)  && (schDates[i]>targetDates[targetIdx]))  ||
                ((!isOnOrAfter) && (schDates[i]>=targetDates[targetIdx])))
            {
                schIdx = i-1;
                break;
            }
        }

        targetRatesLocal[targetIdx] = schRates[schIdx];
    }


    *targetRates = targetRatesLocal;

    status = SUCCESS;
    
  done:

    if (status IS FAILURE)
    {
        FREE_ARRAY(targetRatesLocal);
        GtoErrMsg("%s: Failed.\n", routine);
    }
    return(status); 

}    


/*----------------------------------------------------------------------------
FUNCTION:    DrKapWrapParamInput

CREATED BY:  David Fung -  September 1997

DESCRIPTION: Kapital wrapper convenience routine to read in model parameters
             from either the 4 string lines in the DR wrapper text, or from
             the environment model parameters file.

----------------------------------------------------------------------------*/
GTO_EXPORT(int ) DrKapWrapParamInput(                                                   
    long    nbFactor,                        /* (I) Number of factors       */
    char    overWriteString[3][BUFFER_SIZE], /* (I) string (alpha,beta,rho) */
    char   *fileName,                        /* (I) File name inc extension */
    double *alphas,                          /* (O) size 3 alloc'd on entry */
    double *betas,                           /* (O) size 3 alloc'd on entry */
    double *corrs)                           /* (O) size 3 alloc'd on entry */
{
    int         status = FAILURE;
    static char routine[] = "DrKapWrapParamInput";
    int         i; 
    FILE      *stream = NULL;
    char        buf[BUFFER_SIZE];
    TBoolean    needDefaultValues = FALSE;

    /* array storing entries from the model parameters file */
    /* element i = ith entry in the file */

#undef nbModParamEntries
#define nbModParamEntries 19
    
    double *modParamEntries[nbModParamEntries]; 

    /* initialise modParamEntries */
    for (i=0;i<nbModParamEntries;i++)
    {
        modParamEntries[i] = NULL;
    }


    switch (nbFactor)
    {
    case 1L:

        /***************
         * read alphas *
         ***************/
        if (strstr (overWriteString[0],"nil") IS NULL)
        {
            /*** read overrides ***/
            if (GtoSscanf(overWriteString[0], "%lf \n", 
                          &(alphas[0])) IS FAILURE)
            {   
                GtoErrMsg("%s: Could not find file %s: \
                                       factor weight overwrite required!\n",
                            routine,
                            fileName);
                goto done;
            }
        }
        else
        {
            /*** read from env file ***/
            modParamEntries[1] = &(alphas[0]);
        }

        /***************
         * read betas *
         ***************/
        if (strstr (overWriteString[1],"nil") IS NULL)
        {
            /*** read overrides ***/
            if (GtoSscanf(overWriteString[1], "%lf \n", 
                          &(betas[0])) IS FAILURE)
            {   
                GtoErrMsg("%s: Could not find file %s: \
                                       mean reversion overwrite required!\n",
                            routine,
                            fileName);
                goto done;
            }
        }
        else
        {
            /*** read from env file ***/
            modParamEntries[0] = &(betas[0]);
        }

        
        /* Fill in the unused parameters with N/A values */
        alphas[1] = alphas[2] = -999.;
        betas[1] = betas[2] = -999.;
        corrs[0] = corrs[1] = corrs[2] = -999.;

        break;

    case 2L:
        /***************
         * read alphas *
         ***************/
        if (strstr (overWriteString[0],"nil") IS NULL)
        {
            /*** read overrides ***/
            if (GtoSscanf(overWriteString[0], "%lf %lf \n", 
                          &(alphas[0]), 
                          &(alphas[1])) IS FAILURE)
            {                   
                GtoErrMsg("%s: Could not find file %s: \
                                       factor weight overwrite required!\n",
                            routine,
                            fileName);
                goto done;
            }
        }
        else
        {
            /*** read from env file ***/
            modParamEntries[5] = &(alphas[0]);
            modParamEntries[6] = &(alphas[1]);
        }

        /***************
         * read betas *
         ***************/
        if (strstr (overWriteString[1],"nil") IS NULL)
        {
            /*** read overrides ***/
            if (GtoSscanf(overWriteString[1], "%lf %lf \n", 
                          &(betas[0]), 
                          &(betas[1])) IS FAILURE)
            {   
                GtoErrMsg("%s: Could not find file %s: \
                                       mean reversion overwrite required!\n",
                            routine,
                            fileName);
                goto done;
            }
        }
        else
        {
            /*** read from env file ***/
            modParamEntries[3] = &(betas[0]);
            modParamEntries[4] = &(betas[1]);
        }


        /**************
         * read corrs *
         **************/
        if (strstr (overWriteString[2],"nil") IS NULL)
        {
            /*** read overrides ***/
            if (GtoSscanf(overWriteString[2], "%lf \n", 
                          &(corrs[0])) IS FAILURE)
            {   
                GtoErrMsg("%s: Could not find file %s: \
                                       correlation overwrite required!\n",
                            routine,
                            fileName);
                goto done;
            }
        }
        else
        {
            /*** read from env file ***/
            modParamEntries[7] = &(corrs[0]);
        }


        /* Fill in the unused parameters with N/A values */
        alphas[2] = -999.;
        betas[2] = -999.;
        corrs[1] = corrs[2] = -999.;

        break;

    case 3L:
        /***************
         * read alphas *
         ***************/
        if (strstr (overWriteString[0],"nil") IS NULL)
        {
            /*** read overrides ***/
            if (GtoSscanf(overWriteString[0], "%lf %lf %lf \n", 
                          &(alphas[0]), 
                          &(alphas[1]), 
                          &(alphas[2])) IS FAILURE)
            {                   
                GtoErrMsg("%s: Could not find file %s: \
                                       factor weight overwrite required!\n",
                            routine,
                            fileName);
                goto done;
            }
        }
        else
        {
            /*** read from env file ***/
            modParamEntries[12] = &(alphas[0]);
            modParamEntries[13] = &(alphas[1]);
            modParamEntries[14] = &(alphas[2]);
        }

        /***************
         * read betas *
         ***************/
        if (strstr (overWriteString[1],"nil") IS NULL)
        {
            /*** read overrides ***/
            if (GtoSscanf(overWriteString[1], "%lf %lf %lf \n", 
                          &(betas[0]), 
                          &(betas[1]), 
                          &(betas[2])) IS FAILURE)
            {   
                GtoErrMsg("%s: Could not find file %s: \
                                       mean reversion overwrite required!\n",
                            routine,
                            fileName);
                goto done;
            }
        }
        else
        {
            /*** read from env file ***/
            modParamEntries[9] = &(betas[0]);
            modParamEntries[10] = &(betas[1]);
            modParamEntries[11] = &(betas[2]);
        }


        /**************
         * read corrs *
         **************/
        if (strstr (overWriteString[2],"nil") IS NULL)
        {
            /*** read overrides ***/
            if (GtoSscanf(overWriteString[2], "%lf %lf %lf \n", 
                          &(corrs[0]), 
                          &(corrs[1]), 
                          &(corrs[2])) IS FAILURE)
            {   
                GtoErrMsg("%s: Could not find file %s: \
                                       correlation overwrite required!\n",
                            routine,
                            fileName);
                goto done;
            }
        }
        else
        {
            /*** read from env file ***/
            modParamEntries[15] = &(corrs[0]);
            modParamEntries[16] = &(corrs[1]);
            modParamEntries[17] = &(corrs[2]);
        }

        break;

    default:
        GtoErrMsg("%s: nb of factors (%ld) is not 1,2 or 3.\n", 
                  routine, 
                  nbFactor);
        goto done;
    }

    /* read model parameter entries from the env file */

    /* check of default values are needed */
    needDefaultValues = FALSE;
    for (i=0; ((i<nbModParamEntries) && (!needDefaultValues)) ; i++)
    {
        needDefaultValues = needDefaultValues || (modParamEntries[i] != NULL);
    }

    if (needDefaultValues)
    {

        stream = fopen(fileName, "r");

        if (stream == NULL)
        {
            GtoErrMsg("%s: file %s required but cannot be opened \n",
                      routine,
                      fileName);
            goto done;
        }


        for (i=0; i<nbModParamEntries; i++)
        {
            if (DrFindAndSkipComLine(stream, fileName) IS FAILURE) goto done;
    
            if (modParamEntries[i] IS NULL)
            {
                /* skip this line */            
                if (fgets(buf,BUFFER_SIZE,stream) == NULL) goto done;
            }
            else
            {
                /* read the parameter */
                if (fscanf(stream, "%lf \n", modParamEntries[i]) != 1) 
                {
                    GtoErrMsg("%s: Cannot read model parameter entry (%d).\n", 
                              routine,
                          i);
                goto done;
                }
            }
        }
    }

    /* Check validity of input */
    if (DrKapWrapParamCheck(nbFactor,
                            alphas,
                            betas,
                            corrs) IS FAILURE) goto done;
        
    status = SUCCESS;

done:        
    if (stream ISNT NULL) fclose(stream);

#undef nbModParamEntries

    if (status IS FAILURE)
    {
        GtoErrMsg("%s: Failed.\n", routine);
    }
    return(status); 

}



/*----------------------------------------------------------------------------
FUNCTION:    DrKapWrapParamCheck

CREATED BY:  David Fung -  September 1997

DESCRIPTION: Kapital wrapper convenience routine to verify the model 
             parameters             

----------------------------------------------------------------------------*/
GTO_EXPORT(int ) DrKapWrapParamCheck(
    long      nbFactor,     /* (I) Number of factors     */
    double   *alphas,       /* (I) alphas                */
    double   *betas,        /* (I) mean reversions       */
    double   *corrs)        /* (I) corrleations          */
{
    int status = FAILURE;   /* Error status = FAILURE initially */
    static char routine[] = "DrKapWrapParamCheck";


    #undef  CHECK_IN_RANGE
    #define CHECK_IN_RANGE(num,lo,hi,msg)  \
                if ( ((num)<(lo)) || ((num)>(hi)) ) \
                {                                   \
                    GtoErrMsg(msg);                 \
                    goto done;                      \
                }                                   
            
    switch (nbFactor)
    {
    case 1L:
        CHECK_IN_RANGE(alphas[0],0.0001,10000.0,"Weight #1 out of range!\n");

        CHECK_IN_RANGE(betas[0], 0.0,   10.0,   "Beta #1 out of range!\n");
        break;

    case 2L:
        CHECK_IN_RANGE(alphas[0],0.0001,1000.0,"Weight #1 out of range!\n");
        CHECK_IN_RANGE(alphas[1],0.0001,1000.0,"Weight #2 out of range!\n");

        CHECK_IN_RANGE(betas[0], 0.0,   10.0,  "Beta #1 out of range!\n");
        CHECK_IN_RANGE(betas[1], 0.0,   10.0,  "Beta #2 out of range!\n");

        CHECK_IN_RANGE(corrs[0], -.951, .951,  "Corr out of range!\n");
        break;
    case 3L:
        CHECK_IN_RANGE(alphas[0],0.0001,1000.0,"Weight #1 out of range!\n");
        CHECK_IN_RANGE(alphas[1],0.0001,1000.0,"Weight #2 out of range!\n");
        CHECK_IN_RANGE(alphas[2],0.0001,1000.0,"Weight #3 out of range!\n");

        CHECK_IN_RANGE(betas[0], 0.0,   10.0,  "Beta #1 out of range!\n");
        CHECK_IN_RANGE(betas[1], 0.0,   10.0,  "Beta #2 out of range!\n");
        CHECK_IN_RANGE(betas[2], 0.0,   10.0,  "Beta #3 out of range!\n");

        CHECK_IN_RANGE(corrs[0], -.951, .951, "Corr #1 out of range!\n");
        CHECK_IN_RANGE(corrs[1], -.951, .951, "Corr #2 out of range!\n");
        CHECK_IN_RANGE(corrs[2], -.951, .951, "Corr #3 out of range!\n");
        break;
    }


    status = SUCCESS;
        
done:

    if (status IS FAILURE)
    {
        GtoErrMsg("%s: Failed.\n", routine);
    }

    #undef  CHECK_IN_RANGE    
    return (status);

}


/*----------------------------------------------------------------------------
FUNCTION:    DrExtendTCurve

CREATED BY:  David Fung -  September 1997

DESCRIPTION: Kapital wrapper convenience routine to extend a curve (e.g. zero
             curve) to include dates implied by BaseDate and intervals
             NOTE: baseDate must be >= oriCurve->fBaseDate
----------------------------------------------------------------------------*/

GTO_EXPORT(int ) DrExtendTCurve(
    TCurve         *oriCurve,        /* (I) original curve        */
    TDate           baseDate,        /* (I) start date for intvls */
    long            numIntervals,    /* (I) nb intervals in list  */
    TDateInterval  *intervals,       /* (I) interval list         */
    TCurve        **newCurve)        /* (O) result curve          */
{
    int         status = FAILURE;
    static char routine[] = "DrExtendTCurve";

    TCurve     *newCurveLocal = NULL;

    TDateList  *oriDL = NULL;
    TDateList  *additionalDL = NULL;
    TDateList  *mergedDL = NULL;
    int         i;

    /* make sure base date >= oriCurve->fBaseDate */
    /* to ensure valid interpolation              */
    
    if (baseDate < oriCurve->fBaseDate)
    {
        GtoErrMsg("%s: base Date(%s) < original curve base date (%s) \n",
                  routine,
                  GtoFormatEuroDate(baseDate),
                  GtoFormatEuroDate(oriCurve->fBaseDate));
        goto done;
    }

    /* create combined datelist */
    /* ------------------------ */

    if ((oriDL = GtoNewDateListFromTCurve(oriCurve)) IS NULL) goto done;

    if ((additionalDL = GtoNewEmptyDateList(numIntervals)) IS NULL) goto done;

    for (i=0; i<additionalDL->fNumItems; i++)
    {
        if (GtoDtFwdAny(baseDate,
		                &(intervals[i]),
		                &(additionalDL->fArray[i])) IS FAILURE) goto done;
    }

    mergedDL = GtoMergeDateLists(oriDL,additionalDL);
    if (mergedDL IS NULL) goto done;

    /* create new curve */
    /* ---------------- */

    newCurveLocal = GtoNewTCurveFromDL(mergedDL,
                                       oriCurve->fBaseDate,
                                       oriCurve->fBasis,
                                       oriCurve->fDayCountConv);

    if (newCurveLocal IS NULL) goto done;

    /* fill in the rates */

    for (i=0; i<newCurveLocal->fNumItems; i++)
    {
        if (GtoInterpRate(newCurveLocal->fArray[i].fDate,
                          oriCurve,
                          GTO_LINEAR_INTERP,
                          &(newCurveLocal->fArray[i].fRate)) IS FAILURE) 
            goto done;
    }

    /* pass pointer to output */

    *newCurve = newCurveLocal;

    status = SUCCESS;
        
done:

    GtoFreeDateList(oriDL);
    GtoFreeDateList(additionalDL);
    GtoFreeDateList(mergedDL);

    if (status IS FAILURE)
    {
        GtoFreeTCurve(newCurveLocal);
    }

    return (status);

} /* DrExtendTCurve */


/*f----------------------------------------------------------------------
  FUNCTION:       DRLWrapTDateIntvlToNum
  
  CREATED BY:     Julia Chislenko, October 1997
  
  DESCRIPTION:    Read from a deal file maturity as a TDateInterval
                  and convert it into number of months.

*/   
GTO_EXPORT(int ) DRLWrapTDateIntvlToNum(FILE *fp,
				  char *fileName,
				  char *paramName,  /* for error message */
				  char intvlType,   /* 'M' of 'D' */
				  long *numMonths)  /* (O) */
{
    char   string[80];
    TDateInterval dateIntvl;

    if (DrFindAndSkipComLine(fp, fileName) IS FAILURE) goto done;

    if (fscanf(fp, "%s \n", 
                  string) != 1)                 
    {  
        GtoErrMsg("Cannot read %s.\n", paramName);
        goto done;
    }
    if(GtoStringToDateInterval(string, paramName, &dateIntvl)
       == FAILURE)
	goto done;

    if(dateIntvl.prd_typ != toupper(intvlType))
    {  
        GtoErrMsg("Expected %c interval type for %s.\n", 
		  toupper(intvlType), paramName);
        goto done;
    }

    *numMonths = (long)dateIntvl.prd;

    return SUCCESS;
     
  done:
    return FAILURE;
}

/*f----------------------------------------------------------------------
  FUNCTION:       DRLWrapDayCountToString
  
  CREATED BY:     Julia Chislenko, October 1997
  
  DESCRIPTION:    Read from a deal file London style day count conv
                  and convert it into a string.

*/   
GTO_EXPORT(int ) DRLWrapDayCountToString(FILE *fp,
				  char  *fileName,
				  char *paramName,  /* for error message */
				  char *string)     /* (O) */
{
    char   auxChar;
    long   auxLong;

    if (DrFindAndSkipComLine(fp, fileName) IS FAILURE) goto done;

    if (fscanf(fp, "%c \n", 
                  &auxChar) != 1)                 
    {  
        GtoErrMsg("Cannot read %s.\n", paramName);
        goto done;
    }

    if (DrKapWrapConvDCC(auxChar,&auxLong) IS FAILURE)
        goto done;

    if(DrKapWrapDCCtoString(auxLong,string) IS FAILURE)
        goto done;

    return SUCCESS;
 
  done:
    return FAILURE;
}

/*f----------------------------------------------------------------------
  FUNCTION:       DRLWrapDayCountToString
  
  CREATED BY:     Julia Chislenko, October 1997
  
  DESCRIPTION:    Read from a deal file London style frequency
                  and convert it into a num per year.

*/   
GTO_EXPORT(int ) DRLWrapFrequencyToLong(FILE *fp,
				  char *fileName,
				  char *paramName,     /* for error message */
				  long *frequency)     /* (O) */
{
    char   auxChar;
    long   auxLong;

    if (DrFindAndSkipComLine(fp, fileName) IS FAILURE) goto done;

    if (fscanf(fp, "%c \n", 
                  &auxChar) != 1)                 
    {  
        GtoErrMsg("Cannot read %s.\n", paramName);
        goto done;
    }

    auxChar = (char)toupper(auxChar);
    
    switch (auxChar) 
    {
    case 'P':
	auxLong = 0;
	break;
     case 'A':
	auxLong = 1;
	break;
    case 'S':
	auxLong = 2;
	break;
    case 'Q':
	auxLong = 4;
	break;
    case 'M':
	auxLong = 12;
	break;
    default:
        GtoErrMsg("Invalid input (%c) for %s\n", 
                  auxChar, paramName);
        goto done;
    }

    *frequency = auxLong;

    return SUCCESS;
 
  done:
    return FAILURE;
}
