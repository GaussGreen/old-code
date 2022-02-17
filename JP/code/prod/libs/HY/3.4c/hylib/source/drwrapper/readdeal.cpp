/*------------------------------------------------------------------
C FILE:         readdeal.cpp

CREATED BY:     Neil Yang - Feb 2000

PURPOSE:        DR Wrapper to read deal informations.
---------------------------------------------------------------------- */

#include <ctype.h>
#include <string.h>
#include <math.h>

#include "bastypes.h"
#include "ldate.h"      /* GTO_ACT_... */
#include "strutil.h"    /* GtoStringToUpper */
#include "cerror.h"  
#include "cmemory.h"
#include "macros.h"
#include "tcurve.h"
#include "cfileio.h"
#include "cashflow.h"   /* GtoMakeCFL */
#include "convert.h"    /* GtoMakeDateInterval */
#include "zr2fwd.h"     /* GtoForwardFromZC */
#include "busday.h"     /* GtoDateFwdThenAdjust */

#include "kexception.h"
#include "genfunc.h"
#include "instruments.h"
#if ! defined(_WIN32)
#include "Legacy.h"
#endif
#include "readdeal.h"


/*----------------------------------------------------------------------------
FUNCTION:    DrKapWrapGetCreditDefaultSwap

CREATED BY:  Neil Yang -- Feb 2000

DESCRIPTION: Kapital wrapper function to read the deal data for a credit default swap
			 hycds_t.dat file (in DR format) into the necessary market
             parameters.
             
             Restrictions are applied concerning the format of the data
             file.
----------------------------------------------------------------------------*/
GTO_EXPORT(int) DrKapWrapGetCreditDefaultSwap(
    char            *filename,              /* (I) Usually "resetccs_t.dat"     */
    double          *ppy,                   /* (O) period per year */
	double          *beta,                  /* (O) */
	double          *dps,                   /* (O) */
	double          *x,                     /* (O) */
	double          *lim,                   /* (O) */
	double          *lim1,                  /* (O) */
	double          *lim2,                  /* (O) */
	double          *recovery,              /* (O) */
	double          *notional,              /* (O) */
	TDate           *maturityDate,          /* (O) */
	char            fee[])           /* (O) */
{
    static char routine[]="DrKapWrapGetCreditDefaultSwap";
    int        status = FAILURE;    
        
    FILE        *fp   = NULL;   /* File pointer */
    long        dateDR;
	
    
    /* Open data file   */
    fp = fopen(filename, "r");     
    if (fp IS NULL)
    {
        GtoErrMsg("%s: File %s cannot be opened.\n",  
                  routine,
                  filename);
        goto done;
    }
    
    /* 1 - Skip 2 comment lines and get Maturity Date*/
    if (DrFindAndSkipComLine(fp,filename) IS FAILURE)
        goto done;

    if (DrFindAndSkipComLine(fp,filename) IS FAILURE)
        goto done;

	/* read in debt per share */
    if (fscanf(fp,"%lf\n",dps) != 1)
    {   
        KException("Error reading debt per share.\n");
    }
    
	/* read exponent x */
	if (DrFindAndSkipComLine(fp,filename) IS FAILURE)
        goto done;
	if (fscanf(fp,"%lf\n",x) != 1)
    {   
        KException("Error reading exponent x.\n");
    }

	/* read lim */
	if (DrFindAndSkipComLine(fp,filename) IS FAILURE)
        goto done;
	if (fscanf(fp,"%lf\n",lim) != 1)
    {   
        KException("Error reading lim.\n");
    }

	/* read lim1 */
	if (DrFindAndSkipComLine(fp,filename) IS FAILURE)
        goto done;
	if (fscanf(fp,"%lf\n",lim1) != 1)
    {   
        KException("Error reading lim1.\n");
    }

	/* read lim2 */
	if (DrFindAndSkipComLine(fp,filename) IS FAILURE)
        goto done;
	if (fscanf(fp,"%lf\n",lim2) != 1)
    {   
        KException("Error reading lim2.\n");
    }

	/* read recovery */
	if (DrFindAndSkipComLine(fp,filename) IS FAILURE)
        goto done;
	if (fscanf(fp,"%lf\n",recovery) != 1)
    {   
        KException("Error reading recovery.\n");
    }

	/* read ppy */
	if (DrFindAndSkipComLine(fp,filename) IS FAILURE)
        goto done;
	if (fscanf(fp,"%lf\n",ppy) != 1)
    {   
        KException("Error reading ppy.\n");
    }

	/* read beta */
	if (DrFindAndSkipComLine(fp,filename) IS FAILURE)
        goto done;
	if (fscanf(fp,"%lf\n",beta) != 1)
    {   
        KException("Error reading beta.\n");
    }

	/* read Notional */
	if (DrFindAndSkipComLine(fp,filename) IS FAILURE)
        goto done;
	if (fscanf(fp,"%lf\n",notional) != 1)
    {   
        KException("Error reading notional.\n");
    }

	/* read maturity date */
	if (DrFindAndSkipComLine(fp,filename) IS FAILURE)
        goto done;
	if (fscanf(fp,"%ld\n",&dateDR) != 1)
    {   
        KException("Error reading Mat Date.\n");
    }

    if (DrDateToTDate(dateDR,
                      maturityDate) IS FAILURE)
        goto done;

	/* read fee */
	if (DrFindAndSkipComLine(fp,filename) IS FAILURE)
        goto done;
	if (fscanf(fp,"%s\n",fee) != 1)
    {   
        KException("Error reading fee.\n");
    }

    status = SUCCESS;  
    
done:
  
    if(fp ISNT NULL) fclose(fp);
    
    if (status IS FAILURE)
    {   
        GtoErrMsg("%s: Failed.\n",routine);
    }
        
    return (status);
}
