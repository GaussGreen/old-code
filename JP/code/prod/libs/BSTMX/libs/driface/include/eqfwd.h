/****************************************************************************/
/*      Declaration file for structure of data.                             */
/****************************************************************************/
/*      eqfwd   .h                                                          */
/****************************************************************************/
#ifndef EQFWD_H
#define EQFWD_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cgeneral.h"
#include "bastypes.h"
#include "macros.h"
#include "tcurve.h"
#include "ldate.h"
#include "cerror.h"
#include "drieq.h"

GTO_EXPORT(int)
DrLondonForwardPrice(
        double   spotPrice,       /* (I) spotPrice stock price               */
        int      NbDiv,           /* (I) Number of dividend / forward prices */
        double   *divident,       /* (I) Dividends / Forward prices          */
        long     *divDate,        /* (I) ex-Dividends / forwards dates       */
        long     *payDivDate,     /* (I) Dividend pay dates                  */
        long     *divType,        /* (I) Dividend / forward type             */
        int      NbSettle,        /* (I) Nb of settlement dates(<0 if rolling*/
        TDate    *LastTrading,    /* (I) Last trading date of settl period   */
        TDate    *SettleDate,     /* (I) Corresponding settlement date       */
        char     settleType,      /* (I) Settlement type:'F'ixed or 'R'olling*/
        long     busDayConv,      /* (I) Business day conv                   */
        char     *holidayFile,    /* (I) "NONE" for weekends only
                                   *     "No_Weekends" for no adjustments    */
        TCurve   *zcCurve,        /* (I) Funding Zero Curve                  */
        TDate    *fwdDate,        /* (I) Forward dates to compute fwd price  */
        int      NbFwd,           /* (I) Total number of time points         */
        double   *fwdPrice);      /* (O) Fwd price at eachtime forward date  */

GTO_EXPORT(int)
DrForwardPriceGen(
        double          spotPrice,      /* (I) spot asset price   */
        TEqStatData     *eqStatData,    /* (I) Equity static data */
        TCurve          *zcCoF,         /* (I) Funding zero curve */
        long		busDayConv,	/* (I) Business day conv  */
        char		*holidayFile,	/* (I) Holiday file       */
        long            numFwd,         /* (I) Number of fwd points */
        TDate           *fwdDate,       /* (I) Forward dates to compute fwdPrice
                                         *     fwdDate[0] must be the spotDate
                                         *     where spotPrice is known  */
        double          *fwdPrice);     /* (O) Forward price      */

GTO_EXPORT(int)
DriForwardDiscZero(
        TEqStatData     *eqStatData,    /* (I) Equity static data */
        TCurve          *discZcCurve,   /* (I) Discount zero curve */
        char            *holidayFile,   /* (I) Holiday File       */
        long            numFwd,         /* (I) Number of fwd points */
        TDate           *fwdDate,       /* (I) Forward dates         */
        double          *fwdDiscZero);  /* (O) disc zero at each fwd date  */

#endif
