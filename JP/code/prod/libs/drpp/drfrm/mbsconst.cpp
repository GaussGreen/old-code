/***************************************************************************
 *	SCCS Keyword Information
 *	------------------------
 *	Module name	:  mbsconst.c
 *	Company Name	:  JP Morgan Securities Inc.
 *	Author		:  Davis (Shuenn-Tyan) Lee
 *			   Derivatives Research
 *	Code version    :  1.30
 *	Extracted	:  2/11/97 at 15:13:02
 *	Last Updated	:  2/11/97 at 15:12:59
 ***************************************************************************
 *      Common constants/structures/methods for prepay/pricing code
 *
 *      Copyright 1996 J.P. Morgan & Co. Incorporated. All rights reserved.
 ***************************************************************************/
#include <stdio.h>
#include <stdlib.h>
extern "C" {
#include "cerror.h"
#include "ldate.h"
#include "tcurve.h"
#include "convert.h"
#include "datelist.h"
#include "busday.h"
}
#include "mbsbase.h"
#include "mbsdate.h"
#include "ppbase.h"
#include "mbsconst.h"


/***************************************************************************
 *  PUBLIC METHODS
 ***************************************************************************/

/***************************************************************************
 * MakeMbsDeal();
 * Constructor for TMbsDeal
 * Sets internals of structure to input values
 * If input ptr is NULL, also allocates structure
 * Note: calling code must explicitly construct internal
 * structures (e.g., prepay, mbsio, mbspo) after calling this function
 * Returns SUCCESS/FAILURE
 ***************************************************************************/
int EXPORT
MakeMbsDeal
   (/* UNDERLYING MBS */
    long          mbsType,             /* (I) basic type of MBS (fixed/ARM) */
    char         *mbsSubtype,          /* (I) supplementary info */
    long          mbsAgency,           /* (I) issuing agency 
                                        * (e.g., MBS_AGENCY_FNMA) */
    long          mbsTerm,             /* (I) Either MBS_PP_MBSTERM_30 or
                                        * MBS_PP_MBSTERM_15 */

    TDate         mbsWarmWalaAsOfDate, /* (I) date(month) for which WARM and
                                        * WALA (below) are valid */
    long          mbsWarm,             /* (I) WARM of MBS in as-of month */
    long          mbsWala,             /* (I) WALA of MBS in as-of month */
    /* UNDERLYING MBS CPN BEHAVIOR(S) */
    TMbsIoRule   *mbsNetIoRule,	       /* (I) net cpn behavior of underlying */
    TMbsIoRule   *mbsGrossIoRule,      /* (I) gross cpn behavior of underlying */
    TMbsPoRule   *mbsPoRule,	       /* (I) PO behavior of underlying */
    /* PREPAY INFO */
    TMbsIoRule   *mbsAmortRule,        /* amortization/prepay reset behavior
                                        * of underlying */
    TMbsPrepayAssump *mbsPrepayAssump, /* main block of prepay inputs/assump */
    double       *mbsGrCcSprd,         /* array of MBS_PP_NUM_REFI_TYPES
                                        * spreads (for this mbs type) between
					* gross current coupon and each of all
					* valid refi index rate types */
    /* BASIC DEAL TIMING INFO  */
    long          dealPayIntervalMons, /* length (in months) of regular deal
                                        * accrual period or, equivalently, 
                                        * months between regular cashflow dates;
                                        * typically 1 */
    long          dealMbsAccDelayDays, /* approx. # days (act) between start of
                                        * MBS accrual period and of (later) deal
                                        * regular accrual period (e.g., 19 for
                                        * ARM classic swaps) */
    /* DEAL ACCRUAL/CASHFLOW DATES:
     * We allow two methods of defining these dates:
     *   1) "By rule": User can supply these "general" reset 
     *   parameters--constructor will then build internal lists).... */
    TDate         dealFirstAccStDate,  /* start date of 1st accrual
                                        * period of deal (may be stub);
                                        * 0 if using method 2 */
    TDate         dealRegularAccStDate, /* start date of 1st regular
                                        * accrual period of deal;
                                        * 0 if using method 2 */
    long          dealCashflowDelayDays, /* deal cashflows occur this # days
                                        * after start of NEXT accrual pd;
                                        * typ. 0 for swaps;
                                        * 0 if using method 2 */
    TDate         dealMatDate,	       /* date of last cashflow of deal;
                                        * 0 if using method 2 */
    /*   2) "By list": user can supply lists of cf/accrual dates */
    long          numCashflows,        /* # of dates in next two lists;
                                        * 0 if using method 1 */
    TDate        *dealAccStDates,      /* list of deal accrual start dates;
                                        * NULL if using method 1 */
    TDate        *dealCfDates,         /* list of deal cashflow dates;
                                        * NULL if using method 1 */
    /* LEG A & B OF SWAP */
    TMbsIoRule   *legAIoRule,          /* IO rule for leg A of swap */
    TMbsPoRule   *legAPoRule,          /* PO rule for leg A of swap */
    TMbsIoRule   *legBIoRule,          /* IO rule for leg B of swap */
    TMbsPoRule   *legBPoRule,          /* PO rule for leg B of swap */
    /* PTR TO NEW/EXISTING STRUCTURE */
    TMbsDeal    **mbsDeal)          /* (I/O) structure to alloc/modify */
{
    F_INIT("MakeMbsDeal");
    long     i;
    long    iDate;
    long    iCf;
    TDate   firstMbsCfDateEver;
    TDate  *dealAccEndDates = NULL;

    /* Check inputs
     */
    if (!is_valid_mbs_type(mbsType))
    {
        GtoErrMsg("%s : Invalid MBS type [%d]\n", routine, mbsType);
        goto done;
    }
    if (!is_valid_agency_type(mbsAgency))
    {
        GtoErrMsg("%s : Please check agency type [%d]\n",routine, mbsAgency);
        goto done;
    }
    if (!is_valid_mbs_term(mbsTerm))
    {
        GtoErrMsg("%s : Please check MBS term [%d]\n",routine, mbsTerm);
        goto done;
    }

    /* Prepare an empty structure (may have to de-alloc existing one) */
    if (*mbsDeal ISNT NULL)
    {
        FreeMbsDeal(mbsDeal);
    }
    if ((*mbsDeal = NEW(TMbsDeal)) IS NULL)
    {
        GtoErrMsg("%s : Failed to allocate mbsDeal structure\n", routine);
        goto done;
    }

    /* Just to be safe, null out ptrs for strucs not yet allocated */
    (*mbsDeal)->mbsNetIoRule = NULL;
    (*mbsDeal)->mbsGrossIoRule = NULL;
    (*mbsDeal)->mbsPoRule = NULL;
    (*mbsDeal)->mbsAmortRule = NULL;
    (*mbsDeal)->mbsPrepayAssump = NULL;
    (*mbsDeal)->dealAccStDL = NULL;
    (*mbsDeal)->dealAccEndDL = NULL;
    (*mbsDeal)->dealCfDL = NULL;
    (*mbsDeal)->dealEstRegAccStDL = NULL;
    (*mbsDeal)->dealEstRegAccEndDL = NULL;
    (*mbsDeal)->dealEstRegCfDL = NULL;
    (*mbsDeal)->legAIoRule = NULL;
    (*mbsDeal)->legAPoRule = NULL;
    (*mbsDeal)->legBIoRule = NULL;
    (*mbsDeal)->legBPoRule = NULL;

    /* Set member data 
     */
    (*mbsDeal)->mbsType = mbsType;
    strncpy((*mbsDeal)->mbsSubtype,mbsSubtype,LINESTRLEN-1);
    (*mbsDeal)->mbsAgency = mbsAgency;
    (*mbsDeal)->mbsTerm = mbsTerm;

    /* compute (est.) origination & maturity dates using
     * WARM/WALA supplied */
    (*mbsDeal)->mbsOrigDate = mbsWarmWalaAsOfDate;
    XONFAIL( SetDOM(1,TRUE,&((*mbsDeal)->mbsOrigDate)) );
    XONFAIL( NxtMth((*mbsDeal)->mbsOrigDate,
                    -mbsWala,
                    &((*mbsDeal)->mbsOrigDate)) );
    (*mbsDeal)->mbsMatDate = mbsWarmWalaAsOfDate;
    XONFAIL( SetDOM(1,TRUE,&((*mbsDeal)->mbsMatDate)) );
    XONFAIL( NxtMth((*mbsDeal)->mbsMatDate,
                    mbsWarm,
                    &((*mbsDeal)->mbsMatDate)) );
    /* NB: we alter day-of-month of this mat date to appropriate
     * cashflow day-of-month, below */

    /* Also compute warm+wala (time-invariant) */
    (*mbsDeal)->mbsEffOrigTerm = mbsWarm+mbsWala;

    /* For all MBS, we can assume calendar month accrual */
    (*mbsDeal)->mbsAccrualStDom = 1;

    /* Based on agency/term, set MBS accrual/cashflow day-of-month */
    switch (mbsAgency)
    {
      case MBS_AGENCY_FNMA:
        (*mbsDeal)->mbsCashflowDom = 25;
        break;
      case MBS_AGENCY_FHLMC:
      case MBS_AGENCY_GOLD:
      case MBS_AGENCY_GNMAI:
        (*mbsDeal)->mbsCashflowDom = 15;
        break;
      case MBS_AGENCY_GNMAII:
        (*mbsDeal)->mbsCashflowDom = 20;
        break;
        /* For whole loans, could be anything--user will
         * have to (re-)set this by hand after calling constructor */
      case MBS_AGENCY_WHOLE:
        (*mbsDeal)->mbsCashflowDom = 15; /* dummy value */
        break;
      default:
        GtoErrMsg("%s : Failed to recognize MBS agency %d\n",routine,mbsAgency);
        goto done;
    }
    (*mbsDeal)->mbsPayDelayDays = (*mbsDeal)->mbsCashflowDom
        - (*mbsDeal)->mbsAccrualStDom;
    XONFAIL( SetDOM((*mbsDeal)->mbsCashflowDom,
                    TRUE,
                    &((*mbsDeal)->mbsMatDate)) );


    /* Copy (already constructed!) IO/PO rules */
    if (CopyMbsIoRule(mbsNetIoRule,
                     &((*mbsDeal)->mbsNetIoRule)) ISNT SUCCESS)
    {
        GtoErrMsg("%s : CopyMbsIoRule failed\n",routine);
        goto done;
    }
    if (CopyMbsIoRule(mbsGrossIoRule,
                     &((*mbsDeal)->mbsGrossIoRule)) ISNT SUCCESS)
    {
        GtoErrMsg("%s : CopyMbsIoRule failed\n",routine);
        goto done;
    }
    if (CopyMbsPoRule(mbsPoRule,
                     &((*mbsDeal)->mbsPoRule)) ISNT SUCCESS)
    {
        GtoErrMsg("%s : CopyMbsPoRule failed\n",routine);
        goto done;
    }
    /* Copy (already constructed!) prepay reset rule */
    if (CopyMbsIoRule(mbsAmortRule,
                     &((*mbsDeal)->mbsAmortRule)) ISNT SUCCESS)
    {
        GtoErrMsg("%s : CopyMbsIoRule failed\n",routine);
        goto done;
    }

    /* Copy (already constructed!) prepay info */
    if (CopyMbsPpAssump(mbsPrepayAssump,
                       &((*mbsDeal)->mbsPrepayAssump)) ISNT SUCCESS)
    {
        GtoErrMsg("%s : CopyMbsPpAssump failed\n",routine);
        goto done;
    }
    for (i = 0; i < MBS_PP_NUM_REFI_TYPES; i++)
    {
        (*mbsDeal)->mbsGrCcSprd[i] = mbsGrCcSprd[i];
    }

    /* copy length of pay/accrual period */
    (*mbsDeal)->dealPayIntervalMons = dealPayIntervalMons;
    XONTRUE(dealPayIntervalMons < 1,
            "Invalid deal pay interval (months)");
    /* In current models, we require that deal & MBS accrual start
     * in the same calendar month */
    XONRANGE(dealMbsAccDelayDays,0,31,"dealMbsAccDelayDays");
    (*mbsDeal)->dealMbsAccDelayDays = dealMbsAccDelayDays;


    /* Are we using user-supplied lists? */
    if (numCashflows > 0 ||
        dealAccStDates ISNT NULL ||
        dealCfDates ISNT NULL)
    {
        (*mbsDeal)->numCashflows = numCashflows;

        /* make sure we have all needed info */
        if (numCashflows <= 0)
        {
            GtoErrMsg("%s : Must specify (positive) # "
                      "cashflows in supplied date lists\n",
                routine);
            goto done;
        }
        if (dealAccStDates IS NULL ||
            dealCfDates IS NULL)
        {
            GtoErrMsg("%s : Must supply both cashflow and accrual date lists",
                routine);
            goto done;
        }
        /* check these dates */
        for (iDate = 0; iDate < numCashflows; iDate++)
        {
            /* accrual must begin before cashflow */
            if (dealAccStDates[iDate] >= dealCfDates[iDate])
            {
                GtoErrMsg("%s : Cannot have deal accrual start date "
                          "on/after cashflow\n",
                    routine);
                goto done;
            }
            /* for all but last period, check that accrual end 
             * is on/before cf */
            if (iDate < numCashflows-1)
            {
                if (dealAccStDates[iDate+1] > dealCfDates[iDate])
                {
                    GtoErrMsg("%s : Cannot have deal accrual end AFTER cashflow\n",
                        routine);
                    goto done;
                }
            }
        } /* end iDate loop */
        /* Check that first deal cashflow occurs on/after first-ever
         * cashflow of MBS */
        XONFAIL(GetCfDateFromRegAccrStDate
                ((*mbsDeal)->mbsOrigDate,
                 (*mbsDeal)->mbsPayDelayDays,
                 &firstMbsCfDateEver) );
        if(dealCfDates[0] < firstMbsCfDateEver)
        {
            sprintf(outmesg,
                    "Problem: first deal cashflow (%s) occurs before\n"
                    "first-ever cashflow of this MBS (%s)",
                    GtoFormatDate(dealCfDates[0]),
                    GtoFormatDate(firstMbsCfDateEver));
            XONTRUE(TRUE,outmesg);
        }
        /* Although user didn't supply them, deduce accrual end dates */
        XONTRUE((dealAccEndDates = NEW_ARRAY(TDate,numCashflows)) IS NULL,
                "Failed to alloc temp accrual end date array" );
        for(iCf=0; iCf<numCashflows; iCf++)
        {
            /* for last accrual period, use CF date as accrual end
             * (this assumes deal (not MS) has no pay delay) */
            if(iCf IS numCashflows-1)
            {
                dealAccEndDates[iCf] = dealCfDates[iCf];
            }
            /* otherwise, use start of next accrual period
             * (assumes contiguous accrual periods) */
            else
            {
                dealAccEndDates[iCf] = dealAccStDates[iCf+1];
            }
        }

        /* Create date lists */
        XONTRUE
            (((*mbsDeal)->dealAccStDL = 
              GtoNewDateListFromDates(dealAccStDates,numCashflows)) IS NULL ||
             ((*mbsDeal)->dealAccEndDL = 
              GtoNewDateListFromDates(dealAccEndDates,numCashflows)) IS NULL ||
             ((*mbsDeal)->dealCfDL = 
              GtoNewDateListFromDates(dealCfDates,numCashflows)) IS NULL,
             "Failed to create cashflow/acc date lists");
        /* also create (approx) regular date lists */
        XONTRUE
            (((*mbsDeal)->dealEstRegAccStDL =
              GtoCopyDateList((*mbsDeal)->dealAccStDL)) IS NULL ||
             ((*mbsDeal)->dealEstRegAccEndDL =
              GtoCopyDateList((*mbsDeal)->dealAccEndDL)) IS NULL ||
             ((*mbsDeal)->dealEstRegCfDL =
              GtoCopyDateList((*mbsDeal)->dealCfDL)) IS NULL,
             "Failed to copy/create regular cf/acc date lists" );
        /* Now alter first/last of these dates
         * to correct for (possible) initial/final stubs 
         */
        /* Initial stub (only possible w/more than 1 cashflow) */
        if( numCashflows > 1 )
        {
            XONFAIL(NxtMth
                    ((*mbsDeal)->dealAccEndDL->fArray[0],
                     -(*mbsDeal)->dealPayIntervalMons,
                     &((*mbsDeal)->dealEstRegAccStDL->fArray[0])) );
        }
        /* Final stub: correct both accrual end and cf dates */
        XONFAIL(NxtMth
                ((*mbsDeal)->dealAccStDL->fArray[numCashflows-1],
                 (*mbsDeal)->dealPayIntervalMons,
                 &((*mbsDeal)->dealEstRegAccEndDL->fArray[numCashflows-1])) );
        XONFAIL(NxtMth
                ((*mbsDeal)->dealAccStDL->fArray[numCashflows-1],
                 (*mbsDeal)->dealPayIntervalMons,
                 &((*mbsDeal)->dealEstRegCfDL->fArray[numCashflows-1])) );
    } /* end of method 2: user-supplied lists */
    else
    {
        (*mbsDeal)->numCashflows = 0;

        /* Determine # days between MBS & deal cashflows */
        /* NB: in this case, can compute
         * (*mbsDeal)->dealMbsAccDelayDays
         * directly */
        /* @# do this later... */
/*
        GtoErrMsg("%s : Model currently cannot generate deal cf/accrual dates\n",
            routine);
        goto done;
*/
    } /* end of method 1: generate date lists */


    (*mbsDeal)->doingLegAIo = (legAIoRule ISNT NULL) &&
        !(IS_ALMOST_ZERO(legAIoRule->multiplier));
    (*mbsDeal)->doingLegAPo = (legAPoRule ISNT NULL) &&
        !(IS_ALMOST_ZERO(legAPoRule->multiplier));
    (*mbsDeal)->doingLegBIo = (legBIoRule ISNT NULL) &&
        !(IS_ALMOST_ZERO(legBIoRule->multiplier));
    (*mbsDeal)->doingLegBPo = (legBPoRule ISNT NULL) &&
        !(IS_ALMOST_ZERO(legBPoRule->multiplier));

    /* Copy rules for IO/PO of legs;
     * also null out IO/PO rule if multiplier is zero */
    if( (*mbsDeal)->doingLegAIo )
    {
        if (CopyMbsIoRule(legAIoRule,
                         &((*mbsDeal)->legAIoRule)) ISNT SUCCESS)
        {
            GtoErrMsg("%s : CopyMbsIoRule failed\n",routine);
            goto done;
        }
    }
    if( (*mbsDeal)->doingLegAPo )
    {
        if (CopyMbsPoRule(legAPoRule,
                         &((*mbsDeal)->legAPoRule)) ISNT SUCCESS)
        {
            GtoErrMsg("%s : CopyMbsPoRule failed\n",routine);
            goto done;
        }
    }
    if( (*mbsDeal)->doingLegBIo )
    {
        if (CopyMbsIoRule(legBIoRule,
                         &((*mbsDeal)->legBIoRule)) ISNT SUCCESS)
        {
            GtoErrMsg("%s : CopyMbsIoRule failed\n",routine);
            goto done;
        }
    }
    if( (*mbsDeal)->doingLegBPo )
    {
        if (CopyMbsPoRule(legBPoRule,
                         &((*mbsDeal)->legBPoRule)) ISNT SUCCESS)
        {
            GtoErrMsg("%s : CopyMbsPoRule failed\n",routine);
            goto done;
        }
    }

    status = SUCCESS;
done:
    FREE(dealAccEndDates);
    if (status ISNT SUCCESS)
    {
        FreeMbsDeal(mbsDeal);
    }
    return (status);
}


/***************************************************************************
 * FreeMbsDeal()
 * Destructor for TMbsDeal
 * (also resets pointer to NULL)
 ***************************************************************************/
void EXPORT
FreeMbsDeal
   (TMbsDeal **mbsDeal)         /* (I/O) structure to free */
{
    static char routine[] = "FreeMbsDeal";
    
    if (*mbsDeal ISNT NULL)
    {
        /* free internal dynamic vars */
        FreeMbsPoRule(&(*mbsDeal)->legBPoRule);
        FreeMbsIoRule(&(*mbsDeal)->legBIoRule);
        FreeMbsPoRule(&(*mbsDeal)->legAPoRule);
        FreeMbsIoRule(&(*mbsDeal)->legAIoRule);
        GtoFreeDateList((*mbsDeal)->dealCfDL);
        GtoFreeDateList((*mbsDeal)->dealAccStDL);
        GtoFreeDateList((*mbsDeal)->dealAccEndDL);
        GtoFreeDateList((*mbsDeal)->dealEstRegCfDL);
        GtoFreeDateList((*mbsDeal)->dealEstRegAccStDL);
        GtoFreeDateList((*mbsDeal)->dealEstRegAccEndDL);
        FreeMbsPrepayAssump(&(*mbsDeal)->mbsPrepayAssump);
        FreeMbsIoRule(&(*mbsDeal)->mbsAmortRule);
        FreeMbsPoRule(&(*mbsDeal)->mbsPoRule);
        FreeMbsIoRule(&(*mbsDeal)->mbsGrossIoRule);
        FreeMbsIoRule(&(*mbsDeal)->mbsNetIoRule);

        /* free entire structure */
        FREE(*mbsDeal);
        *mbsDeal = NULL;
    }
}


/***************************************************************************
 * CheckMbsDeal()
 * Checks contents of deal struc for validity
 * Returns SUCCESS/FAILURE
 ***************************************************************************/
int EXPORT
CheckMbsDeal
   (TMbsDeal *mbsDeal)     /* (I) structure to check */
{
#if 0
@# do this later...
    F_INIT("CheckMbsDeal");

    /* First, check ptr itself */
    XONTRUE( mbsDeal IS NULL, "No deal supplied" );
    /* @@ more checking here */
    XONFAIL( CheckMbsIoRule(mbsDeal->mbsIoRule) );
    XONFAIL( CheckMbsPoRule(mbsDeal->mbsPoRule) );

    /* check the legs */
    XONFAIL( CheckMbsIoRule(mbsDeal->legAIoRule) );
    XONFAIL( CheckMbsPoRule(mbsDeal->legAPoRule) );
    if ( (mbsDeal->legBIoRule)->multiplier )
    {
        XONFAIL( CheckMbsIoRule(mbsDeal->legBIoRule) );
        XONFAIL( CheckMbsPoRule(mbsDeal->legBPoRule) );
    }

    F_END;
#endif
    return(SUCCESS);
}

/***************************************************************************
 * PrintMbsDeal()
 * Pint out contents of deal struc for debugging
 * Returns SUCCESS/FAILURE
 ***************************************************************************/
int EXPORT
PrintMbsDeal
   (char        *label,         /* (I) label string telling us what
                                 * this deal applies to */
    long         maxNDates,     /* (I) max # date rows to list (saves space) */
    TMbsDeal *mbsDeal)          /* (I) structure contains MBS deal info. */
{
    F_INIT("PrintMbsDeal");
    long    iCf;
    long    iRefi;
    long    nRows;

    XONTRUE(mbsDeal IS NULL, 
            "Supplied null ptr" );

    GtoErrMsg
     ("------------------------------------------------------------------\n");
    GtoErrMsg("Deal: %s\n",label);
    GtoErrMsg
     ("------------------------------------------------------------------\n");

    GtoErrMsg("mbs type=%ld,  subtype=%s, agency=%ld, term=%ld\n",
       mbsDeal->mbsType, 
       mbsDeal->mbsSubtype,
       mbsDeal->mbsAgency,
       mbsDeal->mbsTerm);
    GtoErrMsg("origDate=%s,  matDate=%s,  effOrigTerm=%ld\n",
       GtoFormatDate(mbsDeal->mbsOrigDate),
       GtoFormatDate(mbsDeal->mbsMatDate),
       mbsDeal->mbsEffOrigTerm);
    GtoErrMsg("U/L MBS: accr st d-o-m=%ld,  CF d-o-m=%ld,  pay delay=%ld days\n",
       mbsDeal->mbsAccrualStDom,
       mbsDeal->mbsCashflowDom,
       mbsDeal->mbsPayDelayDays);
    GtoErrMsg("deal accrual/pay interval = %ld months\n",
       mbsDeal->dealPayIntervalMons);
    GtoErrMsg("# days between start of MBS and deal accrual: %ld days\n",
       mbsDeal->dealMbsAccDelayDays);
    GtoErrMsg("# accrual periods & cashflows: %ld\n",
       mbsDeal->numCashflows);
    GtoErrMsg("                                   ----EstRegularDates-------------\n");
    GtoErrMsg("AccrStDate AccrEnd    CfDate       AccrStDate AccrEnd    CfDate\n");
    nRows = MIN(mbsDeal->numCashflows,maxNDates);
    for(iCf=0; iCf<nRows; iCf++)
    {
        GtoErrMsg("%10s ",
                  (mbsDeal->dealAccStDL IS NULL ? "n/a" :
                   GtoFormatDate(mbsDeal->dealAccStDL->fArray[iCf])));
        GtoErrMsg("%10s ",
                  (mbsDeal->dealAccEndDL IS NULL ? "n/a" :
                   GtoFormatDate(mbsDeal->dealAccEndDL->fArray[iCf])));
        GtoErrMsg("%10s ",
                  (mbsDeal->dealCfDL IS NULL ? "n/a" :
                   GtoFormatDate(mbsDeal->dealCfDL->fArray[iCf])));
        GtoErrMsg("%10s ",
                  (mbsDeal->dealEstRegAccStDL IS NULL ? "n/a" :
                   GtoFormatDate(mbsDeal->dealEstRegAccStDL->fArray[iCf])));
        GtoErrMsg("%10s ",
                  (mbsDeal->dealEstRegAccEndDL IS NULL ? "n/a" :
                   GtoFormatDate(mbsDeal->dealEstRegAccEndDL->fArray[iCf])));
        GtoErrMsg("%10s ",
                  (mbsDeal->dealEstRegCfDL IS NULL ? "n/a" :
                   GtoFormatDate(mbsDeal->dealEstRegCfDL->fArray[iCf])));
        GtoErrMsg("\n");
    } /* iCf */
    if(nRows < mbsDeal->numCashflows)
    {
        GtoErrMsg(" ...\n");
    }

    if(mbsDeal->doingLegAIo)
    {
        XONFAIL(PrintMbsIoRule("Leg A",maxNDates,
                               mbsDeal->legAIoRule));
    }
    if(mbsDeal->doingLegAPo)
    {
        XONFAIL(PrintMbsPoRule("Leg A",
                               mbsDeal->legAPoRule));
    }
    if(mbsDeal->doingLegBIo)
    {
        XONFAIL(PrintMbsIoRule("Leg B",maxNDates,
                               mbsDeal->legBIoRule));
    }
    if(mbsDeal->doingLegBPo)
    {
        XONFAIL(PrintMbsPoRule("Leg B",
                               mbsDeal->legBPoRule));
    }

    XONFAIL(PrintMbsIoRule("Underlying MBS Net Cpn ",
                           maxNDates,
                           mbsDeal->mbsNetIoRule));
    XONFAIL(PrintMbsIoRule("Underlying MBS Gross Cpn",
                           maxNDates,
                           mbsDeal->mbsGrossIoRule));
    XONFAIL(PrintMbsPoRule("Underlying MBS PO Rule",
                           mbsDeal->mbsPoRule));
    XONFAIL(PrintMbsIoRule("Underlying MBS Amort",
                           maxNDates,
                           mbsDeal->mbsAmortRule));

    XONFAIL(PrintMbsPpAssump(FALSE,mbsDeal->mbsPrepayAssump));
    GtoErrMsg("Gross curr cpn spreads: ");
    for(iRefi=0; iRefi<MBS_PP_NUM_REFI_TYPES; iRefi++)
    {
        GtoErrMsg("%.2lfbp ",
           mbsDeal->mbsGrCcSprd[iRefi]*10000.);
    }
    GtoErrMsg("\n");

    F_END;
}


/***************************************************************************
 * MakeMbsIoRule();
 * Constructor for TMbsIoRule
 * Sets internals of structure to input values
 * If input ptr is NULL, also allocates structure
 * Returns SUCCESS/FAILURE
 ***************************************************************************/
int EXPORT
MakeMbsIoRule
   (double        multiplier,          /* (I) all IO cashflows scaled by this */
    /* CPN BEHAVIOR */
    double        cpnIndexMult,        /* (I) multiplier of index rate for cpn;
                                        * use 0.0 for fixed-rate */
    double        cpnRateOrSprd,       /* (I) for fixed cpn: cpn rate;
                                        * for float cpn: sprd to index rate */
    long          cpnIndexFreq,        /* (I) (for flt cpn) pay freq (cpns/yr)
                                        * of index rate (e.g., 2=semiann) */
    long          cpnIndexMatMons,     /* (I) (for flt cpn) maturity (months) of
                                        * index rate (e.g., 12 for 1yr) */
    long          cpnIndexCurve,       /* (I) (for flt cpn) which zero curve 
                                        * to use for cpn (e.g., 1,2) */
    long          cpnIndexDC,          /* (I) (for flt cpn) day count of index
                                        * rate (e.g., GTO_B30_360) */
    long          cpnPayDC,            /* (I) For fixed or flt: day count used 
                                        * in paying cpn (e.g., GTO_B30_360) */
    TBoolean      ruleIsAmortRule,     /* (I) TRUE if IO rule is amort reset info;
                                        * FALSE for normal cpn reset info */
    double        constServSpread,     /* (I) @@ cheat for the gross coupon:
                                        * if the MBS is floating, we often
                                        * approximate curr WAC assuming a
                                        * const serv sprd: wac=currnet+sprd */
    /* FOR FLOATING CPN: KNOWN CPNS */
    long          nKnownCpns,          /* (I) # known coupons 
                                        * (I) (i.e., #values in arrays below) */
    TDate        *knownEffResetDates,  /* (I) Array of effective reset dates
                                        * for each known cpn;
                                        * NB: if this constructor builds list
                                        * of reset dates, these known reset
                                        * dates must be consistent with the
                                        * reset date parameters: first date,
                                        * reset interval, reset d-o-m */
    double       *knownCpns,           /* (I) Array of known coupon
                                        * rates (monthly, decimal) */
    double       *knownResetIndexRates, /* (I) Array of known index rates
                                        * index rates (as of reset date);
                                        * NB: should be raw index rate 
                                        * (e.g., CMT1yr), not including
                                        * cpnIndexMult or cpnRateOrSprd */
    /* FOR FLOATING CPN: ROUNDING */
    double        cpnRounding,         /* (I) reset coupons are rounded to
                                        * nearest multiple of this value;
                                        * use zero for no rounding;
                                        * e.g., GNMA ARM resets involve rounding
                                        * cpns to nearest 1/8th of a point, 
                                        * hence would use 0.00125 */
    /* FLOATING CPN LOOKBACK "RULE" 
     * (always supply this info, even if user supplied reset date list) */
    long          effResetLkbkDays,    /* effective reset date is (approx) this
                                        * many (act) days before start of REGULAR
                                        * accrual period of new cpn 
					* NB: if IO rule describes deal, we mean
					* delay of reset prior to DEAL accrual;
					* if rule is for MBS, then MBS accrual
                                        * NB: if user supplies reset dates
                                        * (in effResetDates[]), this delay
                                        * should be approx. consistent with 
                                        * these reset dates */
    /* FLOATING CPN RESET DATES: 
     * We allow two methods of defining reset dates:
     *   1) "By rule": User can supply these "general" reset 
     *   parameters--constructor will then build internal lists).... */
    TDate         firstEffResetDate,   /* (I) Effective reset date of first
                                        * reset for this IO cpn;
                                        * NB: if user supplies reset dates,
                                        * (meth 2) this param will be 0 */
    long          resetIntervalMons,   /* (I) # months between resets;
                                        * NB: if user supplies reset dates,
                                        * (meth 2) this param will be 0 */
    long          effResetDOM,         /* (I) day-of-month of eff. reset date 
                                        * NB: if user supplies reset dates,
                                        * (meth 2) this param will be 0 */
    /*   2) "By list": user can supply list of reset dates (effResetDates[]) */
    long          numResets,           /* (I) # of resets (& # elements
                                        * in the arrays below);
                                        * use 0 if using meth 1 */
    TDate        *effResetDates,       /* (I) list of effective reset dates;
                                        * NB: all dates should be unique
                                        * (no repeated resets);
                                        * NULL if reset dates are to generated
                                        * automatically (meth 1) */
    TDate        *firstAccStDates,     /* (I) list of first accrual start dates
                                        * corresponding to each reset above;
                                        * NB: if this IO rule describes the
                                        * underlying MBS, these will be
                                        * MBS, not deal, accrual pds */
    TBoolean      firstAccMayBeStub,   /* (I) if TRUE, firstAccStDates[0]
                                        * may be start of stub period */
    /* FLOATING CPN STRIKES (can be NULL):
     *   Depends on method of specifying reset dates:
     *   1) each pointer assumed to point at a single sprd/strike
     *   2) ptr assumed to be an array of numResets sprds/strikes */
    double       *stickyCapSpreads,    /* (I) Cpn spreads of periodic cap for 
                                        * each reset, in decimal (e.g., 0.01) */
    double       *stickyFlrSpreads,    /* (I) Cpn spreads of periodic floor for 
                                        * each reset, in decimal (e.g., -0.01) */
    double       *capStrikes,          /* (I) Normal cap strike for each reset
                                        * (decimal, e.g., 0.11)*/
    double       *flrStrikes,          /* (I) Normal flr strike for each reset
                                        * (decimal, e.g., 0.01)*/
    /* MISC */
    char         *holidayFile,         /* (I) name of GTO holiday file */
    long         *numSettleDays,       /* (I) array of settle days
                                        * for each zero curve */
    TDate         todayDate,           /* (I) (for determining truly historical
                                        * resets) current date */
    /* PTR TO NEW/EXISTING STRUCTURE */
    TMbsIoRule  **mbsIoRule)           /* (I/O) structure to alloc/modify */
{
    F_INIT("MakeMbsIoRule");
    long    iKnown;
    long    iReset;
    long    tmpDom;
    TDate  lookupDate;
    TDateInterval lkbkInterval;
    TBoolean  foundOrigCpnReset;
    static TDate dummyDate = 30000; /* @# */

    /* Check inputs
     */
    /* (checking done below) */

    /* Prepare an empty structure (may have to de-alloc existing one) */
    if (*mbsIoRule ISNT NULL)
    {
        FreeMbsIoRule(mbsIoRule);
    }
    XONTRUE(((*mbsIoRule = NEW(TMbsIoRule)) IS NULL),
            "Failed to allocate mbsIoRule structure");

    /* Before anything, reset internal ptrs to NULL,
     * to indicate that they haven't been allocated yet
     */
    (*mbsIoRule)->cpnIndexRate = NULL;
    (*mbsIoRule)->effResetDL = NULL;
    (*mbsIoRule)->firstAccStDL = NULL;
    (*mbsIoRule)->firstEstRegAccStDL = NULL;
    (*mbsIoRule)->stickyCapSpreads = NULL;
    (*mbsIoRule)->stickyFlrSpreads = NULL;
    (*mbsIoRule)->capStrikes = NULL;      
    (*mbsIoRule)->flrStrikes = NULL;      
    (*mbsIoRule)->isKnownCpn = NULL;
    (*mbsIoRule)->knownCpns = NULL;
    (*mbsIoRule)->knownCpnIndexRates = NULL;
    
    /* Set member data that is always needed, regardless of fixed vs. float
     */
    (*mbsIoRule)->multiplier = multiplier;
    (*mbsIoRule)->cpnPayDC = cpnPayDC;
    (*mbsIoRule)->ruleIsAmortRule = ruleIsAmortRule;
    (*mbsIoRule)->constServSpread = constServSpread;
    strncpy((*mbsIoRule)->holidayFile,holidayFile,MAXPATHLEN-1);
    (*mbsIoRule)->todayDate = todayDate;

    /* For fixed cpn */
    if( IS_ALMOST_ZERO(cpnIndexMult) )
    {
        /* easy--only 2 things to set */
        (*mbsIoRule)->cpnIsFloating = FALSE;
        (*mbsIoRule)->fixedCoupon = cpnRateOrSprd;
        /* zero out float-cpn info (ptrs already null) */
        (*mbsIoRule)->useCpnRounding = FALSE;
        (*mbsIoRule)->cpnRounding = 0.;
        (*mbsIoRule)->effResetLkbkDays = 0;
        (*mbsIoRule)->firstEffResetDate = 0;
        (*mbsIoRule)->resetIntervalMons = 0;
        (*mbsIoRule)->effResetDOM = 0;
        (*mbsIoRule)->numResets = 0;
    }

    else /* floating cpn */
    {
        (*mbsIoRule)->cpnIsFloating = TRUE;
        (*mbsIoRule)->fixedCoupon = 0.; /* not used */
        /* Define index rate */
        XONRANGE(cpnIndexFreq,1,12,"cpn index freq");
        XONRANGE(cpnIndexMatMons,1,120,"index mat mons" );
        XONRANGE(cpnPayDC,1,10,"cpn pay DC");
        XONFAIL(SetSingleFloatArray
                (cpnIndexFreq,    /* (I) cpns per year */
                 cpnIndexMatMons, /* (I) maturity in months */
                 cpnIndexDC,      /* (I) DC convention */
                 cpnIndexCurve,   /* which curve */
                 numSettleDays[cpnIndexCurve],
                 &((*mbsIoRule)->cpnIndexRate)) );
        (*mbsIoRule)->cpnIndexRate->defs[0].weight = cpnIndexMult;
        (*mbsIoRule)->cpnIndexRate->defs[0].spread = cpnRateOrSprd;
        /* rounding */
        XONTRUE(cpnRounding < 0., 
                "Cpn rounding multiple cannot be negative" );
        (*mbsIoRule)->useCpnRounding = ! IS_ALMOST_ZERO(cpnRounding);
        (*mbsIoRule)->cpnRounding = cpnRounding;
        /* always define lookback rule info */
        (*mbsIoRule)->effResetLkbkDays = effResetLkbkDays;
        /* If user supplied a list of resets, copy them */
        if(numResets > 0 &&
           effResetDates ISNT NULL)
        {
            /* first, check that we have at most one reset date = 0
             * (i.e., corresp. to original cpn),
             * and if so, it's the first reset date */
            foundOrigCpnReset = FALSE;
            for(iReset=0; iReset<numResets; iReset++)
            {
                if(effResetDates[iReset] IS 0)
                {
                    XONTRUE(iReset ISNT 0,
                       "Problem: detected 'original cpn reset'\n"
                       "(i.e., reset date=0) for reset that is NOT first reset");
                    XONTRUE(foundOrigCpnReset,
                       "Problem: found more than one 'original cpn reset'\n"
                       "(i.e., reset date=0)");
                    foundOrigCpnReset = TRUE;
                }
            }
            /* reset these "general" reset timing params,
             * since they're not used in this case */
            (*mbsIoRule)->firstEffResetDate = 0;
            (*mbsIoRule)->resetIntervalMons = 0;
            (*mbsIoRule)->effResetDOM = 0;
            /* copy lists */
            XONTRUE(numResets <= 0,
                    "Must supply at least one reset date" );
            (*mbsIoRule)->numResets = numResets;
            XONTRUE(((*mbsIoRule)->effResetDL =
                     GtoNewDateListFromDates(effResetDates,numResets))
                    IS NULL,
                    "Failed to copy reset date list" );
            XONTRUE(((*mbsIoRule)->firstAccStDL =
                     GtoNewDateListFromDates(firstAccStDates,numResets))
                    IS NULL,
                    "Failed to copy 1st acc st date list" );
            /* Now copy acc st dates to est-reg array of these */
            XONTRUE( ((*mbsIoRule)->firstEstRegAccStDL =
                      GtoCopyDateList((*mbsIoRule)->firstAccStDL)) IS NULL,
                    "Failed to copy actual accr. start date list" );
            /* if first accrual start date may be stub,
             * fiddle w/first date to approx reg accr. st date*/
            if( firstAccMayBeStub )
            {
                /* if possible, use 2nd acrual date,
                 * and go back by diff between 2nd/3rd accrual dates */
                if( numResets > 2 )
                {
                    (*mbsIoRule)->firstEstRegAccStDL->fArray[0] =
                        (*mbsIoRule)->firstEstRegAccStDL->fArray[1]
                      -((*mbsIoRule)->firstEstRegAccStDL->fArray[2]
                       -(*mbsIoRule)->firstEstRegAccStDL->fArray[1]);
                    XONFAIL(GetDOM((*mbsIoRule)->firstEstRegAccStDL->fArray[1],
                                   &tmpDom));
                    XONFAIL(ForceToClosestDOM
                            (tmpDom,
                             TRUE,
                             &((*mbsIoRule)->firstEstRegAccStDL->fArray[0])));
                }
                /* otherwise, use lkbk interval & reset date */
                else
                {
                    SET_TDATE_INTERVAL(lkbkInterval,
                                       (*mbsIoRule)->effResetLkbkDays,
                                       'D');
                    XONFAIL(GtoDtFwdAny
                            ((*mbsIoRule)->effResetDL->fArray[0],
                             &lkbkInterval,
                             &((*mbsIoRule)->firstEstRegAccStDL->fArray[0])) );
                }
            }
            XONFAIL( MbsCopyDArray(
                                   stickyCapSpreads,
                                   numResets,
                                   1., /* def value */
                                   &((*mbsIoRule)->stickyCapSpreads)) );
            XONFAIL( MbsCopyDArray(
                                   stickyFlrSpreads,
                                   numResets,
                                   -1., /* def value */
                                   &((*mbsIoRule)->stickyFlrSpreads)) );
            XONFAIL( MbsCopyDArray(
                                   capStrikes,
                                   numResets,
                                   1000., /* def value */
                                   &((*mbsIoRule)->capStrikes)) );
            XONFAIL( MbsCopyDArray(
                                   flrStrikes,
                                   numResets,
                                   0., /* def value */
                                   &((*mbsIoRule)->flrStrikes)) );
        }
        /* otherwise, build them */
        else
        {                
            /* @# For now, we have not included ability to
             * build the reset lists... */
            /* @# scalars dummied for now */
            (*mbsIoRule)->firstEffResetDate = firstEffResetDate;
            (*mbsIoRule)->resetIntervalMons = resetIntervalMons;
            (*mbsIoRule)->effResetDOM = effResetDOM;
            /* @# dummied: 1 reset */
            (*mbsIoRule)->numResets = 1;
            XONTRUE(((*mbsIoRule)->effResetDL =
                     GtoNewDateListFromDates(&firstEffResetDate,
                                             (*mbsIoRule)->numResets))
                    IS NULL,
                    "Failed to copy reset date list" );
            XONTRUE(((*mbsIoRule)->firstAccStDL =
                     GtoNewDateListFromDates(&firstEffResetDate,
                                             (*mbsIoRule)->numResets))
                    IS NULL,
                    "Failed to copy 1st acc st date list" );
            XONTRUE( ((*mbsIoRule)->firstEstRegAccStDL =
                      GtoCopyDateList((*mbsIoRule)->firstAccStDL)) IS NULL,
                    "Failed to copy actual accr. start date list" );

            XONTRUE(((*mbsIoRule)->stickyCapSpreads = 
                     NEW_ARRAY(double,
                               (*mbsIoRule)->numResets)) IS NULL,
                    "Failed to alloc sticky cap sprd array");
            XONTRUE(((*mbsIoRule)->stickyFlrSpreads = 
                     NEW_ARRAY(double,
                               (*mbsIoRule)->numResets)) IS NULL,
                    "Failed to alloc sticky flr sprd array");
            XONTRUE(((*mbsIoRule)->capStrikes = 
                     NEW_ARRAY(double,
                               (*mbsIoRule)->numResets)) IS NULL,
                    "Failed to alloc cap strike array");
            XONTRUE(((*mbsIoRule)->flrStrikes = 
                     NEW_ARRAY(double,
                               (*mbsIoRule)->numResets)) IS NULL,
                    "Failed to alloc flr strike array");


            (*mbsIoRule)->stickyCapSpreads[0] = stickyCapSpreads[0];
            (*mbsIoRule)->stickyFlrSpreads[0] = stickyFlrSpreads[0];
            (*mbsIoRule)->capStrikes[0] = capStrikes[0];
            (*mbsIoRule)->flrStrikes[0] = flrStrikes[0];
        }

        /* Once reset date list ready, copy known cpns/rates into arrays
         * (here, we only take the known values for resets
         * in our list) */
        XONTRUE(((*mbsIoRule)->isKnownCpn = 
                 NEW_ARRAY(TBoolean,(*mbsIoRule)->numResets)) IS NULL,
                "Failed to alloc is-known-cpn array" );
        XONTRUE(((*mbsIoRule)->knownCpns = 
                 NEW_ARRAY(double,(*mbsIoRule)->numResets)) IS NULL,
                "Failed to alloc known cpns array" );
        XONTRUE(((*mbsIoRule)->knownCpnIndexRates =
                 NEW_ARRAY(double,(*mbsIoRule)->numResets)) IS NULL,
                "Failed to alloc known idx rates array" );
        /* reset the "knowns" to -1., to indicate unknown */
        for(iReset=0; iReset<(*mbsIoRule)->numResets; iReset++)
        {
            (*mbsIoRule)->isKnownCpn[iReset] = FALSE;
            (*mbsIoRule)->knownCpns[iReset] = -1.;
            (*mbsIoRule)->knownCpnIndexRates[iReset] = -1.;
        }
	/* verify that user has supplied ptrs to these arrays */
	if( nKnownCpns > 0 )
	{
	    XONTRUE(knownCpns IS NULL,
		    "Must supply array of known cpns to IO rule constructor" );
	    XONTRUE(knownResetIndexRates IS NULL,
		    "Must supply array of known reset index rates to IO rule constructor" );
	}
        /* now copy the known rates supplied by user
         * into these internal arrays */
        for(iKnown=0; iKnown<nKnownCpns; iKnown++)
        {
            /* "known" coupons must be non-negative rates;
             * a value of -1.0 is ignored */
            if(knownCpns[iKnown] >= 0. ||
               knownResetIndexRates[iKnown] >= 0. )
            {
		/* check for -1, and for range */
                XONTRUE(knownCpns[iKnown] < 0.,
                        "For every known coupon, must supply known coupon rate" );
		XONRANGE(knownCpns[iKnown],0.,1., 
			 "known cpn value" );
                if(knownEffResetDates[iKnown] IS 0)
                {
                    /* special case: for initial cpn on instrument,
                     * there is no corresponding index rate, 
                     * so don't check input rate, but instead ignore it
                     * and create fake "index rate" = cpn - margin */
                    knownResetIndexRates[iKnown] = 
                        (knownCpns[iKnown] - cpnRateOrSprd) / cpnIndexMult;
                }
                else
                {
                    /* check for -1, and for range */
                    XONTRUE(knownResetIndexRates[iKnown] < 0.,
                            "For known coupon, must supply known reset index rate" );
                    XONRANGE(knownResetIndexRates[iKnown],0.,1., 
                             "known reset index rate value" );
                }

                /* search for the matching reset date */
                for(iReset=0; iReset<(*mbsIoRule)->numResets; iReset++)
                {
                    if((*mbsIoRule)->effResetDL->fArray[iReset] IS
                       knownEffResetDates[iKnown])
                    {
                        (*mbsIoRule)->isKnownCpn[iReset] = TRUE;
                        (*mbsIoRule)->knownCpns[iReset] =
                            knownCpns[iKnown];
                        /* Important: in structure, we store the
                         * "net" index rate (mult*rawrate+sprd),
                         * to be consistent w/the index rate defn
                         * in cpnIndexRate */
                        (*mbsIoRule)->knownCpnIndexRates[iReset] = 
                            cpnIndexMult * knownResetIndexRates[iKnown] +
                                cpnRateOrSprd;
                        break;
                    }
                }
            }
        }
        /* Now check that every reset with a lookup date
         * on/before today has a known cpn */
        for(iReset=0; iReset<(*mbsIoRule)->numResets; iReset++)        
        {
            XONFAIL(GtoDateFromBusDaysOffset
                    ((*mbsIoRule)->effResetDL->fArray[iReset],
                     -(*mbsIoRule)->cpnIndexRate->defs[0].spotOffset.interval.prd,
                     (*mbsIoRule)->holidayFile,
                     &lookupDate));
            if(lookupDate <= todayDate)
            {
                if( ! (*mbsIoRule)->isKnownCpn[iReset] )
                {
                    sprintf(outmesg,
                         "Reset on %s (lkup=%s) is on/before "
                         "todayDate (%s), yet no known cpn was supplied\n",
                         GtoFormatDate((*mbsIoRule)->effResetDL->fArray[iReset]),
                         GtoFormatDate(lookupDate),
                         GtoFormatDate(todayDate));
                    XONTRUE(TRUE,outmesg);
                }
            }
            else /* lkup after today */
            {
                /* if lkup date after today, quit this loop */
                break;
            }
        }

        /* Check that caps/floors make sense */
        for(iReset=0; iReset<(*mbsIoRule)->numResets; iReset++)        
        {
            XONTRUE((*mbsIoRule)->capStrikes[iReset] < 
                    (*mbsIoRule)->flrStrikes[iReset],
                    "All cap strikes must be at/above floor strikes" );
            XONTRUE((*mbsIoRule)->stickyCapSpreads[iReset] < 0.,
                    "All sticky cap spreads must be >= 0" );
            XONTRUE((*mbsIoRule)->stickyFlrSpreads[iReset] > 0.,
                    "All sticky floor spreads must be <= 0" );
        }
        
    } /* end of floating cpn */


    status = SUCCESS;
done:
    if (status ISNT SUCCESS)
    {
        FreeMbsIoRule(mbsIoRule);
    }
    return (status);
}




/***************************************************************************
 * CopyMbsIoRule();
 * Sets internals of structure to input values
 * If input ptr is NULL, also allocates structure
 * Returns SUCCESS/FAILURE
 ***************************************************************************/
int EXPORT
CopyMbsIoRule
   (TMbsIoRule  *inIoRule,   /* (I) IO cashflow rule struc */
    TMbsIoRule **outIoRule)  /* (O) copy of IO cashflow rule struc */
{
    F_INIT("CopyMbsIoRule");
    long    iReset;
    long    iKnown;
    long   nKnownCpns;
    long   numSettleDays[20];
    TDate  *knownEffResetDates = NULL;
    double *knownCpns = NULL;
    double *knownResetIndexRates = NULL;


    XONTRUE(outIoRule IS NULL,
            "Must supply ptr to ptr for copy of IO rule" );

    /* trivial case: if input rule is null, copy this and quit */
    if(inIoRule IS NULL)
    {
        *outIoRule = NULL;
        status = SUCCESS;
        goto done;
    }

    /* Use the normal constructor, though call it in a way
     * that depends on fixed/float nature of IO rule */
    if( inIoRule->cpnIsFloating )
    {
        /* Set # settle days */
        numSettleDays[inIoRule->cpnIndexRate->curveIndices[0]] =
            inIoRule->cpnIndexRate->defs[0].spotOffset.interval.prd;

        /* For floating, we have to extract an array of known cpns;
         * first, count # known cpns */
        nKnownCpns = 0;
        for(iReset=0; iReset<inIoRule->numResets; iReset++)        
        {
            if(inIoRule->isKnownCpn[iReset])
            {
                nKnownCpns++;
            }
        }
        /* Now allocate and fill these arrays with the (few? none?)
         * known cpns */
        if( nKnownCpns > 0 )
        {
            XONTRUE( (knownEffResetDates =
                      NEW_ARRAY(TDate,nKnownCpns)) IS NULL,
                    "Failed to alloc known eff reset dates" );
            XONTRUE( (knownCpns =
                      NEW_ARRAY(double,nKnownCpns)) IS NULL,
                    "Failed to alloc known cpn array" );
            XONTRUE( (knownResetIndexRates =
                      NEW_ARRAY(double,nKnownCpns)) IS NULL,
                    "Failed to alloc known index rate array" );
            iKnown = 0;
            for(iReset=0; iReset<inIoRule->numResets; iReset++)        
            {
                if(inIoRule->isKnownCpn[iReset])
                {
                    knownEffResetDates[iKnown] = 
                        inIoRule->effResetDL->fArray[iReset];
                    knownCpns[iKnown] =
                        inIoRule->knownCpns[iReset];
                    /* Important: convert cpn reset to raw rate;
                     * since this must be floating cpn, mult should
                     * be non-zero, but check just in case */
                    XONTRUE(IS_ALMOST_ZERO
                            (inIoRule->cpnIndexRate->defs[0].weight),
                            "Cannot have zero index mult for floating cpn" );
                    knownResetIndexRates[iKnown] =
                        (inIoRule->knownCpnIndexRates[iReset] - 
                         inIoRule->cpnIndexRate->defs[0].spread) /
                            inIoRule->cpnIndexRate->defs[0].weight;
                    iKnown++;
                }
            }
        }

        /* now make the copy */
        XONFAIL(MakeMbsIoRule
                (inIoRule->multiplier,
                 inIoRule->cpnIndexRate->defs[0].weight,
                 inIoRule->cpnIndexRate->defs[0].spread,
                 (long) ROUND(12./
                        inIoRule->cpnIndexRate->defs[0].payInterval.prd),
                 /* here we can safely assume mat interval = int # mons */
                 inIoRule->cpnIndexRate->defs[0].matInterval.prd,
                 inIoRule->cpnIndexRate->curveIndices[0],
                 inIoRule->cpnIndexRate->defs[0].dayCountConv,
                 inIoRule->cpnPayDC,
                 inIoRule->ruleIsAmortRule,
                 inIoRule->constServSpread,
                 /* here, use our local list of known cpns */
                 nKnownCpns,           
                 knownEffResetDates,
                 knownCpns,
                 knownResetIndexRates,
                 inIoRule->cpnRounding,
                 inIoRule->effResetLkbkDays,
                 inIoRule->firstEffResetDate,
                 inIoRule->resetIntervalMons,
                 inIoRule->effResetDOM,
                 inIoRule->numResets,
                 inIoRule->effResetDL->fArray,
                 inIoRule->firstAccStDL->fArray,
                 /* if reg acc st doesn't match actual, we have stub */
                 ((inIoRule->firstAccStDL->fArray IS NULL) ? 
                    FALSE :             /* (ignored) */
                   (inIoRule->firstAccStDL->fArray[0] ISNT
                    inIoRule->firstEstRegAccStDL->fArray[0])),
                 inIoRule->stickyCapSpreads,
                 inIoRule->stickyFlrSpreads,
                 inIoRule->capStrikes,
                 inIoRule->flrStrikes,
                 "NONE",               /* simple "holiday file" is enough */
                 numSettleDays,       /* no need to check when copying, */
                 inIoRule->todayDate,
                 outIoRule));
    }
    else                               /* fixed */
    {
        /* for fixed, most inputs not used */
        XONFAIL(MakeMbsIoRule
                (inIoRule->multiplier,
                 0.,                   /* index mult */
                 inIoRule->fixedCoupon,
                 0,                    /* freq */
                 0,                    /* mat */
                 0,                    /* curve */
                 0,                    /* index dc */
                 inIoRule->cpnPayDC,
                 inIoRule->ruleIsAmortRule,
                 inIoRule->constServSpread,
                 0,                    /* # known */
                 NULL,
                 NULL,
                 NULL,
                 0.,                   /* rounding */
                 0,                    /* lkbk days */
                 0,                    /* 1st reset */
                 0,                    /* intervalmons */
                 0,                    /* dom */
                 0,                    /* # resets */
                 NULL,                 /* reset dates */
                 NULL,                 /* first acc st dates */
                 FALSE,                /* 1st acc st stub?--doesn't matter */
                 NULL,                 /* sprds */
                 NULL,                 /* sprds */
                 NULL,                 /* strikes */
                 NULL,                 /* strikes */
                 "NONE",               /* holiday file */
                 numSettleDays, /* ignored */
                 inIoRule->todayDate,
                 outIoRule));
    }


    status = SUCCESS;
done:
    FREE(knownEffResetDates);
    FREE(knownCpns);
    FREE(knownResetIndexRates);

    if (status ISNT SUCCESS)
    {
        FreeMbsIoRule(outIoRule);
    }
    return (status);
}



/***************************************************************************
 * FreeMbsIoRule()
 * Destructor for TMbsIoRule
 * (also resets pointer to NULL)
 ***************************************************************************/
void EXPORT
FreeMbsIoRule
   (TMbsIoRule **mbsIoRule)         /* (I/O) structure to free */
{
    if (*mbsIoRule ISNT NULL)
    {
        /* free internal dynamic vars */
        FREE((*mbsIoRule)->knownCpnIndexRates);
        FREE((*mbsIoRule)->knownCpns);
        FREE((*mbsIoRule)->isKnownCpn);
        FREE((*mbsIoRule)->flrStrikes);
        FREE((*mbsIoRule)->capStrikes);
        FREE((*mbsIoRule)->stickyFlrSpreads);
        FREE((*mbsIoRule)->stickyCapSpreads);
        GtoFreeDateList((*mbsIoRule)->effResetDL);
        GtoFreeDateList((*mbsIoRule)->firstAccStDL);
        GtoFreeDateList((*mbsIoRule)->firstEstRegAccStDL);
        if( (*mbsIoRule)->cpnIndexRate ISNT NULL )
        {
            GtoFloatRateArrayFree( (*mbsIoRule)->cpnIndexRate );
        }
        /* free entire structure */
        FREE(*mbsIoRule);
        *mbsIoRule = NULL;
    }
    return;
}


/***************************************************************************
 * CheckMbsIoRule()
 * Checks contents of io rule struc for validity
 * Returns SUCCESS/FAILURE
 ***************************************************************************/
int EXPORT
CheckMbsIoRule
   (TMbsIoRule *mbsIoRule)     /* (I) structure to check */
{
    F_INIT("CheckMbsIoRule");

    /* @# for now, just check pointer */
    XONTRUE( mbsIoRule IS NULL, "No IO rule supplied" );

    F_END;
}


/***************************************************************************
 * PrintMbsIoRule()
 * Pint out contents of deal struc for debugging
 * Returns SUCCESS/FAILURE
 ***************************************************************************/
int EXPORT
PrintMbsIoRule
   (char        *label,         /* (I) label string telling us what
                                 * this IO rule applies to */
    long         maxNDates,     /* (I) max # date rows to list (saves space) */
    TMbsIoRule  *mbsIoRule)     /* (I) structure contains IO rule info. */
{
    F_INIT("PrintMbsIoRule");
    long   iReset;
    long   nRows;
    char   cpnPayDcStr[36];

    XONTRUE(mbsIoRule IS NULL, 
            "Supplied null ptr" );

    GtoErrMsg("------------------------------------------------------------------\n");
    GtoErrMsg("IO Rule: %s\n",label);
    GtoErrMsg("------------------------------------------------------------------\n");
    XONFAIL(GetDayCountStr(mbsIoRule->cpnPayDC,cpnPayDcStr));
    GtoErrMsg("IO mult=%lf,  payDC=%s,  const serv sprd=%.1lfbp\n%s",
              mbsIoRule->multiplier,
              cpnPayDcStr,
              mbsIoRule->constServSpread*10000.,
              (mbsIoRule->ruleIsAmortRule ? "(is amort rule)\n" : ""));
    /* fixed */
    if(! mbsIoRule->cpnIsFloating)
    {
        GtoErrMsg("FIXED CPN:  cpn=%.3lf%%\n",
           mbsIoRule->fixedCoupon*100.);
    }
    /* floating */
    else                        
    {
        
        GtoErrMsg("FLOATING CPN");
        GtoFloatRateArrayPrint(mbsIoRule->cpnIndexRate,
                               "IO rule index rate");
        if(mbsIoRule->useCpnRounding)
        {
            GtoErrMsg("Rounding used to nearest %.3lf%%\n",
               mbsIoRule->cpnRounding*100.);
        }
        GtoErrMsg("Reset lookback = %ld days\n",
           mbsIoRule->effResetLkbkDays);
        GtoErrMsg("# reset dates = %ld\n",
           mbsIoRule->numResets);
        if(mbsIoRule->numResets > 0)
        {

            GtoErrMsg
                ("RESETS:                                      "
                 "                      ---IsKnown----\n");
            GtoErrMsg
                ("  effReset  1stAccrSt      (reg)   pdCap    "
                 "pdFlr      Cap    Flr    Cpn  Indx+Marg\n");
            nRows = MIN(mbsIoRule->numResets,maxNDates);
            for(iReset=0; iReset<nRows; iReset++)
            {
                sprintf(outmesg,
                     "%10s %10s %10s %6.2lf%% %7.2lf%% %7.2lf%% %5.2lf%%",
                     (mbsIoRule->effResetDL->fArray[iReset] IS 0 ?
                      "    (orig)" :
                      GtoFormatDate(mbsIoRule->effResetDL->fArray[iReset])),
                     GtoFormatDate(mbsIoRule->firstAccStDL->fArray[iReset]),
                     GtoFormatDate(mbsIoRule->firstEstRegAccStDL->fArray[iReset]),
                     (mbsIoRule->stickyCapSpreads IS NULL ? 9.99 :
                      mbsIoRule->stickyCapSpreads[iReset]*100.),
                     (mbsIoRule->stickyFlrSpreads IS NULL ? -9.99 :
                      mbsIoRule->stickyFlrSpreads[iReset]*100.),
                     (mbsIoRule->capStrikes IS NULL ? 99.99 :
                      mbsIoRule->capStrikes[iReset]*100.),
                     (mbsIoRule->flrStrikes IS NULL ? 99.99 :
                      mbsIoRule->flrStrikes[iReset]*100.)
                     );
                if(mbsIoRule->isKnownCpn[iReset])
                {
                    GtoErrMsg("%s %5.2lf%% %5.2lf%%\n",
                              outmesg,
                              mbsIoRule->knownCpns[iReset]*100.,
                              mbsIoRule->knownCpnIndexRates[iReset]*100.
                              );
                }
                else
                {
                    GtoErrMsg("%s\n",outmesg);
                }
            } /* iReset */
            if(nRows < mbsIoRule->numResets)
            {
                GtoErrMsg(" ...\n");
            }
        } /* if any reset rows at all */

    } /* is floating */
              
    F_END;
}



/***************************************************************************
 * MakeMbsPoRule();
 * Constructor for TMbsPoRule
 * Sets internals of structure to input values
 * If input ptr is NULL, also allocates structure
 * Returns SUCCESS/FAILURE
 ***************************************************************************/
int EXPORT
MakeMbsPoRule
   (double        multiplier,        /* (I) all PO cashflows scaled by this */
    TMbsPoRule  **mbsPoRule)         /* (I/O) structure to alloc/modify */
{
    F_INIT("MakeMbsPoRule");

    /* Check inputs
     */
    /* @# */

    /* Prepare an empty structure (may have to de-alloc existing one) */
    if (*mbsPoRule ISNT NULL)
    {
        FreeMbsPoRule(mbsPoRule);
    }
    XONTRUE(((*mbsPoRule = NEW(TMbsPoRule)) IS NULL),
            "Failed to allocate mbsPoRule structure");

    /* Set member data 
     */
    (*mbsPoRule)->multiplier = multiplier;

    status = SUCCESS;
done:
    if (status ISNT SUCCESS)
    {
        FreeMbsPoRule(mbsPoRule);
    }
    return (status);
}


/***************************************************************************
 * CopyMbsPoRule();
 * Sets internals of structure to input values
 * If input ptr is NULL, also allocates structure
 * Returns SUCCESS/FAILURE
 ***************************************************************************/
int EXPORT
CopyMbsPoRule
   (TMbsPoRule  *inPoRule,   /* (I) PO cashflow rule struc */
    TMbsPoRule **outPoRule)  /* (O) copy of PO cashflow rule struc */
{
    int status = FAILURE;
    static char routine[] = "CopyMbsPoRule";
 
    XONTRUE(outPoRule IS NULL,
            "Must supply ptr to ptr for copy of PO rule" );

    /* trivial case: if input rule is null, copy this and quit */
    if(inPoRule IS NULL)
    {
        *outPoRule = NULL;
        status = SUCCESS;
        goto done;
    }


    /* Call normal po rule constructor */
    XONFAIL(MakeMbsPoRule
            (inPoRule->multiplier,
             outPoRule));

    status = SUCCESS;
done:
    if (status ISNT SUCCESS)
    {
        FreeMbsPoRule(outPoRule);
    }
    return (status);
}



/***************************************************************************
 * FreeMbsPoRule()
 * Destructor for TMbsPoRule
 * (also resets pointer to NULL)
 ***************************************************************************/
void EXPORT
FreeMbsPoRule
   (TMbsPoRule **mbsPoRule)         /* (I/O) structure to free */
{
    static char routine[] = "FreeMbsPoRule";
    
    if (*mbsPoRule ISNT NULL)
    {
        /* free dynamics vars */
        /* (none) */
        /* free entire structure */
        FREE(*mbsPoRule);
        *mbsPoRule = NULL;
    }
}


/***************************************************************************
 * CheckMbsPoRule()
 * Checks contents of po rule struc for validity
 * Returns SUCCESS/FAILURE
 ***************************************************************************/
int EXPORT
CheckMbsPoRule
   (TMbsPoRule *mbsPoRule)     /* (I) structure to check */
{
    F_INIT("CheckMbsPoRule");

    /* @# for now, just check pointer */
    XONTRUE( mbsPoRule IS NULL, "No PO rule supplied" );

    F_END;
}


/***************************************************************************
 * PrintMbsPoRule()
 * Pint out contents of deal struc for debugging
 * Returns SUCCESS/FAILURE
 ***************************************************************************/
int EXPORT
PrintMbsPoRule
   (char        *label,         /* (I) label string telling us what
                                 * this PO rule applies to */
    TMbsPoRule  *mbsPoRule)     /* (I) structure contains PO rule info. */
{
    F_INIT("PrintMbsPoRule");

    XONTRUE(mbsPoRule IS NULL, 
            "Supplied null ptr" );


    GtoErrMsg("------------------------------------------------------------------\n");
    GtoErrMsg("PO Rule: %s\n",label);
    GtoErrMsg("------------------------------------------------------------------\n");
    GtoErrMsg("PO mult=%lf\n",mbsPoRule->multiplier);

    F_END;
}


/***************************************************************************
 * MakeMbsPpAssump();
 * Constructor for TMbsPrepayAssump
 * Sets internals of structure to input values
 * If input ptr is NULL, also allocates structure
 * Returns SUCCESS/FAILURE
 ***************************************************************************/
int EXPORT
MakeMbsPpAssump
   (/* BASIC PREPAY ASSUMPTIONS */
    long      prepayModel,       /* (I) prepay model to use */
    TBoolean inclSchedAmort,    /* (I) TRUE=return prepays & sched am,
                                 * FALSE=return just prepays */
    TBoolean inclBullet,        /* (I) TRUE=force total amort on bulletDate;
                                 * does so for ANY prepay model */
    TBoolean useNewWacInfo,     /* (I) TRUE=use new wac info.; default is
                                 * FALSE */
    TDate    bulletDate,        /* (I) date on which to do 100% amort */
    TDate    startDate,         /* (I) first month for which amort needed
                                 * (day-of-month ignored) */
    double   grossCpn,          /* (I) gross coupon rate (in decimal); 
                                 * should be the gross coupon accruing 
                                 * in startDate */ 
    long     wala,              /* (I) WALA (in months) as of start_date */
    long     warm,              /* (I) wgted-avg time to maturity (in months)
                                 * as of start_date */
    long      amortForm,         /* (I) form for output amort
                                 * (e.g., MBS_PP_SPD_PSA, etc.) */
    long      numAmortMons,      /* (I) # months for which to compute amort */
    /* SPECIAL INFO NEEDED FOR SPECIFIC PREPAY MODELS */
    double   grossMargin,       /* (I) gross margin for ARM */
    /* Vector prepay model */
    TCurve  *prepayVector,      /* (I) (for MBS_PP_MODEL_VECTOR only) monthly
                                 * amortization rates (in SMM, decimal), where
                                 * dates should be 1st of reported prepay
                                 * month */
    /* "Std" tweaks for prepay models (valid w/any prepay model);
     * (default) no tweak obtained with 1.0 */
    double   stdSmmSpdTwk,      /* (I) mult. tweak to prepay (SMM) speed */
    double   stdSmmFlrTwk,      /* (I) mult. tweak to floor speed (in SMM)*/
    double   stdSteepTwk,       /* (I) mult. tweak to rationality */
    double   stdIncentPtTwk,    /* (I) mult. tweak to incentive measure */
    /* "FQuinn" tweaks for Mgrp and ARM */
    double   turnoverTwk,       /* (I) additive tweak (in CPR) to turnover rate
                                 * of fully seasoned paper (positive for higher
                                 * turnover prepay rate); use 0.0 for no
                                 * tweak */
    double   incentPtTwk,       /* (I) shift in prepay incentive point--the
                                 * prepay incentive ratio at which refi-based
                                 * prepays begin to kick in is shifted deeper
                                 * in-the-money (lower rates) by this amount;
                                 * use 0.0 for no tweak */
    double   totSteepTwk,       /* (I) multiplicative increase in the
                                 * rationality of the callable region of the
                                 * prepay function is increased by this amount;
                                 * note that this tweak applies both to
                                 * burned-out (base) and hazard components of
                                 * prepay function; use 1.0 for no tweak */
    double   baseSteepTwk,      /* (I) multiplicative increase in the
                                 * rationality of the callable region of the
                                 * prepay function is increased by this amount;
                                 * note that this tweak ONLY applies to the
                                 * burned-out portion of the prepay function;
                                 * use 1.0 for no tweak */
    double   ageRampTwk,        /* (I) Additive increae in month to age */
    TMbsPrepayAssump
          **mbsPrepayAssump)    /* (I/O) structure to alloc/modify */
{
    int status = FAILURE;
    static char routine[] = "MakeMbsPpAssump";

    /* Check inputs
     */
    if (!is_valid_prepay_model(prepayModel))
    {   
        GtoErrMsg("%s : No such prepay model [%d]\n",routine,prepayModel);
        goto done;
    }
    if (!is_valid_prepay_form(amortForm))
    {
        GtoErrMsg("%s : Invalid form [%d] requested for output prepays\n",
            routine,amortForm);
        goto done;
    }
    if (numAmortMons <= 0)
    {
        GtoErrMsg("%s : Must have at least one amortization period [%d]\n",
                  routine,numAmortMons);
        goto done;
    }
    if (numAmortMons >= DEF_MBS_PP_ARM_MAX_PREPAYS)
    {
        GtoErrMsg("%s: Cannot generate more than 400 months of prepays\n",
            routine);
        goto done;
    }
    if (wala < 0 ||
        wala > 480)
    {
        GtoErrMsg("%s: Invalid WALA, out of 0:480 range [%d]\n",routine,wala);
        goto done;
    }
    if (warm + wala > 360)
    {
        GtoErrMsg("%s : WALA(%d) + WARM(%d) larger than 360\n",
                  routine,wala,warm);
        goto done;
    }
    if (prepayModel IS MBS_PP_MODEL_ARM)
    {
        if (grossMargin < 0. ||
            grossMargin > 1.)
        {
            GtoErrMsg("%s : gross margin [%f] must be between 0 and 1\n",
                routine,grossMargin);
            goto done;
        }
    }
    if (grossCpn < 0. ||
        grossCpn > 1.)
    {
        GtoErrMsg("%s : grossCpn [%f] must be between 0 and 1\n",
            routine,grossCpn);
        goto done;
    }


    /* If needed, allocate the structure 
     */
    if (*mbsPrepayAssump ISNT NULL)
    {
        FreeMbsPrepayAssump(mbsPrepayAssump);
    }
    if ((*mbsPrepayAssump = NEW(TMbsPrepayAssump)) IS NULL)
    {
        GtoErrMsg("%s : Failed to allocate mbsPrepayAssump structure\n",
            routine);
        goto done;
    }

    /* Now, in either case, (re)allocate internal dynamic vars */
    /* Actually, here we construct the internal TCurve */
    if (prepayVector ISNT NULL)
    {
        (*mbsPrepayAssump)->prepayVector = GtoCopyCurve(prepayVector);

        if ((*mbsPrepayAssump)->prepayVector IS NULL)
        {
            GtoErrMsg("%s : Failed to allocate internal dynamic vars\n",
                routine);
            goto done;
        }
    }


    /* Set member data 
     */
    (*mbsPrepayAssump)->prepayModel = prepayModel;
    (*mbsPrepayAssump)->inclSchedAmort = inclSchedAmort;
    (*mbsPrepayAssump)->inclBullet = inclBullet;
    (*mbsPrepayAssump)->useNewWacInfo = useNewWacInfo;
    (*mbsPrepayAssump)->bulletDate = bulletDate;
    (*mbsPrepayAssump)->startDate = startDate;
    (*mbsPrepayAssump)->grossCpn = grossCpn;
    (*mbsPrepayAssump)->wala = wala;
    (*mbsPrepayAssump)->warm = warm;
    (*mbsPrepayAssump)->amortForm = amortForm;
    (*mbsPrepayAssump)->numAmortMons = numAmortMons;
    (*mbsPrepayAssump)->grossMargin = grossMargin;
    (*mbsPrepayAssump)->stdSmmSpdTwk = stdSmmSpdTwk;
    (*mbsPrepayAssump)->stdSmmFlrTwk = stdSmmFlrTwk;
    (*mbsPrepayAssump)->stdSteepTwk = stdSteepTwk;
    (*mbsPrepayAssump)->stdIncentPtTwk = stdIncentPtTwk;
    (*mbsPrepayAssump)->turnoverTwk = turnoverTwk;
    (*mbsPrepayAssump)->incentPtTwk = incentPtTwk;
    (*mbsPrepayAssump)->totSteepTwk = totSteepTwk;
    (*mbsPrepayAssump)->baseSteepTwk = baseSteepTwk;
    (*mbsPrepayAssump)->ageRampTwk = ageRampTwk;

    status = SUCCESS;

done:
    if (status ISNT SUCCESS)
    {
        FreeMbsPrepayAssump(mbsPrepayAssump);
    }
    return (status);
}


/***************************************************************************
 * CopyMbsPpAssump();
 * Sets internals of structure to input values
 * If input ptr is NULL, also allocates structure
 * Returns SUCCESS/FAILURE
 ***************************************************************************/
int EXPORT
CopyMbsPpAssump
   (TMbsPrepayAssump  *inPpAssump,   /* (I) input prepay assump struc */
    TMbsPrepayAssump **outPpAssump)  /* (O) copy of prepay assump struc */
{
    F_INIT("CopyMbsPpAssump");


    XONTRUE(outPpAssump IS NULL,
            "Must supply ptr to ptr for copy of pp assump" );
    /* trivial case: if input struc is null, copy this and quit */
    if(inPpAssump IS NULL)
    {
        *outPpAssump = NULL;
        status = SUCCESS;
        goto done;
    }

    /* Use the normal constructor */
    XONFAIL(MakeMbsPpAssump
            (inPpAssump->prepayModel,
             inPpAssump->inclSchedAmort,
             inPpAssump->inclBullet,
             inPpAssump->useNewWacInfo,
             inPpAssump->bulletDate,
             inPpAssump->startDate,
             inPpAssump->grossCpn,
             inPpAssump->wala,
             inPpAssump->warm,
             inPpAssump->amortForm,
             inPpAssump->numAmortMons,
             inPpAssump->grossMargin,
             inPpAssump->prepayVector,
             inPpAssump->stdSmmSpdTwk,
             inPpAssump->stdSmmFlrTwk,
             inPpAssump->stdSteepTwk,
             inPpAssump->stdIncentPtTwk,
             inPpAssump->turnoverTwk,
             inPpAssump->incentPtTwk,
             inPpAssump->totSteepTwk,
             inPpAssump->baseSteepTwk,
             inPpAssump->ageRampTwk,
             outPpAssump));

    status = SUCCESS;
done:
    if (status ISNT SUCCESS)
    {
        FreeMbsPrepayAssump(outPpAssump);
    }
    return (status);
}





/***************************************************************************
 * FreeMbsPrepayAssump()
 * Destructor for TMbsPoRule
 * (also resets pointer to NULL)
 ***************************************************************************/
void EXPORT
FreeMbsPrepayAssump
   (TMbsPrepayAssump **mbsPrepayAssump)         /* (I/O) structure to free */
{
    static char routine[] = "FreeMbsPrepayAssump";
    
    if (*mbsPrepayAssump ISNT NULL)
    {
        /* free internal dynamic vars */
        GtoFreeTCurve( (*mbsPrepayAssump)->prepayVector );
        /* free entire structure */
        FREE(*mbsPrepayAssump);
        *mbsPrepayAssump = NULL;
    }
}


/***************************************************************************
 * CheckMbsPrepayAssump()
 * Checks contents of pp assump struc for validity
 * Returns SUCCESS/FAILURE
 ***************************************************************************/
int EXPORT
CheckMbsPrepayAssump
   (TMbsPrepayAssump *mbsPrepayAssump) /* (I) structure to check */
{
    F_INIT("CheckMbsPrepayAssump");

    /* @# for now, just check pointer */
    XONTRUE( mbsPrepayAssump IS NULL, "No prepay assump supplied" );

    F_END;
}


/***************************************************************************
 * PrintMbsPrepayAssump()
 * Pint out contents of deal struc for debugging
 * Returns SUCCESS/FAILURE
 ***************************************************************************/
int EXPORT
PrintMbsPpAssump
   (TBoolean     printPpModelInfo, /* (I) TRUE to also output info specific to
                                    * prepay-model code */
    TMbsPrepayAssump
          *mbsPrepayAssump)    /* (I) structure contains prepay assump. */
{
    F_INIT("PrintMbsPpAssump");

    XONTRUE(mbsPrepayAssump IS NULL,
            "Failed to supply prepay assump struc" );

    GtoErrMsg
     ("------------------------------------------------------------------\n");
    GtoErrMsg("Prepay Assumptions\n");
    GtoErrMsg
     ("------------------------------------------------------------------\n");
    GtoErrMsg("Prepay Model = %d\n",mbsPrepayAssump->prepayModel);
    GtoErrMsg("Incl sched amort. = %d\n",mbsPrepayAssump->inclSchedAmort);
    GtoErrMsg("Incl bullet = %d\n",mbsPrepayAssump->inclBullet);
    GtoErrMsg("Use New Wac Info = %d\n",mbsPrepayAssump->useNewWacInfo);

    if(printPpModelInfo)
    {
        GtoErrMsg("Bullet date = %d (%s)\n",
                  mbsPrepayAssump->bulletDate,
                  GtoFormatDate(mbsPrepayAssump->bulletDate));
        GtoErrMsg("Start date = %d (%s)\n",
                  mbsPrepayAssump->startDate,
                  GtoFormatDate(mbsPrepayAssump->startDate));
        GtoErrMsg("Gross coupon = %f\n",mbsPrepayAssump->grossCpn);
        GtoErrMsg("WALA = %d\n",mbsPrepayAssump->wala);
        GtoErrMsg("WARM = %d\n",mbsPrepayAssump->warm);
        GtoErrMsg("Amort form = %d\n",mbsPrepayAssump->amortForm);
        GtoErrMsg("# of Amort mths = %d\n",mbsPrepayAssump->numAmortMons);
        GtoErrMsg("Gross margin = %f\n",mbsPrepayAssump->grossMargin);
        GtoErrMsg("Turnover tweak = %f\n",mbsPrepayAssump->turnoverTwk);
        GtoErrMsg("IncentPt tweak = %f\n",mbsPrepayAssump->incentPtTwk);
        GtoErrMsg("TotSteep tweak = %f\n",mbsPrepayAssump->totSteepTwk);
        GtoErrMsg("BaseSteep tweak = %f\n",mbsPrepayAssump->baseSteepTwk);
        GtoErrMsg("AgeRamp tweak = %f\n",mbsPrepayAssump->ageRampTwk);
    }
        
    GtoErrMsg("Std SMM speed twk = %lf\n",mbsPrepayAssump->stdSmmSpdTwk);
    GtoErrMsg("Std SMM floor twk = %lf\n",mbsPrepayAssump->stdSmmFlrTwk);
    GtoErrMsg("Std steep twk = %lf\n",mbsPrepayAssump->stdSteepTwk);
    GtoErrMsg("Std incent pt twk = %lf\n",mbsPrepayAssump->stdIncentPtTwk);
    
    if(mbsPrepayAssump->prepayVector ISNT NULL &&
       mbsPrepayAssump->prepayModel IS MBS_PP_MODEL_VECTOR)
    {
        GtoPrintTCurve(mbsPrepayAssump->prepayVector, "Prepay Vector");
    }

    F_END;
}


/***************************************************************************
 * RefiRateNeeded()
 * From MBS and prepay info, determines type of refi rate needed
 ***************************************************************************/
int EXPORT RefiRateNeeded
   (TMbsDeal         *mbsDeal,           /* (I) deal info */
    TMbsPrepayAssump *mbsPrepayAssump,   /* (I) prepay info */
    long              *refiRateType)      /* (O) type of refi rate needed for
                                          * this MBS and this prepay model
                                          * (e.g., MBS_PP_REFI_TYPE_FH30) */
{
    F_INIT("RefiRateNeeded");

    /* Depending on prepay model & MBS type, set type
     * of refi rate needed--in a sense, these settings
     * depend very much on the nature of our prepay
     * models. For example, MGrp and ARM models currently
     * want to have FH30 rates */
    switch (mbsPrepayAssump->prepayModel)
    {
      case MBS_PP_MODEL_NONE:
        /* no refi rate needed */
        *refiRateType = MBS_PP_REFI_TYPE_NONE;
        break;
      case MBS_PP_MODEL_CONST:
        /* no refi rate needed */
        *refiRateType = MBS_PP_REFI_TYPE_NONE;
        break;
      case MBS_PP_MODEL_VECTOR:
        /* no refi rate needed */
        *refiRateType = MBS_PP_REFI_TYPE_NONE;
        break;
      case MBS_PP_MODEL_FIX_MGRP:
      case MBS_PP_MODEL_ARM:
        switch (mbsDeal->mbsTerm)
        {
          case MBS_PP_MBSTERM_30:
            *refiRateType = MBS_PP_REFI_TYPE_FH30;
            break;
          case MBS_PP_MBSTERM_15:
            *refiRateType = MBS_PP_REFI_TYPE_FH15;
            break;
          default:
            GtoErrMsg("%s : Invalid mbs type\n",routine);
            goto done;
        }
        break;
      default:
        XONTRUE( TRUE, "Unknown prepay model type" );
    }

    F_END;
}


/***************************************************************************
 * MakeMbsRateEnv();
 * Constructor for TMbsRateEnv
 * Sets internals of structure to input values
 * If input ptr is NULL, also allocates structure
 * Returns SUCCESS/FAILURE
 ***************************************************************************/
int EXPORT
MakeMbsRateEnv
   (TDate    tradeDate,          /* (I) date on which market information */
    TCurve  *discZeroCurve,      /* (I) discounting zero curve */
    long      numIndexZeroCurves, /* (I) # of index zero curves */
    TCurve **indexZeroCurves,    /* (I) array of index zero curves */
    long    *numSettleDays,      /* (I) for each zero curve, settle day
                                  * convention (0th value for disc curve,
                                  * then values for index curves) */
    TCurve  *volCurve,           /* (I) Base volatility curve; note */
    TDate   *armResetDates,      /* (I) array of dates for arm reset curve */
    double  *armResetRates,      /* (I) array of rates for dates of curve */
    long      numArmResetPts,     /* (I) # elements in dates & rates array */
    double   armResetFreq,       /* (I) basis. see GtoRateToDiscount */
    long     dayCountConv,       /* (I) see GtoDayCountFraction */
    /* MARKET INFORMATION NEEDED FOR PREPAYS - USER SUPPLIED REFI RATES */
    long     amortIndexType,      /* (I) type of index rate used to drive
                                  * prepays; e.g., MBS_PP_REFI_TYPE_FH30 */
    TDate   amortIndexStartDate, /* (I) month of first amort index rate */
    long     numAmortIndexRates,  /* (I) number of amortIndexRates[] */
    double *amortIndexRates,     /* (I) past & future amort index rates */
    double *fhCommitPts,         /* (I) (only if amort index is FH15/30
                                  * commitment rate) past & future
                                  * monthly-average */
    double  armVsFixedCcSprd,    /* (I) spread between ARM gross current coupon
                                  * and the FHLMC commitment rate (in
                                  * fh_commit[]) as of the prepay startDate;
                                  * spread should be decimal;  e.g., -0.0108);
                                  * this spread is added to fh_commit[] rates
                                  * to obtain future values of the ARM gross
                                  * current coupon rate */
    /* PTR TO NEW/EXISTING STRUCTURE */
    TMbsRateEnv **mbsRateEnv)    /* (I/O) structure to alloc/modify */
{
    F_INIT("MakeMbsRateEnv");
    long  i;
    long  idxCrv;
    long nMons;

    /* Check inputs
     */
    if (armVsFixedCcSprd < -1. ||
        armVsFixedCcSprd > 1.)
    {
        GtoErrMsg("%s : ARM vs. fixed spread [%f] must be between -1 and 1\n",
            routine,armVsFixedCcSprd);
        goto done;
    }

    /* If needed, allocate the structure
     */
    if (*mbsRateEnv ISNT NULL)
    {
        FreeMbsRateEnv(mbsRateEnv);
    }
    if ((*mbsRateEnv = NEW(TMbsRateEnv)) IS NULL)
    {
        GtoErrMsg("%s : Failed to allocate mbsRateEnv structure\n",routine);
        goto done;
    }

    /* copy disc zero curve */
    if (discZeroCurve IS NULL)
    {
        (*mbsRateEnv)->discZeroCurve = NULL;
    }
    else
    {
        XONTRUE( ((*mbsRateEnv)->discZeroCurve = 
                  GtoCopyCurve(discZeroCurve)) IS NULL,
                "Failed to allocate indexZeroCurve" );
    }

    /* copy index zero curves */
    XONTRUE(numIndexZeroCurves < 0,
            "Cannot have negative numIndexZeroCurves" );
    (*mbsRateEnv)->numIndexZeroCurves = numIndexZeroCurves;
    if( numIndexZeroCurves > 0 )
    {
        XONTRUE(((*mbsRateEnv)->indexZeroCurves = 
                 NEW_ARRAY(TCurve*,
                           numIndexZeroCurves)) IS NULL,
                "Failed to alloc array of index zc ptrs" );
        for(idxCrv=0; idxCrv<numIndexZeroCurves; idxCrv++)
        {
            if (indexZeroCurves[idxCrv] IS NULL)
            {
                (*mbsRateEnv)->indexZeroCurves[idxCrv] = NULL;
            }
            else
            {
                XONTRUE( ((*mbsRateEnv)->indexZeroCurves[idxCrv] = 
                          GtoCopyCurve(indexZeroCurves[idxCrv])) IS NULL,
                        "Failed to allocate indexZeroCurve" );
            }
        }
    }
    else
    {
        (*mbsRateEnv)->indexZeroCurves = NULL;
    }

    /* copy array of # settle days */
    if (numSettleDays ISNT NULL)
    {
        XONTRUE( ((*mbsRateEnv)->numSettleDays = 
                  NEW_ARRAY(long,1+numIndexZeroCurves)) IS NULL,
                "Failed to alloc settle days array" );
        for (idxCrv = 0; idxCrv <= numIndexZeroCurves; idxCrv++)
        {
            (*mbsRateEnv)->numSettleDays[idxCrv] = numSettleDays[idxCrv];
        }
    }

    /* copy vol curve */
    if (volCurve IS NULL)
    {
        (*mbsRateEnv)->volCurve = volCurve;
    }
    else
    {
        XONTRUE( ((*mbsRateEnv)->volCurve = GtoCopyCurve(volCurve)) IS NULL,
                "Failed to allocate volCurve");
    }


    (*mbsRateEnv)->armResetIndex = NULL;
    if (numArmResetPts > 0)
    {
        (*mbsRateEnv)->armResetIndex = GtoMakeTCurve(tradeDate,
                                                     armResetDates,
                                                     armResetRates,
                                                     numArmResetPts,
                                                     armResetFreq,
                                                     dayCountConv);
    
        if ((*mbsRateEnv)->armResetIndex IS NULL)
        {
            GtoErrMsg("%s : Failed to allocate armResetIndex\n", routine);
            goto done;
        }

        /* In this model, resets should be yearly */
        for (i = 1; i < ((*mbsRateEnv)->armResetIndex)->fNumItems; i++)
        {
            if (MonsDiff(((*mbsRateEnv)->armResetIndex)->fArray[i-1].fDate,
                         ((*mbsRateEnv)->armResetIndex)->fArray[i].fDate,
                        &nMons) ISNT SUCCESS)
            {
                GtoErrMsg("%s : MonsDiff failed\n",routine);
                goto done;
            }
            if (nMons ISNT 12)
            {
                GtoErrMsg("%s: reset index is not yearly in (%s) and (%s)\n",
                    routine,
                    GtoFormatDate(((*mbsRateEnv)->armResetIndex)->fArray[i-1].fDate),
                    GtoFormatDate(((*mbsRateEnv)->armResetIndex)->fArray[i].fDate));
                goto done;
            }
        }
    }
    /* @# emergency bug fix 9/3/96:
     * if reset rate TCurve is empty,
     * supply it with at least one rate on maturity date, to
     * keep code from crashing
     */
    else if (numArmResetPts IS 0)
    {
        numArmResetPts = 1;
        armResetDates = NEW_ARRAY(TDate, numArmResetPts);
        armResetRates = NEW_ARRAY(double, numArmResetPts);
        armResetRates[0] = 0.10;  /*  dummy rate--doesn't matter */
        /* date DOES matter: make it 10 years after trade date */
        NxtMth(tradeDate,120,&(armResetDates[0]));

        (*mbsRateEnv)->armResetIndex = GtoMakeTCurve(tradeDate,
                                                     armResetDates,
                                                     armResetRates,
                                                     numArmResetPts,
                                                     armResetFreq,
                                                     dayCountConv);
        FREE_ARRAY(armResetDates);
        FREE_ARRAY(armResetRates);
    }

    /* Set member data
     */
    (*mbsRateEnv)->tradeDate = tradeDate;
    (*mbsRateEnv)->amortIndexType = amortIndexType;
    (*mbsRateEnv)->amortIndexStartDate = amortIndexStartDate;
    (*mbsRateEnv)->numAmortIndexRates = numAmortIndexRates;
    (*mbsRateEnv)->armVsFixedCcSprd = armVsFixedCcSprd;
    if (numAmortIndexRates > 0)
    {
        (*mbsRateEnv)->amortIndexRates = NEW_ARRAY(double,numAmortIndexRates);
        (*mbsRateEnv)->fhCommitPts = NEW_ARRAY(double,numAmortIndexRates);
        for (i = 0; i < numAmortIndexRates; i++)
        {
            (*mbsRateEnv)->amortIndexRates[i] = amortIndexRates[i];
            (*mbsRateEnv)->fhCommitPts[i] = fhCommitPts[i];
        }
    }

    status = SUCCESS;

done:
    if (status ISNT SUCCESS)
    {
        FreeMbsRateEnv(mbsRateEnv);
    }
    return (status);
}


/***************************************************************************
 * FreeMbsRateEnv()
 * Destructor for TMbsRateEnv
 * (also resets pointer to NULL)
 ***************************************************************************/
void EXPORT
FreeMbsRateEnv
   (TMbsRateEnv **mbsRateEnv)         /* (I/O) structure to free */
{
    static char routine[] = "FreeMbsRateEnv";
    long   idxCrv;

    if (*mbsRateEnv ISNT NULL)
    {
        /* free internal dynamic vars */
        GtoFreeTCurve( (*mbsRateEnv)->discZeroCurve );
        for(idxCrv=0; idxCrv<(*mbsRateEnv)->numIndexZeroCurves; idxCrv++)
        {
            GtoFreeTCurve( (*mbsRateEnv)->indexZeroCurves[idxCrv] );
        }
        FREE_ARRAY((*mbsRateEnv)->indexZeroCurves);
        FREE_ARRAY((*mbsRateEnv)->numSettleDays);
        GtoFreeTCurve( (*mbsRateEnv)->volCurve );
        GtoFreeTCurve( (*mbsRateEnv)->armResetIndex );
        FREE_ARRAY((*mbsRateEnv)->amortIndexRates);
        FREE_ARRAY((*mbsRateEnv)->fhCommitPts);
        FREE_ARRAY((*mbsRateEnv)->histAmortIndexDates);
#ifndef HARDCODE_V
        FREE_ARRAY((*mbsRateEnv)->histAmortIndexRates);
        FREE_ARRAY((*mbsRateEnv)->histFhCommitPts);
#endif
        /* free entire structure */
        FREE(*mbsRateEnv);
        *mbsRateEnv = NULL;
    }
}


/***************************************************************************
 * CheckMbsRateEnv()
 * Checks contents of rate env. struc for validity
 * Returns SUCCESS/FAILURE
 ***************************************************************************/
int EXPORT
CheckMbsRateEnv
   (TMbsRateEnv *mbsRateEnv) /* (I) structure to check */
{
    F_INIT("CheckMbsRateEnv");

    /* @# for now, just check pointer */
    XONTRUE( mbsRateEnv IS NULL, "No rate environment supplied" );

    F_END;
}

/***************************************************************************
 * PrintMbsRateEnv()
 * Pint out contents of deal struc for debugging
 * Returns SUCCESS/FAILURE
 ***************************************************************************/
int EXPORT
PrintMbsRateEnv
   (TBoolean     printPrepayInfo, /* (I) TRUE to also output info specific to
                                   * prepay-model code */
    TMbsRateEnv *mbsRateEnv)    /* (I) structure contains MBS rate Env. info. */
{
    F_INIT("PrintMbsRateEnv");
    long  i;
    long  idxCrv;
    char curveLbl[48];

    XONTRUE(mbsRateEnv IS NULL,
            "Failed to provide struc" );

    GtoErrMsg
     ("------------------------------------------------------------------\n");
    GtoErrMsg("Rate Environment\n");
    GtoErrMsg
     ("------------------------------------------------------------------\n");
    GtoErrMsg("Trade date = %s\n",
        GtoFormatDate(mbsRateEnv->tradeDate));
    if(printPrepayInfo)
    {
        GtoErrMsg("Amort index type = %d\n",mbsRateEnv->amortIndexType);
        GtoErrMsg("Amort index start date = %d (%s)\n",
                  mbsRateEnv->amortIndexStartDate,
                  GtoFormatDate(mbsRateEnv->amortIndexStartDate));
        GtoErrMsg("Arm vs. fixed CcSprd = %f\n",mbsRateEnv->armVsFixedCcSprd);
        GtoErrMsg("# of Amort index rates = %d\n",mbsRateEnv->numAmortIndexRates);
        for (i = 0; i < mbsRateEnv->numAmortIndexRates; i++)
        {
            GtoErrMsg("I = %d Amort index rate = %f Pts = %f\n",
                      i,mbsRateEnv->amortIndexRates[i],mbsRateEnv->fhCommitPts[i]);
        }
        if (mbsRateEnv->armResetIndex ISNT NULL)
        {
            GtoPrintTCurve(mbsRateEnv->armResetIndex, "Arm reset index");
        }
    }
    GtoErrMsg("# index zero curves: %ld\n",
       mbsRateEnv->numIndexZeroCurves);
    GtoErrMsg("# Settle days: %ld (disc)",
       mbsRateEnv->numSettleDays[0]);
    for(idxCrv=1; idxCrv<=mbsRateEnv->numIndexZeroCurves; idxCrv++)
    {
        GtoErrMsg(" %ld(indx %ld)",
                  mbsRateEnv->numSettleDays[idxCrv],
                  idxCrv);
    }
    GtoErrMsg("\n");
    if(mbsRateEnv->discZeroCurve ISNT NULL)
    {
        GtoPrintTCurve(mbsRateEnv->discZeroCurve,
                       "Discount Zero Curve");
    }
    for(idxCrv=1; idxCrv<=mbsRateEnv->numIndexZeroCurves; idxCrv++)
    {
        sprintf(curveLbl,"Index Zero Curve %ld",idxCrv);
        if(mbsRateEnv->indexZeroCurves ISNT NULL)
        {
            GtoPrintTCurve(mbsRateEnv->indexZeroCurves[idxCrv-1], 
                           curveLbl);
        }
    }
    if (mbsRateEnv->volCurve ISNT NULL)
    {
        GtoPrintTCurve(mbsRateEnv->volCurve, "Vol. Curve");
    }

    F_END;
}


/***************************************************************************
 * MakeMbsPricingAssump();
 * Constructor for TMbsPricingAssump
 * Sets internals of structure to input values
 * If input ptr is NULL, also allocates structure
 * Returns SUCCESS/FAILURE
 ***************************************************************************/
int EXPORT
MakeMbsPricingAssump
   (TDate           valueDate,      /* (I) date for which to compute PV */
    TBoolean        detFirstCfUsingMbsAcc, 
                                    /* (I) affects how first cashflow to price
                                     * is determined from valueDate:
                                     * TRUE: use mbs accrual periods
                                     * FALSE: use deal/swap accrual periods */
    TDate           pvBalanceDate,  /* (I) price to be computed rel. to MBS 
                                     * balance on this date */
    double          oas,            /* (I) option-adjusted spread */
    TBoolean        computeOasFromPv, /* FALSE ("normal") for PV from OAS;
                                       * TRUE to compute OAS from PV */
    double          pvForOas,       /* PV to use if computing OAS from PV */
    /* MANNER OF PRICING */
    TBoolean        priceBareOptions, /* (I) FALSE=compute deal PV 
                                       * w/embed. options; TRUE=compute PV 
                                       * of some/all options only*/
    char           *optionsToPrice, /* (I) (when pricing bare options) 
                                     * 2-char string:
                                     * 1st char:
                                     *        GTO_MODEL_CAP      'C'
                                     *        GTO_MODEL_FLOOR    'F'
                                     *        GTO_MODEL_COLLAR   'L'
                                     * 2nd char: dep on context
                                     * (e.g., is 'N','S','B' in ARM model)*/
    /* TYPE OF MODEL */
    long            irModel,        /* (I) int rate model to use 
                                     * (e.g., MBS_IRMODEL_TREE1FAC) */
    /* TREE MODELS */
    long            minPpy,         /* (I) pds per year (early) */
    long            minPpy2,        /* (I) later pds per year */
    TDate           ppySwitchDate,  /* (I) early/later transition date */
    long            maxNumSvDims,   /* (I) # values in numStateVarGridPts */
    long           *numStateVarGridPts, /* (I) # sample pts for each sv dim */
    /* MISC */
    char           *holidayFile,    /* (I) for holiday/date calc */
    char           *dataFileDir,    /* (I) dir for suppl. data files */
    /* PTR TO NEW/EXISTING STRUCTURE */
    TMbsPricingAssump
                  **mbsPricingAssump) /* (I/O) structure to alloc/modify */
{
    F_INIT("MakeMbsPricingAssump");
    long   iSv;

    /* Check inputs
     */

    /* If needed, allocate the structure
     */
    if (*mbsPricingAssump ISNT NULL)
    {
        FreeMbsPricingAssump(mbsPricingAssump);
    }
    if ((*mbsPricingAssump = NEW(TMbsPricingAssump)) IS NULL)
    {
        GtoErrMsg("%s : Failed to allocate mbsPricingAssump structure\n",routine);
        goto done;
    }
 
    /* Set member data
     */
    (*mbsPricingAssump)->valueDate = valueDate;
    (*mbsPricingAssump)->detFirstCfUsingMbsAcc = detFirstCfUsingMbsAcc;
    (*mbsPricingAssump)->pvBalanceDate = pvBalanceDate;
    /* force balance date to 1st of month */
    XONFAIL( SetDOM(1,TRUE,&((*mbsPricingAssump)->pvBalanceDate)) );
    (*mbsPricingAssump)->oas = oas;
    (*mbsPricingAssump)->computeOasFromPv = computeOasFromPv;
    (*mbsPricingAssump)->pvForOas = pvForOas;
    /* manner of pricing */
    (*mbsPricingAssump)->priceBareOptions = priceBareOptions;
    (*mbsPricingAssump)->optionsToPrice[0] = toupper(optionsToPrice[0]);
    (*mbsPricingAssump)->optionsToPrice[1] = toupper(optionsToPrice[1]);
    (*mbsPricingAssump)->optionsToPrice[2] = '\0';
    /* type of int rate model */
    (*mbsPricingAssump)->irModel = irModel;
    /* info for tree models */
    (*mbsPricingAssump)->minPpy = minPpy;
    (*mbsPricingAssump)->minPpy2 = minPpy2;
    (*mbsPricingAssump)->ppySwitchDate = ppySwitchDate;
    /* copy state var grid resolutions */
    XONRANGE(maxNumSvDims,1,20,"maxNumStateVarDims");
    (*mbsPricingAssump)->maxNumSvDims = maxNumSvDims;
    XONTRUE(((*mbsPricingAssump)->numStateVarGridPts = 
             NEW_ARRAY(long,(*mbsPricingAssump)->maxNumSvDims)) IS NULL,
            "Failed to alloc grid pts array" );
    for(iSv=0; iSv<(*mbsPricingAssump)->maxNumSvDims; iSv++)
    {
        XONRANGE(numStateVarGridPts[iSv],2,1000,"numStateVarGridPts");
        (*mbsPricingAssump)->numStateVarGridPts[iSv] = numStateVarGridPts[iSv];
    }
    /* misc */
    strncpy((*mbsPricingAssump)->holidayFile,holidayFile,MAXPATHLEN-1);
    strncpy((*mbsPricingAssump)->dataFileDir,dataFileDir,MAXPATHLEN-1);

    status = SUCCESS;
done:
    if (status ISNT SUCCESS)
    {
        FreeMbsPricingAssump(mbsPricingAssump);
    }
    return (status);
}

 

/***************************************************************************
 * FreeMbsPricingAssump()
 * Destructor for TMbsPricingAssump
 * (also resets pointer to NULL)
 ***************************************************************************/
void EXPORT
FreeMbsPricingAssump
   (TMbsPricingAssump **mbsPricingAssump)   /* (I/O) structure to free */
{
    static char routine[] = "FreeMbsPricingAssump";
 
    if (*mbsPricingAssump ISNT NULL)
    {
        /* free internal dynamic vars */
        FREE((*mbsPricingAssump)->numStateVarGridPts);
        /* free entire structure */
        FREE(*mbsPricingAssump);
        *mbsPricingAssump = NULL;
    }
}


/***************************************************************************
 * CheckMbsPricingAssump()
 * Checks contents of pricing assump. struc for validity
 * Returns SUCCESS/FAILURE
 ***************************************************************************/
int EXPORT
CheckMbsPricingAssump
   (TMbsPricingAssump *mbsPricingAssump) /* (I) structure to check */
{
    F_INIT("CheckMbsPricingAssump");

    /* @# for now, just check pointer */
    XONTRUE( mbsPricingAssump IS NULL, "No pricing assumptions supplied" );

    F_END;
}


/***************************************************************************
 * PrintMbsPricingAssump()
 * Pint out contents of deal struc for debugging
 * Returns SUCCESS/FAILURE
 ***************************************************************************/
int EXPORT
PrintMbsPricingAssump
   (TMbsPricingAssump  *mbsPricingAssump) /* (I) structure contains princing
                                           * assumptions  */
{
    F_INIT("PrintMbsPricingAssump");
    long  iSv;

    XONTRUE(mbsPricingAssump IS NULL,
            "Failed to supply struc" );


    GtoErrMsg
     ("------------------------------------------------------------------\n");
    GtoErrMsg("Pricing Assumptions\n");
    GtoErrMsg
     ("------------------------------------------------------------------\n");
    GtoErrMsg("Value date = %s, PV balance date = %s\n",
        GtoFormatDate(mbsPricingAssump->valueDate),
        GtoFormatDate(mbsPricingAssump->pvBalanceDate));
    GtoErrMsg("OAS = %f\n",mbsPricingAssump->oas);
    if(mbsPricingAssump->priceBareOptions)
    {
        GtoErrMsg("Bare options being priced: (%c, %c)\n",
           mbsPricingAssump->optionsToPrice[0],
           mbsPricingAssump->optionsToPrice[1]);
    }
    GtoErrMsg("IR model=%ld\n",mbsPricingAssump->irModel);
    GtoErrMsg("Min pds/yr = %ld",mbsPricingAssump->minPpy);
    if(mbsPricingAssump->ppySwitchDate ISNT 0)
    {
        GtoErrMsg(" (%ld after %s)",
           mbsPricingAssump->minPpy2,
           GtoFormatDate(mbsPricingAssump->ppySwitchDate));
    }
    GtoErrMsg("\n");
    if(mbsPricingAssump->maxNumSvDims > 0)
    {
        GtoErrMsg("%ld state var dimens, sampled with:",
           mbsPricingAssump->maxNumSvDims);
        for(iSv=0; iSv<mbsPricingAssump->maxNumSvDims; iSv++)
        {
            GtoErrMsg(" %ld",mbsPricingAssump->numStateVarGridPts[iSv]);
        }
        GtoErrMsg("\n");
    }
    GtoErrMsg("Holiday file: %s\n",
       mbsPricingAssump->holidayFile);
    GtoErrMsg("Data files in: %s\n",
       mbsPricingAssump->dataFileDir);
    
    status = SUCCESS;
done:
    return (status);
}

 
/***************************************************************************
 * MakeMbsPricingOutput
 * Constructor for pricing output struc 
 * Note: allocates only; does not fill member data
 ***************************************************************************/
int EXPORT
MakeMbsPricingOutput
   (TMbsPricingOutput    **mbsPricingOutput) /* (I/O) structure to alloc/mod */
{
    F_INIT("MakeMbsPricingOutput");


    /* Toss any existing structure
     */
    if (*mbsPricingOutput ISNT NULL)
    {
        FreeMbsPricingOutput(mbsPricingOutput);
    }
    /* (Re) allocate structure */
    XONTRUE( (*mbsPricingOutput = NEW(TMbsPricingOutput)) IS NULL,
            "Failed to alloc pricing output struc" );
    
    /* Init member data */
    (*mbsPricingOutput)->pv = 0.;
    (*mbsPricingOutput)->oasFromPv = 0.;

    F_END;
}


/***************************************************************************
 * FreeMbsPricingOutput()
 * Destructor for TMbsPricingOutput
 * (also resets pointer to NULL)
 ***************************************************************************/
void EXPORT
FreeMbsPricingOutput
   (TMbsPricingOutput **mbsPricingOutput)     /* (I/O) structure to free */
{
    static char routine[] = "FreeMbsPricingOutput";
 
    if (*mbsPricingOutput ISNT NULL)
    {
        /* free lonernal dynamic vars */
	/* (none for now) */
        /* free entire structure */
        FREE(*mbsPricingOutput);
        *mbsPricingOutput = NULL;
    }
}
