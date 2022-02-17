/*
***************************************************************************
** rate.c
**
** Simple rate functions.
** Based on cfinanci.c from ALIB.
***************************************************************************
*/
#include "irx/rate.h"

#include <ctype.h>
#include <math.h>

#include <irx/error.h>
#include <irx/dateutils.h>
#include <irx/macros.h>

IrxTRateType IRX_CONTINUOUS_RATE = 
{
    IRX_RT_CONTINUOUS,
    0
};

IrxTRateType IRX_ANNUAL_RATE =
{
    IRX_RT_FREQUENCY,
    1
};

IrxTRateType IRX_SIMPLE_RATE = 
{
    IRX_RT_SIMPLE,
    0
};

IrxTRateType IRX_DISCOUNT_RATE =
{
    IRX_RT_DISCOUNT,
    0
};

IrxTRateType IRX_COMPOUND_RATE(int frequency)
{
    IrxTRateType compoundRate;

    compoundRate.type      = IRX_RT_FREQUENCY;
    compoundRate.frequency = frequency;

    return compoundRate;
}

IrxTRateType* irxRateTypeFromString (char* str)
{
    static char routine[] = "irxRateTypeFromString";

    IrxTRateType *rateType = NULL;

    REQUIRE (str != NULL);

    switch (toupper(str[0]))
    {
    case 'S':
        rateType = irxRateTypeMake(IRX_RT_SIMPLE, 0);
        break;
    case 'C':
        rateType = irxRateTypeMake(IRX_RT_CONTINUOUS, 0);
        break;
    case 'D':
        rateType = irxRateTypeMake(IRX_RT_DISCOUNT, 0);
        break;
    case 'F':
    {
        int frequency;
        if (sscanf(&str[1], "%d", &frequency) != 1)
        {
            irxError ("%s: Could not parse frequency following F in (%s)\n",
                      routine, str);
            goto RETURN; /* failure */
        }
        rateType = irxRateTypeMake(IRX_RT_FREQUENCY, frequency);
    }
    default:
        irxError ("%s: Unknown rate type %c\n", routine, toupper(str[0]));
        goto RETURN; /* failure */
    }

 RETURN:

    if (rateType == NULL)
        irxErrorFailure (routine);

    return rateType;
}


IrxTRateType* irxRateTypeFromLong (long freq)
{
    static char routine[] = "irxRateTypeFromLong";
    
    IrxTRateType* rateType = NULL;

    REQUIRE (freq >= 0);

    if (freq == 0)
    {
        rateType = irxRateTypeMake(IRX_RT_SIMPLE, 0);
    }
    else if (freq > 365)
    {
        rateType = irxRateTypeMake(IRX_RT_CONTINUOUS, 0);
    }
    else
    {
        rateType = irxRateTypeMake(IRX_RT_FREQUENCY, (int)freq);
    }

 RETURN:

    if (rateType == NULL)
        irxErrorFailure (routine);

    return rateType;
}


int irxRateTypeValidate (IrxTRateType *rateType)
{
    static char routine[] = "irxRateTypeValidate";
    int         status    = FAILURE;

    REQUIRE (rateType != NULL);

    switch (rateType->type)
    {
    case IRX_RT_FREQUENCY:
        switch (rateType->frequency)
        {
        case 1:
        case 2:
        case 4:
        case 12:
        case 365:
            break;
        default:
            irxError ("%s: Invalid compounding frequency %d - should be "
                      "1,2,4,12 or 365\n", routine, rateType->frequency);
            return FAILURE;
        }
        break;
    case IRX_RT_SIMPLE:
    case IRX_RT_DISCOUNT:
    case IRX_RT_CONTINUOUS:
        REQUIRE(rateType->frequency == 0);
        break;
    default:
        PROGRAM_BUG();
        goto RETURN; /* failure */
    }

    status = SUCCESS;

 RETURN:
    return status;
}


int irxDiscountToRate
(double           discount,         /* (I) Discount factor */
 IrxTDate         startDate,        /* (I) Start date */
 IrxTDate         endDate,          /* (I) End date */
 IrxTDayCountConv rateDayCountConv, /* (I) Day count convention for rate */
 IrxTRateType     rateBasis,        /* (I) Basis for the rate */
 double          *rate)             /* (O) output rate */
{
    static char routine[] = "irxDiscountToRate" ;
    int status = FAILURE;               /* Until proven successful */
    
    double      rateYF ;
    
    REQUIRE (discount > 0.0);
    REQUIRE (startDate != endDate);
    REQUIRE (rate != NULL);

    /* get year fractions for rate */
    if (irxDayCountFraction(startDate, endDate, rateDayCountConv,
                            &rateYF) != SUCCESS)
        goto RETURN;

    if (irxDiscountToRateYearFrac(discount,
                                  rateYF,
                                  rateBasis,
                                  rate) != SUCCESS)
        goto RETURN;

    status = SUCCESS;

 RETURN:
    if (status != SUCCESS)
        irxErrorFailure (routine);

    return (status);
}


int irxRateToDiscount
(double         rate,             /* (I) Rate */
 IrxTDate          startDate,        /* (I) Start date */
 IrxTDate          endDate,          /* (I) End date */
 IrxTDayCountConv  rateDayCountConv, /* (I) Day count convention for rate */
 IrxTRateType      rateBasis,        /* (I) Basis for the rate */
 double        *discount)         /* (O) Discount factor */
{
    int status = FAILURE;               /* Until proven successful */
    static char routine[] = "irxRateToDiscount" ;
    double      rateYF ;

    /* get year fractions for rate */
    if (irxDayCountFraction(startDate, endDate, rateDayCountConv,
                            &rateYF) != SUCCESS)
        goto RETURN;
    
    if (irxRateToDiscountYearFrac(rate, 
                                  rateYF, 
                                  rateBasis,
                                  discount) != SUCCESS)
        goto RETURN;

    status = SUCCESS;
    
 RETURN:
    
    if (status != SUCCESS)
        irxErrorFailure (routine);

    return status;
}

int irxDiscountToRateYearFrac(
    double       discount,                  /* (I) Discount factor */
    double       yearFraction,              /* (I) Year fraction */ 
    IrxTRateType basis,                     /* (I) Basis for the rate */
    double      *rate)                     /* (O) Output rate */
{
    static char routine[] = "irxDiscountToRateYearFrac";
    int         status    = FAILURE;

    REQUIRE (discount > 0.0);
    REQUIRE (IS_NOT_ZERO(yearFraction));
    REQUIRE (rate != NULL);

    switch (basis.type)
    {
    case IRX_RT_SIMPLE:
        *rate = (1.0/discount - 1.0)/yearFraction ;
        break;

    case IRX_RT_DISCOUNT:
        *rate = (1.0 - discount) / yearFraction;
        break;

    case IRX_RT_CONTINUOUS:
        *rate = (-log(discount)/yearFraction) ;
        break;

    case IRX_RT_FREQUENCY:
    {
        double freq = (double)basis.frequency;
        /* We prefer to be stodgy and not do the log/exp thing
         * here for performance.
         */
        *rate = freq * (pow(discount, -1.0/(freq*yearFraction)) - 1.0);
        break;
    }
    default:
        PROGRAM_BUG();
        goto RETURN; /* failure */
    }

    status = SUCCESS;

 RETURN:

    if (status != SUCCESS)
        irxErrorFailure(routine);

    return status;
}

int irxRateToDiscountYearFrac(
    double  rate,           /* (I) The rate */
    double  yearFraction,   /* (I) Year fraction */ 
    IrxTRateType basis,          /* (I) Basis for the rate */
    double *discount )      /* (O) Output discount rate */
{
    static char routine[] = "irxRateToDiscountYearFrac";
    int         status    = FAILURE;

    REQUIRE (discount != NULL);

    switch (basis.type)
    {
    case IRX_RT_SIMPLE:
    {
        double denom = 1.0 + rate * yearFraction;
        if (denom <= 0.0)
        {
            irxError("%s: Invalid simple interest rate:%f\n", routine, rate);
            *discount = 0.0;
            goto RETURN; /* failure */
        }
        *discount = 1.0 / denom;
        break;
    }
    case IRX_RT_DISCOUNT:
        *discount = 1.0 - rate * yearFraction;
        if (*discount <= 0.0)
        {
            irxError("%s: Invalid discount rate:%f\n", routine, rate);
            goto RETURN; /* failure */
        }
        break;

    case IRX_RT_CONTINUOUS:
        *discount = exp(-rate*yearFraction);
        break;

    case IRX_RT_FREQUENCY:
    {
        double tmp;
        double freq = (double)basis.frequency;
        REQUIRE(basis.frequency >= 1);
        REQUIRE(basis.frequency <= 365);

        tmp = 1.0 + rate / freq;
        if (tmp <= 0.0)
        {
            irxError("%s: Bad rate: %f.\n", routine, rate);
            goto RETURN; /* failure */
        }
        /* We prefer to be stodgy and not do the log/exp thing
         * here for performance.
         */
        *discount = pow( tmp, -freq*yearFraction);
        break;
    }
    default:
        PROGRAM_BUG();
        goto RETURN; /* failure */
    }

    status = SUCCESS;

 RETURN:

    if (status != SUCCESS)
        irxErrorFailure(routine);

    return status;
}

int irxRateValid(
    const char  *routine,           /* (I) Routine name to print */
    double rate,              /* (I) Rate to validate */
    IrxTDate  startDate,         /* (I) Starting date */
    IrxTDate  endDate,           /* (I) Ending date */
    IrxTDayCountConv  rateDayCountConv,  /* (I) Day count convention */
    IrxTRateType rateBasis)         /* (I) Compounding basis */
{
    int status = FAILURE;

    double yearFraction;

    switch (rateBasis.type)
    {
    case IRX_RT_SIMPLE:
    case IRX_RT_DISCOUNT:
        /* Only need year fraction in these cases.
         */
        if (irxDayCountFraction(startDate, 
                                endDate, 
                                rateDayCountConv,
                                &yearFraction) != SUCCESS)
            goto RETURN; /* failure */
        break;
        
    default:
        yearFraction = 1.0;
    }

    if (irxRateValidYearFrac(routine,
                             rate,
                             yearFraction,
                             rateBasis) != SUCCESS)
        goto RETURN; /* failure */

    status = SUCCESS;

 RETURN:

    return (status);
}


int irxRateValidYearFrac(
    const char  *routine,           /* (I) Routine name to print */
    double rate,              /* (I) Rate to validate */
    double yearFraction,      /* (I) Fraction of year */
    IrxTRateType basis)             /* (I) Compounding basis */
{
    int status = FAILURE;

    switch (basis.type)
    {
    case IRX_RT_SIMPLE:
        if (rate * yearFraction <= -1.0)
        {
            irxError("%s: Simple Rate (%f) * Year Fraction (%f) must "
                      "be > -1.0.\n",
                     routine, rate, yearFraction);
            goto RETURN;
        }
        break;
        
    case IRX_RT_DISCOUNT:
        if (rate * yearFraction >= 1.0)
        {
            irxError("%s: Discount Rate (%f) * Year Fraction (%f) must "
                     "be < 1.0.\n",
                     routine, rate, yearFraction);
            goto RETURN;
        }
        break;
        
    case IRX_RT_CONTINUOUS:
        /* Any rate is valid */
        break;

    case IRX_RT_FREQUENCY:
    {
        double freq = (double)(basis.frequency);
        if (rate <= -freq)
        {
            irxError("%s: Rate (%f) must be greater than (%f).\n",
                      routine, rate, -freq);
            goto RETURN;
        }
        break;
    }
    default:
        PROGRAM_BUG();
        goto RETURN; /* failure */
    }


    status = SUCCESS;

 RETURN:
    return (status);
}

