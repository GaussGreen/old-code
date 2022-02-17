/***************************************************************************
 *      SCCS Keyword Information
 *      ------------------------
 *      Module name     :  ppconst.c
 *      Company Name    :  JP Morgan Securities Inc.
 *      Author          :  Davis (Shuenn-Tyan) Lee
 *                         Derivatives Research
 *      Code version    :  1.7
 *      Extracted       :  5/19/97 at 13:54:17
 *      Last Updated    :  5/19/97 at 13:54:07
 ***************************************************************************
 *      Older constructors/destructors used by Primus interface
 *
 *      Copyright 1996 J.P. Morgan & Co. Incorporated. All rights reserved.
 ***************************************************************************/
#include <stdio.h>
#include <stdlib.h>
extern "C" {
#include "cerror.h"
#include "tcurve.h"
}

#include "mbsbase.h"
#include "ppconst.h"
#include "mbsconst.h"
#include "ppbase.h"
/*static double DEF_MBS_PP_MGRP_GNMA_LT_AGES[] = {0.06, 1.0625, 0.085, \
                                                0.0100, 0.0100, 0.0938}; 
static double DEF_MBS_PP_MGRP_FNMA_LT_AGES[] = {0.06, 1.0625, 0.085, \
                                                0.1340, 0.99, 0.0521};	 */
/***************************************************************************
 *  InitMbsPrepayDefaults
 *  This function has to call once before running MbsAmort function which
 *  will setup all necessary GLOBAL default numbers to mbsPrepayDefaults'
 *  data structure.
 *  "Contingency" function to allow user to modify params of ARM/MGRP FIXED
 *  model
 *
 *  !! WARNING: this function should be used with caution--please
 *  contact Forrest Quinn, Davis Lee or Bob Lenk before using !!
 *
 *  INPUTS: if indicated parameter is NOT to be altered, use NULL;
 *          otherwise, the default model parameter(s) will be changed
 *          to the specified value, and will affect all subsequent calls
 *          to MbsAmort().
 *  FUNCTION RETURNS: 0=success, nonzero=failure
 ***************************************************************************/
int EXPORT
InitMbsPrepayDefaults
   (long prepayModel,
    long mbsAgency,
    long amortIndexType,
    /* Following inputs are user supplied to alter global default numbers;
     * for both ARM model and MGRP FIXED model
     */
    double *seasonality,        /* (I) ptr to new seasonality array
                                 * (12 values) */
    /* Following inputs are user supplied to alter global default numbers;
     * for ARM model only
     */
    long     maxLag,             /* (I) new value for max lag (months) */
    double  absoluteRateEffect, /* (I) absolute rate effect coeff. */
    double  histRefi30Idx,      /* (I) new value for historic value of
                                 * refi_index rate (FH30 commit rate) */
    double  histRefi15Idx,      /* (I) new value for historic value of
                                 * refi_index rate (FH15 commit rate) */
    double  histCmt10Idx,       /* (I) new value for historic value of
                                   refi_index rate (CMT10 rate) */
    double *armLagWgts,         /* (I) ptr to new array of lag wgts (expected
                                 * to contain max_lag values) */
    double *armBaseLogistic,    /* (I) ptr to new array of logistic params
                                 * (currently 4 values) for base prepay
                                 * function logistic */
    double *armSeasLogistic,    /* (I) ptr to new array of logistic params
                                 * (currently 4 values) for seasoning ramp
                                 * logistic */
    /* Following inputs are user supplied to alter global default numbers;
     * for MGRP FIXED model only
     */
    long    numGroups,           /* (I) the number of groups used in
                                 * multi-group */
    double critRat,             /* (I) shift refinance incentive according
                                 * to critRat */
    double gpAddlSmm,           /* (I) Additional smm for a group that is
                                 * triggered */
    double *scurveLogistic,     /* (I) ptr to new array of logistic params
                                 * (currently 6 values; max, min, width,
                                 * inflec,skew,disc decay) for base prepay
                                 * function logistic */
    double *seasLogistic,       /* (I) ptr to new array of logistic params
                                 * (currently 4 values) for seasoning ramp
                                 * logistic */
    double *groupPs,            /* (I) is an array with the six parameters
                                 * for the group distribution*/
    double *ltAges,             /* (I) praameters used in adjusting the
                                 * refinance ratio */
    double *mgrpLagWgts,        /* (I) Array of wgts to use w/prev refi
                                 * rates */
    TMbsPrepayDefaults **mbsPrepayDefaults)
{
    int status = FAILURE;
    long i;
    static char routine[] = "InitMbsPrepayDefaults";

    /* If needed, allocate the structure
     */
    if (*mbsPrepayDefaults IS NULL)
    {
        if ((*mbsPrepayDefaults = NEW(TMbsPrepayDefaults)) IS NULL)
        {
            GtoErrMsg("%s : Failed to allocate mbsPrepayDefaults structure\n",
                routine);
            goto done;
        }
    }

    (*mbsPrepayDefaults)->desiredRefiPoints = DEF_MBS_PP_DESIRED_REFI_POINTS;

    /* Determine IO multiple */
    switch (amortIndexType)
    {
    case MBS_PP_REFI_TYPE_FH30:
        (*mbsPrepayDefaults)->ioMultiplier = hist_conv30_io_multiple;
        break;
    case MBS_PP_REFI_TYPE_FH15:
        (*mbsPrepayDefaults)->ioMultiplier = hist_conv15_io_multiple;
        break;
    default:
        GtoErrMsg("%s: Invalid refi rate type: %d\n",routine, amortIndexType);
        goto done;
    }

    if (prepayModel == MBS_PP_MODEL_ARM)
    {
        (*mbsPrepayDefaults)->maxLag = DEF_MBS_PP_ARM_MAX_LAG;
        for (i = 0; i < DEF_MBS_PP_ARM_MAX_NUM_LAGMONS; i++)
        {
            (*mbsPrepayDefaults)->armLagWgts[i] = DEF_MBS_PP_ARM_LAGWGTS[i];
        }

        (*mbsPrepayDefaults)->smmLogistRight = DEF_MBS_PP_ARM_SMM_LOGIST_RIGHT;
        (*mbsPrepayDefaults)->smmLogistLeft = DEF_MBS_PP_ARM_SMM_LOGIST_LEFT;
        (*mbsPrepayDefaults)->smmLogistWidth = DEF_MBS_PP_ARM_SMM_LOGIST_WIDTH;
        (*mbsPrepayDefaults)->smmLogistInflec = DEF_MBS_PP_ARM_SMM_LOGIST_INFLEC;

        (*mbsPrepayDefaults)->seasLogistRight = DEF_MBS_PP_ARM_SEAS_LOGIST_RIGHT;
        (*mbsPrepayDefaults)->seasLogistLeft = DEF_MBS_PP_ARM_SEAS_LOGIST_LEFT;
        (*mbsPrepayDefaults)->seasLogistWidth  = DEF_MBS_PP_ARM_SEAS_LOGIST_WIDTH;
        (*mbsPrepayDefaults)->seasLogistInflec = DEF_MBS_PP_ARM_SEAS_LOGIST_INFLEC;

        (*mbsPrepayDefaults)->histFh30Idx = DEF_MBS_PP_ARM_HIST_FH30_IDX;
        (*mbsPrepayDefaults)->histFh15Idx = DEF_MBS_PP_ARM_HIST_FH15_IDX;
        (*mbsPrepayDefaults)->histCmt10Idx = DEF_MBS_PP_ARM_HIST_CMT10_IDX;

        (*mbsPrepayDefaults)->refiIncShift = DEF_MBS_PP_ARM_INC_SHIFT;
        (*mbsPrepayDefaults)->newWacRefiIncShift = DEF_MBS_PP_ARM_NEW_WAC_INC_SHIFT;

        (*mbsPrepayDefaults)->absoluteRateEffect = DEF_MBS_PP_ARM_ABS_RATE_EFF;
        for (i = 0; i < DEF_MBS_PP_NUM_SEASONALITY; i++)
        {
            (*mbsPrepayDefaults)->seasonality[i] = DEF_MBS_PP_ARM_SEASONALITY[i];
        }
    }
    else if (prepayModel == MBS_PP_MODEL_FIX_MGRP)
    {
        (*mbsPrepayDefaults)->numGroups = DEF_MBS_PP_MGRP_NUM_GROUPPS;
        (*mbsPrepayDefaults)->critRat = DEF_MBS_PP_MGRP_CRITRAT;
        if (mbsAgency == MBS_AGENCY_GNMAI ||
            mbsAgency == MBS_AGENCY_GNMAII )
        {
            (*mbsPrepayDefaults)->gpAddlSmm = DEF_MBS_PP_MGRP_GNMA_GPADDLSMM;
            for (i = 0; i < DEF_MBS_PP_NUM_SEASONALITY; i++)
            {
                (*mbsPrepayDefaults)->seasonality[i] = DEF_MBS_PP_MGRP_GNMA_SEASONALITY[i];
            }
            for (i = 0; i < DEF_MBS_PP_MGRP_NUM_SCURVE_LOGIST; i++)
            {
                (*mbsPrepayDefaults)->scurveLogist[i] = DEF_MBS_PP_MGRP_GNMA_SCURVE_LOGIST[i];
            }
            for (i = 0; i < DEF_MBS_PP_MGRP_NUM_SEAS_LOGIST; i++)
            {
                (*mbsPrepayDefaults)->seasLogist[i] = DEF_MBS_PP_MGRP_GNMA_SEAS_LOGIST[i];
            }
            for (i = 0; i < DEF_MBS_PP_MGRP_NUM_GROUPPS; i++)
            {
                (*mbsPrepayDefaults)->groupPs[i] = DEF_MBS_PP_MGRP_GNMA_GROUPPS[i];
            }
            for (i = 0; i < DEF_MBS_PP_MGRP_NUM_LTAGES; i++)
            {
                (*mbsPrepayDefaults)->ltAges[i] = DEF_MBS_PP_MGRP_GNMA_LT_AGES[i];
            }
            for (i = 0; i < DEF_MBS_PP_MGRP_NUM_WTLAGS; i++)
            {
                (*mbsPrepayDefaults)->mgrpLagWgts[i] = DEF_MBS_PP_MGRP_GNMA_WTLAGS[i];
            }
        }
        else if (mbsAgency == MBS_AGENCY_FNMA ||
                 mbsAgency == MBS_AGENCY_FHLMC ||
                 mbsAgency == MBS_AGENCY_GOLD ||
                 mbsAgency == MBS_AGENCY_WHOLE)
        {
            (*mbsPrepayDefaults)->gpAddlSmm = DEF_MBS_PP_MGRP_FNMA_GPADDLSMM;
            for (i = 0; i < DEF_MBS_PP_NUM_SEASONALITY; i++)
            {
                (*mbsPrepayDefaults)->seasonality[i] = DEF_MBS_PP_MGRP_FNMA_SEASONALITY[i];
            }
            for (i = 0; i < DEF_MBS_PP_MGRP_NUM_SCURVE_LOGIST; i++)
            {
                (*mbsPrepayDefaults)->scurveLogist[i] = DEF_MBS_PP_MGRP_FNMA_SCURVE_LOGIST[i];
            }
            for (i = 0; i < DEF_MBS_PP_MGRP_NUM_SEAS_LOGIST; i++)
            {
                (*mbsPrepayDefaults)->seasLogist[i] = DEF_MBS_PP_MGRP_FNMA_SEAS_LOGIST[i];
            }
            for (i = 0; i < DEF_MBS_PP_MGRP_NUM_GROUPPS; i++)
            {
                (*mbsPrepayDefaults)->groupPs[i] = DEF_MBS_PP_MGRP_FNMA_GROUPPS[i];
            }
            for (i = 0; i < DEF_MBS_PP_MGRP_NUM_LTAGES; i++)
            {
                (*mbsPrepayDefaults)->ltAges[i] = DEF_MBS_PP_MGRP_FNMA_LT_AGES[i];
            }
            for (i = 0; i < DEF_MBS_PP_MGRP_NUM_WTLAGS; i++)
            {
                (*mbsPrepayDefaults)->mgrpLagWgts[i] = DEF_MBS_PP_MGRP_FNMA_WTLAGS[i];
            }
        }
    }
    else
    {
        GtoErrMsg("%s : No such prepay model [%d]\n",routine,prepayModel);
        goto done;
    }

    /* Only modify if user actually supplied something */
    if (seasonality)
    {
        for (i = 0; i < DEF_MBS_PP_NUM_SEASONALITY; i++)
        {
            (*mbsPrepayDefaults)->seasonality[i] = seasonality[i];
        }
    }
    if (prepayModel == MBS_PP_MODEL_ARM)
    {
        if (maxLag)
        {
            if (maxLag > DEF_MBS_PP_ARM_MAX_NUM_LAGMONS)
            {
                GtoErrMsg("%s: max_lag too big--at most 12 months\n",routine);
                goto done;
            }
            (*mbsPrepayDefaults)->maxLag = maxLag;
        }
        if (armLagWgts)
        {
            if (!maxLag)
            {
                GtoErrMsg("%s: To modify lagwgts, must also supply maxLag\n",
                    routine);
                goto done;
            }
            for (i = 0; i < maxLag; i++)
            {
                (*mbsPrepayDefaults)->armLagWgts[i] = armLagWgts[i];
            }
        }
        if (absoluteRateEffect)
        {
            (*mbsPrepayDefaults)->absoluteRateEffect = absoluteRateEffect;
        }
        if (histRefi30Idx)
        {
            (*mbsPrepayDefaults)->histFh30Idx = histRefi30Idx;
        }
        if (histRefi15Idx)
        {
            (*mbsPrepayDefaults)->histFh15Idx = histRefi15Idx;
        }
        if (histCmt10Idx)
        {
            (*mbsPrepayDefaults)->histCmt10Idx = histCmt10Idx;
        }
        if (armBaseLogistic)
        {
            (*mbsPrepayDefaults)->smmLogistRight = armBaseLogistic[0];
            (*mbsPrepayDefaults)->smmLogistLeft = armBaseLogistic[1];
            (*mbsPrepayDefaults)->smmLogistWidth = armBaseLogistic[2];
            (*mbsPrepayDefaults)->smmLogistInflec = armBaseLogistic[3];
        }
        if (armSeasLogistic)
        {
            (*mbsPrepayDefaults)->seasLogistRight = armSeasLogistic[0];
            (*mbsPrepayDefaults)->seasLogistLeft = armSeasLogistic[1];
            (*mbsPrepayDefaults)->seasLogistWidth = armSeasLogistic[2];
            (*mbsPrepayDefaults)->seasLogistInflec = armSeasLogistic[3];
        }
    }
    else if (prepayModel == MBS_PP_MODEL_FIX_MGRP)
    {
        if (critRat)
        {
            (*mbsPrepayDefaults)->critRat = critRat;
        }
        if (numGroups)
        {
            (*mbsPrepayDefaults)->numGroups = numGroups;
        }
        if (gpAddlSmm)
        {
            (*mbsPrepayDefaults)->gpAddlSmm = gpAddlSmm;
        }
        if (scurveLogistic)
        {
            for (i = 0; i < DEF_MBS_PP_MGRP_NUM_SCURVE_LOGIST; i++)
            {
                (*mbsPrepayDefaults)->scurveLogist[i] = scurveLogistic[i];
            }
        }
        if (seasLogistic)
        {
            for (i = 0; i < DEF_MBS_PP_MGRP_NUM_SEAS_LOGIST; i++)
            {
                (*mbsPrepayDefaults)->seasLogist[i] = seasLogistic[i];
            }
        }
        if (groupPs)
        {
            for (i = 0; i < DEF_MBS_PP_MGRP_NUM_GROUPPS; i++)
            {
                (*mbsPrepayDefaults)->groupPs[i] = groupPs[i];
            }
        }
        if (ltAges)
        {
            for (i = 0; i < DEF_MBS_PP_MGRP_NUM_LTAGES; i++)
            {
                (*mbsPrepayDefaults)->ltAges[i] = ltAges[i];
            }
        }
        if (mgrpLagWgts)
        {
            for (i = 0; i < DEF_MBS_PP_MGRP_NUM_WTLAGS; i++)
            {
                (*mbsPrepayDefaults)->mgrpLagWgts[i] = mgrpLagWgts[i];
            }
        }
    }
    
    status = SUCCESS;
 
done:
    if (status IS FAILURE)
    {
        GtoErrMsg("%s : Failed\n",routine);
    }
    return(status);
}


/***************************************************************************
 * FreeMbsPrepayDefaults()
 * Destructor for TMbsPrepayDefaults
 * (also resets pointer to NULL)
 ***************************************************************************/
void EXPORT
FreeMbsPrepayDefaults
   (TMbsPrepayDefaults **mbsPrepayDefaults)
{
    static char routine[] = "FreeMbsPrepayDefaults";

    if (*mbsPrepayDefaults ISNT NULL)
    {
        /* free internal dynamic vars */
        /* @@ none for now */
        /* free entire structure */
        FREE(*mbsPrepayDefaults);
        *mbsPrepayDefaults = NULL;
    }
}


/***************************************************************************
 * PrintMbsPrepayDefaults()
 * Pint out contents of deal struc for debugging
 * Returns SUCCESS/FAILURE
 ***************************************************************************/
int EXPORT
PrintMbsPrepayDefaults
   (long prepayModel,
    TMbsPrepayDefaults *mbsPrepayDefaults) /* (I) structure contains default
                                            * numbers */
{
    long i;

    GtoErrMsg("\n*********************************************\n");
    GtoErrMsg("***      FOR CHECKING/REFERENCE ONLY      ***\n");
    GtoErrMsg("*** Print contents of MBS Prepay Defaults ***\n");
    GtoErrMsg("*********************************************\n");
    for (i = 0; i < DEF_MBS_PP_NUM_SEASONALITY; i++)
    {
        GtoErrMsg("I = %d Seasonality = %f\n",
            i,mbsPrepayDefaults->seasonality[i]);
    }
    if (prepayModel == MBS_PP_MODEL_ARM)
    {
        GtoErrMsg("Max. lag = %d\n",mbsPrepayDefaults->maxLag);
        for (i = 0; i < mbsPrepayDefaults->maxLag; i++)
        {
            GtoErrMsg("I = %d Lagwght = %f\n",
                i,mbsPrepayDefaults->armLagWgts[i]);
        }
        GtoErrMsg("Abs. rate effect = %f\n",
            mbsPrepayDefaults->absoluteRateEffect);
        GtoErrMsg("Hist fh30 idx. = %f\n",mbsPrepayDefaults->histFh30Idx);
        GtoErrMsg("Hist fh30 idx. = %f\n",mbsPrepayDefaults->histFh15Idx);
        GtoErrMsg("Hist cmt10 idx. = %f\n",mbsPrepayDefaults->histCmt10Idx);
        GtoErrMsg("SMM logist right = %f\n",mbsPrepayDefaults->smmLogistRight);
        GtoErrMsg("SMM logist left = %f\n",mbsPrepayDefaults->smmLogistLeft);
        GtoErrMsg("SMM logist width = %f\n",mbsPrepayDefaults->smmLogistWidth);
        GtoErrMsg("SMM logist inflec = %f\n",
            mbsPrepayDefaults->smmLogistInflec);
        GtoErrMsg("SEAS logist right = %f\n",
            mbsPrepayDefaults->seasLogistRight);
        GtoErrMsg("SEAS logist left = %f\n",mbsPrepayDefaults->seasLogistLeft);
        GtoErrMsg("SEAS logist width = %f\n",
            mbsPrepayDefaults->seasLogistWidth);
        GtoErrMsg("SEAS logist inflec = %f\n",
            mbsPrepayDefaults->seasLogistInflec);
    }
    else if (prepayModel == MBS_PP_MODEL_FIX_MGRP)
    {
        GtoErrMsg("CritRat = %f\n",mbsPrepayDefaults->critRat);
        GtoErrMsg("# of groups = %d\n",mbsPrepayDefaults->numGroups);
        GtoErrMsg("Group Add SMM = %f\n",mbsPrepayDefaults->gpAddlSmm);
        for (i = 0; i < DEF_MBS_PP_MGRP_NUM_SCURVE_LOGIST; i++)
        {
            GtoErrMsg("I = %d Scurve logist = %f\n",
                i,mbsPrepayDefaults->scurveLogist[i]);
        }
        for (i = 0; i < DEF_MBS_PP_MGRP_NUM_SEAS_LOGIST; i++)
        {
            GtoErrMsg("I = %d Seas logist = %f\n",
                i,mbsPrepayDefaults->seasLogist[i]);
        }
        for (i = 0; i < DEF_MBS_PP_MGRP_NUM_GROUPPS; i++)
        {
            GtoErrMsg("I = %d Group distribution = %f\n",
                i,mbsPrepayDefaults->groupPs[i]);
        }
        for (i = 0; i < DEF_MBS_PP_MGRP_NUM_LTAGES; i++)
        {
            GtoErrMsg("I = %d longer term aging parameters = %f\n",
                i,mbsPrepayDefaults->ltAges[i]);
        }
        for (i = 0; i < DEF_MBS_PP_MGRP_NUM_WTLAGS; i++)
        {
            GtoErrMsg("I = %d MRGP lagwght = %f\n",
                i,mbsPrepayDefaults->mgrpLagWgts[i]);
        }
    }

    return(SUCCESS);
}


/***************************************************************************
 * make_mbs_pp_modelinfo_const();
 * Constructor for MBS_PP_MODELINFO_CONST
 * Sets internals of structure to input values
 * If input ptr is NULL, also allocates structure
 * Returns SUCCESS/FAILURE
 ***************************************************************************/
int
make_mbs_pp_modelinfo_const
   (double  const_prepay,       /* (I) constant prepay rate to use */
    long     const_prepay_form,  /* (I) form of constant prepay rate;
                                   MBS_PP_SPD_SMM, etc. */
    /* PTR TO NEW/EXISTING STRUCTURE */
    MBS_PP_MODELINFO_CONST 
          **model_data)         /* (I/O) structure to alloc/modify */
{
    int status = FAILURE;
    static char routine[] = "make_mbs_pp_modelinfo_const";


    /* Check inputs
     */
    if (!is_valid_prepay_form(const_prepay_form))
    {
        GtoErrMsg("%s : Invalid form for const prepays: %d\n",
                  routine, const_prepay_form);
        goto done;
    }

    /* If needed, allocate the structure 
     */
    if (*model_data IS NULL)
    {
        *model_data = 
            (MBS_PP_MODELINFO_CONST*) malloc(sizeof(MBS_PP_MODELINFO_CONST));
        if (*model_data IS NULL)
        {
            GtoErrMsg("%s : Failed to allocate model data structure\n",
                      routine);
            goto done;
        }
    }

    /* Set member data 
     */
    (*model_data)->prepayModel = MBS_PP_MODEL_CONST;
    (*model_data)->const_prepay = const_prepay;
    (*model_data)->const_prepay_form = const_prepay_form;

    status = SUCCESS;

done:
    if (status ISNT SUCCESS)
    {
        free_mbs_pp_modelinfo_const(model_data);
    }
    return (status);
}


/***************************************************************************
 * free_mbs_pp_modelinfo_const()
 * Destructor for MBS_PP_MODELINFO_CONST
 * (also resets pointer to NULL)
 ***************************************************************************/
void
free_mbs_pp_modelinfo_const
   (MBS_PP_MODELINFO_CONST **model_data)     /* (I/O) structure to free */
{
    static char routine[] = "free_mbs_pp_modelinfo_const";
    
    if (model_data ISNT NULL)
    {
        /* free internal dynamic vars */
        /* (none for now) */
        /* free entire structure */
        free(*model_data);
        *model_data = NULL;
    }
}

/***************************************************************************
 * make_mbs_pp_modelinfo_vector();
 * Constructor for MBS_PP_MODELINFO_VECTOR
 * Sets internals of structure to input values
 * If input ptr is NULL, also allocates structure
 * Returns SUCCESS/FAILURE
 ***************************************************************************/
int
make_mbs_pp_modelinfo_vector
   (long      vector_prepay_form,   /* (I) form of prepay speeds in vector table;
                                      e.g., MBS_PP_SPD_SMM */
    long      vector_num_ages,      /* (I) # age categories (rows) in prepay 
                                      vector table */
    long      vector_num_ratemoves, /* (I) # rate move columns in prepay vector
                                      table */
    long     *vector_ages,          /* (I) ages (in months) of rows of vector
                                      table */
    double  *vector_ratemoves,     /* (I) rate move sizes (decimal, e.g., +100bp
                                      as 0.01) */
    double **vector_prepays,       /* (I) table of prepay speeds, indexed as:
                                      vector_prepays[i_age][j_ratemove]  */
    long      vector_prepay_lag,    /* (I) # months lag between rate changes
                                      and effect on prepays */
    /* PTR TO NEW/EXISTING STRUCTURE */
    MBS_PP_MODELINFO_VECTOR 
           **model_data)           /* (I/O) structure to alloc/modify */
{
    int status = FAILURE;
    static char routine[] = "make_mbs_pp_modelinfo_vector";
    long iage,
        irate;

    /* Check inputs
     */
    if (!is_valid_prepay_form(vector_prepay_form))
    {
        GtoErrMsg("%s : Invalid form for vector prepays: %d\n",
                  routine, vector_prepay_form);
        goto done;
    }
    if (vector_num_ages < 1)
    {
        GtoErrMsg("%s : Must have at least one age row in prepay table;\n"
                  " user specified %d\n",
                  routine, vector_num_ages);
        goto done;
    }
    if (vector_num_ratemoves < 1)
    {
        GtoErrMsg("%s : Must have at least one ratemove column "
                  "in prepay table;\n"
                  " user specified %d\n",
                  routine, vector_num_ratemoves);
        goto done;
    }

    /* If needed, allocate the structure 
     */
    if (*model_data IS NULL)
    {
        *model_data = 
            (MBS_PP_MODELINFO_VECTOR*) malloc(sizeof(MBS_PP_MODELINFO_VECTOR));
        if (*model_data IS NULL)
        {
            GtoErrMsg("%s : Failed to allocate model data structure\n",
                      routine);
            goto done;
        }
    }
    /* For an existing structure, de-allocate internal dynamic vars */
    else
    {
        mbs_free_ivector( (*model_data)->vector_ages, 
                         0, (*model_data)->vector_num_ages-1);
        mbs_free_dvector( (*model_data)->vector_ratemoves, 
                         0, (*model_data)->vector_num_ratemoves-1);
        mbs_free_dmatrix( (*model_data)->vector_prepays,
                         0,(*model_data)->vector_num_ages-1,
                         0,(*model_data)->vector_num_ratemoves-1);
    }

    /* Now, in either case, (re)allocate internal dynamic vars */
    if ((mbs_ivector(0,vector_num_ages-1, 
                    &((*model_data)->vector_ages)) ISNT SUCCESS ) ||
        (mbs_dvector(0,vector_num_ratemoves-1, 
                    &((*model_data)->vector_ratemoves)) ISNT SUCCESS ) ||
        (mbs_dmatrix(0,vector_num_ages-1,
                     0,vector_num_ratemoves-1,
                    &((*model_data)->vector_prepays)) ISNT SUCCESS ))
    {
        GtoErrMsg("%s : Failed to allocate internal dynamic vars\n",
                  routine);
        goto done;
    }

    /* Set member data 
     */
    (*model_data)->prepayModel = MBS_PP_MODEL_VECTOR;
    (*model_data)->vector_prepay_form = vector_prepay_form;
    (*model_data)->vector_num_ages = vector_num_ages;
    (*model_data)->vector_num_ratemoves = vector_num_ratemoves;
    for (iage=0; iage<vector_num_ages; iage++)
    {
        (*model_data)->vector_ages[iage] = vector_ages[iage];
    }
    for (irate=0; irate<vector_num_ratemoves; irate++)
    {
        (*model_data)->vector_ratemoves[irate] = vector_ratemoves[irate];
    }
    for (iage=0; iage<vector_num_ages; iage++)
    {
        for (irate=0; irate<vector_num_ratemoves; irate++)
        {
            (*model_data)->vector_prepays[iage][irate] =
                vector_prepays[iage][irate];
        }
    }
    (*model_data)->vector_prepay_lag = vector_prepay_lag;

    status = SUCCESS;

done:
    if (status ISNT SUCCESS)
    {
        free_mbs_pp_modelinfo_vector(model_data);
    }
    return (status);
}


/***************************************************************************
 * free_mbs_pp_modelinfo_vector()
 * Destructor for MBS_PP_MODELINFO_VECTOR
 * (also resets pointer to NULL)
 ***************************************************************************/
void
free_mbs_pp_modelinfo_vector
   (MBS_PP_MODELINFO_VECTOR **model_data)      /* (I/O) structure to free */
{
    static char routine[] = "free_mbs_pp_modelinfo_vector";
    
    if (model_data ISNT NULL)
    {
        /* free internal dynamic vars */
        mbs_free_ivector( (*model_data)->vector_ages, 
                         0, (*model_data)->vector_num_ages-1);
        mbs_free_dvector( (*model_data)->vector_ratemoves, 
                         0, (*model_data)->vector_num_ratemoves-1);
        mbs_free_dmatrix( (*model_data)->vector_prepays,
                         0,(*model_data)->vector_num_ages-1,
                         0,(*model_data)->vector_num_ratemoves-1);

        /* free entire structure */
        free(*model_data);
        *model_data = NULL;
    }
}
