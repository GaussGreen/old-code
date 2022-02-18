/**********************************************************************
 SRT_F_SETGRFNPARAMS.C
 **********************************************************************/

#include "srt_h_all.h"

Err srt_f_set_default_GrfnParams(SrtGrfnParam* param)
{
    Err err = NULL;

    if (!param)
        return serror("Null param pointer passed in default_GrfnParams");

    memset(param, 0, sizeof(SrtGrfnParam));

    param->num_MCarlo_paths   = 5000;
    param->max_time_per_slice = 14.0 / 365.0;
    param->min_nodes_per_path = 50;
    param->force_mc           = SRT_NO;
    param->jumping            = SRT_NO;
    param->sample_type        = RANDOM_GAUSS_DYNAMIC;
    param->difference_scheme  = EULER;
    param->mc_renormalize     = SRT_YES;

    param->lsm         = SRT_FOR;
    param->minim       = SRT_NO;
    param->colmininf   = 0;
    param->colminsup   = 0;
    param->colmincible = 0;
    param->minmaxtime  = 0.0;
    param->minfreedom  = 0;

    param->rand_seed          = RANDINIT;
    param->width_phi          = 5;
    param->imm_exit           = SRT_NO;
    param->skip_evaluation    = SRT_NO;
    param->insert_vol_dates   = SRT_NO;
    param->end_of_day_flg     = SRT_NO;
    param->end_of_day_fixings = SRT_NO;
    param->end_of_day_payment = SRT_NO;

    param->prin_steps = SRT_NO;
    param->prin_deal  = SRT_NO;
    param->prin_ok    = SRT_NO;

    param->spline             = SRT_NO;
    param->dbg                = SRT_NO;
    param->max_nodes_per_path = 500;
    param->insert_tenor_dates = SRT_NO;

    param->step_num     = 500;
    param->imp_mdl_type = NONE;
    param->imp_type     = NONE;

    param->node_dim         = 1;
    param->closed_form_type = REAL_CLSDFRM;

    param->calib      = SRT_NO;
    param->first_time = SRT_YES;

    param->measure              = SPOT_MEAS;
    param->exfrontier           = 100;
    param->recpay               = SRT_PAYER;
    param->exfrontierfirst_time = SRT_NO;
    param->exfrontierlast_time  = SRT_NO;
    param->exfrontiercounter    = 0;
    param->exfrontiernumpoints  = 2;
    param->gammaexfrontier      = 0.1;
    param->MaxNumDisc           = 0;
    param->MinNumDisc           = 0;
    return err;
}

Err srt_f_set_GrfnParams(
    int numParams, String* paramStrings, String* valueStrings, SrtGrfnParam* param)
{
    Err            err = NULL;
    int            tmpInt, i;
    long           tmpLng;
    double         tmpDbl;
    SrtMCSamType   mcsam;
    SrtMCDfSchType mcdfs;
    SrtMdlType     mdl;
    SrtMdlDim      mdl_dim;
    SrtMeasureType meas;

    if (err = srt_f_set_default_GrfnParams(param))
        return err;

    for (i = 0; i < numParams; i++)
    {
        strupper(paramStrings[i]);
        strip_white_space(paramStrings[i]);

        /*1*/
        if ((!strcmp(paramStrings[i], "NUMPATH")) || (!strcmp(paramStrings[i], "PATHNUM")))
        {
            if (sscanf(valueStrings[i], "%ld", &tmpLng) == 1)
            {
                if (tmpLng > 0)
                {
                    param->num_MCarlo_paths = tmpLng;
                }
                else
                {
                    return serror("NUMPATH must be > 0");
                }
            }
            else
            {
                return serror("Failed to set NUMPATH");
            }
        }
        /*2*/
        else if (!strcmp(paramStrings[i], "MAXTIME"))
        {
            if (sscanf(valueStrings[i], "%lf", &tmpDbl) == 1)
            {
                if (tmpDbl > 0.0)
                {
                    param->max_time_per_slice = tmpDbl;
                }
                else
                {
                    return serror("MAXTIME must be > 0");
                }
            }
            else
            {
                return serror("Failed to set MAXTIME");
            }
        }
        /*3*/
        else if ((!strcmp(paramStrings[i], "MINNODE")) || (!strcmp(paramStrings[i], "MINSTEP")))
        {
            if (sscanf(valueStrings[i], "%ld", &tmpLng) == 1)
            {
                if (tmpLng > 0)
                {
                    param->min_nodes_per_path = tmpLng;
                }
                else
                {
                    return serror("MINNODE must be real > 0");
                }
            }
            else
            {
                return serror("Failed to set MINNODE");
            }
        }
        /*4*/
        else if (!strcmp(paramStrings[i], "DIFFSCHEME"))
        {
            err = srt_f_interp_mcdiffscheme(valueStrings[i], &mcdfs);
            if (err)
            {
                return err;
            }
            else
            {
                param->difference_scheme = mcdfs;
            }
        }
        /*5*/
        else if (!strcmp(paramStrings[i], "SAMPLETYPE"))
        {
            err = srt_f_interp_mcsample(valueStrings[i], &mcsam);
            if (err)
            {
                return err;
            }
            else
            {
                param->sample_type = mcsam;
            }
        }
        /*6*/
        else if (!strcmp(paramStrings[i], "RENORM"))
        {
            if (!strcmp(valueStrings[i], "YES"))
            {
                param->mc_renormalize = SRT_YES;
            }
            else if (!strcmp(valueStrings[i], "NO"))
            {
                param->mc_renormalize = SRT_NO;
            }
            else
            {
                return serror("RENORM must be YES or NO");
            }
        }
        /*7*/
        else if (!strcmp(paramStrings[i], "INSVOL"))
        {
            if (!strcmp(valueStrings[i], "YES"))
            {
                param->insert_vol_dates = SRT_YES;
            }
            else if (!strcmp(valueStrings[i], "NO"))
            {
                param->insert_vol_dates = SRT_NO;
            }
            else
            {
                return serror("INSVOL must be YES or NO");
            }
        }
        /*8*/
        else if (!strcmp(paramStrings[i], "INSTENOR"))
        {
            if (!strcmp(valueStrings[i], "YES"))
            {
                param->insert_tenor_dates = SRT_YES;
            }
            else if (!strcmp(valueStrings[i], "NO"))
            {
                param->insert_tenor_dates = SRT_NO;
            }
            else
            {
                return serror("INSVOL must be YES or NO");
            }
        }
        /*9*/
        else if (!strcmp(paramStrings[i], "FORCEMC"))
        {
            if (!strcmp(valueStrings[i], "YES"))
            {
                param->force_mc = SRT_YES;
            }
            else if (!strcmp(valueStrings[i], "NO"))
            {
                param->force_mc = SRT_NO;
            }
            else
            {
                return serror("FORCEMC must be YES or NO");
            }
        }
        /*10*/
        else if (
            (!strcmp(paramStrings[i], "WIDTHPHI")) || (!strcmp(paramStrings[i], "NUMPHI")) ||
            (!strcmp(paramStrings[i], "NUMPHIS")))
        {
            if (sscanf(valueStrings[i], "%d", &tmpInt) == 1)
            {
                if (tmpInt > 0)
                {
                    param->width_phi = tmpInt;
                }
                else
                {
                    return serror("WIDTHPHI must be integer > 0");
                }
            }
            else
            {
                return serror("Failed to set WIDTHPHI");
            }
        }
        /*11*/
        else if ((!strcmp(paramStrings[i], "DEBUG")) || (!strcmp(paramStrings[i], "DBG")))
        {
            if (!strcmp(valueStrings[i], "YES"))
            {
                param->dbg = SRT_YES;
            }
            else if (!strcmp(valueStrings[i], "NO"))
            {
                param->dbg = SRT_NO;
            }
            else
            {
                return serror("DEBUG must be YES or NO");
            }
        }
        /*12*/
        else if ((!strcmp(paramStrings[i], "PRINOK")) || (!strcmp(paramStrings[i], "PRINDEAL")))

        {
            if (!strcmp(valueStrings[i], "YES"))
            {
                param->prin_ok = SRT_YES;
            }
            else if (!strcmp(valueStrings[i], "NO"))
            {
                param->prin_ok = SRT_NO;
            }
            else
            {
                return serror("PRINDEAL must be YES or NO");
            }
        }
        /*13*/
        else if (!strcmp(paramStrings[i], "PRINSTEPS"))
        {
            if (!strcmp(valueStrings[i], "YES"))
            {
                param->prin_steps = SRT_YES;
            }
            else if (!strcmp(valueStrings[i], "NO"))
            {
                param->prin_steps = SRT_NO;
            }
            else
            {
                return serror("PRINSTEPS must be YES or NO");
            }
        }
        /*14*/
        else if ((!strcmp(paramStrings[i], "MAXNODE")) || (!strcmp(paramStrings[i], "MAXSTEPS")))
        {
            if (sscanf(valueStrings[i], "%ld", &tmpLng) == 1)
            {
                if (tmpLng > 0)
                {
                    param->max_nodes_per_path = tmpLng;
                }
                else
                {
                    return serror("MAXNODE must be real > 0");
                }
            }
            else
            {
                return serror("Failed to set MAXNODE");
            }
        }
        /*15*/
        else if (!strcmp(paramStrings[i], "SPLINE"))

        {
            if (!strcmp(valueStrings[i], "YES"))
            {
                param->spline = SRT_YES;
            }
            else if (!strcmp(valueStrings[i], "NO"))
            {
                param->spline = SRT_NO;
            }
            else
            {
                return serror("SPLINE must be YES or NO");
            }
        }
        /*16: for compatibility with old grfn only */
        else if ((!strcmp(paramStrings[i], "MODEL")) || (!strcmp(paramStrings[i], "MDL")))
        {
            err = srt_f_interp_model(valueStrings[i], &mdl, &mdl_dim);
            if (err)
            {
                return err;
            }
            else
            {
                param->imp_mdl_type = mdl;
            }
        }
        /*17: to skip the grfn evaluation */
        else if (
            (!strcmp(paramStrings[i], "SKIPEVAL")) ||
            (!strcmp(paramStrings[i], "SKIPEVALUATION")) ||
            (!strcmp(paramStrings[i], "SKIP_EVALUATION")))
        {
            if (!strcmp(valueStrings[i], "YES"))
            {
                param->skip_evaluation = SRT_YES;
            }
            else if (!strcmp(valueStrings[i], "NO"))
            {
                param->skip_evaluation = SRT_NO;
            }
            else
            {
                return serror("SKIPEVAL must be YES or NO");
            }
        }
        /*18: for an imminent exit */
        else if ((!strcmp(paramStrings[i], "IMMEXIT")) || (!strcmp(paramStrings[i], "EXIT")))
        {
            if (!strcmp(valueStrings[i], "YES"))
            {
                param->imm_exit = SRT_YES;
            }
            else if (!strcmp(valueStrings[i], "NO"))
            {
                param->imm_exit = SRT_NO;
            }
            else
            {
                return serror("IMMEXIT must be YES or NO");
            }
        }
        /* 19 : end of day flag */
        else if (
            (!strcmp(paramStrings[i], "ENDOFDAY")) || (!strcmp(paramStrings[i], "END_OF_DAY")) ||
            (!strcmp(paramStrings[i], "END_OF_DAY_FLAG")))
        {
            if (!strcmp(valueStrings[i], "YES"))
            {
                param->end_of_day_flg = SRT_YES;
            }
            else if (!strcmp(valueStrings[i], "NO"))
            {
                param->end_of_day_flg = SRT_NO;
            }
            else
            {
                return serror("ENDOFDAY must be YES or NO");
            }
        }
        /* 20 : end of day fixings ( if yes: will have to get them from database */
        else if (
            (!strcmp(paramStrings[i], "FIXINGS")) ||
            (!strcmp(paramStrings[i], "ENDOFDAYFIXINGS")) ||
            (!strcmp(paramStrings[i], "USETODAYFIXINGS")))
        {
            if (!strcmp(valueStrings[i], "YES"))
            {
                param->end_of_day_fixings = SRT_YES;
            }
            else if (!strcmp(valueStrings[i], "NO"))
            {
                param->end_of_day_fixings = SRT_NO;
            }
            else
            {
                return serror("ENDOFDAY must be YES or NO");
            }
        }
        /* 21 : initial seed for random number generation */
        else if (
            (!strcmp(paramStrings[i], "SEED")) || (!strcmp(paramStrings[i], "RANDSEED")) ||
            (!strcmp(paramStrings[i], "RAND_SEED")))
        {
            if (sscanf(valueStrings[i], "%ld", &tmpLng) == 1)
            {
                if (tmpLng < 0)
                {
                    param->rand_seed = tmpLng;
                }
                else
                {
                    return serror("RANDSEED must be < 0");
                }
            }
            else
            {
                return serror("Failed to set RANDSEED (try -123456789");
            }
        }
        /* 22: closed form type : Real One or Quick */
        else if (
            (!strcmp(paramStrings[i], "CLOSEDFORM")) || (!strcmp(paramStrings[i], "CLSDFRM")) ||
            (!strcmp(paramStrings[i], "CLOSED_FORM")))
        {
            if (!strcmp(valueStrings[i], "REAL"))
            {
                param->closed_form_type = REAL_CLSDFRM;
            }
            else if (!strcmp(valueStrings[i], "QUICK"))
            {
                param->closed_form_type = QUICK_CLSDFRM;
            }
            else
            {
                return serror("CLOSEDFORM must be REAL or QUICK");
            }
        }
        /*23: Jumping numeraire method */
        else if (
            (!strcmp(paramStrings[i], "JUMPING")) ||
            (!strcmp(paramStrings[i], "JUMPING_NUMERAIRE")) ||
            (!strcmp(paramStrings[i], "JUMPINGNUMERAIRE")) ||
            (!strcmp(paramStrings[i], "JUMPNUM")) || (!strcmp(paramStrings[i], "JUMP_NUM")))
        {
            if (!strcmp(valueStrings[i], "YES"))
            {
                param->jumping = SRT_YES;
            }
            else if (!strcmp(valueStrings[i], "NO"))
            {
                param->jumping = SRT_NO;
            }
            else
            {
                return serror("JUMPING must be YES or NO");
            }
        }
        /*24: Computation Measure*/
        else if (
            (!strcmp(paramStrings[i], "MEASURE")) || (!strcmp(paramStrings[i], "MEAS")) ||
            (!strcmp(paramStrings[i], "COMPUTATION_MEASURE")) ||
            (!strcmp(paramStrings[i], "COMPUTATION_MEAS")) ||
            (!strcmp(paramStrings[i], "COMPUTATIONMEASURE")) ||
            (!strcmp(paramStrings[i], "COMPUTATIONMEAS")))
        {
            err = srt_f_interp_measuretype(valueStrings[i], &meas);
            if (err)
            {
                return err;
            }
            else
            {
                param->measure = meas;
            }
        }
        /* 25: LSMMC */
        else if ((!strcmp(paramStrings[i], "LSM")))
        {
            if (!strcmp(valueStrings[i], "FORWARD"))
            {
                param->lsm = SRT_FOR;
            }
            else if (!strcmp(valueStrings[i], "BACKWARD"))
            {
                param->lsm = SRT_BACK;
            }
            else if (!strcmp(valueStrings[i], "FORWARD-BACKWARD"))
            {
                param->lsm = SRT_FORBACK;
            }
            else
            {
                return serror("LSM must be FORWARD, BACKWARD or FORWARD-BACKWARD");
            }
        }
        /* 26: MINIM */
        else if ((!strcmp(paramStrings[i], "MINIM")))
        {
            if (!strcmp(valueStrings[i], "YES"))
            {
                param->minim = SRT_YES;
            }
            else if (!strcmp(valueStrings[i], "NO"))
            {
                param->minim = SRT_NO;
            }
            else
            {
                return serror("MINIM must be YES or NO");
            }
        }
        /* 27 : colmininf */
        else if ((!strcmp(paramStrings[i], "COLMININF")) || (!strcmp(paramStrings[i], "COLLSMINF")))
        {
            if (sscanf(valueStrings[i], "%ld", &tmpLng) == 1)
            {
                if (tmpLng >= 0)
                {
                    param->colmininf = tmpLng;
                }
                else
                {
                    return serror("COLMININF must between 0 and colmax");
                }
            }
        }
        /* 28 : colminsup */
        else if ((!strcmp(paramStrings[i], "COLMINSUP")) || (!strcmp(paramStrings[i], "COLLSMSUP")))
        {
            if (sscanf(valueStrings[i], "%ld", &tmpLng) == 1)
            {
                if (tmpLng >= 0)
                {
                    param->colminsup = tmpLng;
                }
                else
                {
                    return serror("COLMINSUP must between 0 and colmax");
                }
            }
        }
        /* 29 : colmincible */
        else if (!strcmp(paramStrings[i], "COLMINCIBLE"))
        {
            if (sscanf(valueStrings[i], "%ld", &tmpLng) == 1)
            {
                if (tmpLng >= 0)
                {
                    param->colmincible = tmpLng;
                }
                else
                {
                    return serror("COLMINCIBLE must between 0 and colmax");
                }
            }
        }
        /* 30 : minmaxtime */
        else if (!strcmp(paramStrings[i], "MINMAXTIME"))
        {
            if (sscanf(valueStrings[i], "%lf", &tmpDbl) == 1)
            {
                if (tmpDbl > 0.0)
                {
                    param->minmaxtime = tmpDbl;
                }
                else
                {
                    return serror("MINMAXTIME must be > 0");
                }
            }
            else
            {
                return serror("Failed to set MINMAXTIME");
            }
        }
        /* 31 : colmincible */
        else if (
            (!strcmp(paramStrings[i], "COLMINFREEDOM")) || (!strcmp(paramStrings[i], "COLLSMDOM")))
        {
            if (sscanf(valueStrings[i], "%ld", &tmpLng) == 1)
            {
                if (tmpLng >= 0)
                {
                    param->minfreedom = tmpLng;
                }
                else
                {
                    return serror("COLMINFREEDOM must between 0 and colmax");
                }
            }
        }

        /* 32 : ExFrontier */
        else if ((!strcmp(paramStrings[i], "EXFRONTIER")))
        {
            if (sscanf(valueStrings[i], "%ld", &tmpLng) == 1)
            {
                if (tmpLng >= 0)
                {
                    param->exfrontier = tmpLng;
                }
                else
                {
                    return serror("EXFRONTIER must be a nonnegative integer");
                }
            }
        }

        /* 33 : ExFrontierNumPoints */
        else if ((!strcmp(paramStrings[i], "EXFRONTIERNUMPOINTS")))
        {
            if (sscanf(valueStrings[i], "%ld", &tmpLng) == 1)
            {
                if (tmpLng >= 0)
                {
                    param->exfrontiernumpoints = tmpLng;
                }
                else
                {
                    return serror("EXFRONTIERNUMPOINTS must be a nonnegative integer");
                }
            }
        }

        /* 34: RecPay */
        else if ((!strcmp(paramStrings[i], "RECPAY")))
        {
            if (!strcmp(valueStrings[i], "REC"))
            {
                param->recpay = SRT_RECEIVER;
            }
            else if (!strcmp(valueStrings[i], "PAY"))
            {
                param->recpay = SRT_PAYER;
            }

            else
            {
                return serror("RECPAY must be REC or PAY");
            }
        }

        /* 35 : GammaExFrontier */
        else if ((!strcmp(paramStrings[i], "GAMMAEXFRONTIER")))
        {
            if (sscanf(valueStrings[i], "%lf", &tmpDbl) == 1)
            {
                if (tmpDbl >= 0)
                {
                    param->gammaexfrontier = tmpDbl;
                }
                else
                {
                    return serror("GAMMAEXFRONTIER must be a nonnegative integer");
                }
            }
        }

        /* 36 : MaxNumDisc */
        else if ((!strcmp(paramStrings[i], "MAXNUMDISC")))
        {
            if (sscanf(valueStrings[i], "%ld", &tmpLng) == 1)
            {
                if (tmpLng >= 0)
                {
                    param->MaxNumDisc = tmpLng;
                }
                else
                {
                    return serror(" MaxNumDisc must be a non negative integer");
                }
            }
        }

        /* 37 : MinNumDisc */
        else if ((!strcmp(paramStrings[i], "MINNUMDISC")))
        {
            if (sscanf(valueStrings[i], "%ld", &tmpLng) == 1)
            {
                if (tmpLng >= 0)
                {
                    param->MinNumDisc = tmpLng;
                }
                else
                {
                    return serror(" MinNumDisc must be a non negative integer");
                }
            }
        }
        else if ((!strcmp(paramStrings[i], "MINNUMPDEMESH")))
        {
            if (sscanf(valueStrings[i], "%ld", &tmpLng) == 1)
            {
                if (tmpLng >= 0)
                {
                    param->min_pde_num_mesh = tmpLng;
                }
                else
                {
                    return serror(" min_pde_num_mesh must be a non negative integer");
                }
            }
        }
        else if (
            (!strcmp(paramStrings[i], "PAY")) || (!strcmp(paramStrings[i], "ENDOFDAYPAY")) ||
            (!strcmp(paramStrings[i], "USETODAYPAY")))
        {
            if (!strcmp(valueStrings[i], "YES"))
            {
                param->end_of_day_payment = SRT_YES;
            }
            else if (!strcmp(valueStrings[i], "NO"))
            {
                param->end_of_day_payment = SRT_NO;
            }
            else
            {
                return serror("ENDOFDAYPAY must be YES or NO");
            }
        }

        else
        {
            return serror("Model Parameter not recognised: %s", paramStrings[i]);
        }
    } /* end loop over numParams */
    /* 38 : PAYMENTFIXINGS */

    return NULL;
}
