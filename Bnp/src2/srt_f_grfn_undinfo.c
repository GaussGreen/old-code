/* ===============================================================================

   FILENAME	:	srt_f_und_info.c

   PURPOSE :   the structure to store and the functions to work with the global
               underlyings info for a multi-underlyings evaluation in Grfn:
               SrtUndInfo

   =============================================================================== */

#include "srt_h_all.h"
#include "srt_h_grfn_undinfo.h"

/* --------------------------------------------------------------------------------- */
/* Set all the SrtUndInfo flags from the list of undrlyings stored in the SrtUndInfo */
Err analyze_und_info(SrtUndInfo* und_info)
{
    int           i;
    char          domestic[SRTBUFSZ];
    String        und_name;
    String*       used_und_names = NULL;
    SrtUndPtr     und            = NULL;
    Err           err            = NULL;
    SrtCorrLstPtr cls            = NULL;
    SrtMdlType    mdl_type;
    SrtMdlDim     mdl_dim;

    /* there should at least be one underlying returned from the call, which
       is the domestic underlying */

    /* Default values for the flags */

    /* number of brownian motions to model */
    und_info->no_of_brownians = 0;

    /* initialises to NO the use of a two factor model (default) */
    und_info->two_factor_model = SRT_NO;

    /* initialises to NO the use of a stochastic volatility model (default) */
    und_info->use_stochastic_vol = SRT_NO;

    /* assume that a jumping discretisation is used, until we know otherwise */
    und_info->jumping = SRT_YES;

    /* domestic underlying interest rate is assumed to be first input */
    strcpy(domestic, und_info->und_data[0].und_name);

    /* initialise correlation Term Structure (needed to prevent bugs)*/
    und_info->corr_ts = NULL;

    /* by default we use the correlation saved in the static corr list */
    und_info->use_corr_ts = SRT_YES;

    /* Set up Underlying specific information */

    for (i = 0; i < und_info->no_of_underlyings; i++)
    {
        /* Get the ith underlying name in the SrtUndInfo */
        und_name = und_info->und_data[i].und_name;

        /* Check that the underlying is defined */
        if ((und = lookup_und(und_name)) == NULL)
            return (serror("Underlying %s not initialize", und_name));

        /* See if stochastic process is required for underlying */
        if (err = get_underlying_mdltype(und, &mdl_type))
            return err;

        if (mdl_type == NONE) /* not stochastic */
        {
            und_info->und_data[i].stochastic = SRT_NO;
        }

        else /* stochastic */
        {
            /* Set flag for stochastic underlying */
            und_info->und_data[i].stochastic = SRT_YES;

            if (err = get_underlying_mdldim(und, &mdl_dim))
                return err;

            /* Check for two factor cases */
            if (mdl_dim == TWO_FAC)
            {
                /* Two factor model is not implemented for multiple underlyings */
                if (und_info->no_of_underlyings == 1)
                {
                    /* set flag : there is a two factor model*/
                    und_info->two_factor_model = SRT_YES;

                    /* two factor model requires two brownians */
                    und_info->no_of_brownians = 2;
                }
                else
                {
                    /*
                   return(serror("Two Factor not available in Multi Asset"));
                    */

                    und_info->no_of_brownians += 2;
                    /*und_info->two_factor_model = SRT_YES;*/
                }
            } /* END if mdl_dim == STOCH_VOL  */

            else
                /* Check for stochastic volatility */
                if (mdl_type == LGM_STOCH_VOL || mdl_type == CHEY_STOCH_VOL ||
                    mdl_type == CHEY_BETA_STOCH_VOL || mdl_type == EQ_STOCH_RATES_SRVGS)
            {
                /* Stochastic volatility model not implemented in Multi-Underlying Exept in
                 * EQ_STOCH_RATES_SRVGS (fudge) */
                if ((!(mdl_type == EQ_STOCH_RATES_SRVGS)) && (!(mdl_type == LGM_STOCH_VOL)))
                {
                    if (und_info->no_of_underlyings > 1)
                        return (serror("Stochastic Volatility Model not available in Multi Asset"));
                    else
                    {
                        /* Set flag : we use stochastic vol */
                        und_info->use_stochastic_vol = SRT_YES;

                        /* A stochastic vol model requires two brownians */
                        und_info->no_of_brownians = 2;
                    }
                }
                else if (mdl_type == EQ_STOCH_RATES_SRVGS)
                {
                    und_info->no_of_brownians += 2;
                }
                else if (mdl_type == LGM_STOCH_VOL)
                {
                    und_info->no_of_brownians += 2;
                    und_info->use_stochastic_vol = SRT_YES;
                }

            } /* END if mdl_dim == TWO_FAC */

            else
            /* ONE FACTOR MODEL */
            {
                /* one factor model requires one brownian */
                und_info->no_of_brownians++;
            }

        } /* END if (mdl_type !=NONE)   <==> stochastic */

        /* Cheyette model does not exist in jumping numeraire format */
        if ((mdl_type == CHEY) || (mdl_type == CHEY_BETA))
            und_info->jumping = SRT_NO;
        else
            /* BetaEta model does not exist in jumping numeraire format */
            if (mdl_type == ETABETA)
            und_info->jumping = SRT_NO;
        else
            /* OVE && REGIS && JL: FX_STOCH_RATES model does not allow jumping numeraire YET */
            if (mdl_type == FX_STOCH_RATES)
            und_info->jumping = SRT_YES;
        else
            /* FXBETADLM model doesn't use the corr matrix */
            if (mdl_type == FX_BETADLM)
            und_info->use_corr_ts = SRT_NO;

        /* STOCH_VOL models do not exist in jumping numeraire format */
        if (und_info->use_stochastic_vol == SRT_YES)
            und_info->jumping = SRT_NO;
        /* EQ_STOCH_VOL model no jumping numeraire */
        if ((mdl_type == EQ_STOCH_RATES) || (mdl_type == EQ_STOCH_RATES_SRVGS))
            und_info->jumping = SRT_NO;

        /* NEW MODEL TYPE */
        if (is_model_New_type(mdl_type))
            und_info->jumping = SRT_NO;
        /* C.Godart Black_Und does not use corr matrix */
        if (!strcmp(und->underl_lbl, "BLACK_UND"))
            und_info->use_corr_ts = SRT_NO;

        /* SABR_UND does not use corr matrix */
        if (!strcmp(und->underl_lbl, "SABR_UND"))
            und_info->use_corr_ts = SRT_NO;

    } /* END   	for (i=0 ; i<und_info->no_of_underlyings ; i++) loop */

    /* Builds and attaches the correlation term structure for these underlyings */
    if (und_info->no_of_brownians > 1 && und_info->use_corr_ts == SRT_YES)
    {
        if (und_info->two_factor_model == SRT_YES || und_info->use_stochastic_vol == SRT_YES)
        {
            /* Attach a NULL SrtCorrLst to the und_info */
            und_info->corr_ts = NULL;
        } /* END if (two_factor_model == SRT_YES || stoch_vol == SRT_YES ) */
        else
        {
            /* Just stacks the underlying names in one single vector */
            used_und_names = svector(0, und_info->no_of_underlyings - 1);
            for (i = 0; i < und_info->no_of_underlyings; i++)
                used_und_names[i] = und_info->und_data[i].und_name;

            /* Creates and fills the local SrtCorrLst with relevent
                    underlying correlation matrix*/
            err = srt_f_make_deal_corrlist(
                used_und_names,
                und_info->no_of_underlyings,
                "Local deal correlation list",
                srt_f_GetTheCorrelationList(),
                &cls);
            if (err)
            {
                free_svector(used_und_names, 0, und_info->no_of_underlyings - 1);
                return err;
            }

            /* Attach the local SrtCorrLst to the und_info */
            und_info->corr_ts = cls;

            /*  Frees the svector of used_und_names: do not need it anymore */
            free_svector(used_und_names, 0, und_info->no_of_underlyings - 1);

        } /* END if (und_info->two_factor_model != 1) */

    } /* END if no_of_brownians >1 */

    /* Return a success message */
    return NULL;

} /* END Err analyze_und_info(SrtUndInfo *und_info)


/* -------------------------------------------------------------------------- */

/* Get the index in SrtUndInfo corresponding an underlying name */
Err get_index_from_und_info(SrtUndInfo* und_info, char* name, int* index)
{
    int i;

    for (i = 0; i < und_info->no_of_underlyings; i++)
    {
        if (strcmp(und_info->und_data[i].und_name, name) == 0)
        {
            *index = i;
            return NULL;
        }
    }

    return serror("Couldn't find an underlying with such name");

} /* END Err get_index_from_und_info(...) */

/* ------------------------------------------------------------------------- */

/* Add some underlyings which are not already in the SrtUndInfo list but are
   required for the computation (dom and for IR_UND for an FX_STOCH_RATES...) */
Err add_more_underlyings_to_und_info(SrtUndInfo* und_info)
{
    int               i, j, k;
    SrtErr            err;
    SrtMdlType        mdl_type;
    SrtUndPtr         und;
    String            tempdomccy, tempforccy, und_ccy;
    String            dom_ccy_str, und_name_str;
    SrtCorrLstPtr     corrlist;
    SrtCorrLstVal*    corrval;
    SrtUnderlyingType und_type;
    SRT_Boolean       found = SRT_NO;

    /* Get the domestic currency */
    und = lookup_und(und_info->und_data[0].und_name);
    if (und == NULL)
    {
        return serror("Underlying %s not initialise", und_info->und_data[0].und_name);
    }
    dom_ccy_str = get_underlying_ccy(und);

    /* Check if we are going to need the Corr TS
    if (err = get_underlying_mdltype(und, &mdl_type))
            return err;

    if (mdl_type == FX_BETADLM || !strcmp(und->underl_lbl,"BLACK_UND"))
    {
            und_info->use_corr_ts &= SRT_NO;
    }
    else
    {
            und_info->use_corr_ts &= SRT_YES;
    }
    */
    /* Check if we are going to need the Corr TS */
    /* C. Godart: at this point we want to check that all underlyings agree
            about using the static correlation matrix */
    und_info->use_corr_ts = SRT_YES;
    for (i = 0; i < und_info->no_of_underlyings; i++)
    {
        und = lookup_und(und_info->und_data[i].und_name);
        if (und == NULL)
        {
            return serror("Underlying %s not initialise", und_info->und_data[i].und_name);
        }

        /* Check if we are going to need the Corr TS */
        if (err = get_underlying_mdltype(und, &mdl_type))
            return err;

        if (mdl_type == FX_BETADLM || !strcmp(und->underl_lbl, "BLACK_UND") ||
            !strcmp(und->underl_lbl, "SABR_UND"))
        {
            und_info->use_corr_ts &= SRT_NO;
        }
        else
        {
            und_info->use_corr_ts &= SRT_YES;
        }
    }

    /* Loop on the different underlying to see if we need to add some new ones */
    for (i = 0; i < und_info->no_of_underlyings; i++)
    {
        und = lookup_und(und_info->und_data[i].und_name);
        if (und == NULL)
        {
            return serror("Underlying %s not initialise", und_info->und_data[i].und_name);
        }

        /* Get the underlying type */
        und_type = get_underlying_type(und);

        /* Get the model type */
        if (err = get_underlying_mdltype(und, &mdl_type))
            return err;

        /* EQ_STOCH_RATES add the IR underlying */
        if ((mdl_type == EQ_STOCH_RATES) || (mdl_type == EQ_STOCH_RATES_SRVGS))
        {
            und_name_str = get_discname_from_underlying(und);
            for (j = 0; j < und_info->no_of_underlyings; j++)
            {
                if (strcmp(und_name_str, und_info->und_data[j].und_name) == 0)
                    break;
            }
            /* The domestic interest rate is not in the list: add it */
            if (j >= und_info->no_of_underlyings)
            {
                strcpy(und_info->und_data[und_info->no_of_underlyings].und_name, und_name_str);
                und_info->no_of_underlyings++;
            }

            /* If EQ_UND is P&L underlying, change NUMERAIRE index to domestic IR_UND for discount
             */
            if (i == 0)
                und_info->numeraire_index = j;
        } /* end of EQ_STOCH_RATES */

        /* FX_UND with stochastic rates : dom and for IR_UND need to be there */
        else if ((und_type == FOREX_UND) && (mdl_type == FX_LGMSV))
        {
            /* Check for the domestic interest rate model */
            und_name_str = get_domname_from_fxund(und);
            for (j = 0; j < und_info->no_of_underlyings; j++)
            {
                if (strcmp(und_name_str, und_info->und_data[j].und_name) == 0)
                    break;
            }

            /* The domestic interest rate is not in the list: add it */
            if (j >= und_info->no_of_underlyings)
            {
                strcpy(und_info->und_data[und_info->no_of_underlyings].und_name, und_name_str);
                und_info->no_of_underlyings++;
            }

            /* If FX_UND is P&L underlying, change NUMERAIRE index to domestic IR_UND for discount
             */
            if (i == 0)
                und_info->numeraire_index = j;

            /* Check for the foreign interest rate  underlying */
            und_name_str = get_forname_from_fxund(und);
            for (j = 0; j < und_info->no_of_underlyings; j++)
            {
                if (strcmp(und_name_str, und_info->und_data[j].und_name) == 0)
                    break;
            }
            /* The foreign interest rate is not in the list: add it */
            if (j >= und_info->no_of_underlyings)
            {
                strcpy(und_info->und_data[und_info->no_of_underlyings].und_name, und_name_str);
                und_info->no_of_underlyings++;
            }

        } /* END Fx stoch rates */

        /* FX_UND with stochastic rates : dom and for IR_UND need to be there */
        else if ((und_type == FOREX_UND) && (mdl_type == FX_BETADLM))
        {
            /* Check for the domestic interest rate model */
            und_name_str = get_domname_from_fxund(und);
            for (j = 0; j < und_info->no_of_underlyings; j++)
            {
                if (strcmp(und_name_str, und_info->und_data[j].und_name) == 0)
                    break;
            }

            /* The domestic interest rate is not in the list: add it */
            if (j >= und_info->no_of_underlyings)
            {
                strcpy(und_info->und_data[und_info->no_of_underlyings].und_name, und_name_str);
                und_info->no_of_underlyings++;
            }

            /* If FX_UND is P&L underlying, change NUMERAIRE index to domestic IR_UND for discount
             */
            if (i == 0)
                und_info->numeraire_index = j;

            /* Check for the foreign interest rate  underlying */
            und_name_str = get_forname_from_fxund(und);
            for (j = 0; j < und_info->no_of_underlyings; j++)
            {
                if (strcmp(und_name_str, und_info->und_data[j].und_name) == 0)
                    break;
            }
            /* The foreign interest rate is not in the list: add it */
            if (j >= und_info->no_of_underlyings)
            {
                strcpy(und_info->und_data[und_info->no_of_underlyings].und_name, und_name_str);
                und_info->no_of_underlyings++;
            }

        } /* END Fx Beta DLM stoch rates */

        /* FX_UND with stochastic rates : dom and for IR_UND need to be there */
        if ((und_type == FOREX_UND) && (mdl_type == FX_STOCH_RATES))
        {
            /* Check for the domestic interest rate model */
            und_name_str = get_domname_from_fxund(und);
            for (j = 0; j < und_info->no_of_underlyings; j++)
            {
                if (strcmp(und_name_str, und_info->und_data[j].und_name) == 0)
                    break;
            }

            /* The domestic interest rate is not in the list: add it */
            if (j >= und_info->no_of_underlyings)
            {
                strcpy(und_info->und_data[und_info->no_of_underlyings].und_name, und_name_str);
                und_info->no_of_underlyings++;
            }

            /* If FX_UND is P&L underlying, change NUMERAIRE index to domestic IR_UND for discount
             */
            if (i == 0)
                und_info->numeraire_index = j;

            /* Check for the foreign interest rate  underlying */
            und_name_str = get_forname_from_fxund(und);
            for (j = 0; j < und_info->no_of_underlyings; j++)
            {
                if (strcmp(und_name_str, und_info->und_data[j].und_name) == 0)
                    break;
            }
            /* The foreign interest rate is not in the list: add it */
            if (j >= und_info->no_of_underlyings)
            {
                strcpy(und_info->und_data[und_info->no_of_underlyings].und_name, und_name_str);
                und_info->no_of_underlyings++;
            }

        } /* END Fx stoch rates */

        else

            /* Check for different currencies (QUANTO) : need the FX to be defined */
            if (((und_type == FOREX_UND) || (und_type == INTEREST_RATE_UND) ||
                 (und_type == EQUITY_UND)) &&
                und_info->use_corr_ts == SRT_YES)
        {
            /* Get the currency of the underlying to check for QUANTO adjustment */
            und_ccy = get_underlying_ccy(und);

            /* Check if the und currency is the same; if not: quanto adjustment */
            if (strcmp(und_ccy, dom_ccy_str) != 0)
            {
                found = SRT_NO;

                /* Extract the global correlation list ( for the underlying names) */
                corrlist = srt_f_GetTheCorrelationList();
                if (!corrlist)
                    return serror("Correaltion list improperly initialised...");
                if (!corrlist->head)
                    return serror("Correaltion list improperly initialised...");
                corrval = (SrtCorrLstVal*)corrlist->head->element->val.pval;

                /* Go through the correlation list to find an FX und with both currencies */
                for (k = 0; k < corrval->nund; k++)
                {
                    /* Get the underlying associated to the k-th name in the correlation list */
                    und = lookup_und(corrval->und_names[k]);
                    if (und == NULL)
                    {
                        return serror("Underlying %s not initialise", corrval->und_names[k]);
                    }

                    /* Get the underlying type */
                    und_type = get_underlying_type(und);

                    /* If it is an FX underlying : this might be of interest */
                    if (und_type == FOREX_UND)
                    {
                        /* Get the domestic and foreign currencies from the FX underlying */
                        err = get_fx_underlying_currencies(und, &tempdomccy, &tempforccy);

                        /* Look for the FX with the same or opposite convention */
                        if (((strcmp(tempdomccy, dom_ccy_str) == 0) &&
                             (strcmp(tempforccy, und_ccy) == 0)) ||
                            ((strcmp(tempforccy, dom_ccy_str) == 0) &&
                             (strcmp(tempdomccy, und_ccy) == 0)))
                        {
                            found = SRT_YES;

                            /* Check if this underlying belongs to the Grfn underlying list */
                            for (j = 0; j < und_info->no_of_underlyings; j++)
                            {
                                if (strcmp(und_info->und_data[j].und_name, corrval->und_names[k]) ==
                                    0)
                                    break;
                            }

                            /* OVE: the FX underlying should not be added to the list of underlyings
                               to be diffused if it is not there
                                                                            /* It is not in the
                               list, we add it ( it will have to be diffused...) if (j >=
                               und_info->no_of_underlyings)
                                                                            {
                                                                                    strcpy(und_info->und_data[und_info->no_of_underlyings].und_name,
                               corrval->und_names[k]); und_info->no_of_underlyings++;
                                                                            }
                            */

                        } /* End of a forex with the right currency exists */

                    } /* End of a forex exists in the correlation matrix */

                } /* End of loop on correlation matrix underlyings */

                /* We have not found a forex underlying for the quanto adjustment */
                if (found == SRT_NO)
                {
                    return serror("Need a FOREX underlying [%s/%s]", dom_ccy_str, und_ccy);
                }

            } /* End of currency different from domestic one (QUANTO) */

        } /* End of different underlying */

    } /* End of the loop on all underlying in und_info */

    /* Return a success message */
    return NULL;

} /* END Err add_more_underlyings_to_und_info(...) */

/* ---------------------------------------------------------------------------------- */
