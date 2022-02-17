/* -------------------------------------------------------------------------------
   FILENAME:      srt_f_tre_core.c
   
   MAIN FUNCTION: TreeCore

   PURPOSE:       The main funcion to prcieusing a tree: sends to the right tree
   ------------------------------------------------------------------------------- */

#include "srt_h_all.h"
#include "srt_h_vegatreche2dr.h"
#include "srt_h_cheybetatree.h"
#include "srt_h_cheybeta2ftree.h"
/*
#include "srt_h_newbeta2ftree.h" 
#include "srt_h_chey1equ_tree.h"
#include "srt_h_ammc_cheybeta.h"
#include "srt_h_ammc_lgm2f.h"
*/
#include "srt_h_lgm2ftree_minim.h"
#include "srt_h_lsmcsv.h"

Err TreeCore
	(
	SrtUndPtr         und, 
	SrtGrfnParam     *grfnparam, 
	SrtStpPtr         step, 
	GrfnDeal         *gd,
	EvalEventFct      evalcf, 
	void             *iolist,
	SrtUndInfo       *und_info
	)
{
Err 			err = NULL;
SrtIOStruct* 	io_request = iolist;
String          und_name;
SrtMdlDim 		mdl_dim;
SrtMdlType		mdl_type;

	if (und_info->no_of_underlyings == 1 )
	{
	/* Check underlying is the domestic one (underlying number 0 in this case) */
		und_name = get_underlying_name(und);
		if (strcmp(und_name , und_info->und_data[0].und_name) )
		{
			return serror("Market in TreeCore is != from domestic one !");
		}
		err = get_underlying_mdltype(und,&mdl_type);
		if (err)
		{
			return err;
		}
		err = get_underlying_mdldim(und,&mdl_dim);
		if (err)
		{
			return err;
		}
		
	/* -> 1 FACTOR and 2 FACTOR CHEYETTE TREES <- */

		if ( mdl_type == CHEY ) 
		{
			switch ( mdl_dim )
			{
			case ONE_FAC:
				err = srt_f_vegatreche2dr(und,grfnparam,step,gd,evalcf,io_request,und_info);
				if (err)
				{
					return err;
				}
				break;
			case TWO_FAC:
				err = srt_f_cheybeta2dtree( und,
					                        grfnparam,
										    step,
										    gd,
										    evalcf,
										    io_request,
										    und_info);

				if (err)
				{
					return (err) ;
				}
				break;
			default:
	 			return serror("Sorry, you will have to wait to have this tree");
			break;
			}
		} /* END  if (mdl_type == CHEY ) */

	/* -> 1 FACTOR and 2 FACTOR CHEY BETA TREE <- */

		else 
		if(mdl_type == CHEY_BETA)
		{
			switch ( mdl_dim )
			{
			case ONE_FAC:				
				if ((grfnparam->lsm == SRT_BACK) || (grfnparam->lsm == SRT_FORBACK))
				{					                        
					/*
					err = srt_f_ammccheybeta (grfnparam, step, gd, evalcf, io_request, und_info); 
					if (err) return err;
					*/
				}
				else
				{
					err = srt_f_CheyBetaTree (und, grfnparam, step, gd, evalcf, io_request, und_info); 
					if (err)   return err;
				}
				break;
			case TWO_FAC:
				if ((grfnparam->lsm == SRT_BACK) || (grfnparam->lsm == SRT_FORBACK))
				{					                        
					/*
					err = srt_f_newbeta2ftree( und, grfnparam, step, gd, evalcf, io_request, und_info); 
					*/
				}
				else
				{
					err = srt_f_cheybeta2dtree( und, grfnparam, step, gd, evalcf, io_request, und_info); 
				}

				if (err)
				{
					return (err) ;
				}
				break;
			default:
	 			return serror("Sorry, you will have to wait to have this tree");
				break;
			}
		} /* END if mdl_type == CHEY_BETA */   

	/* -------------- 1 FACTOR and 2 FACTOR LGM TREES -------------------- */

		else 
		if( mdl_type == LGM)
		{
			switch ( mdl_dim )
			{
			case ONE_FAC:
				err = srt_f_vegashifttrelgm1d(und,grfnparam,step,gd,evalcf,io_request,und_info);
				if (err) 
				{
					return err;
				}
		/*  we try to get the prices for the shift in the yield curve            
			err = srt_f_yctrelgm1d(m,grfnparam,step,gd,evalcf,iolist); 
			not yet implemented */
				break;
			case TWO_FAC:
				if ((grfnparam->lsm == SRT_BACK) || (grfnparam->lsm == SRT_FORBACK))
				{					                        
					/*
					err = srt_f_ammclgm2f(grfnparam, step, gd, evalcf, io_request, und_info); 
					if (err) return err;
					*/
				}
				else
				{
					if (grfnparam->minim == SRT_YES)
					{						
						err = srt_f_lgm2dtree_proj(und,grfnparam,step,gd,evalcf,io_request,und_info);
						if (err) return err;						
					}
					else 
					{
						err = srt_f_lgm2dtree(und,grfnparam,step,gd,evalcf,io_request,und_info);
						if (err) return err;
					}
				}
				break;
			default:
	 			return serror("Sorry, the LGM model required does not exist");
			break;
			}
		} /* END if mdl_type == LGM */   

		else 
		if ( mdl_type == BDT)
		{
			switch (mdl_dim)
			{
				case ONE_FAC:
					err = srt_f_trebdt1f(und,grfnparam,step,gd,evalcf,io_request,und_info);
					if (err) 
					{
						return err;
					}
					break;
				case TWO_FAC:
				/*	Modified AS: BDT 2F Tree is unavailable for the moment, The file is out of the project...  */
				/* err = srt_f_trebdt2f(und,mdl,step,gd,evalcf,io_request,und_info); 
				if (err) 
				{
					return err;
				} */
					return serror ("BDT 2F Tree is not available");
				break;
			}
		} /* END if mdl_type == BDT */   

		else 
	
	/* -------------------- ONE FACTOR BLACK SCHOLES LIKE TREE ------------------------- */
		if ((mdl_type == BLACK_SCHOLES))
		{
			err = srt_f_logtree(und,grfnparam,step,gd,evalcf,io_request,und_info);
			if (err) 
			{
				return err;
			}
		} /* END if mdl_type == BLACK_SCHOLES */

		else 
		if (mdl_type == NORMAL_BS)
		{
			return serror ("Normal Tree not yet implemented");
		}

        /*------------------ CHEY_STOCH_VOL ---------------*/
		else
		if( mdl_type == CHEY_BETA_STOCH_VOL)
		{
			switch ( mdl_dim )
			{
			case ONE_FAC:
				err = srt_f_csv_lsm(grfnparam, step, gd, evalcf, io_request, und_info); 
			    if (err) return err;
				break;
			default:
	 			return serror("Sorry, the  model required does not exist");
			break;
			}
		} /* END if mdl_type == CHEY_STOCH_VOL */   

		else 
			return serror("Model requested has not been implemented");

	} /* END if no_of_underlyings == 1 */    
	else
	/* -------------------- TWO UNDERLYINGS TREE ------------------------ */
	
	if (und_info->no_of_underlyings == 2 )
	{
		SrtMdlDim 		mdl_dim1, mdl_dim2;
		SrtMdlType		mdl_type1, mdl_type2;
		SrtUndPtr       und1, und2;

		und1 = lookup_und(und_info->und_data[0].und_name);
		err = get_underlying_mdldim(und1,&mdl_dim1 );
		if (err)  return err;
		err = get_underlying_mdltype(und1,&mdl_type1 );
		if (err)  return err;

		und2 = lookup_und(und_info->und_data[1].und_name);
		err = get_underlying_mdldim(und2,&mdl_dim2 );
		if (err)  return err;
		err = get_underlying_mdltype(und2,&mdl_type2 );
		if (err)  return err;

		if ( (mdl_dim1 == ONE_FAC) && (mdl_dim2 == ONE_FAC) && (mdl_type2 != EQ_STOCH_RATES_SRVGS ))
		{
			if ( ((mdl_type1 == LGM) && ((mdl_type2 == LGM) || (mdl_type2 == BLACK_SCHOLES)
											|| (mdl_type2 == EQ_STOCH_RATES)))
				|| ((mdl_type1 == BLACK_SCHOLES) && ( (mdl_type2 == LGM) || (mdl_type2 == EQ_STOCH_RATES) )) )
			{

				err = srt_f_twoundtree(grfnparam, step, gd, evalcf, io_request,und_info);
				if (err)  return err;
			}
			else if ( (mdl_type1 == CHEY_BETA) && (mdl_type2 == BLACK_SCHOLES) )
			{
				/*
				err = srt_f_cheybeta1f1etree(grfnparam, step, gd, evalcf, io_request,und_info);
				*/
				if (err)  return err;
			}
		}
		else
		{
				return serror("Cannot build a two factor tree with these models");
		}
	}
	else
	{
		return serror("Cannot build a tree with more than two underlyings");
	}

	return err;
}               

