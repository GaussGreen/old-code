/*******************************************************************************
**                                                                          
**		Copyright (c) 1993 PARIBAS Capital Markets Group             
**                                                                          
********************************************************************************
**									    
**	SYSTEM:		SRT	SORT, Fixed Income 2020 Addins	            
**      SUB_SYSTEM:  	OPT	Option Tools				    
**	                                                                    
**	MODULE NAME:	OPTION_TOOLS					    
**    									    
**	PURPOSE:	Standard 2020 Functions used in Option Tools	    
**                                                                          
**	AUTHORS:	Guillaume AMBLARD, Remy KLAMMERS		    
**									    
**	DATE:		8th October, 1992				    
**                                                                          
**	VERSION:        XX                                                  
**	                  						    
**	DESCRIPTION:	XX                                                  
**                                     					    
**	FUNCTIONS USED:	XXX_X_XXXXXXXXX                                     
**		Must include all imported function call made by the module  
**                     							    
**	PARAMETERS:  	<not applicable>                                    
**	        							    
**	RETURNS:                                                            
**	                                                                    
**	DATA ACCESSED:	                                                    
**		Tables, files, and global variables accessed. 		    
**                                                                          
********************************************************************************
**			Amendment History				    
********************************************************************************
**	AMENDED BY:	Ken N LINTON	                                    
**	                                                                    
**	DATE:	  	18th March, 1994                                    
**	                                                                    
**	VERSION:  	<not applicable>                            	    
**	                                                                    
**	REASON:         Restructuring for re-use             		    
**	                                                                    
**	REQUEST NO:  	<not applicable>	                            
**                                                                          
**	DESCRIPTION: 	<not applicable>	                            
*******************************************************************************
**	AMENDED BY:	Julia Matsumoto	                                    
**	DATE:	  	18th July, 1994                                    
**	REASON:         Life Options              		    
**	DESCRIPTION: 	<not applicable>	                            
*******************************************************************************
**	AMENDED BY:	Julia Matsumoto	                                    
**	DATE:	  	18th April, 1995                                    
**	REASON:         Add Proba_tools, trivariate function  		    
**	DESCRIPTION: 	<not applicable>	       
******************************************************************************/
                                                                            


/* ==========================================================================
   declarations for all functions
   ==========================================================================*/

/* ========================================================================== 

				PART 1
			STANDARD OPTIONS

  ========================================================================== */


int    	black_scholes_addin(
			int 	argc, 
			UF_ARG 	*argv, 
			UF_ARG 	*ret
			); 

int     implied_vol_bs_addin(
			int 	argc, 
			UF_ARG 	*argv, 
			UF_ARG 	*ret);

int     implied_strike_bs_addin(
			int 	argc, 
			UF_ARG 	*argv, 
			UF_ARG 	*ret);

int     bs_d1_d2_addin(
			int 	argc, 
			UF_ARG 	*argv, 
			UF_ARG 	*ret);

int     bs_normal_addin(
			int 	argc, 
			UF_ARG 	*argv, 
			UF_ARG 	*ret
			);
 
int 	amer_binomial_addin(
			int 	argc, 
			UF_ARG 	*argv, 
			UF_ARG 	*ret);

int 	impvol_amer_addin(
			int 	argc, 
			UF_ARG 	*argv, 
			UF_ARG 	*ret);


#define STANDARD_OPT_UFDESC \
{\
	"black_scholes",\
	"DDDDDWVWVW",\
	"D",\
	"@black_scholes( fwd, strike, vol, mat, disc, c/p, [grk], [lgn_nrm])",\
	black_scholes_addin\
}\
,\
{\
	"implied_vol_bs",\
	"DDDDDWVW",\
	"D",\
	"@implied_vol_bs( premium, fwd, strike, mat, disc, c/p,[lgn_nrm])",\
	implied_vol_bs_addin\
},\
{\
	"implied_strike_bs",\
	"DDDDDWVW",\
	"D",\
	"@implied_strike_bs( premium, fwd, vol, mat, disc, c/p, [lgn_nrm])",\
	implied_strike_bs_addin\
},\
{\
	"bs_d1_d2",\
	"DDDDI",\
	"D",\
	"@bs_d1_d2( fwd, strike, vol, mat, info)",\
	bs_d1_d2_addin\
},\
{\
  	"amer",\
	"DDDDDDWVIVS",\
	"D",\
	"@amer( fwd, spot, strike, vol, mat, disc, c/p, [nsteps], [grk])",\
	amer_binomial_addin\
},\
{\
	"impvol_amer",\
	"DDDDDDWVI",\
	"D",\
	"@impvol_amer( premium, fwd, spot, strike, mat, disc, c/p, [nsteps])",\
	impvol_amer_addin\
}


/* ========================================================================== 

				PART 2
			   BARRIER OPTIONS

  ========================================================================== */
    
int     extinguish_addin(
			int 	argc, 
			UF_ARG 	*argv, 
			UF_ARG 	*ret); 

int     part_extinguish_addin(
			int 	argc, 
			UF_ARG 	*argv, 
			UF_ARG 	*ret); 


int    	fwd_extinguish_addin(
			int 	argc, 
			UF_ARG 	*argv, 
			UF_ARG 	*ret);


int     lightable_addin(
			int 	argc, 
			UF_ARG 	*argv, 
			UF_ARG 	*ret);

int     part_lightable_addin(
			int 	argc, 
			UF_ARG 	*argv, 
			UF_ARG 	*ret); 


int    	fwd_lightable_addin(
			int 	argc, 
			UF_ARG 	*argv, 
			UF_ARG 	*ret);


int     dbarrier_gauss_addin(
			int 	argc, 
			UF_ARG *argv, 
			UF_ARG *ret
			);

int    	part_dbarrier_addin(
			int 	argc, 
			UF_ARG 	*argv, 
			UF_ARG 	*ret);


#define EXTING_OPT_UFDESC \
{\
  	"extinguish",\
	"DDDDDDDWWVS",\
	"D",\
	"@extinguish( fwd, spot, strike, barrier, vol, mat, disc, c/p, d/u, [grk])",\
	extinguish_addin\
},\
{\
  	"part_extinguish",\
	"DDDDDDDDDDWWVS",\
	"D",\
	"@part_extinguish(fwd1,fwd2,spot,strike,barrier,vol1,vol2,mat1,mat2,disc,c/p,d/u,[grk])",\
	part_extinguish_addin\
},\
{\
  	"fwd_extinguish",\
	"DDDDDDDDDWWVS",\
	"D",\
	"@fwd_extinguish(fwd1,fwd2,strike,barrier,vol1,vol2,mat1,mat2,disc,c/p,d/u,[grk])",\
	fwd_extinguish_addin\
},\
{\
  	"lightable",\
	"DDDDDDDWWVS",\
	"D",\
	"@lightable( fwd, spot, strike, barrier, vol, mat, disc, c/p, d/u, [grk])",\
	lightable_addin\
},\
{\
	"part_lightable",\
	"DDDDDDDDDWWVS",\
	"D",\
	"@part_lightable(fwd1,fwd2,spot,strike,barrier,vol1,vol2,mat1, mat2, df, c/p, d_u, [grk])",\
	part_lightable_addin\
},\
{\
  	"fwd_lightable",\
	"DDDDDDDDDWWVS",\
	"D",\
	"@fwd_lightable(fwd1,fwd2,strike,barrier,vol1,vol2,mat1,mat2,disc,c/p,d/u,[grk])",\
	fwd_lightable_addin\
},\
{\
	"dbarrier",\
	"DDDDDDDDWVIVS",\
	"D",\
	"@dbarrier(fwd,spot,strike,bar_do,bar_up,vol,mat,disc,cp,[nterms],[grk])",\
        dbarrier_gauss_addin\
},\
{\
	"part_dbarrier",\
	"DDDDDDDDDDDWVIVS",\
	"D",\
	"@part_dbarrier(fwd1,fwd2,spot,strike,bar_do,bar_up,vol1,vol2,mat1,mat2,df,cp,[nterms],[grk])",\
        part_dbarrier_addin\
}


/* ========================================================================== 

				PART 3
			 TWO UNDERLYING OPTIONS

  ========================================================================== */

int 	spread_numer_addin(
			int 	argc, 
			UF_ARG 	*argv, 
			UF_ARG 	*ret);

int     two_factor_tree_amer_addin(
			int 	argc, 
			UF_ARG 	*argv, 
			UF_ARG 	*ret); 

int     nono_barrier_addin(
			int 	argc, 
			UF_ARG 	*argv, 
			UF_ARG 	*ret); 

int 	scud_out_addin(
			int 	argc, 
			UF_ARG 	*argv, 
			UF_ARG 	*ret); 

int 	scud_in_addin(
			int 	argc, 
			UF_ARG 	*argv, 
			UF_ARG 	*ret); 
   
int 	multiple_addin(
			int 	argc, 
			UF_ARG 	*argv, 
			UF_ARG 	*ret);



#define TWO_UND_OPT_UFDESC \
{\
	"spread_numer",\
	"DDDDDDDDDDWVW",\
	"D",\
	"@spread_numer( fwdx, fwdy, K, a, b, sigx, sigy, rho, mat, disc, cp, [grk])",\
	spread_numer_addin\
},\
{\
	"two_factor_tree_american",\
	"DDDDDDDDDDDDWIVS",\
	"D",\
	"@two_factor_tree_american(fwd1,fwd2,spot1,spot2,strike,a_inp,b_inp,vol1,vol2,rho,mat,disc,c/p,steps,[grk])",\
	two_factor_tree_amer_addin\
},\
{\
	"scud_out",\
	"DDDDDDDDDDDWWVS",\
	"D",\
	"@scud_out(fwdx,fwdy,spoty,strike,barrier,sigx,sigy,rho,mat1,mat2,disc,c/p,scud_type,[grk])",\
	scud_out_addin\
},\
{\
	"scud_in",\
	"DDDDDDDDDDDWWVS",\
	"D",\
	"@scud_in(fwdx,fwdy,spoty,strike,barrier,sigx,sigy,rho,mat1,mat2,disc,c/p,scud_type,[grk])",\
	scud_in_addin\
},\
{\
	"nono_barrier",\
	"DDDDDDDDDDDDDDIWVS",\
	"D",\
	"@nono_barrier(fx,fy,sx,sy,K,Bup_x,Bdwn_x,Bup_y,Bdwn_y,sigx,sigy,rho,mat,df,steps,c/p,[grk])",\
	nono_barrier_addin,\
},\
{\
	"multiple_opt",\
	"DDDDDDDDDWWVS",\
	"D",\
	"@multiple_opt(fwdX,strX,fwdY,strY,volX,volY,corr,mat,disc,c/p_X,c/p_Y, [grk])",\
	multiple_addin\
}


/* ========================================================================== 

				PART 4
			PATH DEPENDENT OPTIONS

  ========================================================================== */

int     ratchet_addin(	
			int 	argc, 
			UF_ARG 	*argv, 
			UF_ARG 	*ret); 

int     inc_strike_opt_ext_smth_addin(
			int 	argc, 
			UF_ARG 	*argv, 
			UF_ARG 	*ret);

int     ladder_addin(
			int argc, 
			UF_ARG *argv, 
			UF_ARG *ret); 

int     one_touch_addin(
			int argc, 
			UF_ARG *argv, 
			UF_ARG *ret);
 
int     look_back_addin(
			int argc, 
			UF_ARG *argv, 
			UF_ARG *ret);             

int		partial_look_back_addin(
			int argc, 
			UF_ARG *argv, 
			UF_ARG *ret);

int     extrema_addin(
			int argc, 
			UF_ARG *argv, 
			UF_ARG *ret);             

int     life_opt_addin(
			int 	argc, 
			UF_ARG *argv, 
			UF_ARG *ret
			);

int     fwd_life_opt_addin(
			int 	argc, 
			UF_ARG *argv, 
			UF_ARG *ret
			);

int 	life_out_addin(
			int 	argc, 
			UF_ARG 	*argv, 
			UF_ARG 	*ret);


#define PATH_DEP_OPT_UFDESC \
{\
	"ratchet",\
	"DDDDDDDDDWGGVS",\
	"D",\
	"@ratchet( fwd1, fwd2, spot, hist_ext, vol1, vol2, mat1, mat2, disc, c/p, rungs[], strikes[], [grk])",\
	ratchet_addin\
},\
{\
  	"iko_ext",\
	"DDDDDDDDDDWGGVIVS",\
	"D",\
	"@iko_ext( fwd1, fwd2, spot, barrier, hist_ext, vol1, vol2, mat1, mat2, disc, c/p,  rungs[], strikes[], [n_steps], [greek])",\
	inc_strike_opt_ext_smth_addin\
},\
{\
	"ladder",\
	"DDDDDDDWWGVS",\
	"D",\
	"@ladder( fwd, spot, strike, hist_ext, vol, mat, disc, c/p, g/b, rungs[], [grk])",\
	ladder_addin\
},\
{\
	"one_touch",\
	"DDDDDDDWVS",\
	"D",\
	"@one_touch( fwd, spot, strike, barrier, vol, mat, disc, c/p, [grk])",\
	one_touch_addin\
},\
{\
	"look_back",\
	"DDDDDDWVS",\
	"D",\
	"@look_back( fwd, spot, hist_extr, vol, mat, disc, c/p, [grk])",\
	look_back_addin\
},\
{\
	"partial_look_back",\
	"DDDDDDDDDWVS",\
	"D",\
	"@partial_look_back( fwd1, fwd2, spot, hist_extr, vol1, vol2, mat1, mat2, disc, c/p, [grk])",\
	partial_look_back_addin\
},\
{\
	"extrema",\
	"DDDDDDDWVS",\
	"D",\
	"@extrema( fwd, spot, strike, hist_extr, vol, mat, disc, c/p, [grk])",\
	extrema_addin\
},\
{\
	"life",\
        "DDDDDDWVS",\
        "D",\
        "@life( fwd, spot, barrier, vol, mat, disc, d/u, [grk])",\
        life_opt_addin\
},\
{\
	"fwd_life",\
        "DDDDDDDDDWVS",\
        "D",\
        "@fwd_life( fwd1, fwd2, spot, barrier, vol1, vol2, mat1, mat2, disc1, disc2, d/u, [grk])",\
        fwd_life_opt_addin\
},\
{\
	"life_out",\
	"DDDDDDDVIVS",\
	"D",\
	"@life_out( fwd, spot, b_life, b_out, vol, mat, disc, [nterms], [grk])",\
	life_out_addin,\
}

/* ===========================================================================

				PART 5
			 DIGITAL OPTIONS

  ============================================================================*/

int	euro_digital_addin(	
			int 	argc, 
			UF_ARG 	*argv, 
			UF_ARG 	*ret
			);

int	amer_digital_addin(		
			int 	argc, 
			UF_ARG 	*argv, 
			UF_ARG 	*ret
			);           

int	fwd_amer_digital_addin(		
			int 	argc, 
			UF_ARG 	*argv, 
			UF_ARG 	*ret
			);           

int	euro_slalom_digit_addin(	
			int 	argc, 
			UF_ARG 	*argv, 
			UF_ARG 	*ret
			);
                                      
int	am_slalom_digit_addin(	
			int 	argc, 
			UF_ARG 	*argv, 
			UF_ARG 	*ret
			);
                                      
int     slalom_addin(
			int 	argc, 
			UF_ARG *argv, 
			UF_ARG *ret
			);

int	touch_pay_amedig_addin(	
			int 	argc, 
			UF_ARG 	*argv, 
			UF_ARG 	*ret
			);

int	first_second_addin(	
			int 	argc, 
			UF_ARG 	*argv, 
			UF_ARG 	*ret
			);



#define DIGITAL_OPT_UFDESC \
{\
	"euro_digital",\
	"DDDDDSVS",\
	"D",\
	"@euro_digital( fwd, barrier, vol, mat, disc, b/a, [grk])",\
	euro_digital_addin\
},\
{\
	"amer_digital",\
	"DDDDDDIVIVS",\
	"D",\
	"@amer_digital( fwd, spot, barrier, vol, mat, disc, b/a, [min/max],[grk])",\
	amer_digital_addin\
},\
{\
	"fwd_am_digit",\
	"DDDDDDDDII",\
	"D",\
	"@fwd_am_digit( fwd1, fwd2, barrier, vol1, vol2, mat1, mat2, disc, b/a, [min/max], [grk])",\
	fwd_amer_digital_addin\
},\
{\
	"touch_pay_amedigit",\
	"DDDDDDDDDSZS",\
	"D",\
	"@touch_pay_amedigit( fwd1, fwd2, barrier, vol1, vol2, mat1, mat2, disc1, disc2, b/a, [grk])",\
	touch_pay_amedig_addin\
},\
{\
	"first_second_digital",\
	"DDDDDDDII",\
	"D",\
	"@first_second_digital( fwd, spot, first_bar, secnd_bar, vol, mat, disc, [nterms], [grk])",\
	first_second_addin\
},\
{\
  	"euro_slalom_digit",\
	"DDDDDDDDSVS",\
	"D",\
	"@euro_slalom_digit( fwd1, fwd2, barrier, vol1, vol2, mat1, mat2, disc, b/a, [grk])",\
	euro_slalom_digit_addin\
},\
{\
	"amer_slalom_digit",\
	"DDDDDDDDSVS",\
	"D",\
	"@amer_slalom_digit( fwd1, fwd2, strike, vol1, vol2, mat1, mat2, disc, b/a, [grk])",\
	am_slalom_digit_addin\
},\
{\
	"slalom",\
	"DDDDDIGGIVI",\
	"D",\
	"@slalom( fwd, spot, bar, vol, disc, today, nb_dates, date[], up_do[], step_num, [dbg])",\
        slalom_addin\
}

/* ========================================================================== 

				PART 6
		  DISTRIBUTIONS and PROBABILITIES

  ========================================================================== */

int     normal_addin(
			int 	argc, 
			UF_ARG 	*argv, 
			UF_ARG *ret);    

int     bivariate_addin(
			int 	argc, 
			UF_ARG 	*argv, 
			UF_ARG 	*ret);

int     trivariate_addin(
			int 	argc, 
			UF_ARG 	*argv, 
			UF_ARG 	*ret);

int 	stay_inside_addin(
			int 	argc, 
			UF_ARG 	*argv, 
			UF_ARG 	*ret
			);

int 	stay_inside_max_addin(
			int 	argc, 
			UF_ARG 	*argv, 
			UF_ARG 	*ret
			);

int 	stay_inside_fwd_addin(
			int 	argc, 
			UF_ARG 	*argv, 
			UF_ARG 	*ret
			);
  
int 	stay_inside_fwd_max_addin(
			int 	argc, 
			UF_ARG 	*argv, 
			UF_ARG 	*ret
			);

int 	hit_one_barrier_addin(
			int 	argc, 
			UF_ARG 	*argv, 
			UF_ARG 	*ret
			);                                    



#define PROBA_OPT_UFDESC \
{\
 	"normal",\
	"D",\
	"D",\
	"@normal( d)",\
	normal_addin\
},\
{\
   	"bivariate",\
	"DDD",\
	"D",\
	"@bivariate( x, y, rho)",\
	bivariate_addin\
},\
{\
   	"trivariate",\
	"DDDDDD",\
	"D",\
	"@trivariate( a1, a2, a3, ro12, ro23, ro31)",\
	trivariate_addin\
},\
{\
	"stay_inside",\
	"DDDDDDVI",\
	"D",\
	"@stay_inside( spot, fwd, bar_low, bar_up, vol, mat, [flag, nb_terms]",\
	stay_inside_addin\
},\
{\
	"stay_inside_max",\
	"DDDDDVIVI",\
	"D",\
	"@stay_inside_max( spot, fwd, bar_gap, vol, mat, [flag, nb_terms]",\
	stay_inside_max_addin\
},\
{\
	"stay_inside_fwd",\
	"DDDDDDDDDVIVI",\
	"D",\
	"@stay_inside_fwd( spot, fwd1, fwd2, bar_low, bar_up, vol0, vol, t1, t2, [flag,n_terms]",\
	stay_inside_fwd_addin\
},\
{\
	"stay_inside_fwd_max",\
	"DDDDDDDVIVI",\
	"D",\
	"@stay_inside fwd_max( spot, fwd1, fwd2, bar_gap, vol, t1, t2, [flag, nb_terms]",\
	stay_inside_fwd_max_addin\
},\
{\
	"hit_one_barrier",\
	"DDDDDDVIVI",\
	"D",\
	"@hit_one barrier( spot, fwd, bar_yes, bar_no, vol, mat, [flag, nb_terms]",\
	hit_one_barrier_addin\
}

/* ========================================================================== 

				  PART 7
		  		DOC ADDIN

  ========================================================================== */

int	optdoc_addin(
			int argc,
			UF_ARG *argv,
			UF_ARG *ret);



/* ========================================================================= */
/* ========================================================================= 

				  PART 8
		  		MESSY ADDIN

  ========================================================================== */

int 	lightable_swap_addin(
			int argc, 
			UF_ARG *argv, 
			UF_ARG *ret);

int     two_factor_tree_addin(
			int 	argc, 
			UF_ARG 	*argv, 
			UF_ARG 	*ret); 

int 	mltspd_addin	(
			int 	argc, 
			UF_ARG 	*argv, 
			UF_ARG 	*ret);

int 	dbscud_addin(
			int 	argc, 
			UF_ARG 	*argv, 
			UF_ARG 	*ret); 

int 	triscud_addin(
			int 	argc, 
			UF_ARG 	*argv, 
			UF_ARG 	*ret);    
   
int 	thrust_addin(
			int 	argc, 
			UF_ARG 	*argv, 
			UF_ARG 	*ret);    

int     rebate_addin(
			int argc, 
			UF_ARG *argv, 
			UF_ARG *ret);


#define MESSY_OPT_UFDESC \
{\
	"lightable_swap",\
	"IDDDDDDDDIIDI",\
	"D",\
	"@lightable_swap( num_path,fwdx,sigx,fwdy,sigy,rho,spread,t, df,compd,nfp,coupon,lt_gt)",\
       	lightable_swap_addin\
},\
{\
  	"rebate",\
	"DDDDDDDWVS",\
	"D",\
	"@rebate( fwd, spot, barrier, rebate, vol, mat, disc, c/p, [grk])",\
	rebate_addin\
},\
{\
 	"two_factor_tree",\
	"DDDDDDDDDDWIVS",\
	"D",\
	"two_factor_tree(fwd1,fwd2,strike,a_inp,b_inp,vol1,vol2,rho,mat,disc,c/p,steps,[grk])",\
	two_factor_tree_addin\
},\
{\
	"multi_sprd",\
	"DDDDDDDDDDDDDDWDVS",\
	"D",\
	"@multi_sprd(fx,fy,fz,sx,sy,sz,rxy,rxz,ryz,b_y,c_z,K,mat,df, c/p,path,[grk])",\
	mltspd_addin\
},\
{\
	"dbscud",\
	"DDDDDDDDDDDDWII",\
	"D",\
	"@dbscud(fx,fy,sx,sy,K,barx,bary,sigx,sigy,rho,mat,disc,c/p,d/u_x,d/u_y)",\
	dbscud_addin\
},\
{\
	"triscud",\
	"DDDDDDDDDDDDDWIIVS",\
	"D",\
	"@triscud(fx,fy,sx,sy,K_x,BU_x,BD_x,B_y,sigx,sigy,rho,mat,df,cp_x,d/u_y,nterms,[grk])",\
	triscud_addin\
},\
{\
	"thrust",\
	"DDDDDDDDDDDDDDDDWIIVS",\
	"D",\
	"@thrust(fx,fy,fz,sy,sz,K_x,B_yz,sigx,sigy,sigz,rxy,ryz,rxz,mat,mate,df,cp,du,el,[grk])",\
	thrust_addin\
}




/* ========================================================================= */
/* ========================================================================== 
   Library info structure
   ==========================================================================*/

static	UF_INFO	lib_info =
{
   	"OPTION_TOOLS: Option Pricing Library -- Version SRT006",
	__DATE__ ,
	__TIME__
};


/* ==========================================================================
   Function descriptor table
   ==========================================================================*/

static UF_DESC desc_info[] = 
{
	{	"optdbg",
		"I",
		"S",
		"@optdbg(int)", 
		debug_addin
	},
	STANDARD_OPT_UFDESC,               
	EXTING_OPT_UFDESC,
	TWO_UND_OPT_UFDESC,
	PATH_DEP_OPT_UFDESC,
	DIGITAL_OPT_UFDESC,
	PROBA_OPT_UFDESC,
	MESSY_OPT_UFDESC,
	{
		"optdoc",			
		"WVI",		
		"S",
		"@optdoc( index, [arglist])",
	  	optdoc_addin
	},


  "","","","",0

} ;



/* =========================================================================

    	


   ========================================================================= */







