/* ====================================================================================
   FILENAME:  srt_h_mc_evolve.h

   PURPOSE:   description of all the functions used to evolve IR or Stock/FX models
              in a Monte-Carlo type discretisation

   DESCRIPTION: all the functions have the same type of name, and the same inputs
   ==================================================================================== */

Err monte_CHE_1f_Euler_evolve(
				double 		*rndm_vec,   
				SrtStpPtr 	top, 
				SrtSample 	*sam, 
				int 		index);     
/*----------------------------------------------------------------------------*/

Err monte_CHE_1f_Milshtein_evolve(
				double      *rndm_vec,  
				SrtStpPtr   top, 
				SrtSample   *sam,
				int         index);     
/*----------------------------------------------------------------------------*/

Err monte_CHE_2f_Euler_evolve(
				double      **rndm_mat,   
				SrtStpPtr   top,  
				SrtSample   *sam, 
				int         index ) ;    

/*----------------------------------------------------------------------------*/

Err monte_LGM_1f_Jumping_evolve(
				double      *rndm_vec,  
				SrtStpPtr   top, 
				SrtSample   *sam,
				int         index) ;    
/*----------------------------------------------------------------------------*/

Err monte_LGM_1f_Euler_evolve(
				double      *rndm_vec,  
				SrtStpPtr   top, 
				SrtSample   *sam,
				int         index );     
/*----------------------------------------------------------------------------*/

Err monte_NEWLGM_1f_Euler_evolve(
				double      *rndm_vec,  
				SrtStpPtr   top, 
				SrtSample   *sam,
				int         index ); 

/*----------------------------------------------------------------------------*/

Err monte_vasicek_1f_euler_evolve(double      *rndm_vec,  
								  SrtStpPtr   top, 
  	                			  SrtSample   *sam,
								  int         und_index);

/*----------------------------------------------------------------------------*/

Err monte_NEWCHEYBETA_1f_Euler_evolve(
				double      *rndm_vec,  
				SrtStpPtr   top, 
				SrtSample   *sam,
				int         index ); 
    
/*----------------------------------------------------------------------------*/

Err monte_LGM_2f_evolve(
				double      **rndm_mat,  
				SrtStpPtr   top, 
				SrtSample   *sam,
				int         index) ;    

/*  ------------------------------------------------------------------------------ */
Err monte_BETAETA_1f_evolve(
				double      *rndm_vec,  
				SrtStpPtr   top, 
				SrtSample   *sam,
				int         index)     ;
/*----------------------------------------------------------------------------*/

Err monte_LGM_1f_STOCHVOL_evolve(
				double       **rndm_mat,  
				SrtStpPtr    top, 
				SrtSample    *sam,
				double       *Pt_by_Zt,
				int          index ) ; 
/*----------------------------------------------------------------------------*/
Err monte_CHE_1f_STOCHVOL_evolve(
				double 		**rndm_mat,   
				SrtStpPtr 	top, 
				SrtSample 	*sam, 
				double 		*Pt_by_Zt,
				int 		index ) ;

/*----------------------------------------------------------------------------*/
Err monte_CHEYBETA_1f_STOCHVOL_evolve(
				double 		**rndm_mat,   
				SrtStpPtr 	top, 
				SrtSample 	*sam, 
				double 		*Pt_by_Zt,
				int 		index) ;
/*----------------------------------------------------------------------------*/

Err monte_BLACKSCHOLES_evolve(
				double      *rndm_vec, 
				SrtStpPtr   top, 
				SrtSample   *sam,
				int         index ) ;    

/*----------------------------------------------------------------------------*/

Err monte_NORMALBS_evolve(
				double      *rndm_vec, 
				SrtStpPtr   top, 
				SrtSample   *sam,
				int         index ) ;    

/*----------------------------------------------------------------------------*/

Err monte_FX_STOCH_RATES_evolve(double       *rndm_vec,
								SrtStpPtr    top,
								SrtSample    *sam,
								int          index,
								int			 dom_index,
								int			 for_index);

Err monte_EQ_STOCH_RATES_evolve(double       *rndm_vec,
								SrtStpPtr    top,
								SrtSample    *sam,
								int          index,
								int			 dom_index);

Err monte_EQ_STOCH_RATES_SRVGS(double		*rndm_vec_spot,
							   double		*rndm_vec_vol,
							   SrtStpPtr		top,
							   SrtSample     *sam,
							   int			und_index,
							   int			dom_index);

Err monte_FX_STOCH_RATES_Jumping_evolve(
						double  *rndm_vect,	  
						SrtStpPtr top,
						SrtSample *sam,
						int und_index,
						int dom_index,
						int for_index);



/*----------------------------------------------------------------------------*/
/* set sam[t].numeraire at each time step to 1 for a deterministic 
   underlying.  */
Err monte_DETERMINISTIC_evolve(
				SrtStpPtr   top, 
				SrtSample   *sam);


/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Err monte_CHEYBETA_1f_Euler_evolve(
					double 		*rndm_vec,   
					SrtStpPtr 	top, 
					SrtSample 	*sam, 
					int 		index);

/*----------------------------------------------------------------------------*/
Err monte_CHEYBETA_2f_Euler_evolve(
				double       **rndm_mat,   
				SrtStpPtr    top, 
				SrtSample    *sam, 
				int          index);

/*----------------------------------------------------------------------------*/
Err monte_MIXEDBETA_2f_Euler_evolve(
				double       **rndm_mat,   
				SrtStpPtr    top, 
				SrtSample    *sam, 
				int          index);

/*----------------------------------------------------------------------------*/
