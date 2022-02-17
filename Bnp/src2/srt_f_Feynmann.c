
#include "SRT_H_ALL.H>
#include "srt_h_BrownianMotion.H>
#include "srt_h_FeynmannKac.H>
#include "math.h"

#define PDE_THETA 0.5
#define LIM_IN_STDEV	6 
#define ONE_HALF 0.5
#define MAX_TIME 0.0416666666666667

Err srt_f_feynmann_kac_pde(long			    num_mesh,
						 long			    num_time_step,
						 double             dX,
						 double				time_step,
						 double				*X,
						 double				**conv_mat, /* CONVECTION MATRIX */
						 double				**diff_mat, /* DIFFUSION MATRIX */
						 SrtPdePDEType		pde_type,
						 SrtPdeBoundaryCond	pde_bound_cond,
						 double				*pde_bound,
						 double				*term_payoff, /* FINAL MATURITY PAYOFF */
						 double				**amer_payoff, /* AMERICAN PAYOFF */

						 double				**div_mat,
					  
						 double				**U)
{
	Err			err = NULL;
	long		i,j;
	double		*diag_A = NULL, *diag_B = NULL,*diag_C = NULL,*vect_R = NULL,*vect_U = NULL,
				*PrevPrePtr = NULL,	PDE_Param[4];

	if(pde_bound_cond == DIRICHLET) /* FIRST CASE: THE BOUDARY CONDITIONS ARE DIRICHLET ONES */
	{
		/* ALLOCATIONS OF THE TRIDAG STRUCTURES */
		diag_A = dvector(1,num_mesh-1);
		diag_B = dvector(1,num_mesh-1);
		diag_C = dvector(1,num_mesh-1);
		vect_R = dvector(1,num_mesh-1);
		vect_U = dvector(1,num_mesh-1);

		PrevPrePtr = dvector(1,num_mesh-1);

		PDE_Param[0] = ONE_HALF*PDE_THETA*time_step/(dX*dX);
		PDE_Param[1] = ONE_HALF*(PDE_THETA-1.0)*time_step/(dX*dX);
		
		PDE_Param[2] = ONE_HALF*PDE_THETA*time_step/dX;
		PDE_Param[3] = ONE_HALF*(PDE_THETA-1.0)*time_step/dX;
	
		/*PAYOFF FUNCTION AT THE LAST STEP*/
		for(j = 1; j < num_mesh; j++) PrevPrePtr[j] = term_payoff[j];
		
		/* START OF THE forward ALGORITHM */
		for(i = 1; i <= num_time_step; i++) /* LOOP ON THE TIME STEP */
		{
			for(j = 1; j < num_mesh; j++) /* LOOP ON THE SPACE POINTS */
			{
				diag_A[j] = - PDE_Param[2]*conv_mat[i][j] + PDE_Param[0]*diff_mat[i][j]*diff_mat[i][j];
				diag_B[j] = -2*PDE_Param[0]*diff_mat[i][j]*diff_mat[i][j] - 1.0;
				diag_C[j] = PDE_Param[2]*conv_mat[i][j] + PDE_Param[0]*diff_mat[i][j]*diff_mat[i][j];			
			
				/* COMPUTE THE RIGHT HAND VECTOR */
				if( j == 1) /* BOUNDARY CONDITION TO APPLY */
				{
					
					vect_R[j] = 0;					
					vect_R[j] = (PDE_Param[1]*diff_mat[i-1][j]*diff_mat[i-1][j] + PDE_Param[3]*conv_mat[i-1][j])*PrevPrePtr[j+1]
								-( 1.0 + 2.0*PDE_Param[1]*diff_mat[i-1][j]*diff_mat[i-1][j])*PrevPrePtr[j]
								+ (PDE_Param[1]*diff_mat[i-1][j]*diff_mat[i-1][j]-PDE_Param[3]*conv_mat[i-1][j])*pde_bound[0]
								- diag_A[j]*pde_bound[0]
								
								+(1-PDE_THETA)*time_step*div_mat[i-1][j]
								+PDE_THETA*time_step*div_mat[i][j];
								
				
				} 
				else if( j == (num_mesh-1)) /* BOUNDARY CONDITION TO APPLY */
				{
					
					vect_R[j] = 0;
					vect_R[j] = (PDE_Param[1]*diff_mat[i-1][j]*diff_mat[i-1][j] + PDE_Param[3]*conv_mat[i-1][j])*pde_bound[1]
								-( 1.0 + 2.0*PDE_Param[1]*diff_mat[i-1][j]*diff_mat[i-1][j])*PrevPrePtr[j]
								+ (PDE_Param[1]*diff_mat[i-1][j]*diff_mat[i-1][j]-PDE_Param[3]*conv_mat[i-1][j])*PrevPrePtr[j-1]
								- diag_C[j]*pde_bound[1]

								+(1-PDE_THETA)*time_step*div_mat[i-1][j]
								+PDE_THETA*time_step*div_mat[i][j];

				}
				else /* NO BOUNDARY CONDITION TO APPLY */
				{
					
					vect_R[j] = 0;
					vect_R[j] = (PDE_Param[1]*diff_mat[i-1][j]*diff_mat[i-1][j] + PDE_Param[3]*conv_mat[i-1][j])*PrevPrePtr[j+1]
								-( 1.0 + 2.0*PDE_Param[1]*diff_mat[i-1][j]*diff_mat[i-1][j])*PrevPrePtr[j]
								+ (PDE_Param[1]*diff_mat[i-1][j]*diff_mat[i-1][j]-PDE_Param[3]*conv_mat[i-1][j])*PrevPrePtr[j-1]

								+(1-PDE_THETA)*time_step*div_mat[i-1][j]
								+PDE_THETA*time_step*div_mat[i][j];
			
				}
			
			} /* LOOP ON THE SPACE POINTS */

			err = tridag(diag_A,diag_B,diag_C,vect_R,vect_U,num_mesh-1);
			if(err) return err;

			for(j = 1; j < num_mesh ; j++)
			{
				if(pde_type == LINEAR) PrevPrePtr[j] = vect_U[j];
				else  PrevPrePtr[j] =  max(vect_U[j],amer_payoff[i][j]);
				
			}

		}/* LOOP ON THE TIME STEP */
	
	} /* END OF FIRST CASE */

	for(j = 1; j < num_mesh; j++) (*U)[j] = vect_U[j];

	if(diag_A) free_dvector(diag_A,1,num_mesh-1); diag_A = NULL;
	if(diag_B) free_dvector(diag_B,1,num_mesh-1); diag_B = NULL;
	if(diag_C) free_dvector(diag_C,1,num_mesh-1); diag_C = NULL;
	if(vect_R) free_dvector(vect_R,1,num_mesh-1); vect_R = NULL;
	if(vect_U) free_dvector(vect_U,1,num_mesh-1); vect_U = NULL;
	if(PrevPrePtr) free_dvector(PrevPrePtr,1,num_mesh-1); PrevPrePtr = NULL;

	return (err);

}

/* FUNCTION TO TEST THE FEYNMANN KAC SOLVER */

Err srt_f_black_scholes_pde(Date		eval_date,
							Date		exp_date,
							double      spot,
							double      forward,
							double      equity_strike,
							double      equity_vol,
							double		min_node,
							double      min_mesh,
							
							/*OUTPUT */
							double	*bs_pv)
{
	Err					err  = NULL;
	long				i,j, k, num_mesh, num_time_step;
	double				YrtoExpDate,X0,Xmin,Xmax,dXmax,time_step = MAX_TIME,dX, 
						*X = NULL, *U = NULL,**diff_mat = NULL, **conv_mat = NULL,**div_mat = NULL,
						*term_payoff = NULL, pde_bound[2],drift;
	SrtPdePDEType		pde_type = LINEAR;
	SrtPdeBoundaryCond	pde_bound_cond = DIRICHLET;

	YrtoExpDate = (double)(exp_date-eval_date)*YEARS_IN_DAY;

	X0 = log(spot);
	Xmin = X0 - ONE_HALF*equity_vol*equity_vol*YrtoExpDate - LIM_IN_STDEV*equity_vol*sqrt(YrtoExpDate);
	Xmax = X0 - ONE_HALF*equity_vol*equity_vol*YrtoExpDate + LIM_IN_STDEV*equity_vol*sqrt(YrtoExpDate);
	dXmax = equity_vol*equity_vol/fabs(-ONE_HALF*equity_vol*equity_vol);

	drift = log(forward/spot)/YrtoExpDate;

	time_step = min(time_step,YrtoExpDate/min_node);
	num_time_step = ((long) (YrtoExpDate/ time_step));

	/* COMPUTE THE SPACE STEP */
	dX = (Xmax-Xmin)/min_mesh;
	dX = min(dX,dXmax);

	/* ADAPT dX SO THAT X0 IS A MESH POINT */
	k = (long)((X0-Xmin)/dX) + 1; 
	dX = (X0-Xmin)/k;
	num_mesh = (long)((Xmax-Xmin)/dX);
	
	/* ALLOCATE AND FILL X, U, conv_mat, diff_mat */
	X			   = dvector(1,num_mesh-1);
	term_payoff	   = dvector(1,num_mesh-1);
	U			   = dvector(1,num_mesh-1);
	diff_mat	   = dmatrix(0,num_time_step,1,num_mesh);
	conv_mat	   = dmatrix(0,num_time_step,1,num_mesh);
	div_mat		   = dmatrix(0,num_time_step,1,num_mesh);
	
	for(j = 1; j < num_mesh; j++) 
	{
		X[j] = Xmin + j*dX;
		term_payoff[j] = max(equity_strike-exp(X[j]),0);

		for(i = 0; i <= num_time_step; i++)
		{
			conv_mat[i][j] = drift-ONE_HALF*equity_vol*equity_vol;
			diff_mat[i][j] = equity_vol;
			div_mat[i][j] = 0.0;
		}
	}

	pde_bound[0] = max(equity_strike-exp(Xmin),0);
	pde_bound[1] = max(equity_strike-exp(Xmax),0);
	
	err = srt_f_feynmann_kac_pde(num_mesh,
							  num_time_step,
							  dX,
							  time_step,
							  X,
							  conv_mat,
							  diff_mat,
							  pde_type,
							  pde_bound_cond,
							  pde_bound,
							  term_payoff,
							  NULL,
							  div_mat,
							  &U);
	if(err) return err;

	(*bs_pv) = U[k];

	if(X) free_dvector(X,1,num_mesh-1); X = NULL;
	if(U) free_dvector(U,1,num_mesh-1); U = NULL;
	if(diff_mat) free_dmatrix(diff_mat,0,num_time_step,1,num_mesh); diff_mat = NULL;
	if(conv_mat) free_dmatrix(conv_mat,0,num_time_step,1,num_mesh); conv_mat = NULL;
	if(div_mat) free_dmatrix(div_mat,0,num_time_step,1,num_mesh); div_mat = NULL;
	if(term_payoff) free_dvector(term_payoff,1,num_mesh-1); term_payoff = NULL;

	return err;


}


#undef PDE_THETA 
#undef LIM_IN_STDEV	
#undef ONE_HALF 
#undef MAX_TIME
