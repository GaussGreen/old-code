
#ifndef Fx3FTreeHelpH
#define	Fx3FTreeHelpH

void fill_fwd_var(	
					long		nstp,
					double		*time,
					double		*date,
					double		*sig_dom,		
					double		*sig_for, 
					double		*sig_fx,
					double		dom_lam, 
					double		for_lam, 
					double		corr_dom_for,
					double		corr_dom_fx,
					double		corr_for_fx,
					char		*dom_yc,
					char		*for_yc,
					/*	To be allocated by caller, filled by the function */
					double		*dom_ifr,
					double		*dom_fwd,
					double		*dom_var,
					double		*for_ifr,
					double		*for_fwd,
					double		*for_var,
					double		*fx_fwd,
					double		*fx_var);

void fill_fwd_var_corr(	
					long		nstp,
					double		*time,
					double		*date,
					double		*sig_dom,		
					double		*sig_for, 
					double		*sig_fx,
					double		dom_lam, 
					double		for_lam, 
					double		*corr_dom_for,
					double		*corr_dom_fx,
					double		*corr_for_fx,
					char		*dom_yc,
					char		*for_yc,
					/*	To be allocated by caller, filled by the function */
					double		*dom_ifr,
					double		*dom_fwd,
					double		*dom_var,
					double		*for_ifr,
					double		*for_fwd,
					double		*for_var,
					double		*fx_fwd,
					double		*fx_var);

#endif
