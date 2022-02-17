
#ifndef Fx3FMultiGrfnH
#define	Fx3FMutliGrfnH

#include "grf_h_mdlcomm.h"

typedef struct
{
	int						num_und;		/* number of underlyings */
	int						nb_dates;		/* number of dates */
	int						*type;			/* 0 for domestic LGM, 1 for foreign LGM, 2 for FX */
	int						*dom_forex;		/* index of the corresponding domestic to the forex i*/
	int						*for_forex;		/* index of the corresponding foreign to the forex i*/
	int						*fx_index;		/* index of the fx underlying FOREIGN LGM / DOMESTIC LGM */
	int						dom_index;		/* index of the domestic underlying */
	char					**yc;
	double					*lambda;
	double					*spot_fx;
	double					*sig_dates;
	long					nb_sig_dates;
	double					**sig_curve;
	double					**phi;
	double					**fwd;
	double					**std;
	double					**beta;
	double					*dom_bond_pay;
	double					*dom_beta_pay;
	double					***correlations;	/* correlation between underlyings i and j for each of nb_sig_dates */
	double					***covariance;	

}	link_und, *LINK_UND;



typedef struct
{
	GRFNCOMMSTRUCT			global;
	FIRSTMktAtT				*local;
	
	long					*num_df;		

	LINK_UND				link;
	double					**df_tms;
	long					**df_dts;
	double					**dff;
	double					**gam;
	double					**gam2;

	int						*do_und;
	
} grfn_parm_multi_mc, *GRFNPARM_MULTIMC;



Err allocate_link_und(
						LINK_UND	link
					 );

Err free_link_und	(
					LINK_UND link
					);

Err grfn_payoff_multi_3dfx_mc (
								/* Event */
								double		evt_date,
								double		evt_time,
								void		*func_parm, 
								/* Market data */
								LINK_UND	link,
								double		*sv,
								/* Results */
								int			num_col,
								double		*res);
#endif
