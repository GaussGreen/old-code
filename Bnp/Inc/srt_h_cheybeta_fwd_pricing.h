#ifndef SRT_H_CHEYBETA_FWD_PRICING_H
#define SRT_H_CHEYBETA_FWD_PRICING_H

/*	Swaption data struct	*/

/*	Data to be stored for each option to be priced	*/
typedef struct
{
	int				ncpn;
	
	/*	Times from today	*/
	double			exp_time;
	double			*pay_time;

	double			*pay_cpn;
	double			*pay_dff;
	
	/*	Reconstruction formula parameters, the only ones to be updated
		when model parameters change	*/
	double			*pay_gamma;
	double			*pay_0_50_gamma_sqr;

	/*	Pricing result	*/
	double			prem;
}	SWAPTION_DATA;

/*	Global data on all swaptions and expiries	*/
typedef struct
{
	/*	Number of swaptions	*/
	int				nswapt;
	/*	Swaption data	*/
	SWAPTION_DATA	*data;

	/*	Number of instruments	*/
	int				ninst;
	/*	Number of swaptions per instrument	*/
	int				*ndata_i;
	/*	Index of swaptions per instrument	*/
	int				**data_i; 

	/*	Number of expiries	*/
	int				nex;
	/*	Expiry times	*/
	double			*ex;
	/*	Number of swaptions per expiry	*/
	int				*ndata_e;
	/*	Index of swaptions per expiry	*/
	int				**data_e;
}	SWAPTIONS_DATA;

Err add_new_swaption (
/*	Underlyin info	*/
CHEYBETA_MDL					*pmdl,
char							*yc_name,
long							today,
int								spot_lag,
/*	Index of the instrument the option belongs to	*/
int								inst_idx,
/*	Option info	*/
long							expiry,
long							start,
long							end,
char							*freq,
char							*basis,
double							strike,
double							bond_strike,
char							*pay_rec,
char							*refrate,
SWAPTIONS_DATA					*g);

/*	Init swaption data	*/

Err chey_beta_initmem_swaption_data(
int								ninst,
SWAPTIONS_DATA					*g );

Err chey_beta_fwd_init_swaption_data(
CHEYBETA_MDL					*pmdl,
int								ninst,
long							*start,
long							*end,
char							**freq_str,
char							**basis_str,
double							*strike,
double							*bond_strike,
char							**type_str,
char							**pay_rec_str,
char							**refrate_str,
SWAPTIONS_DATA					*g);

/*	Update gammas in swaption data as model parameters change	*/

Err chey_beta_fwd_update_swaption_data_for_model(
CHEYBETA_MDL					*pmdl,
SWAPTIONS_DATA					*g);

/*	Free globals	*/

Err chey_beta_fwd_free_swaption_data (
SWAPTIONS_DATA					*g);

/*	Do price	*/

Err chey_beta_fwd_price(
CHEYBETA_MDL					*pmdl,
int								nt,
int								nx,
int								nphi,
SWAPTIONS_DATA					*g,
/*	Allocated inside, to be freed by caller	*/
int								*ninst,
double							**inst_prem);

#endif