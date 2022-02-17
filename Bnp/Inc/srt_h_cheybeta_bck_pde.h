#ifndef SRT_H_CHEYBETA_BCK_PDE_H
#define SRT_H_CHEYBETA_BCK_PDE_H
 
/*	Backward PDE Functions	*/

/*	Temp variable structure	*/
typedef struct
{
	double						*coef;
	double						**val;
	double						*mu;
	double						*var;
	double						*r;
} CHEYBETA_BCK_PDE_TEMP;

/*	Allocate a few temp vectors	*/
void srt_f_chey_beta_bck_pde_alloc_temp(
CHEYBETA_BCK_PDE_TEMP			*tmp,
int								nx,
int								nphi,
int								num_prod);

/*	Free temp vectors	*/
void srt_f_chey_beta_bck_pde_free_temp(
CHEYBETA_BCK_PDE_TEMP			*tmp);

/*	Go backward one step	*/
void srt_f_chey_beta_bck_pde_one_step(
CNPDE_TEMP						*tmppde,
CHEYBETA_BCK_PDE_TEMP			*tmpcb,
/*	Model	*/
CHEYBETA_MDL					*mdl,
/*	Grid	*/
CHEYBETA_GRID					*grid,
/*	Current time step idx	*/
int								t1,
/*	Next time step idx	*/
int								t2);

/*	Go backward n steps	*/
void srt_f_chey_beta_bck_pde_n_steps(
CHEYBETA_BCK_PDE_TEMP			*tmpcb,
/*	Model	*/
CHEYBETA_MDL					*mdl,
/*	Grid	*/
CHEYBETA_GRID					*grid,
/*	Current time step idx	*/
int								*ti,
/*	Number of steps backward	*/
int								n);

#endif