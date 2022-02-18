#ifndef SRT_H_CHEYBETA_FWD_PDE_H
#define SRT_H_CHEYBETA_FWD_PDE_H

/*	Forward PDE Functions	*/

/*	Temp variable structure	*/
typedef struct
{
    double* coef;
    double* ad;
    double* mu;
    double* var;
    double* r;
} CHEYBETA_FWD_PDE_TEMP;

/*	Allocate a few temp vectors	*/
void srt_f_chey_beta_fwd_pde_alloc_temp(CHEYBETA_FWD_PDE_TEMP* tmp, int nx, int nphi);

/*	Free temp vectors	*/
void srt_f_chey_beta_fwd_pde_free_temp(CHEYBETA_FWD_PDE_TEMP* tmp);

/*	Get Arrow-Debreu prices @ first step	*/
void srt_f_chey_beta_fwd_pde_first_step(
    /*	Model	*/
    CHEYBETA_MDL* mdl,
    /*	Grid	*/
    CHEYBETA_GRID* grid);

/*	Go forward one step	*/
void srt_f_chey_beta_fwd_pde_one_step(
    CNPDE_TEMP*            tmppde,
    CHEYBETA_FWD_PDE_TEMP* tmpcb,
    /*	Model	*/
    CHEYBETA_MDL* mdl,
    /*	Grid	*/
    CHEYBETA_GRID* grid,
    /*	Current time step idx	*/
    int t1,
    /*	Next time step idx	*/
    int t2);

/*	Go forward n steps	*/
void srt_f_chey_beta_fwd_pde_n_steps(
    CHEYBETA_FWD_PDE_TEMP* tmpcb,
    /*	Model	*/
    CHEYBETA_MDL* mdl,
    /*	Grid	*/
    CHEYBETA_GRID* grid,
    /*	Current time step idx	*/
    int* ti,
    /*	Number of steps forward	*/
    int n);

#endif