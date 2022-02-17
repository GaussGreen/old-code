#ifndef SRT_H_CHEYBETA_MDL_DEF_H
#define SRT_H_CHEYBETA_MDL_DEF_H

/*	Model definition	*/

typedef struct {
  /*	Today	*/
  double today_date;
  /*	Name of the yield curve	*/
  char *yc_name;
  /*	Number of sigma points	*/
  int num_sigma;
  /*	Sigma times	*/
  double *sigma_times;
  /*	Sigmas	*/
  double *sigma;
  /*	Number of tau points	*/
  int num_lambda;
  /*	Lambda times	*/
  double *lambda_times;
  /*	Lambdas	*/
  double *lambda;
  /*	Number of beta points	*/
  int num_beta;
  /*	Beta times	*/
  double *beta_times;
  /*	Betas	*/
  double *beta;
  /*	Vol adjustment: vol = (Rt + vol_adj)^beta	*/
  double vol_adj;
} CHEYBETA_MDL;

/*	Initialise      , i.e. set pointers to NULL	*/
void chey_beta_mdl_init(CHEYBETA_MDL *mdl);

/*	Free	*/
void chey_beta_mdl_free(CHEYBETA_MDL *mdl);

/*	Build from description	*/
void chey_beta_mdl_build(CHEYBETA_MDL *mdl,
                         /*	Name of the yield curve	*/
                         char *yc_name,
                         /*	Number of sigma points	*/
                         int num_sigma,
                         /*	Sigma times	*/
                         double *sigma_times,
                         /*	Sigmas	*/
                         double *sigma,
                         /*	Number of tau points	*/
                         int num_lambda,
                         /*	Lambda times	*/
                         double *lambda_times,
                         /*	Lambdas	*/
                         double *lambda,
                         /*	Beta	*/
                         double beta);

/*	Build from underlying	*/
void chey_beta_mdl_build_from_und(CHEYBETA_MDL *mdl, SrtUndPtr und);

/*
        Diffusion parameters between two given dates:

                sigma      , lambda      , beta      , dt
*/

typedef struct {
  double sigma;
  double lambda;
  double beta;
  double dt;
} CHEYBETA_PARAM;

/*	Fill the diffusion parameters structure at a given date	*/
void chey_beta_mdl_param(
    /*	Model	*/
    CHEYBETA_MDL *mdl,
    /*	Times	*/
    double t1, double t2,
    /*	Result	*/
    CHEYBETA_PARAM *param);

/*	Get normal instantaneous variance at t	*/
double chey_beta_mdl_norm_var_at_t(
    /*	Param	*/
    CHEYBETA_PARAM *param,
    /*	Statevars	*/
    double x, double phi,
    /*	Forward rate	*/
    double f,
    /*	Max var	*/
    double maxvar,
    /*	Min var	*/
    double minvar);

/*	Get drift      , actually Et2 [ Xt2 - Xt1 / Ft1 ]	*/
double chey_beta_mdl_drift(
    /*	Model	*/
    CHEYBETA_MDL *mdl,
    /*	Param	*/
    CHEYBETA_PARAM *param,
    /*	Statevars	*/
    double x, double phi);

/*	Get var      , actually Vt2 [ Xt2 - Xt1 / Ft1 ]	*/
double chey_beta_mdl_var(
    /*	Model	*/
    CHEYBETA_MDL *mdl,
    /*	Param	*/
    CHEYBETA_PARAM *param,
    /*	Statevars	*/
    double x, double phi,
    /*	Norm var	*/
    double norm_var);

/*	Get forward phi	*/
double chey_beta_mdl_phi(
    /*	Model	*/
    CHEYBETA_MDL *mdl,
    /*	Param	*/
    CHEYBETA_PARAM *param,
    /*	Statevars	*/
    double x, double phi,
    /*	Output from chey_beta_mdl_var	*/
    double var);

#endif