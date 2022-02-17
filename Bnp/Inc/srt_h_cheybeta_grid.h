#ifndef SRT_H_CHEYBETA_GRID_H
#define SRT_H_CHEYBETA_GRID_H

/*	Grid definition	*/

typedef struct {
  /*	Time grid	*/
  int nt;
  double *t;

  /*	Max var / Min var info	*/
  double *maxvar;
  double *minvar;

  /*	Forward info	*/
  /*	Forward rate from ti to ti+1 (last = 0)	*/
  double *fwdrt;

  /*	X grid	*/
  int nx;
  double *x;

  /*	Phi grid	*/
  /*	Maximum number of phis	*/
  int maxnphi;
  /*	Number of phis at time step	*/
  int *nphi;
  /*	Phis at time step	*/
  double **phi;

  /*	For the forward PDE	*/

  /*	AD grids  , idx1 = phi  , idx2 = x	*/

  /*	AD prices @ step t1	*/
  double **adt1;
  /*	AD prices @ step t2	*/
  double **adt2;

  /*	For the backward PDE	*/

  int num_prod;

  /*	Value grids  , idx1 = phi  , idx2 = x  , idx3 = product id	*/

  /*	Prices @ step t1	*/
  double ***vt1;
  /*	Prices @ step t2	*/
  double ***vt2;

} CHEYBETA_GRID;

/*	Initialise  , i.e. set pointers to NULL	*/
void chey_beta_grid_init(CHEYBETA_GRID *grid);

/*	Free	*/
void chey_beta_grid_free(CHEYBETA_GRID *grid);

/*	Construct time grid	*/
void chey_beta_grid_time(
    /*	Model	*/
    CHEYBETA_MDL *mdl,
    /*	Vector of required stop times	*/
    int num_req_times, double *req_times,
    /*	Number of steps required	*/
    int nt,
    /*	Grid	*/
    CHEYBETA_GRID *grid);

/*	Construct phi and x grids	*/
void chey_beta_grid_phi_x(
    /*	Model	*/
    CHEYBETA_MDL *mdl,
    /*	Grid	*/
    CHEYBETA_GRID *grid,
    /*	Num x	*/
    int nx,
    /*	Num phi	*/
    int nphi);

/*	For the forward PDE	*/

/*	Allocate AD grids	*/
void chey_beta_grid_ad(
    /*	Grid	*/
    CHEYBETA_GRID *grid);

/*	Copy density	*/
void chey_beta_density_at_t(
    /*	Grid	*/
    CHEYBETA_GRID *grid,
    /*	Time step index	*/
    int ti,
    /*	X grid	*/
    int *nx, double **x,
    /*	Phi grid	*/
    int *nphi, double **phi,
    /*	AD grid	idx1 = phi  , idx2 = x	*/
    double ***ad);

/*	Free density	*/
void chey_beta_free_density_at_t(
    /*	X grid	*/
    int *nx, double **x,
    /*	Phi grid	*/
    int *nphi, double **phi,
    /*	AD grid	idx1 = phi  , idx2 = x	*/
    double ***ad);

/*	For the backward PDE	*/

/*	Allocate value grids	*/
void chey_beta_grid_prod_val(
    /*	Grid	*/
    CHEYBETA_GRID *grid,
    /*	Number of products to be valued	*/
    int num_prod);

/*	Get	product values	*/
void chey_beta_prod_val_at_t(
    /*	Grid	*/
    CHEYBETA_GRID *grid,
    /*	Time step index	*/
    int ti,
    /*	X value required	*/
    double x,
    /*	Phi value required	*/
    double phi,
    /*	Product values	*/
    int *num_prod,
    /*	Allocated inside  , to be freed by caller	*/
    double **prod_val);

#endif