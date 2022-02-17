/* Workspace for adaptive integrators */
typedef struct
{
    size_t limit;
    size_t size;
    size_t nrmax;
    size_t i;
    size_t maximum_level;
    double *alist;
    double *blist;
    double *rlist;
    double *elist;
    size_t *order;
    size_t *level;
}
gsl_integration_workspace;

gsl_integration_workspace * gsl_integration_workspace_alloc (const size_t n);
void gsl_integration_workspace_free (gsl_integration_workspace * w);

typedef int (*PTR_FUNCTION)(double, double *, void *);

/**********************************************************/
/* INTEGRATION ON A FINITE INTERVAL                       */
/**********************************************************/
int gsl_integration_qag (   const PTR_FUNCTION f,
			                double a,
                            double b,
			                double epsabs,
                            double epsrel,
                            size_t limit,
			                int key,
			                gsl_integration_workspace * workspace,
			                double *result,
                            double *abserr,
                            void *param);

/**********************************************************/
/* INTEGRATION ON AN INFINITE INTERVAL                    */
/**********************************************************/
int gsl_integration_qagi (  PTR_FUNCTION f,
			                double epsabs,
                            double epsrel,
                            size_t limit,
			                gsl_integration_workspace * workspace,
			                double *result,
                            double *abserr,
                            void *param);

/**********************************************************/
/* INTEGRATION ON A SEMIFINITE INTERVAL [a,+infinity[     */
/**********************************************************/
int gsl_integration_qagiu ( PTR_FUNCTION f,
			                double a,
			                double epsabs,
                            double epsrel,
                            size_t limit,
			                gsl_integration_workspace * workspace,
			                double *result,
                            double *abserr,
                            void *param);

/**********************************************************/
/* INTEGRATION ON A SEMIFINITE INTERVAL ]-infinity,b]     */
/**********************************************************/
int gsl_integration_qagil ( PTR_FUNCTION f,
			                double b,
			                double epsabs,
                            double epsrel,
                            size_t limit,
			                gsl_integration_workspace * workspace,
			                double *result,
                            double *abserr,
                            void *param);

/**********************************************************/
/* INTEGRATION ON A SEMIFINITE INTERVAL                       */
/**********************************************************/
int gsl_integration_qags (  const PTR_FUNCTION f,
			                double a,
                            double b,
			                double epsabs,
                            double epsrel,
                            size_t limit,
			                gsl_integration_workspace * workspace,
			                double *result,
                            double *abserr,
                            void *param);




