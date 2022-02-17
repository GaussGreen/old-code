
typedef struct
{
    int kmax;
    int nbPoints;
    double dxsav;
    double *xp;
    double **yp;
} DEBUGINFO;


int rkqs(   double y[],
            double dydx[],
            int n,
            double *x,
            double htry,
            double eps,
	        double yscal[],
            double *hdid,
            double *hnext,
	        int (*derivs)(  double,
                            double [],
                            double [],
                            void *),
            void *param);

int odeint( double ystart[],
            int nvar,
            double x1,
            double x2,
            double eps,
            double h1,
	        double hmin,
            int maxNbPoints,
            int *nok,
            int *nbad,
	        int (*derivs)(  double,
                            double [],
                            double [],
                            void *),
	        int (*rkqs)(    double [],
                            double [],
                            int,
                            double *,
                            double,
                            double,
                            double [],
	                        double *,
                            double *,
                            int (*)(    double,
                                        double [],
                                        double [],
                                        void *),
                            void *param),
            void *param,
            DEBUGINFO *debugInfo);

