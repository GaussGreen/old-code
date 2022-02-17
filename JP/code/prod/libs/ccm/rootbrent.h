/* $Header$ */
typedef int (*TObjectFunc)(double, void*, double *);

int RootFindBrent(           
   TObjectFunc funcd,                   /* (I) Function to call */
   void       *data,                    /* (I) Data to pass into funcd */
   double      boundLo,                 /* (I) Lower bound on legal X */
   double      boundHi,                 /* (I) Upper bound on legal X */
   int         numIterations,           /* (I) Maximum number of iterations */
   double      guess,                   /* (I) Initial guess */
   double      initialXStep,            /* (I) Size of step in x */
   double      initialFDeriv,           /* (I) Initial derivative or 0*/
   double      xacc,                    /* (I) X accuracy tolerance */
   double      facc,                    /* (I) Function accuracy tolerance */
   double      *solution);              /* (O) Root found */


