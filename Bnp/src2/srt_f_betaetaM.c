#include "math.h"
#include "num_h_simpso.h"
#include "srt_h_all.h"
#include "srt_h_betaetaM.h"
#include "srt_h_ts_irm.h"
/*
        This file contains the code to compute the M matrix for the beta-eta
   model. The routine make_M_cube should be called from init_ir_und in grfn.

          Glossary:
          1) make_M_matrix()
          2) I_nu  , K_nu()
          4) integrand_1()
          5) integrand_2()
          6) integrand_3()
          7) integrand_4()
          8) print_M()
          9) calc_bigM()

*/

#define EPS_SIMPSO (10e-6)

/* Global static variables for use with Simpson's rule */
static double STATIC_BETA;
static double STATIC_ETA;
static double STATIC_X_K;
static double STATIC_LAMBDA;
static double STATIC_DELTAU;
static double STATIC_START;
static double STATIC_NUMOFPER;

/*<%%STA******************************************************************
  FUNCNAME        :make_M_matrix
  AUTHOR          :A.Berner
  DESCRIPTION     :Driver routine to compute the M matrix for the beta-eta model

  CALLS: calc_bigM

******************************************************************************/

Err make_M_matrix(TermStruct *l) {

  IrmTermStructVal *c; /* c  = pval of current list element	*/
  SrtLst *lc;          /* lc = current list element			*/
  Err err = NULL;
  double beta; /* scaling quantity for distribution	*/
  double eta;  /* power for distribution formulation	*/
  double **M;
  double x;         /* state variable						*/
  double delta_tau; /*										*/
  double lambda;    /*										*/
  int s, theta;     /* M(s  ,theta)     */
  double nu;        /* nu of bessel function				*/
  double a, b;      /* integration range					*/
  double y;
  double log(double);

  lc = l->head; /* first element of the list */
  c = (IrmTermStructVal *)lc->element->val.pval;
  /* Check beta  , eta arguments */
  if (c->beta <= 0.)
    return serror("input for beta invalid");
  if (c->eta < 0. || c->eta > 1.)
    return serror("input for eta invalid");

  /*
  Compute M[s][theta] for 0<s<SMAX  , 0<theta<THETAMAX.
  Hardwire bet and x for computation of the M matrix.
  We can do this because of the way the matrix is stored  ,
  Since s and theta are a combo of lambda  ,beta  , etc. we can choose the
  easiest combo of lambda  ,beta  , etc. to compute and store M(s  ,theta).
  So we set  beta=1  , x=0 which makes s dependent only on delta_tau
  and theta dependent only on lambda. */

  eta = c->eta;
  beta = 1.;
  x = 0.;

  M = dmatrix(0, SMAX - 1, 0, THETAMAX - 1);
  if (M == NULL)
    return serror("Error on M memory allocation");

  /* Compute M[s][theta] based on eta */
  s = 0;
  while (s < SMAX) {
    theta = 0;
    while (theta < THETAMAX) {

      /* With beta = 1  , x = 0 -> delta_tau = s/SMAX  , lambda = theta/THETAMAX
       */
      delta_tau = (double)s / SMAX;
      lambda = (double)theta / THETAMAX;

      if ((s == 0) || (theta == 0)) /* zero soln on boundary     */
        M[s][theta] = 0.;
      else if (eta == 0) /* Analytic soln for eta=.0  */
        M[s][theta] = (.5 * delta_tau * lambda * lambda);
      else if (eta == .5) /* Analytic soln for eta=.5  */
        M[s][theta] =
            (.5 * delta_tau * lambda * lambda) / (1. + .5 * delta_tau * lambda);
      else {
        /* Use simpson's rule to perform integration over density for rest of
         * cases */

        /* Compute range of integration for the two cases eta<1  , eta==1  */
        if (eta < 1.) /* integration range needs to stay positive */
        {
          /* Compute range of integration out to 5 S.D.'s */
          y = (pow(fabs(1. + beta * x), 1 - eta)) / (beta * (1. - eta));
          a = y - 5. * sqrt(sqrt(2. * delta_tau));
          if (a <= 0.) { /* check if distribution goes negative and reset close
                            to zero. */
            a = .000000001;
          }
          /* Need to make a better choice on b */
          b = y + 5. * sqrt(sqrt(2. * delta_tau));
          b = 25.;
        } else /* ETA == 1 */
        {
          /* Need to make a better choice here */
          a = -25.;
          b = 25.;
        }

        /* Compute M based on case 1-4 from eta */
        if (0 <= eta && eta < .5)
          /* Case 1:  No barrier b.c.'s */
          M[s][theta] = log(sm_qsimp(integrand_1, a, b, EPS_SIMPSO * .0000001));
        else if (.5 <= eta && eta < 1) { /*  Case 2 Natural b.c.'s */
          nu = 1. / (2 * (1. - eta));
          y = pow(fabs(1 + beta * x), 1 - nu) / (beta * (1 - eta));
          M[s][theta] = log(sm_qsimp(integrand_3, a, b, EPS_SIMPSO) +
                            exp(-lambda * (1. / beta - x) - gammln(nu)) *
                                (gammq(nu, (y * y) / (2. * delta_tau)) -
                                 pow((y * y) / (2. * delta_tau), nu - 1.) *
                                     exp(-(y * y) / (2. * delta_tau))));
        } else
            /* Case 4 */
            if (eta == 1)
          M[s][theta] = log(sm_qsimp(integrand_4, a, b, EPS_SIMPSO));
      }

      /*						err = calc_bigM(M  ,
         beta  , eta  , lambda  , delta_tau  , x  , s_index  , theta_index); if
         (err) return err;
      */
      theta++;
    }
    s++;
  }

  c->M_beta_eta = M; /* Put M into termstructure */
  c->is_M_allocated_here = SRT_YES;

  while (lc->next != NULL) /* At least 2 elements in the list */
  {
    /*
    No need to create extra copies of M for each node  , just have each
    pointer point to the first one
    */
    lc = lc->next;
    c = (IrmTermStructVal *)lc->element->val.pval; /* <=> [i]   */
    c->M_beta_eta = M;
    c->is_M_allocated_here = SRT_NO;
  }
  return err;
}

/* Integrands for cases 1-4 */

static double integrand_1(double y_bar) {

  /* NB: Case 1a now contains both 1a and 1b so the 1b integrand
     is no longer needed */

  /* Computes the integrand for Case 1: 0<=STATIC_ETA<.5  .5 <= nu < 1
     No barrier b.c.'s
   */

  double ans;
  double nu;
  double y;
  double xb_yb;
  double dxb_yb;
  double F;
  double F_xb;

  double xt_yb;
  double i_nu, i_minus_nu, k_nu;
  Err err = NULL;

  nu = 1. / (2. * (1. - STATIC_ETA));
  y = (pow(fabs(1. + STATIC_BETA * STATIC_X_K), 1 - STATIC_ETA)) /
      (STATIC_BETA * (1. - STATIC_ETA));
  xb_yb = (pow(fabs(STATIC_BETA * (1 - STATIC_ETA) * y_bar), 2 * nu) - 1.) /
          STATIC_BETA;
  dxb_yb = 2 * nu * pow((1. - STATIC_ETA) * STATIC_BETA * y_bar, 2 * nu - 1);
  F = exp(-STATIC_LAMBDA * (xb_yb - STATIC_X_K));
  F_xb = -STATIC_LAMBDA * F;
  xt_yb = (-pow(fabs(STATIC_BETA * (1 - STATIC_ETA) * y_bar), 2 * nu) - 1.) /
          STATIC_BETA;

  err = I_nu(-(1 - nu), y * y_bar / STATIC_DELTAU, &i_nu);
  err = I_nu(1 - nu, y * y_bar / STATIC_DELTAU, &i_minus_nu);
  err = K_nu(1 - nu, y * y_bar / STATIC_DELTAU, &k_nu);

  ans =
      /* 1a */
      (pow(y * y / STATIC_DELTAU, nu - 1) *
       (STATIC_LAMBDA *
            pow(STATIC_BETA * (1 - STATIC_ETA) * y_bar, 2. * nu - 1.) +
        y_bar / STATIC_DELTAU) *
       exp(-STATIC_LAMBDA * (xb_yb - STATIC_X_K)) *
       pow(y_bar * y / STATIC_DELTAU, 1 - nu) * (.5 * i_nu + .5 * i_minus_nu) *
       /* exp(-y*y_bar/STATIC_DELTAU)*  taken out due to e^{-x} from N.R. I_nu
        */
       exp(-((y_bar - y) * (y_bar - y)) /
           (2 *
            STATIC_DELTAU))) - /* minus sign comes from integration by parts */
                               /* 1b */
      ((pow(y * y / STATIC_DELTAU, nu - 1)) * (sin(SRT_PI * nu) / SRT_PI) *
       (-STATIC_LAMBDA *
            pow(STATIC_BETA * (1. - STATIC_ETA) * y_bar, 2 * nu - 1) +
        y_bar / STATIC_DELTAU) *
       exp(-STATIC_LAMBDA * (xt_yb - STATIC_X_K)) *
       pow(y * y_bar / STATIC_DELTAU, 1 - nu) * k_nu *
       exp(-(y_bar + y) * (y_bar + y) / (2 * STATIC_DELTAU)));

  return ans;
}

static double integrand_2(double y_bar) {

  /* Computes the integrand for Case 2: 0<=STATIC_ETA<.5  .5 <= nu < 1
     with reflecting b.c.'s  This routine is currently not used
     nor hooked into the code */

  double y;
  double xb_yb;
  double dxb_yb;
  double nu;
  double ans; /* value of integrand */
  double F;
  double F_xb;
  double i_nu;
  Err err = NULL;

  nu = 1 / (2 * (1 - STATIC_ETA));
  y = (pow(fabs(1 + STATIC_BETA * STATIC_X_K), 1 - STATIC_ETA)) /
      (STATIC_BETA * (1 - STATIC_ETA));
  xb_yb = (pow(fabs(STATIC_BETA * (1 - STATIC_ETA) * y_bar), 2 * nu) - 1) /
          STATIC_BETA;
  dxb_yb = pow(STATIC_BETA * (1 - STATIC_ETA) * y_bar, 2 * nu - 1);
  F = exp(-STATIC_LAMBDA * (xb_yb - STATIC_X_K));
  F_xb = -STATIC_LAMBDA * F;

  err = I_nu(1. - nu, y * y_bar / STATIC_DELTAU, &i_nu);

  ans = ((y_bar / STATIC_DELTAU) * F - F_xb * dxb_yb) * pow(y_bar / y, 1 - nu) *
        i_nu *
        /*  exp(y*y_bar/STATIC_DELTAU)*    */ /* taken out due to NR scaling */
        exp(-(y_bar - y) * (y_bar - y) / (2 * STATIC_DELTAU));

  return ans;
}

static double integrand_3(double y_bar) {

  /* Computes the integrand for Case 3: .5 <=eta< 1   nu >= 1 */

  double y;
  double xb_yb;
  double dxb_yb;
  double nu;
  double ans; /* value of integrand */
  double F;
  double F_xb;
  double i_nu;
  Err err = NULL;

  nu = 1 / (2 * (1 - STATIC_ETA));
  y = (pow(fabs(1 + STATIC_BETA * STATIC_X_K), 1 - STATIC_ETA)) /
      (STATIC_BETA * (1 - STATIC_ETA));
  xb_yb = (pow(fabs(STATIC_BETA * (1 - STATIC_ETA) * y_bar), 2 * nu) - 1) /
          STATIC_BETA;
  dxb_yb = pow(STATIC_BETA * (1 - STATIC_ETA) * y_bar, 2 * nu - 1);
  F = exp(-STATIC_LAMBDA * (xb_yb - STATIC_X_K));
  F_xb = -STATIC_LAMBDA * F;

  err = I_nu(nu - 1., y * y_bar / STATIC_DELTAU, &i_nu);

  ans = ((y_bar / STATIC_DELTAU) * F - F_xb * dxb_yb) * pow(y_bar / y, 1 - nu) *
        i_nu *
        /*	  exp(-y*y_bar/STATIC_DELTAU)*    taken out due to NR scaling */
        exp(-(y_bar - y) * (y_bar - y) / (2 * STATIC_DELTAU));

  return ans;
}

static double integrand_4(double y_bar) {

  /* Computes the integrand for Case 4: eta = 1 (bdy at - infty) */

  double y;
  double xb_yb;
  double ans; /* value of integrand */
  double F;

  y = log(1. + STATIC_BETA * STATIC_X_K);
  xb_yb = (exp(STATIC_BETA * y_bar) - 1.) / STATIC_BETA;
  F = exp(-STATIC_LAMBDA * (xb_yb - STATIC_X_K));

  ans = F *
        exp(-(y_bar - y + .5 * STATIC_BETA * STATIC_DELTAU) *
            (y_bar - y + .5 * STATIC_BETA * STATIC_DELTAU) /
            (2 * STATIC_DELTAU)) /
        sqrt(2 * SRT_PI * STATIC_DELTAU);

  return ans;
}

#if 0
Err calc_bigM(double **M  , double beta  , double eta  , double lambda  , double deltau  , double x_k  , 
			   int i  , int j)
{
 
	double nu;
	double log(double);
	double a  ,b;
	double y;


		/* Assign parameter inputs to the global static variables */
		STATIC_BETA   = beta;
		STATIC_ETA    = eta;
		STATIC_LAMBDA = lambda;
		STATIC_DELTAU = deltau;
		STATIC_X_K    = x_k;

    	nu = 1./(2*(1.- STATIC_ETA));

		/* Compute range of integration out to 5 S.D.'s 
		   a = left bound ; b = right bound */
	    if (STATIC_ETA != 1)
		{/* exponent < 1 so range stays positive */
			y = (pow(fabs(1.+STATIC_BETA*STATIC_X_K)  ,1-STATIC_ETA))/(STATIC_BETA*(1.-STATIC_ETA)); 
			a = y - 5.*sqrt(sqrt(2.*STATIC_DELTAU));
			if (a <= 0.) 
			{/* check if distribution goes negative and reset close to zero. */
				a = .000000001;	
			}
			/* Need to make a better choice on b */
			b = y + 5.*sqrt(sqrt(2.*STATIC_DELTAU));
	
			b = 25.;
		}
		else
		{/* STATIC_ETA == 1 */
			/* Need to make a better choice here */
			a =-25.;		
			b = 25.;
		}
		/* Determine case 1-4 from eta  , nu */ 

		
	
		if (0 <= STATIC_ETA && STATIC_ETA < .5)
		{/* Case 1:  No barrier b.c.'s */
			/* No barrier */
		
			M[i][j] = log(sm_qsimp(integrand_1  , a  ,b  ,EPS_SIMPSO*.0000001));
			
		}
		else 
		
		if (0 <= STATIC_ETA && STATIC_ETA < .5 )
		{/* Case 2: Reflecting b.c.'s
		            Need to differentiate this from Case 1 somehow  , 
			        since this case is currently not used   */
			M[i][j]= log(sm_qsimp(integrand_2  , a  ,b  , EPS_SIMPSO));
		}
		else
		/* Case 3 */
		if (.5 <= STATIC_ETA && STATIC_ETA < 1 )
		{ /* Natural b.c.'s */
			y = pow(fabs(1+STATIC_BETA*x_k)  ,1-nu)/(STATIC_BETA*(1-STATIC_ETA));
			M[i][j]= 
				log(sm_qsimp(integrand_3  , a  ,b  ,EPS_SIMPSO)
				+  
				exp(-STATIC_LAMBDA * (1./STATIC_BETA - STATIC_X_K)-gammln(nu)) *
				 (gammq(nu  ,(y*y)/(2.*STATIC_DELTAU)) - pow((y*y)/(2.*STATIC_DELTAU)  ,nu-1.)*
				 exp(-(y*y)/(2.*STATIC_DELTAU)))
				 );
		}
		else
		/* Case 4 */
		if (STATIC_ETA == 1)
		{
			M[i][j]= 
				log(sm_qsimp(integrand_4  , a  ,b  , EPS_SIMPSO));

		}

	return 0;

}

#endif

double print_M(double beta, double eta, double lambda, double deltau,
               double x_k) {
  /* Function used for debugging M only */

#if 0

double   **M;
double   x[100]  ,y[100];
Err      err=NULL;
int      i  ,j;
int      N;
double   y_b;
int      flag;
double   tmp1;

double   lambda_t  , lambda_T;
double   df;
double   price;

	lambda_t = beta;
	lambda_T = eta ;
	df       = lambda;


	/* print out the integrand for 1a&b */
	
	/* determine choice
	     0: value of distribution at y_bar(input using eta so this is only
											valid for the eta=0 case
		 1: value of M
	*/

	flag = (int) x_k2;
	x_k2 = 0.;

	beta   = beta2;
	eta    = eta2; 
	lambda = lambda2;
	deltau = deltau2;
	x_k    = x_k2;
	
	if (flag == 0)
	{/* print out distribution */
		y_b = eta;   
		eta = 0.;      /* change here to vary value of eta */
		val =  integrand_1(y_b); 
	}
	else
	if (flag == 1)
	{/* compute a value of M */
		i = 0; j = 0;
		M = dmatrix(0  ,1  ,0  ,1);   
		err = calc_bigM(M  , beta  , eta  , lambda  , deltau  , x_k  , i  , j);
		val = M[i][j];
		free_matrix(M  ,0  ,1  ,0  ,1);
	}
	return val;

#endif

#if 0 
	/* test the lq_inter_2d routine */

	i = 0; j = 0;
	N = 100;

	 M = dmatrix(0  ,N-1  ,0  ,N-1);   

	for (i=0;i<N;i++)
		for (j=0;j<N;j++)
			M[i][j] = (double) (.01*i)*(.1*j)*(.1*j);

	for (i=0;i<N;i++)
		x[i] = (double) .01*i;


	for (j=0;j<N;j++)
		y[j] = (double) .1*j;

	val = lq_inter_2d(beta2  ,eta2  ,x  ,y  ,M  ,100  ,100);

	return (val);

#endif

  return (0.);
}

#undef EPS_SIMPSO
