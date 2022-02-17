#include <nag.h>	
#include "ICMKernel\util\icm_RootFinderND.h"
#include "ICMKernel\optim\uncoptim.h"
	
#include <nags.h>		
#include <nagg01.h>		
#include <nage04.h>	
#include <nagg05.h>	

# include <math.h>
# include <time.h>
#include "ICMKernel\random\icm_randomnag.h"
#include "ARMKernel\util\rand-gen.h"
#include "gpnumlib\levmarq.h"
#include "glob\icm_calibrator.h"

//pswarm ---------------------
#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <signal.h>
#include <setjmp.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>


//#include "pswarm_main.h"
#include "ICMKernel\optim\pswarm.h"


using namespace std;
//void* OptimTools::GansoLibContext = NULL;

void NAG_CALL FuncConfun(long n, long ncnlin, long needc[], double x[],
                         double conf[], double cjac[],
                         Nag_Comm* comm)
{
}

double  OptimTools::NagMultiOptimisator(void* context,
							  NAG_E04JBC_FUN f,
							  vector<double>& X,
							  vector<double>& bound_inf_X,
							  vector<double>& bound_sup_X,
							  double tol,
							  int maxiter,
							  double step_max,
							  double linesearch)
{
	double step_max2=step_max;
    double objf;
    Nag_BoundType bound;
   
    NagError fail ; 
    INIT_FAIL(fail); // 

       // Initialise options structure 
    Nag_E04_Opt options ;
    ICM_NAGSILENT(options) ;
    options.init_state=Nag_Init_None; 
    options.optim_tol=tol;
    options.max_iter=maxiter;
	options.linesearch_tol=linesearch;

	if (linesearch>1.)
		options.linesearch_tol = 0.999;
	else if (linesearch<0.)
		options.linesearch_tol = 0.;

	if (step_max2<tol) {step_max2=tol;}
	options.step_max=step_max2;
   
    Nag_Comm comm;
    comm.p = (Pointer)(context) ;

	bound = Nag_Bounds;

	vector<double> gvals;
	gvals.resize(X.size());
	
	//Opt routine
// FIXMEFRED: mig.vc8 (28/05/2007 10:53:37):cast
	e04jbc(X.size(), f, bound, &(*bound_inf_X.begin()), &(*bound_sup_X.begin()), &(*X.begin()), &objf,
		        &(*gvals.begin()), &options, &comm, &fail);

	ICM_NAGFREE(options); 
		
	return objf;
}


double  OptimTools::GansoMultiOptimisator(void* context,
							  OptFunction f,
							  std::vector<double>& X,
							  std::vector<double>& bound_inf_X,
							  std::vector<double>& bound_sup_X,
							  double tol,
							  int maxiter,
							  double step_max,
							  double linesearch)
{
/*
	double val = 0.;
	OptimTools::GansoLibContext = context;

	Ganso G;
	int retcode= G.MinimizeECAM(X.size(),X.begin(),&val,f,0,0,maxiter,
					NULL,NULL,NULL,NULL, bound_inf_X.begin(), bound_sup_X.begin(), NULL,5000,20);

 */
	return 0.;

}

double  OptimTools::UNCMultiOptimisator(void* context,
							  NAG_E04JBC_FUN f,
							  std::vector<double>& X,
							  std::vector<double>& bound_inf_X,
							  std::vector<double>& bound_sup_X,
							  double tol,
							  int maxiter,
							  double step_max,
							  double linesearch)
{
	double step_max2=step_max;
    double objf=0.;
    Nag_BoundType bound;
   
    NagError fail ; 
    INIT_FAIL(fail); // 

       // Initialise options structure 
    Nag_E04_Opt options ;
    ICM_NAGSILENT(options) ;
    options.init_state=Nag_Init_None; 
    options.optim_tol=tol;
    options.max_iter=maxiter;
	options.linesearch_tol=linesearch;

	if (linesearch>1.)
		options.linesearch_tol = 0.999;
	else if (linesearch<0.)
		options.linesearch_tol = 0.;

	if (step_max2<tol) {step_max2=tol;}
	options.step_max=step_max2;
   
    Nag_Comm comm;
    comm.p = (Pointer)(context) ;

	bound = Nag_Bounds;

	vector<double> gvals;
	gvals.resize(X.size());
	
	Opt O(X.size(),0);

// FIXMEFRED: mig.vc8 (28/05/2007 10:55:41):cast
	memcpy(O.x,&(*X.begin()),sizeof(double)*X.size());
	O.functionNag = f;
	O.Nagcomm = &comm;

	//Opt routine
	O.PRoptimize(0);

// FIXMEFRED: mig.vc8 (28/05/2007 10:55:54):cast
	memcpy(&(*X.begin()),O.x,sizeof(double)*X.size());

	ICM_NAGFREE(options); 
		
	return 0.;
}

void GradientLM(double *x,double *jac, int m, int n, void *adata)
{
	double f;
	OptimTools::LMfunction(x,&f,m,n,adata);

	double f0 = f;

	for (int i=0; i<n; i++)
	{	
		x[i] += GRAD_EPSILON;
		OptimTools::LMfunction(x,&f,m,n,adata);
		jac[i] = (f - f0)/GRAD_EPSILON;
		x[i] -= GRAD_EPSILON;
	}

//	*f = f0;
}

double  OptimTools::LMMultiOptimisator(void* context,
							  LMFunction f,
							  std::vector<double>& X,
							  std::vector<double>& bound_inf_X,
							  std::vector<double>& bound_sup_X,
							  int nbFuncs,
							  double tol,
							  int maxiter,
							  double step_max,
							  double linesearch)
{

	class ICM_LMCalibration : public ARM::ARM_LEVMARQFunc
	{
	private:
		void*	context;
		LMFunction* fun;

	public:

		ICM_LMCalibration(void* ctxt,LMFunction* f)
		{context = ctxt;fun=f;}

		~ICM_LMCalibration()
		{}

		void operator()(double p[], double hx[], int m, int n, void * adata = NULL) const
		{
			(*fun)(p,hx,m,n,adata);
		}
	};

	ICM_LMCalibration LMfunc(context,&f);
	ARM::ARM_GP_Vector X0(X.size(),X.begin());
	ARM::ARM_GP_Vector X0_INF(X.size(),bound_inf_X.begin());
	ARM::ARM_GP_Vector X0_SUP(X.size(),bound_sup_X.begin());

/*
	vector<double> foncs(nbFuncs);
	LMfunc(X.begin(),foncs.begin(), X.size(),nbFuncs,NULL);
	ARM::ARM_GP_Vector FX(nbFuncs,foncs.begin());
*/
	ARM::ARM_GP_Vector FX(nbFuncs);

	double opts[LM_OPTS_SZ], info[LM_INFO_SZ];

	opts[0]=LM_INIT_MU; opts[1]=1E-15; opts[2]=1E-15; opts[3]=1E-20;
	opts[4]=LM_DIFF_DELTA; // relevant only if the finite difference jacobian version is used 

	int status = ARM::LEVMARQMin(LMfunc,NULL,X0,FX,&X0_INF,&X0_SUP,info,maxiter,opts);
	
	double out = 0.;
	
	for (int i=0;i<nbFuncs;i++)
	{out+=ABS(FX[i]);}

	return (out);
}

// ------------------------------------------------------------------------
// optimizator BFGS NAG
// ------------------------------------------------------------------------

double  OptimTools::NagMultiOptimisator_BFGS(void* context,
										NAG_E04UNC_FUN func,
										std::vector<double>& X,
										std::vector<double>& bound_inf_X,
										std::vector<double>& bound_sup_X,
										int nbFuncs,
										double precision,
										int maxIter)
{
	int nbParams = X.size();

	double* x = &(*X.begin());
	double* LB = &(*bound_inf_X.begin());
	double* UB = &(*bound_sup_X.begin());

    long m, n, nclin, ncnlin, i, tda, tdfjac, nmax;
    double *xl=NULL, *xu=NULL;
    double IFAIL=0;

    m      = nbFuncs;
    n      = nbParams;
    ncnlin = 0;
    nclin  = 0;

    nmax   = n+nclin+ncnlin;

    static NagError fail, fail2;

    double *a = NULL, *y, *fjac, objf, *f;

    tda    = n;
    tdfjac = n;

    y = new double[m];
    f = new double[m];

    fjac = new double[m*n];

    for (i = 0; i < m; i++)
    { y[i] = 0.0; }

    xl = new double[nmax];
    xu = new double[nmax];

    for (i = 0; i < n; i++)
    {
        xl[i] = LB[i];
        xu[i] = UB[i];
    }

	Nag_E04_Opt options ;

	ICM_NAGSILENT(options); 
	options.max_iter = maxIter ;
	options.optim_tol = precision ;

	options.con_deriv = FALSE ;

	fail.print = TRUE;

    Nag_Comm comm;
    comm.p = (Pointer)(context) ;

    try
    {
        e04unc(m, n, nclin, ncnlin, a, tda, xl, xu,
               y, func,FuncConfun /*NULL*/, x, &objf, f, fjac,
               tdfjac, &options, &comm, &fail);
		ICM_NAGFREE(options);
    }
    catch(Exception& expt)
    {
        delete y;
        delete f;
        delete fjac;
        delete xl;
        delete xu;
     
        throw expt;
    }
    catch(...)
    {
        delete y;
        delete f;
        delete fjac;

        delete xl;
        delete xu;

        throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
             "Unrecognized Failure in : ARM_NAgOptimizationCall-->e04unc");
    }

	for (i=0;i<m;i++)
	{IFAIL += f[i]*f[i];}

	IFAIL = sqrt(IFAIL);

    delete y;
    delete f;
    delete fjac;

    delete xl;
    delete xu;

    if ( IFAIL == 1 )
    {IFAIL = RET_OK;}

    return(IFAIL);
}


void OptimTools::GenerateX0(ARM_MMTGenerator& RandGen,
										std::vector<double>& X0,
										std::vector<double>& bound_inf_X,
										std::vector<double>& bound_sup_X)
{
	double value=0.; 

	for (int i=0;i<X0.size();i++)
	{
	value = RandGen.Generate();
	X0[i] = bound_inf_X[i] + value * ( bound_sup_X[i] - bound_inf_X[i] ); 
	}
}

// ------------------------------------------------------------------------
// Generic optimizator
// ------------------------------------------------------------------------

double  OptimTools::MultiOptimisator(qOPTIMIZE_TYPE& type,
							  void* context,
							  void* f,
							  std::vector<double>& X,
							  std::vector<double>& bound_inf_X,
							  std::vector<double>& bound_sup_X,
							  int nbFuncs,
							  double tol,
							  int maxiter,
							  double step_max,
							  double linesearch)
{

	//OptimTools::GansoLibContext = context;
	double value = 0.;

	switch (type)
	{

		case qOPT_NAG_NOLIN_BC_BFGS :
			{
			value =  NagMultiOptimisator_BFGS(context,(NAG_E04UNC_FUN) f,
							  X,bound_inf_X,bound_sup_X,nbFuncs,tol,maxiter);
			break;
			}
		case qOPT_LEVMAR :
			{
			value =  LMMultiOptimisator(context,(LMFunction) f,
							  X,bound_inf_X,bound_sup_X,nbFuncs,
							  tol,maxiter,step_max,linesearch);
			break;
			}
		case qOPT_GLOBAL_UNCON :
			{
			value =  UNCMultiOptimisator(context,(NAG_E04JBC_FUN) f,
							  X,bound_inf_X,bound_sup_X,tol,
							  maxiter,step_max,linesearch);
			break;
			}
		case qOPT_GANSO :
			{
			value =  GansoMultiOptimisator(context,(OptFunction) f,
							  X,bound_inf_X,bound_sup_X,tol,
							  maxiter,step_max,linesearch);
			break;
			}
		case qOPT_PSWARM :
			{
			value =  PSwarmMultiOptimisator(context,(PSWGENFunction) f,
							  X,bound_inf_X,bound_sup_X,nbFuncs,tol,
							  maxiter,step_max,linesearch);
			break;
			}
		case qOPT_NAG_NOLIN_BC :  
		default :
			{
			value =  NagMultiOptimisator(context,(NAG_E04JBC_FUN) f,X,
							  bound_inf_X,bound_sup_X,tol,
							  maxiter,step_max,linesearch);
			break;
			}
	}

	return value;
}

// ------------------------------------------------------------------------
// Global alea search
// ------------------------------------------------------------------------

double OptimTools::AleaOptimSearch(qOPTIMIZE_TYPE& type,
							  void* context,
							  void* f,
							  std::vector<double>& X,
							  std::vector<double>& bound_inf_X,
							  std::vector<double>& bound_sup_X,
							  int nbFuncs,
							  long nbsimuls,
							  double maxtime,
							  double tol,
							  int maxiter,
							  double step_max,
							  double linesearch,
							  double TRACE)
{
	long i=0,j=0;
	long t=0;
	std::vector<double> X0 = X;
	std::vector<double> bound_inf_X0=bound_inf_X;
	std::vector<double> bound_sup_X0=bound_sup_X;
	double error =1.e20,erraux =0.;
	ARM_MMTGenerator	RandGen;
	time_t start_time, cur_time;
	time(&start_time);
	double coef=0.5;

	FILE *stream = NULL;

	if (TRACE)
	{stream = fopen("c:\\temp\\trace_optim.txt", "w+");}

	do
	{

		if (type != qOPT_PSWARM)
			OptimTools::GenerateX0(RandGen,X0,bound_inf_X,bound_sup_X);
		else
		{
		/*	X0=X;
			for (j=0;j<X.size();j++)
			{bound_inf_X0=(1.-VB_SGN(X0[j])*coef)*X0[j];
			bound_sup_X0=(1.+VB_SGN(X0[j])*coef)*X0[j];}
		*/
		}


		erraux = OptimTools::MultiOptimisator(type,context,f,X0,
							  bound_inf_X0,bound_sup_X0,nbFuncs,tol,maxiter,step_max,linesearch);
		if (erraux<error)
		{
			error = erraux;
			X = X0;

			if (TRACE){
			fprintf(stream,"\n");
			fprintf(stream,"N°SIM:%i\t, ERROR: %.1lf,\t",i,error);
			for (int y=0;y<X.size();y++)
			{fprintf(stream,"X[%i]= %1.5lf, ",y,X[y]);}
			}
		}

		i++;
		time(&cur_time);
	}
	while ((i<nbsimuls) && ((cur_time - start_time) < (maxtime*60.)) && (error>tol));

	if (TRACE)
	{fclose(stream);}

	return (error);
}

// ------------------------------------------------------------------------
// Ps Warm
// ------------------------------------------------------------------------

#ifdef MPE
/* for MPE */
int ComputeID_begin, ComputeID_end, SendID_begin, SendID_end, RecvID_begin,
  RecvID_end;
#endif

typedef struct { char *msg; } Exit_code;
extern "C" struct Stats stats;


static Exit_code exit_codes[] = {
  {"Normal exit"},                /* 0 */
  {"Abnormal exit"},              /* 1 */
  {"Failed to allocate memory"},  /* 2 */
  {"No simple bound and no initial guess"},  /* 3 */
};

char proc[]="00";
char resultfile[]="tipo00.xi2";
char datafile[]="dados00.xi2";


double objfun(int, double *, double *, double *);

#ifdef MPI
void MPI_objfun_deamon(int n, double *, double *, int);
#endif

extern "C" int PSwarm(double*,int, PSWGENFunction,double *, double *, double **, double *, double *);
extern "C" struct Options opt;

extern "C" void save_cache_file(int n, int m);
extern "C" void load_cache_file(int n, int m);
extern int read_cesam_cache(int n, double *x, double *age, 
			    double *teff, double *lum, double *r);
extern "C" void write_cesam_cache(int n, double *x, double *age,
			      double *teff, double *lum, double *r);

#ifndef AMPL

//void set_problem(double *x, double *lb, double *ub);
void set_problem_dimension(int *n);

#endif

//extern "C" double *scale;


#ifdef AMPL
static double objsign;
static fint NERROR = -1;
#define asl cur_ASL

char pswarm_version[]="PSwarm v1.0";

/* This struct member names must be in alphabetic order,
   for binary search */

keyword keywds[] = {
  KW("cognitial"   , pswarm_opt_d, (Char*)&opt.mu,
     "Cognitial parameter"),
  KW("fweight"   , pswarm_opt_d, (Char*)&opt.fweight,
     "Final weight (inercial parameter)"),
  KW("iweight"   , pswarm_opt_d, (Char*)&opt.iweight,
     "Initial weight (inercial parameter)"),
  KW("maxf", pswarm_opt_i, (Char*)&opt.maxf,
	 "Maximum number of function evaluations times problem dimension"),
  KW("maxit", pswarm_opt_i, (Char*)&opt.maxiter,
	 "Maximum number of iterations times problem dimension"),
  KW("size"        , pswarm_opt_i, (Char*)&opt.s,
     "Swarm size"),
  KW("social"        , pswarm_opt_d, (Char*)&opt.nu,
     "Social parameter"),
};


struct Option_Info Oinfo = { "Particle Swarm", "PSwarm", "pswarm_options",
			     keywds, nkeywds, 1, pswarm_version, 0, NULL};



/**********************************************
Set options. String type
**********************************************/
char *pswarm_opt_s(Option_Info *oi, keyword *kw, char *value)
{
  char *s;


  /* never echo options */
  oi->option_echo &= ~ASL_OI_echo; 


  if(!strncmp("method", kw->name, 6)){
    s=value;
    while(*s!=' ' && *s!=0)
      s++;
    if(s<=value)
      return value;

    if(!strncmp("disc_hett", value, 9)){
      *(int *)kw->info=0; /*DISC_METHOD;*/
      printf("Discretization method selected Hettich version\n");
      return s;
    }

    /* unknown method */
    return value;
  }


  /* not implemented option */
  return value;
}


/**********************************************
Set options. Integer type
**********************************************/
char *pswarm_opt_i(Option_Info *oi, keyword *kw, char *value)
{
  long optval;
  char *s;

  /* never echo options */
  oi->option_echo &= ~ASL_OI_echo;

  optval=strtol(value, &s, 10);
  if(s > value){
    /* existing integer number */
    *(int *)kw->info=(int)optval;
    printf("\nDefault option %s=%d changed\n", kw->name, *(int *)kw->info);
    return s;
  }

return value;
}


/**********************************************
Set options. Double type
**********************************************/
char *pswarm_opt_d(Option_Info *oi, keyword *kw, char *value)
{
  double optval;
  char *s;

  /* never echo options */
  oi->option_echo &= ~ASL_OI_echo;  

  optval=strtod(value, &s);

  if(s > value){
    /* existing double number */
    *(double *)kw->info=optval;
    printf("\nDefault option %s=%.6f changed\n", kw->name, *(double *)kw->info);
    return s;
  }

return value;
}

#endif /* AMPL */


static jmp_buf Jb;

void catchfpe(int n)
{
#ifdef AMPL
  report_where(asl);
#endif /* AMPL */
  printf("\nFloating point error.\n");
  fflush(stdout);
  longjmp(Jb,1);
}


double OptimTools::PSwarmMultiOptimisator(void* &context,
							  PSWGENFunction fu,
							  std::vector<double>& X,
							  std::vector<double>& bound_inf_X,
							  std::vector<double>& bound_sup_X,
							  int nbFuncs,
							  double tol,
							  int in_maxiter,
							  double step_max,
							  double linesearch)
{

  int exit_code, i;
  double *sol=NULL;
  double *f=NULL;
  double *lb, *ub;
  double minublb;

  opt.maxiter = in_maxiter;
  //opt.maxf = tol;

  int n=X.size();
  double *X0;
  
  if (!setjmp(Jb)){
    signal(SIGFPE, catchfpe);
  
  /* lower and upper bounds on variables */

  lb=(double*)malloc(n*sizeof(double));    
  ub=(double*)malloc(n*sizeof(double));
    
// FIXMEFRED: mig.vc8 (28/05/2007 15:07:46):cast
  memcpy(lb,&(*bound_inf_X.begin()),sizeof(double)*n);
  memcpy(ub,&(*bound_sup_X.begin()),sizeof(double)*n);

  X0=(double*)malloc(n*sizeof(double));
  memcpy(X0,&(*X.begin()),sizeof(double)*n);

  _context_cal* p = (_context_cal*)context ;
  p->m_scale = OptimTools::scale_=(double*)malloc(n*sizeof(double));
 
  minublb=ub[0]-lb[0];
  for(i=1;i<n;i++){
     if(minublb>ub[i]-lb[i])
	minublb=ub[i]-lb[i];
   }

  for(i=0;i<n;i++){
    OptimTools::scale_[i]=(ub[i]-lb[i])/minublb;
    ub[i]/=OptimTools::scale_[i];
    lb[i]/=OptimTools::scale_[i];
    if(X0)
	X0[i]/=OptimTools::scale_[i];
  }

  /*
  printf("Variables scaled by:\n");
  for(i=0;i<n;i++)
    printf("scale[%d]=%lf\n", i, scale[i]);
  */
    
  sprintf(proc, "00");
  sprintf(datafile, "dados00.ini");
  sprintf(resultfile, "tipo00.xi2");

  exit_code=PSwarm(OptimTools::scale_,n, (PSWGENFunction)fu, lb, ub, &sol, f, NULL); /*NULL); /*X0);*/
  //printf("\n%s\n", exit_codes[exit_code].msg);

  memcpy(&(*X.begin()),sol,sizeof(double)*n);
  }

  vector<double> XFin=X;
  for (i=0;i<X.size();i++)
  {XFin[i]=X[i]/OptimTools::scale_[i];}

  double final = fu(n,&(*XFin.begin()),lb,ub);

  if (lb) free(lb);
  if (ub) free(ub);
  if (OptimTools::scale_) free(OptimTools::scale_);

  return (final);
}


/*
double objfun(int n, double *y, double *lb, double *ub)
{
 
  int i;

  double *x;
  double fx;

  x=(double*)malloc(n*sizeof(double));
  
  stats.objfunctions++;

  for(i=0;i<n;i++){
	  if(y[i]<lb[i] || y[i]>ub[i]){

      return 1e2;
	  }
    x[i]=y[i]*scale[i];
  }

  fx=pow(x[0],2.0);

  free(x);
  return fx;
}
*/


void set_problem_dimension(int *n)
{
		/* USER PROBLEM DIMENSION */

  *n=6;
}




