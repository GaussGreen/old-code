#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<math.h>
#include<memory.h>

#ifdef MPI
#include "mpi.h"
#endif
#ifdef MPE
#include "mpe.h"
#endif
#include "pattern.h"
#include "pswarm.h"


#ifdef MPE
/* for MPE */
extern int ComputeID_begin, ComputeID_end, SendID_begin, SendID_end, RecvID_begin,
  RecvID_end;
#endif



void print_pop(int, int, int , struct swarm *);
void init_pop(int, int, struct swarm *, double *, double *, int, double *);
double projection(double, double, double);
void matlab_write_pop(int n, int gbest, int s, struct swarm *pop, int iter );
void print_best(double *scale,int n, int gbest, int s, struct swarm *pop, int iter );

extern void pollstep(int n, int pi, double (*objf)(), double *lb, double *ub, struct poll_vector **last_sucess);
extern void init_pattern(int);
extern void clean_pattern();

#if SYS_RANDOM!=1
static long rand_seed;

#define SHUFFLE 256		/* size of random array */
#define DBL_MIN 2.2250738585072014e-308

static long rand_seed;
static double randflt(long *);
static double resettable_randflt(long *rand_seed, int reset);
#endif




struct swarm pop;
struct Stats stats;

/* default options */
struct Options opt = { 
	42,      /* swarm size */ 
	0.5,     /* cognitial parameter */
	0.5,     /* social parameter */
	0.5,     /* maximum velocity factor */
	2000,    /* maximum of iterations */
	2000,    /* maximum of function evaluations */
	0.9,     /* initial weight */
	0.4,     /* final weight */
	0.5,     /* max norm 2 for gradient */
	10,      /* bound limit */
	1.0e-5,  /* tolerance for stopping criteria */
	2.0,     /* initial delta */
	2.0,     /* factor for initial delta */
	2,       /* increase delta by a factor of */
	0.5,     /* decrease delta by a factor of */
	N2
};

typedef double (*PSWGENFunction_)(int, double *, double *, double *);

/* Pattern Swarm algorithm */
int PSwarm(double *scale,int n, PSWGENFunction_ objf, double *lb, double *ub, double **sol,
	     double *f, double *x)
{
  int i,j, iter, gbest, success, actives, iterunsuc=0, process;
  double *maxv, maxnormv, normtmp, weight, mindelta, normtmp2;
  char *buff;
  time_t tt;
  FILE *fp;
  static struct poll_vector *last_success=NULL;

#ifdef MPI
  MPI_Status status;
  int received, particle, MPI_numprocs, action;
#endif
  
  /* Initial time */
  time(&tt);
  buff=ctime(&tt);
  printf("\nInitial time: %s\n",buff);

#ifdef MPI
  MPI_Comm_size(MPI_COMM_WORLD, &MPI_numprocs);
#endif
  
  
#if SYS_RANDOM==1
  /* seed random number with time */
  srand((unsigned int)tt);
#else
  rand_seed=(long)abs(tt);
  resettable_randflt (&rand_seed, 1); /* initialize random number generator */
  /* seed random number with time */
  srand((unsigned int)tt);
#endif
  
  
  /* checkup on variables bounds */
  for(i=0;i<n;i++){
    if((lb[i]<=-Inf && !x) || (ub[i]>=Inf && !x)){
      printf("Not all variables have finite bound and no initial guess given\n");
      printf("All variables must have finite simple bounds or an initial guess should be provided\n");
      return EXIT_INITIAL;
    }
  }
  
  
  /* allocate memory for swarm */
  pop.x     =malloc(opt.s*n*sizeof(double));
  pop.v     =malloc(opt.s*n*sizeof(double));
  pop.y     =malloc(opt.s*n*sizeof(double));
  pop.fx    =malloc(opt.s * sizeof(double));
  pop.fy    =malloc(opt.s * sizeof(double));
  pop.active=malloc(opt.s * sizeof(int));
  
  /* allocate memory for maximum velocity allowed */
  maxv=malloc(n*sizeof(double));
  
  if(!pop.x || !pop.v || !pop.y || !pop.fx || !pop.fy || !pop.active || ! maxv)
    return EXIT_MEM;
  
  
  /* initialize maximum velocity allowed and compute delta*/
  mindelta=Inf;
  for(j=0;j<n;j++){
    if(lb[j]>-Inf && ub[j]<Inf){
      if(mindelta>(ub[j]-lb[j]))
	mindelta=(ub[j]-lb[j]);
      maxv[j]=(ub[j]-lb[j])*opt.maxvfactor;
    } else {
      maxv[j]=Inf;
    }
  }
  
  if(mindelta>= Inf || mindelta<2*opt.tol)
    opt.delta=2*opt.tol;
  else opt.delta=mindelta/opt.fdelta;
  
  
  printf("Delta for pattern search: %f\n", opt.delta);
  
  
  /* initialize population */
  if(x){
    printf("Initial guess provided, including in initial population\n\n");
    init_pop(n, opt.s, &pop, lb, ub, 1, x);
  } else {
    init_pop(n, opt.s, &pop, lb, ub, 0, NULL);
  }
  
  actives=opt.s;
  
  
  
  
  iter=0;
  stats.pollsteps=0;
  stats.sucpollsteps=0;
  stats.objfunctions=0;
  gbest=0; /* global best */
  maxnormv=Inf; /* don't stop in first iteration */
  
  init_pattern(n); /* initialize pattern search */
  
  /* Main cycle */
  while(iter < opt.maxiter && stats.objfunctions < opt.maxf){
    
	  print_best(scale,n, gbest, opt.s, &pop, iter);

    if(maxnormv<opt.tol && pop.delta<opt.tol){
      printf("\n\nStopping due to velocity and tolerance\n\n");
      break;
    }
    
    if(actives<=1 && pop.delta<opt.tol){
      printf("\n\nStopping due to single particle and tolerance\n\n");
      break;
    }
    
    iter++;
    
    
    /* inertia factor is a linear interpolation from iweight to fweight */
    weight = opt.iweight - (opt.iweight-opt.fweight)*((double)(iter))/((double)opt.maxiter);
    
    success=0; /* controls if gbest was updated with success */

#ifdef MPI    
    /* send jobs for computing objective functions */
    printf("Sending jobs in pswarm.c\n");
    
    i=0;
    action=2; /* compute objective value */
    received=0;
    while(!pop.active[i])i++; /* first active particle */
#ifdef MPE
    MPE_Log_event(SendID_begin, 0, NULL);
#endif
    for(process=1;process<MPI_numprocs && i<opt.s;process++){
      MPI_Send(&action, 1, MPI_INT, process, 99, MPI_COMM_WORLD);
      MPI_Send(&i, 1, MPI_INT, process, 99, MPI_COMM_WORLD);
      MPI_Send(&pop.x[i*n], n, MPI_DOUBLE, process, 99, MPI_COMM_WORLD);
      i++;
      while(!pop.active[i])i++; /* next active particle */
    }
#ifdef MPE
    MPE_Log_event(SendID_end, 0, NULL);
#endif


    while(i<opt.s){ /* we have more particle to assign */
      /* receive one objective function */
#ifdef MPE
      MPE_Log_event(RecvID_begin, 0, NULL);
#endif
      MPI_Recv(&particle, 1, MPI_INT, MPI_ANY_SOURCE, 99, MPI_COMM_WORLD, &status);
      process=status.MPI_SOURCE;
      MPI_Recv(&pop.fx[particle], 1, MPI_DOUBLE, process, 99, MPI_COMM_WORLD, &status);
#ifdef MPE
      MPE_Log_event(RecvID_end, 0, NULL);
#endif
      received++;
      stats.objfunctions++;
      /* Sent another one to the same processor */
#ifdef MPE
      MPE_Log_event(SendID_begin, 0, NULL);
#endif
      MPI_Send(&action, 1, MPI_INT, process, 99, MPI_COMM_WORLD);
      MPI_Send(&i, 1, MPI_INT, process, 99, MPI_COMM_WORLD);
      MPI_Send(&pop.x[i*n], n, MPI_DOUBLE, process, 99, MPI_COMM_WORLD);
      i++;
#ifdef MPE
      MPE_Log_event(SendID_end, 0, NULL);
#endif
      while(!pop.active[i])i++; /* next active particle */
    }

#ifdef MPE
    MPE_Log_event(RecvID_begin, 0, NULL);
#endif
    while(received<actives){ /* wait for the remaining ones */
      MPI_Recv(&particle, 1, MPI_INT, MPI_ANY_SOURCE, 99, MPI_COMM_WORLD, &status);
      process=status.MPI_SOURCE;
      MPI_Recv(&pop.fx[particle], 1, MPI_DOUBLE, status.MPI_SOURCE, 99, MPI_COMM_WORLD, &status);
      received++;
      stats.objfunctions++;
    }
#ifdef MPE
    MPE_Log_event(RecvID_end, 0, NULL);

    MPE_Log_event(ComputeID_begin, 0, NULL);
#endif
    
    printf("Done sending jobs in pswarm.c\n");
    
    /* All jobs sent */
#else /* MPI */
    /* no MPI, so just compute all the objective function values */

    for(i=0;i<opt.s;i++){
      if(pop.active[i]){
	pop.fx[i]=objf(n, &pop.x[i*n], lb, ub);
      }
    }

#endif

    process=0;
    for(i=0;i<opt.s;i++){
      if(pop.active[i]){
	if(pop.fy[i]>pop.fx[i]){      /* progress obtained */
	  pop.fy[i]=pop.fx[i];           /* Search step */
	  memcpy(&pop.y[i*n], &pop.x[i*n], n*sizeof(double));
	  
	  /* check if a new leader is available or if a progress was
	     obtained on the leader */
	  if(pop.fy[gbest]>pop.fy[i] || gbest==i){
	    gbest=i; /* global best indice */
	    success=1; /* success for leader obtained */
	    last_success=NULL; /* reset successful direction on poll step */
	  }
	}
      }
    }

    
    
    if(!success){ /* no success for the gbest particle in one generation, so performe a poll step */
      if(pop.delta>=opt.tol){

#ifdef MPE
	MPE_Log_event(ComputeID_end, 0, NULL); /* we log in pattern.c */
#endif
	/* performe a poll step, update y and delta */
	pollstep(n, gbest, objf, lb, ub, &last_success);

#ifdef MPE
	MPE_Log_event(ComputeID_begin, 0, NULL);
#endif

	stats.pollsteps++;
	iterunsuc=0;
      } else {
	iterunsuc++;
	//printf("Consecutive unsuccesseful iterations %d\n", iterunsuc);
      }
    } else { /* success for the gbest particle */
      iterunsuc=0;
      //printf("Success in Swarm iteration\n");
      /* increase delta */
      if(pop.delta<opt.delta){
	pop.delta*=opt.idelta;
	//	printf("Increasing delta in search step\n");
      }
      /* allow at least one more iteration */
      if(pop.delta<opt.tol)
	pop.delta=2*opt.tol;
    }
    
    
    
    for(i=0;i<opt.s;i++){ /* for each particle */
      
      if(pop.active[i]){ /* active particle */
	/* update velocity */
	for(j=0;j<n;j++){
	  pop.v[i*n+j] =
#if SYS_RANDOM==1
	    projection(weight*pop.v[i*n+j]+
		       opt.mu*(rand()/(RAND_MAX+1.0))*(pop.y[i*n+j]-pop.x[i*n+j])
		       +opt.nu*(rand()/(RAND_MAX+1.0))*(pop.y[gbest*n+j]-pop.x[i*n+j]),
		       -maxv[j],maxv[j]);
#else
	  projection(weight*pop.v[i*n+j]+
		     opt.mu*(randflt (&rand_seed))*(pop.y[i*n+j]-pop.x[i*n+j])
		     +opt.nu*(randflt (&rand_seed))*(pop.y[gbest*n+j]-pop.x[i*n+j]),
		     -maxv[j],maxv[j]);
#endif
	  
	}
	
	
	/* update particle and check bound limits */
	for(j=0;j<n;j++){
	  pop.x[i*n+j]=projection(pop.x[i*n+j]+pop.v[i*n+j],lb[j],ub[j]);
	}
      }
    }
    
    
    /* check for all norm velocities to zero */
    
    /* first for gbest */
    for(j=0;j<n;j++)
      normtmp+=pow(pop.v[gbest*n+j],2.0);
    maxnormv=sqrt(normtmp);
    
    
    /* remove particle close to gbest and compute maximum velocity */
    actives=0;
    for(i=0;i<opt.s;i++){ /* for each particle */
      if(pop.active[i] && i!=gbest){ /* active particle and not the gbest */
	normtmp=0.0; normtmp2=0.0;
	for(j=0;j<n;j++){
	  normtmp+=pow(pop.y[i*n+j]-pop.y[gbest*n+j],2.0);
	  normtmp2+=pow(pop.v[i*n+j],2.0);
	}
	normtmp=sqrt(normtmp); normtmp2=sqrt(normtmp2);
	if(normtmp<opt.delta && normtmp2<opt.delta){ //(fabs((double)(iter-(opt.maxiter)/100.0)))){
	  pop.active[i]--; /* remove particle */
	  // printf("Particle %d inactive iter=%d\n", i, iter);
	} else { /* particle not removed, so account for maxnormv */
	  if(maxnormv<normtmp2)
	    maxnormv=normtmp2;
	}
      }
      if(pop.active[i])
	actives++; /* account in actives */
    }
    
    //    printf("Maximum velocity norm: %f\n", maxnormv);
    
    //printf("%d;%.20f\n",stats.objfunctions,pop.fy[gbest]);

#ifdef MPE
    MPE_Log_event(ComputeID_end, 0, NULL);
#endif
  }
  
  if (iter >= opt.maxiter){
    printf("\n\nStopping due to maximum number of iterations reached\n\n");
  }
  
  if(stats.objfunctions >= opt.maxf){
    printf("\n\nStopping due to maximum number of function evaluations reached\n\n");
  }
  
  print_pop(n, gbest, opt.s, &pop);
  print_best(scale,n, gbest, opt.s, &pop, iter);
  
  // matlab_write_pop(n, gbest, opt.s, &pop, iter);
  
  printf("maxnormv=%.20f\n",maxnormv);
  
  printf("delta=%.20f\n",pop.delta);
  
  printf("%d iterations\n", iter);
  
  printf("%d function evaluations\n", stats.objfunctions);
  
  printf("%d poll steps performed\n", stats.pollsteps);
  
  printf("%d poll steps performed with success\n", stats.sucpollsteps);
  
  printf("%d & %d & %d & %d & %.4f\n",iter, stats.objfunctions, stats.pollsteps, stats.sucpollsteps, pop.fy[gbest]);
  
  fp = fopen("./results", "a+");
  if(fp){
    fprintf(fp, " & %d & %d & %d & %d & %.2f & %.6f\\\\\n", n, stats.objfunctions,
	    iter, stats.pollsteps,
	    (double)(100*(double)stats.sucpollsteps/(double)stats.pollsteps),
	    pop.fy[gbest]);
    fclose(fp);
  }
  
  clean_pattern();
  
  time(&tt);
  buff=ctime(&tt);
  printf("Final time: %s\n", buff);
  
  /* returning the best of the population */
  *sol=malloc(n*sizeof(double));

  if(!(*sol))
    return EXIT_MEM;
  
  for (i=0;i<n;i++)
  {
	  (*sol)[i]=scale[i]*pop.y[gbest*n+i];
  }
  //memcpy(*sol,&pop.y[gbest*n],n*sizeof(double));
  
  /* free allocated memory */
  free(pop.x);
  free(pop.y);
  free(pop.v);
  free(pop.fx);
  free(pop.fy);
  free(pop.active);
  free(maxv);
  
  return EXIT_OK;
  
}



void init_pop(int n, int s, struct swarm *pop, double *lb, double *ub,
          int ninitials, double *initials)
{
  int i, j;
  double normtmp=10.0; /* should never be used by default */

  /*   The user can provide an initial guess to include in the initial swarm
       A reset in the population can also proposed a fixed number of point
     to be in the next swarm */

  /* Do simple check in the simple bound limits */

  if(ninitials && initials){
	  for(i=0;i<ninitials;i++){
		  //printf("Init %d: %f\n", i, objfun(n,&initials[i*n]));
		  for(j=0;j<n;j++)
			  pop->x[i*n+j]=projection(initials[i*n+j],lb[j],ub[j]);
		  pop->fy[i]=+Inf*10; /* in first iteration, y will be set */
		  pop->active[i]=1;   /* chances to be near the gbest */
	  }
  }
  
  if(ninitials){
    /* compute standard deviation of first particle */
    normtmp=0.0;
    for(j=0;j<n;j++)
      normtmp+=pow(initials[j],2.0);
    if(normtmp<10)
      normtmp=opt.blim;
  }
      
  
  for(i=ninitials;i<s;i++){ /* for all remaining particle   */
	  for(j=0;j<n;j++){   /* for all dimensions */
		  if(lb[j]>-Inf && ub[j]<Inf)
#if SYS_RANDOM==1
			  pop->x[i*n+j]=(rand()/(RAND_MAX+1.0))*(ub[j]-lb[j])+lb[j];
#else
			  pop->x[i*n+j]=(randflt (&rand_seed))*(ub[j]-lb[j])+lb[j];
#endif
		  else {
			  if(lb[j]<=-Inf && ub[j]>=Inf){ /* both limits infinite */
#if SYS_RANDOM==1
				  pop->x[i*n+j]=2*(rand()/(RAND_MAX+1.0)-0.5)*normtmp+initials[j];
#else
				  pop->x[i*n+j]=2*(randflt (&rand_seed)-0.5)*normtmp+initials[j];
#endif
				  printf("x[%d]=%f\n",j,pop->x[i*n+j]);
			  } else {
				  if(lb[j]<=-Inf){ /* lower infinite and upper finite */
#if SYS_RANDOM==1
					  pop->x[i*n+j]=2*(rand()/(RAND_MAX+1.0)-0.5)*(ub[j]-initials[j])+initials[j];
#else
					  pop->x[i*n+j]=2*(randflt (&rand_seed)-0.5)*(ub[j]-initials[j])+initials[j];
#endif
				  } else { /* upper infinite and lower finite */
#if SYS_RANDOM==1
					  pop->x[i*n+j]=2*(rand()/(RAND_MAX+1.0)-0.5)*(initials[j]-lb[j])+initials[j];
#else
					  pop->x[i*n+j]=2*(randflt (&rand_seed)-0.5)*(initials[j]-lb[j])+initials[j];
#endif
				  }
			  }
    
      }
    }
    pop->fy[i]=+Inf*10; /* in first iteration, y will be set */
    pop->active[i]=1;   /* chances to be near the gbest */
  }

  pop->delta=opt.delta;
  /* initialize  velocity */
  memset(pop->v, 0, s*n*sizeof(double));
}


double projection(double xi, double lbi, double ubi)
{
  if(xi<lbi)
    return lbi;
  if(xi>ubi)
    return ubi;
  return xi;
}


/* print the best of each particle in a population */
void print_pop(int n, int gbest, int s, struct swarm *pop)
{
  int i, j, inactive;

  printf("Printing the best so far for each particle\n");

  inactive=0;
  for(i=0;i<s;i++){ /* for each particle */
	  if(pop->active[i]){ /* active particle */
		printf("x(%d)=[", i);
		for(j=0;j<n-1;j++) /* for each coordinate */
			printf("%.4f,", pop->x[i*n+j]);
		printf("%.4f];\n", pop->x[i*n+n-1]);
		printf("y(%d)=[", i);
		for(j=0;j<n-1;j++) /* for each coordinate */
			printf("%.4f,", pop->y[i*n+j]);
		printf("%.4f];\n", pop->y[i*n+n-1]);
		printf("v(%d)=[", i);
		for(j=0;j<n-1;j++) /* for each coordinate */
			printf("%.4f,", pop->v[i*n+j]);
		printf("%.4f];\n", pop->v[i*n+n-1]);
		printf("f(%d)=%.20f\n", i, pop->fy[i]);
	  } else {
		  inactive++;
	  }
  }

  printf("%d inactive particles\n", inactive);
    
  printf("\n The very best\n");
  printf("p(%d)=[", gbest);
  for(j=0;j<n-1;j++) /* for each coordinate */
    printf("%.10f,", pop->y[gbest*n+j]);
  printf("%.10f];\n", pop->y[gbest*n+n-1]);
  printf("f(%d)=%.10f\n", gbest, pop->fy[gbest]);

  
}

/* write the best of each particle in a population */
void matlab_write_pop(int n, int gbest, int s, struct swarm *pop, int iter)
{
  int i;
  FILE *mfile;

  if(n!=2)return;

  if(iter==1)mfile=fopen("pop.m", "w");
  else mfile=fopen("pop.m", "a");

  if(!mfile)return;

  fprintf(mfile,"xa1=[");
  for(i=0;i<s;i++){ /* for each particle */
	  if(pop->active[i])
		  fprintf(mfile,"%.2f,", pop->y[i*n]);
  }
  fprintf(mfile,"];");

  fprintf(mfile,"xa2=[");
  for(i=0;i<s;i++){ /* for each particle */
	  if(pop->active[i])
		  fprintf(mfile,"%.2f,", pop->y[i*n+1]);
  }
  fprintf(mfile,"];");
  fprintf(mfile,"hold off;\nir2;\nhold on;\nplot(xa1,xa2,'or');\n");
//fprintf(mfile,"plot(xa1,xa2,'ob');\n");

//  fprintf(mfile,"xb1=[");
//  for(i=0;i<s;i++){ /* for each particle */
//	  if(!pop->active[i])
//		  fprintf(mfile,"%.2f,", pop->y[i*n]);
//  }
//  fprintf(mfile,"];");

//  fprintf(mfile,"xb2=[");
//  for(i=0;i<s;i++){ /* for each particle */
//	  if(!pop->active[i])
//		  fprintf(mfile,"%.2f,", pop->y[i*n+1]);
//  }
//  fprintf(mfile,"];");
//  fprintf(mfile,"hold off;\nir2;\nhold on;\nplot(xa1,xa2,'or');\nplot(xb1,xb2,'ob');\n");
//  fprintf(mfile,"plot(xb1,xb2,'or');\n");
  fprintf(mfile,"title('iter=%d, best fx=%.4f, pollsteps=%d, suc=%d, delta=%.8f nfx=%d');\npause;\n", iter, pop->fy[gbest],stats.pollsteps, stats.sucpollsteps, pop->delta, stats.objfunctions);

  fclose(mfile);

}


/* write the best swarm particle */
void print_best(double *scale,int n, int gbest, int s, struct swarm *pop, int iter)
{
//  double *scale=NULL;

  int i, nactive;
  FILE *file;


  if(iter==1)file=fopen("results.txt", "w");
  else file=fopen("results.txt", "a");

  nactive=0;
  for(i=0;i<s;i++)
    if(pop->active[i])nactive++;

  if(!file)return;

  fprintf(file,"x=[");
  for(i=0;i<n;i++){ /* for each dimension */
	  fprintf(file,"%.8f,", scale[i]*pop->y[gbest*n+i]);
  }
  fprintf(file,"]  fx=%lf\n", pop->fy[gbest]);

  fprintf(file,"Nobj=%d  Npoll=%d  Nsucpoll=%d Active=%d\n", stats.objfunctions, stats.pollsteps, stats.sucpollsteps, nactive);

  fclose(file);

}




#if SYS_RANDOM!=1


#define LONG_INT long

/*
  The next two functions, myrand and randflt, were copied from
  user.c in ASA.
*/

#define MULT ((LONG_INT) 25173)
#define MOD ((LONG_INT) 65536)
#define INCR ((LONG_INT) 13849)
#define FMOD ((double) 65536.0)

/***********************************************************************
* double myrand - returns random number between 0 and 1
*	This routine returns the random number generator between 0 and 1
***********************************************************************/

static double
myrand (LONG_INT * rand_seed)
{
#if TRUE			/* (change to FALSE for alternative RNG) */
  *rand_seed = (LONG_INT) ((MULT * (*rand_seed) + INCR) % MOD);
  return ((double) (*rand_seed) / FMOD);
#else
  /* See "Random Number Generators: Good Ones Are Hard To Find,"
     Park & Miller, CACM 31 (10) (October 1988) pp. 1192-1201.
     ***********************************************************
     THIS IMPLEMENTATION REQUIRES AT LEAST 32 BIT INTEGERS
     *********************************************************** */
#define _A_MULTIPLIER  16807L
#define _M_MODULUS     2147483647L	/* (2**31)-1 */
#define _Q_QUOTIENT    127773L	/* 2147483647 / 16807 */
#define _R_REMAINDER   2836L	/* 2147483647 % 16807 */
  long lo;
  long hi;
  long test;

  hi = *rand_seed / _Q_QUOTIENT;
  lo = *rand_seed % _Q_QUOTIENT;
  test = _A_MULTIPLIER * lo - _R_REMAINDER * hi;
  if (test > 0)
    {
      *rand_seed = test;
    }
  else
    {
      *rand_seed = test + _M_MODULUS;
    }
  return ((double) *rand_seed / _M_MODULUS);
#endif /* alternative RNG */
}

/***********************************************************************
* double randflt
***********************************************************************/

static double
randflt (LONG_INT *rand_seed)
{
  return (resettable_randflt (rand_seed, 0));
}

/***********************************************************************
* double resettable_randflt
***********************************************************************/
static double
resettable_randflt (LONG_INT * rand_seed, int reset)

  /* shuffles random numbers in random_array[SHUFFLE] array */

{

  /* This RNG is a modified algorithm of that presented in
   * %A K. Binder
   * %A D. Stauffer
   * %T A simple introduction to Monte Carlo simulations and some
   *    specialized topics
   * %B Applications of the Monte Carlo Method in statistical physics
   * %E K. Binder
   * %I Springer-Verlag
   * %C Berlin
   * %D 1985
   * %P 1-36
   * where it is stated that such algorithms have been found to be
   * quite satisfactory in many statistical physics applications. */

  double rranf;
  unsigned kranf;
  int n;
  static int randflt_initial_flag = 0;
  LONG_INT initial_seed;
  static double random_array[SHUFFLE];	/* random variables */

  if (*rand_seed < 0)
    *rand_seed = -*rand_seed;

  if ((randflt_initial_flag == 0) || reset)
    {
      initial_seed = *rand_seed;

      for (n = 0; n < SHUFFLE; ++n)
	random_array[n] = myrand (&initial_seed);

      randflt_initial_flag = 1;

      for (n = 0; n < 1000; ++n)	/* warm up random generator */
	rranf = randflt (&initial_seed);

      rranf = randflt (rand_seed);

      return (rranf);
    }

  kranf = (unsigned) (myrand (rand_seed) * SHUFFLE) % SHUFFLE;
  rranf = *(random_array + kranf);
  *(random_array + kranf) = myrand (rand_seed);

  return (rranf);
}
#endif
