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


extern struct swarm pop;
extern struct Options opt;
extern struct Stats stats;


struct poll_vector *D=NULL; /* linked list of poll vectors */
struct poll_vector *last_D=NULL;
struct poll_vector *TC=NULL; /* extra poll vectors for tangent cone */


/* performe a poll step for y[pi] position */
void pollstep(int n, int pi, double (*objf)(), double *lb, double *ub, struct poll_vector **last_success)
{
  int i;
  double *poll_point, fx, minfx;
  struct poll_vector *tmp, *minvector;

#ifdef MPI
  MPI_Status status;
  int process, MPI_numprocs, action, received, directions;
  struct poll_vector **vectors;
#endif



  poll_point=malloc(n*sizeof(double));

  if(!poll_point){
    printf("Unable to alocate memory in pattern.c\n");
    exit(1);
  }


  /* performe a poll step for each poll vector in D and TC */
  
  /* connect both linked lists */
  if(last_D)
    last_D->next=TC;

#ifdef MPI
  MPI_Comm_size(MPI_COMM_WORLD, &MPI_numprocs);

  vectors=malloc(MPI_numprocs*sizeof(struct poll_vector *));

  if(!vectors){
    printf("Unable to alocate memory in pattern.c\n");
    exit(1);
  }

  printf("Sending jobs in pattern.c\n");

#ifdef MPE
  MPE_Log_event(SendID_begin, 0, NULL);
#endif

  directions=0; /* how many poll_point do we have ? */
  action=1;
  tmp=D;
  for(process=1;process<MPI_numprocs && tmp;process++){
    for(i=0;i<n;i++)
      poll_point[i]=pop.y[pi*n+i]+pop.delta*tmp->vector[i];

    MPI_Send(&action, 1, MPI_INT, process, 99, MPI_COMM_WORLD);
    MPI_Send(poll_point, n, MPI_DOUBLE, process, 99, MPI_COMM_WORLD);
    vectors[process]=tmp; /* keep direction that each process has */
    tmp=tmp->next;
    directions++;
  }

#ifdef MPE
  MPE_Log_event(SendID_end, 0, NULL);
#endif

  received=0;
  minvector=NULL;
  minfx=Inf;
  while(tmp){ /* we have more trial points to assign */
#ifdef MPE
    MPE_Log_event(RecvID_begin, 0, NULL);
#endif
    MPI_Recv(&fx, 1, MPI_DOUBLE, MPI_ANY_SOURCE, 99, MPI_COMM_WORLD, &status);
    process=status.MPI_SOURCE;
    
#ifdef MPE
    MPE_Log_event(RecvID_end, 0, NULL);
#endif
    
    if(minfx>fx){
      minfx=fx;
      minvector=vectors[process]; /* we can break here in an oportunistic version */
    }
    received++;
    stats.objfunctions++;
    /* Sent another one to the same processor */
#ifdef MPE
    MPE_Log_event(SendID_begin, 0, NULL);
#endif
    MPI_Send(&action, 1, MPI_INT, process, 99, MPI_COMM_WORLD);
    for(i=0;i<n;i++)
      poll_point[i]=pop.y[pi*n+i]+pop.delta*tmp->vector[i];
    
    MPI_Send(poll_point, n, MPI_DOUBLE, process, 99, MPI_COMM_WORLD);
    vectors[process]=tmp;
#ifdef MPE
    MPE_Log_event(SendID_end, 0, NULL);
#endif
    tmp=tmp->next;
    directions++;
  }

#ifdef MPE
  MPE_Log_event(RecvID_begin, 0, NULL);
#endif
  while(received<directions){ /* wait for the remaining ones */
    MPI_Recv(&fx, 1, MPI_DOUBLE, MPI_ANY_SOURCE, 99, MPI_COMM_WORLD, &status);
    process=status.MPI_SOURCE;
    received++;
    if(minfx>fx){
      minfx=fx;
      minvector=vectors[process];
    }
    stats.objfunctions++;
  }
#ifdef MPE
  MPE_Log_event(RecvID_end, 0, NULL);
  
  MPE_Log_event(ComputeID_begin, 0, NULL);
#endif
  
  printf("Done reading jobs in pattern.c\n");
  
  
#else /* MPI */

  tmp=D;
  minvector=NULL;
  minfx=Inf;
  while(tmp){
    for(i=0;i<n;i++)
      poll_point[i]=pop.y[pi*n+i]+pop.delta*tmp->vector[i];
    fx=objf(n, poll_point, lb, ub);
    if(minfx>fx){
      minfx=fx;
      minvector=tmp; /* if oportunistic then break if fx < leader fx */
    }
    tmp=tmp->next;
  }
#endif

  
  if(pop.fy[pi]>minfx){ /* a successful poll point */
    stats.sucpollsteps++;
    for(i=0;i<n;i++)
      pop.y[pi*n+i]=pop.y[pi*n+i]+pop.delta*minvector->vector[i];
    pop.fy[pi]=minfx;
    if(*last_success==minvector){ /* last success vector twice, increase delta */
      pop.delta*=opt.idelta;
      //printf("Increasing delta in poll step %f\n", pop.delta);
    } else { /* last success vector is different */
      *last_success=minvector;
    }
  } else {
    pop.delta*=opt.ddelta;
    //printf("Decreasing delta\n");
    *last_success=NULL;
  }
  
  
  /* free variables */
  free(poll_point);
  
  /* disconnect lists */
  if(last_D)
    last_D->next=NULL;
  
}


void init_pattern(int n)
{
  
  init_D(n);
  //	print_D(n);
  
}

void clean_pattern()
{
  
  clean_D();
  
  
}


void clean_D()
{
  struct poll_vector *tmp1, *tmp2;
  
  tmp1=D;
  while(tmp1){
    tmp2=tmp1->next;
    free(tmp1->vector);
    free(tmp1);
    tmp1=tmp2;
  }
  
  D=NULL;
}




void init_D(int n)
{
  int i;
  struct poll_vector *tmp1, *tmp2;
  
  
  if(D) /* D already initialized */
    return;
  
  switch(opt.pollbasis){
  default:
    printf("\n Poll basis order not defined\n");
    printf("\n Using I -I order\n");
  case N2: /* I -I */
    D=malloc(sizeof(struct poll_vector));
    D->next=NULL;
    D->vector=malloc(n*sizeof(double));
    memset(D->vector, 0, n*sizeof(double));
    D->vector[0]=+1.0;
    tmp2=D;
    for(i=1;i<2*n;i++){
      tmp1=malloc(sizeof(struct poll_vector));
      tmp2->next=tmp1;
      tmp1->vector=malloc(n*sizeof(double));
      memset(tmp1->vector, 0, n*sizeof(double)); /* clear memory */
      if(i<n)
	tmp1->vector[i]=+1.0;
      else
	tmp1->vector[i-n]=-1.0;
      tmp1->next=NULL ;
      tmp2=tmp1;
    }
    tmp1=malloc(sizeof(struct poll_vector));
    tmp2->next=tmp1;
    tmp1->vector=malloc(n*sizeof(double));
    for(i=0;i<n;i++)
      tmp1->vector[i]=+1.0;
    tmp1->next=NULL;
    tmp2=tmp1;
    tmp1=malloc(sizeof(struct poll_vector));
    tmp2->next=tmp1;
    tmp1->vector=malloc(n*sizeof(double));
    for(i=0;i<n;i++)
      tmp1->vector[i]=-1.0;
    tmp1->next=NULL;
    tmp2=tmp1;
    last_D=tmp2;
  }
}


void print_D(n)
{
  struct poll_vector *tmp;
  
  tmp=D;
  while(tmp){
    print_poll_vector(n, tmp->vector);
    tmp=tmp->next;
  }
  
}


void print_poll_vector(int n, double *vector)
{
  int i;
  
  
  if(!vector)
    return;
  
  printf("D=(%.2f, ", vector[0]);
  for(i=1;i<n;i++)
    printf("%.2f ", vector[i]);
  printf(")\n");
  
}
