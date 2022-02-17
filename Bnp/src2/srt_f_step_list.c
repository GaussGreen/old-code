/* ============================================================================  
   FILENAME:       srt_f_step_list.c
   PURPOSE:        A step contains the information for events at that time and 
                   for the time slice that follows.  It is used for discretizing a stochastic
                   differential equation in time.

                   Basically, a SrtStrpPtr is a linked list of time steps
  ============================================================================ */
#include "math.h"
#include "srt_h_all.h"
#include "grf_h_public.h"

SrtStpPtr node_alloc()
/* allocates a node in the list */
{
    SrtStpPtr a;
    a=(SrtStpPtr)srt_calloc(1,sizeof(SrtStp));
    if (a==NULL){
	fprintf(stderr,"NODES FAILED TO ALLOCATE\n");
	exit(1);
    }
    a->index=(-1);
    a->date=(-1);
    a->time=(-1.0);
    a->ddate=(-1.0);
    a->delta_t=(-1.0);
    a->sqrt_delta_t=(-1.0);
    return(a) ;
}

SrtStpPtr free_node(SrtStpPtr a)
/* frees a node in the list */
/* WARNING nodes should only be freed from the top down.  This function */
/* returns a pointer to the next node.  If it returns NULL then the list */
/* is gone */
/* We terminate tminf[i] with 0 in the allocation process, see srtstptminfalloc() */

{
   SrtStpPtr b;
   int i = 0 ;
   b=a->next;
   if(a->e)
	srt_free(a->e);
   if(a->tminf)
   {
      while (a->tminf[i]!=0)
      {
        srt_free(a->tminf[i]);
        i++ ;
      }
      srt_free(a->tminf) ;
   }
   if(a->trinf)
	srt_free(a->trinf);

/* This free might be redundant with the free of the SrtCorrList,
	as a->correl... only points to the correl of the SrtCorrLst */

   if ((a->und_num !=0) && (a->correl != NULL))
       	free_dmatrix(a->correl, 0, a->und_num-1, 0, a->und_num-1);
   a->correl = NULL;	

   srt_free(a);
   return(b);
}

SrtStpPtr rem_node(SrtStpPtr a)
/* frees a, and updates a->next and a->prev to take this into account */
/* does not redo indexing of list */
/* returns a pointer to the next node.  If it returns NULL then the list */
/* is gone */
{
   SrtStpPtr b,c;
   b=a->next;
   c=a->prev;
   if(b)b->prev=c;
   if(c)c->next=b;
   if(c && b)
   {
     c->delta_t = b->time - c->time;
     c->sqrt_delta_t = sqrt(c->delta_t);
   }
   return free_node(a);
}

SrtStpPtr rem_eventless_nodes(SrtStpPtr a, Date today)
/* removes all nodes from list that have NULL event field, except today  */
/* returns NULL if this is all the nodes, else top of what remains */
{
  SrtStpPtr b=NULL;
  a = gototop(a);
  while(a)
  {
     if(!a->e && a->date != today)
     {
       a=rem_node(a);
     }
     else
     {
       if(!b)b=a;
       a=a->next;
     }
  }
  return b;
}

SrtStpPtr gototop(SrtStpPtr a)
/* goes to the top of the list  */
{
   SrtStpPtr b;
   for(b=a;b->prev!=NULL;b=b->prev);
   return(b);
}

SrtStpPtr gotobot(SrtStpPtr a)
/* goes to the bottom of the list  */
{
   SrtStpPtr b;
   for(b=a;b->next!=NULL;b=b->next);
   return(b);
}

/* ----------------------------------------------------------------- */
SrtStpPtr free_list(SrtStpPtr a)
/* frees the list of which a is a part */
{
   for(a=gototop(a);a!=NULL;a=free_node(a));
   return(NULL);
}

/* ----------------------------------------------------------------- */

SrtStpPtr insert_after(SrtStpPtr a)
/* Inserts a new node after the current one.
   If NULL is passed creates a new list. P
   Points to new node */
{
   SrtStpPtr b;
   if(a==NULL)
   {
       a=node_alloc();
       a->prev=NULL;
       a->next=NULL;
   }
   else 
   {
       b=node_alloc();
       b->prev=a;
       b->next=a->next;
       if(a->next!=NULL)
       {
	   (a->next)->prev=b;
       }
       a->next=b;
       a=b;
   }
   return(a);
}

/* ----------------------------------------------------------------- */

SrtStpPtr insert_before(SrtStpPtr a)
/* Inserts a new node before the current one.  
   If NULL is passed creates a new list.
   Points to old node */
{
   SrtStpPtr b;
   if(a==NULL)
   {
       a=node_alloc();
       a->prev=NULL;
       a->next=NULL;
   }
   else 
   {
       b=node_alloc();
       b->next=a;
       b->prev=a->prev;
       if(a->prev!=NULL)
       {
	   (a->prev)->next=b;
       }
       a->prev=b;
   }
   return(a);
}

/* ---------------------------------------------------------- */

long create_index(SrtStpPtr a)
/* creates an index for a.  useful for debugging.  returns the top index*/
{
   long i;
   SrtStpPtr b;
   for(i=0,b=gototop(a);b->next!=NULL;b->index=i,b=b->next,i++);
   b->index=i;
   return (i);
}

SrtStpPtr step_list_to_step_array(SrtStpPtr a)
/* used to test speed difference between arrays and lists (neglegible) */
/* puts the list a in contiguous cells; list structure is maintained */
/* this is rather dangerous to use if freeing is not watched over carefully*/
/* a is not freed; and should not be until after the list
  returned by this function is also freed */
{
   long i,lastindex = create_index(a);
   SrtStpPtr b,new;

   new = (SrtStpPtr)srt_malloc((lastindex+1)*sizeof(SrtStp));
   for(b = gototop(a),i=0;i<=lastindex;i++){
      new[i] = *b;
      b = b->next;
   }
   for(i=0;i<lastindex;i++){
      new[i].next = &new[i+1];
      new[i+1].prev = &new[i];
   }
   new[0].prev = NULL;
   new[lastindex].next = NULL;
   return new;
}


/** 
  calloc space sz to each element in stp, either
  stp -> tminf if which == 0 or
  stp -> trinf if which == 1
**/
Err srtstpalloc (SrtStpPtr stp, size_t sz, int which, int index)
{
	void *vptr;

	stp = gototop (stp);
	while (stp)
	{
		vptr = srt_calloc (1, sz);
		if(!vptr)
			return serror("SrtStpPtr alloc failed");
    
		if (which == 1)
		{
			stp -> trinf = vptr;
		}

		else
		{
			stp -> tminf[index] = vptr;
		}
		
		stp = stp -> next;
	}
	return NULL;
}

/* Allocates space for tminf (time info) = market information at each step of stp structure */
/* We need to allocate one extra tminf and terminate it with 0, otherwise free_node would not work */
/* Size corrsponds to the number of underlyings used: it will allocates space for the **tminf vector */
Err srtstptminfalloc (SrtStpPtr stp, int size)
{
	void *vptr;

	stp = gototop(stp);
	while (stp)
	{
		vptr = srt_calloc (size + 1, sizeof (void*));
		if (!vptr) 
			return serror("SrtStpPtr tminf alloc failed");
		stp -> tminf = vptr;
		stp -> tminf[size] = 0;
		stp = stp -> next;
	}
	return NULL;
}     

/* Allocate space for the correlation matrix in SrtStpPtr... */
/* und_num corresponds to the number of underlyings used:
	it will allocates space for the **correl matrix */
Err srtstpcorrelalloc(SrtStpPtr stp, int und_num)
{
  double **correl ;

  stp = gototop(stp) ;
  while (stp)
  {
    if (und_num == 0 )
	stp->correl = NULL;
    else
    {
	correl = dmatrix(0, und_num-1, 0, und_num-1)  ;
    	if (!correl) 
		return serror("SrtStpPtr correlation alloc fail") ;
    }
    stp->correl = correl ;
    stp = stp->next ;
  }
  return NULL ;
}     

/* ======================================================================== */

/* Do not need to initialise the new one, 
   but have to call srtstptminfalloc after for the tminf allocation*/

SrtStpPtr copy_srtstpptr_to_srtstpptr(SrtStpPtr stpr)
{
   SrtStpPtr cur_stpr, new_stp_cur=NULL, new_stpr = NULL;

   cur_stpr = gototop(stpr);
   while (cur_stpr)
   {
	new_stpr = insert_after(new_stpr);
   	new_stpr->index = cur_stpr->index;
   	new_stpr->date  = cur_stpr->date;
   	new_stpr->time = cur_stpr->time;
   	new_stpr->ddate = cur_stpr->ddate;
   	new_stpr->delta_t = cur_stpr->delta_t;
   	new_stpr->sqrt_delta_t = cur_stpr->sqrt_delta_t;
	cur_stpr = cur_stpr->next;
   }
   new_stpr = gototop(new_stpr);
   return(new_stpr);
}	

