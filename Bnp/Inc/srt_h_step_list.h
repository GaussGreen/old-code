#ifndef SRT_H_STEP_LIST_H
#define SRT_H_STEP_LIST_H


/* -----------------------------------------------------------------------
  TYPE            :SrtStp,SrtStpPtr 
  AUTHOR          :E.Auld
  DESCRIPTION     :A step contains the information for events at that time and 
for the time slice that follows.  It is used for discretizing a stochastic
differential equation in time.
  DEFINITION      :

  ------------------------------------------------------------------------- */


typedef void    *SrtEvent;

typedef struct SrtStp{
/* Linked list definition: pointers to previous and next */
  struct   SrtStp    *prev;
  struct   SrtStp    *next;

 /* Time step information : index and date, time,... */
  long                index;                 /* INDEX OF NODE */
  long                date;                  /* INTEGER DATE AT START NODE */
  double              time;                  /* TIME IN YEARS FROM TODAY TO NODE */   
  double              ddate;                 /* DOUBLE DATE AT START NODE DTOL(date) */   
  double              delta_t;               /* TIME IN YEARS FROM THIS NODE TO NEXT */   
  double              sqrt_delta_t;          /* SQUARE ROOT OF delta_t */   

/* Market / model / deal specific information */
  SrtEvent            e;          /* financial event */
  void              **tminf;        /* mkt inf : one *tminf for each underlying */
  void               *trinf;         /* tree inf */
  
/* Underlyings correlation information */
  double            **correl;		/* The correlation matrix for the underlyings*/
  double            **coeff;		/* The linear coefficients used to correlate underlyings */
  int 	              und_num;       /* The number of underlying involved*/

} SrtStp,*SrtStpPtr;

/* ------------------------------------------------------------------------------ */
SrtStpPtr node_alloc();
/* allocates a node in the list */

/* ------------------------------------------------------------------------------ */
SrtStpPtr free_node(SrtStpPtr a);
/* frees a node in the list */
/* WARNING nodes should only be freed from the top down.  This function */
/* returns a pointer to the next node.  If it returns NULL then the list */
/* is gone */

/* ------------------------------------------------------------------------------ */
SrtStpPtr rem_node(SrtStpPtr a);
/* frees a, and updates a->next and a->prev to take this into account */
/* does not redo indexing of list */
/* returns a pointer to the next node.  If it returns NULL then the list */
/* is gone */

/* ------------------------------------------------------------------------------ */
SrtStpPtr rem_eventless_nodes(SrtStpPtr a, Date today);
/* removes all nodes from list that have NULL event field, except today  */
/* returns NULL if this is all the nodes, else top of what remains */
/* does not redo indexing of list */

/* ------------------------------------------------------------------------------ */
SrtStpPtr gototop(SrtStpPtr a);
/* goes to the top of the list  */

/* ------------------------------------------------------------------------------ */
SrtStpPtr gotobot(SrtStpPtr a);
/* goes to the bottom of the list  */

/* ------------------------------------------------------------------------------ */
SrtStpPtr free_list(SrtStpPtr a);
/* frees the list of which a is a part */

/* ------------------------------------------------------------------------------ */
SrtStpPtr insert_after(SrtStpPtr a);
/* Inserts a new node after the current one.  If NULL is passed creates a new*/
/* list. points to new node */

/* ------------------------------------------------------------------------------ */
SrtStpPtr insert_before(SrtStpPtr a);
/* Inserts a new node before the current one.  If NULL is passed creates a new*/
/* list. points to old node */

/* ------------------------------------------------------------------------------ */
long create_index(SrtStpPtr a);
/* creates an index for a.  useful for debugging.  returns the top index*/

/* ------------------------------------------------------------------------------ */
SrtStpPtr step_list_to_step_array(SrtStpPtr a);

/* ------------------------------------------------------------------------------ */

Err srtstpalloc(SrtStpPtr stp, size_t sz, int which, int index);
/* ------------------------------------------------------------------------------ */

Err srtstptminfalloc(SrtStpPtr stp, int size);

#endif
