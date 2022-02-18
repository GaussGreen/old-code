#include "grf_h_all.h"
#include "stdio.h"
#include "stdlib.h"

/*<%%STA>-----------------------------------------------------------------
  DESCRIPTION     :Utility functions for dealing with a COMLL_PTR linked list

<%%END>---------------------------------------------------------------------*/

/*<%%STA>-----------------------------------------------------------------
  FUNCNAME        :COMLL_PTR comll_node_alloc()
  DESCRIPTION     :returns pointer to new COMLL_STR, or NULL if it fails.
Everything in the node will be set to zero except the index, which will
be set to -1.
<%%END>---------------------------------------------------------------------*/

COMLL_PTR comll_node_alloc()
{
    COMLL_PTR a = NULL;

    a = (COMLL_PTR)srt_calloc(1, sizeof(COMLL_STR));
    if (a == NULL)
    {
        return NULL;
    }
    a->index = (-1);

    return (a);
}

/*<%%STA>-----------------------------------------------------------------
  FUNCNAME        :COMLL_PTR comll_free_node(COMLL_PTR a)
  DESCRIPTION     :frees a node in the list
WARNING nodes should only be freed from the top down.  This function
returns a pointer to the next node.  If it returns NULL then the list
is gone. Frees the cvg and dfind arrays if they are non NULL.
<%%END>---------------------------------------------------------------------*/

COMLL_PTR comll_free_node(COMLL_PTR a)
{
    COMLL_PTR b = NULL;
    int       i = 0;

    b = a->next;

    if (a->ptr)
        srt_free(a->ptr);

    if (a->cvg)
        srt_free(a->cvg);

    if (a->cvgfloat)
        srt_free(a->cvgfloat);

    if (a->dfind != NULL)
        srt_free(a->dfind);

    if (a->dfindfloat != NULL)
        srt_free(a->dfindfloat);

    if (a->spread != NULL)
        srt_free(a->spread);

    if (a->dstore != NULL)
        srt_free(a->dstore);

    if (a->sstore != NULL)
        srt_free(a->sstore);

    srt_free(a);

    return (b);
}

/*<%%STA>-----------------------------------------------------------------
  FUNCNAME        :COMLL_PTR comll_gototop(COMLL_PTR a)
  DESCRIPTION     :returns top of list, or NULL.
<%%END>---------------------------------------------------------------------*/
COMLL_PTR comll_gototop(COMLL_PTR a)
{
    COMLL_PTR b = NULL;

    if (!a)
        return NULL;

    for (b = a; b->prev != NULL; b = b->prev)
        ;

    return (b);
}

/*<%%STA>-----------------------------------------------------------------
  FUNCNAME        :COMLL_PTR comll_gotobot(COMLL_PTR a)
  DESCRIPTION     :returns bottom of list, or NULL.
<%%END>---------------------------------------------------------------------*/
COMLL_PTR comll_gotobot(COMLL_PTR a)
{
    COMLL_PTR b = NULL;

    if (!a)
        return NULL;

    for (b = a; b->next != NULL; b = b->next)
        ;

    return (b);
}

/*<%%STA>-----------------------------------------------------------------
  FUNCNAME        :COMLL_PTR comll_free_list(COMLL_PTR a)
  DESCRIPTION     :if(a) frees a using comll_gototop() and comll_free_node(),
<%%END>---------------------------------------------------------------------*/
COMLL_PTR comll_free_list(COMLL_PTR a)
{
    for (a = comll_gototop(a); a != NULL; a = comll_free_node(a))
        ;

    return (NULL);
}

/*<%%STA>-----------------------------------------------------------------
  FUNCNAME        :COMLL_PTR comll_insert_after(COMLL_PTR a)
  DESCRIPTION     :Inserts a new node after the current one.  If NULL is
passed creates a new list. Returns pointer to new node.  Returns NULL if fails.
<%%END>---------------------------------------------------------------------*/
COMLL_PTR comll_insert_after(COMLL_PTR a)
{
    COMLL_PTR b = NULL;

    if (a == NULL)
    {
        a = comll_node_alloc();
        if (!a)
            return NULL;

        a->prev = NULL;
        a->next = NULL;
    }
    else
    {
        b = comll_node_alloc();
        if (!b)
            return NULL;

        b->prev = a;
        b->next = a->next;

        if (a->next != NULL)
        {
            (a->next)->prev = b;
        }

        a->next = b;
        a       = b;
    }

    return (a);
}

/*<%%STA>-----------------------------------------------------------------
  FUNCNAME        :COMLL_PTR comll_insert_before(COMLL_PTR a)
  DESCRIPTION     :Inserts a new node before the current one.  If NULL is
passed creates a new list. Returns pointer to new node.  Returns NULL if fails.
<%%END>---------------------------------------------------------------------*/
COMLL_PTR comll_insert_before(COMLL_PTR a)
{
    COMLL_PTR b = NULL;

    if (a == NULL)
    {
        b = comll_node_alloc();
        if (!b)
            return NULL;

        b->prev = NULL;
        b->next = NULL;
    }
    else
    {
        b = comll_node_alloc();
        if (!b)
            return NULL;

        b->next = a;
        b->prev = a->prev;

        if (a->prev != NULL)
        {
            (a->prev)->next = b;
        }

        a->prev = b;
    }

    return (b);
}

/*<%%STA>-----------------------------------------------------------------
  FUNCNAME        :long comll_create_index(COMLL_PTR a)
  DESCRIPTION     :creates an index for the list --- index of top element is 0,
then 1, 2,...; returns index of last element. Returns -1 if a is NULL.
<%%END>---------------------------------------------------------------------*/
long comll_create_index(COMLL_PTR a)
{
    long      i = 0;
    COMLL_PTR b = NULL;

    if (!a)
        return -1;

    for (i = 0, b = comll_gototop(a); b->next != NULL; b->index = i, b = b->next, i++)
        ;

    b->index = i;

    return (i);
}

/*<%%STA>-----------------------------------------------------------------
  FUNCNAME        :COMLL_PTR comll_join(COMLL_PTR a, COMLL_PTR b)
  DESCRIPTION     :joins together two lists, so that a->b; returns a unless a
is NULL
<%%END>---------------------------------------------------------------------*/
COMLL_PTR comll_join(COMLL_PTR a, COMLL_PTR b)
{
    COMLL_PTR enda = NULL;

    if (!a)
        return b;
    if (!b)
        return a;

    b = comll_gototop(b);

    enda       = comll_gotobot(a);
    enda->next = b;
    b->prev    = enda;

    return a;
}

/*<%%STA>-----------------------------------------------------------------
  FUNCNAME        :void comll_cut(COMLL_PTR a)
  DESCRIPTION     :frees list after a.
<%%END>---------------------------------------------------------------------*/
void comll_cut(COMLL_PTR a)
{
    COMLL_PTR b = NULL;

    if (!a || !a->next)
        return;

    b       = a->next;
    a->next = NULL;
    b->prev = NULL;
    comll_free_list(b);
}

/*<%%STA>-----------------------------------------------------------------
  FUNCNAME        :int comll_atom(COMLL_PTR a,COMLLType t)
  DESCRIPTION     :returns true if a is a single node of type t
<%%END>---------------------------------------------------------------------*/
int comll_atom(COMLL_PTR a, COMLLType t)
{
    return (a->next == NULL && a->prev == NULL && a->type == t);
}

/*
assumes that an index has been created in b. Returns ptr to
  node that contains index ind, or NULL if there is no such.
*/
static COMLL_PTR comll_find_index(COMLL_PTR b, long ind)
{
    while (b)
    {
        if (b->index > ind)
            b = b->prev;
        else if (b->index < ind)
            b = b->next;
        else
            break;
    }

    return b;
}

/** b is assumed to be a copy of c, but lacking the jmp pointers ***/
static void comll_copy_jmp(COMLL_PTR b, COMLL_PTR c)
{
    comll_create_index(b);
    comll_create_index(c);

    while (b)
    {
        if (c->jmp)
        {
            b->jmp = comll_find_index(b, c->jmp->index);
        }

        b = b->next;
        c = c->next;
    }
}

/*<%%STA>-----------------------------------------------------------------
  FUNCNAME        :COMLL_PTR comll_copy(COMLL_PTR a)
  DESCRIPTION     :copy a comll(even jumps, and allocates new cvg arrays,etc.)
<%%END>---------------------------------------------------------------------*/
/****************************************************************************
COMLL_PTR comll_copy(COMLL_PTR a)
{
COMLL_PTR b=NULL,c = NULL;

        if (!a)
                return NULL;

        c = comll_gototop(a);

        while(c)
        {
                b = comll_insert_after(b);
                memcpy(b,c,sizeof(COMLL_STR));

                if(b->dfindlen)
                {
                        b->dfind = (long *)srt_calloc(b->dfindlen,sizeof(int));
                        b->cvg = (double *)srt_calloc(b->dfindlen,sizeof(double));
                        memcpy(b->dfind,c->dfind,b->dfindlen*sizeof(int));
                        memcpy(b->cvg,c->cvg,b->dfindlen*sizeof(double));
                }
                c = c->next;
        }

        c = comll_gototop(a);
        b = comll_gototop(b);
        comll_copy_jmp(b,c);

        return b;
}
******************************************************************************/
