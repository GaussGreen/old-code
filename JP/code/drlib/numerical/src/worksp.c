#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif

#if defined(USE_ALIB_MALLOC)

/*---------------------------------------------------
 * imsl_alib_malloc
 * 
 * Calls malloc but initializes the allocated memory 
 * block to zero via a call to memset(). 
 * 
 * Jim M.
 */ 

#if defined(ANSI)
void * imsl_alib_malloc (long nbytes)
#else
void * imsl_alib_malloc (nbytes)
       long nbytes; 
#endif
{
    char   * p_ret = NULL; /* return value */ 
    size_t   n_bytes = 0;
    if (nbytes > 0L) 
    {
        n_bytes = (size_t) nbytes; 
        p_ret = (char *)malloc( n_bytes ); 
        if ( p_ret == NULL)
             return (NULL); /* failed */

        memset (p_ret, 0, n_bytes ); 
        return ((void *)p_ret); 
    }
    else 
    { 
        return (NULL); 
    }
}

/*---------------------------------------------------
 * imsl_alib_free
 * 
 * Calls free but checks if the pointer is NULL before
 * making the call. 
 * 
 * Jim M.
 */ 

#if defined(ANSI) 
void   imsl_alib_free (void * p_block)
#else 
void   imsl_alib_free (p_block)
     void * p_block; 
#endif
{ 
     if (p_block != NULL) 
       free( p_block ); 
}

#else  /* either USE_IMSL_MALLOC or imsl_malloc ... arre defined via macros. */

                /* define functions only if not defined via a macro */
#ifndef imsl_malloc

static int active_allocs = 0;
static int total_allocs = 0;
static int active_allocs_size = 0;
static int max_active_allocs_size = 1<<30;
static int extra = 4*sizeof(int);
static int begin_marker = 123457;
static int end_marker = 7543218;

#if defined(ANSI)
void *imsl_malloc (int nbytes)
#else
void *imsl_malloc(nbytes)
    int	    nbytes;
#endif
{
    char    *p;
    int    size, k;
    int    *ip;
    char    *p_block; 

    active_allocs++;
    total_allocs++;
    active_allocs_size += nbytes;
    if (active_allocs_size <= max_active_allocs_size) 
    {
      size = nbytes + ((4-(nbytes&0x3))&0x3);	/* force size to mult of 4 */
      p  = malloc(size+extra);
      ip = (int*)p;
    } 
    else
    {
      p = NULL;
    }
    if (p == NULL) 
    {
      active_allocs--;
      active_allocs_size -= nbytes;
      return NULL;
    }

    ip[0] = nbytes;
    ip[1] = size;
    ip[2] = begin_marker;
    ip[3+size/4] = end_marker;

    /*
    ** Zero out the allocated block. 
    */
    p_block = (char*)(&ip[3]);
    for (k=0; k < size; k++, p_block++) 
    {
      *p_block = 0; 
    }

    return (p+3*sizeof(int));
}


#if defined(ANSI)
void imsl_free (void *p)
#else
void imsl_free(p)
    void *p;
#endif
{
    int	    nbytes;
    int	    *ip = (int *)((char* )p-3*sizeof(int));
    int	    size;

    if (p == NULL) return;

    active_allocs--;
    if (ip[2] != begin_marker) {
/*	(5, 2, "Error in workspace.  Probably overwrote past the start of the allocation.");
*/
        imsl_ermes(IMSL_TERMINAL, IMSL_WORKSPACE_ERROR);
	return;
    }
    nbytes = ip[0];
    size   = ip[1];
    if (ip[3+size/4] != end_marker) {
/*	(5, 2, "Error in workspace.  Probably overwrote past the end of the allocation.");
*/
        imsl_ermes(IMSL_TERMINAL, IMSL_WORKSPACE_ERROR);
	return;
    }
    active_allocs_size -= nbytes;
    if (p != NULL) free((char*)ip);
}


void set_max_worksp(nbytes)
    int	    nbytes;
{
    if (max_active_allocs_size < active_allocs_size) {
	imsl_e1sti(1, active_allocs_size);
	imsl_e1sti(2, nbytes);
/*	(5, 3, "The maximum allocation size = %(i2) is less than the current active allocation size.");
*/
        imsl_ermes(IMSL_TERMINAL, IMSL_ERROR_ACTIVE_ALLOC_SIZE);
    } else
	max_active_allocs_size = nbytes;
}

void worksp_report()
{
    fprintf(stdout,
	"\n%s%10d\n%s%10d\n%s%10d\n%s%10d\n",
	"Current number of active allocations       ", active_allocs,
	"Total number of allocations                ", total_allocs,
	"Size of current active allocations (bytes) ", active_allocs_size,
	"Maximum size of active allocations allowed ", max_active_allocs_size);
}
#endif


/* USE_ALIB_MALLOC */
#endif


