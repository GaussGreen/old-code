/* ==========================================================================
   FILE_NAME:	utalloc.c

   PURPOSE:     for debugging memory allocation or tracing it
   ========================================================================== */

#include "utallhdr.h"

/**************************************************

        Debugging functions for allocation
        these functions allocate memory  , and also
        print to a file the address just allocated.
        Free prints to a similar file. The idea is to compare
        the two files to check for memory leakage.

***************************************************/

EXTERND char srt_dbg_alloc_buffer[200];

void *srt_dbg_calloc(size_t nobj, size_t size) {
  void *mem = calloc(nobj, size);
  FILE *f;
  if (!mem) {
    f = fopen("alloc.out", "a");
    if (f) {
      fprintf(f, "srt_calloc failed size %d\n", (int)size);
      fclose(f);
    }
    return mem;
  }
  f = fopen("alloc.out", "a");
  if (f) {
    fprintf(f, "%s %lu\n", srt_dbg_alloc_buffer, (unsigned long)mem);
    fclose(f);
  }
  return mem;
}
void *srt_dbg_malloc(size_t size) {
  void *mem = malloc(size);
  FILE *f;
  if (!mem) {
    f = fopen("alloc.out", "a");
    if (f) {
      fprintf(f, "srt_malloc failed size %d\n", (int)size);
      fclose(f);
    }
  }
  f = fopen("alloc.out", "a");
  if (f) {
    fprintf(f, "%s %lu\n", srt_dbg_alloc_buffer, (unsigned long)mem);
    fclose(f);
  }
  return mem;
}

void srt_dbg_free(void *mem) {
  FILE *f;
  f = fopen("free.out", "a");
  if (f) {
    fprintf(f, "%s %lu\n", srt_dbg_alloc_buffer, (unsigned long)mem);
    fclose(f);
  }
  free(mem);
  mem = NULL;
}

void srt_dbg_alloc_init(void) {
  FILE *f;
  f = fopen("free.out", "w");
  if (f)
    fclose(f);
  f = fopen("alloc.out", "w");
  if (f)
    fclose(f);
}
void srt_dbg_alloc_end(void) {
  /*	FILE *f1  ,*f2;

          f1 = fopen("free.out"  ,"r");
          f2 = fopen("alloc.out"  ,"r");
          if(!f1 || !f2){
                  if(f1)fclose(f1);
                  if(f2)fclose(f2);
                  return;
          }
  */
}

void free_and_zero(void **p) {
  if (*p)
    free(*p);
  *p = 0;
}
