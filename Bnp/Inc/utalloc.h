#ifndef UTALLOC_H
#define UTALLOC_H

void srt_dbg_alloc_init(void);
void srt_dbg_alloc_end(void);
void *srt_dbg_calloc(size_t nobj, size_t size);
void srt_dbg_free(void *mem);
void *srt_dbg_malloc(size_t size);
void free_and_zero(void **p);

EXTERN char srt_dbg_alloc_buffer[];

#ifndef SRT_DBG_ALLOC
#define srt_calloc calloc
#define srt_malloc malloc
#define srt_free(x) (free(x), x = NULL)
#define srt_alloc_init
#define srt_alloc_end
#else
#define srt_calloc(x, y)                                                       \
  (sprintf(srt_dbg_alloc_buffer, "cal x y"), srt_dbg_calloc(x, y))
#define srt_malloc(x)                                                          \
  (sprintf(srt_dbg_alloc_buffer, "mal x"), srt_dbg_malloc(x))
#define srt_free(x) (sprintf(srt_dbg_alloc_buffer, "fre x"), srt_dbg_free(x))
#define srt_alloc_init srt_dbg_alloc_init
#define srt_alloc_end srt_dbg_alloc_end
#endif

#endif
