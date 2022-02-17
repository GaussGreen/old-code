#ifndef NAGM01
#define NAGM01
#ifdef __cplusplus
extern "C"
{
#endif


/* <nagm01.h>
 *
 * Copyright 1996 Numerical Algorithms Group
 *
 * Include file for NAG C Library m01 Chapter
 *
 * Mark 4, revised, 1996.
 * Mark 5 revised. IER-2165 (Feb 1998).
 */

#define SIZE_T_MAX 2147483647   /* maximum size_t value */
#define SIZE_T_BITS 32            /* no. of bits in size_t */
#define PTRDIFF_T_MAX (SIZE_T_MAX/2) /* maximum ptrdiff_t value */
#define MAX_LENGTH (PTRDIFF_T_MAX/sizeof(size_t)) /* maximum size of an array */
#define ELEMENT_SPACE 256         /* max. inline scratch element */

#define NAG_MCHAP_ERROR_BUF_LEN 20

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL m01cac(double NAG_HUGE vec[], size_t n, Nag_SortOrder order, NagError NAG_HUGE *fail);
#else
extern void m01cac();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL m01csc(Pointer vec, size_t n, size_t size, ptrdiff_t stride,
            NAG_M01_FUN compare,
            Nag_SortOrder order, NagError NAG_HUGE *fail);
#else
extern void m01csc();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL m01ctc(Pointer vec, size_t n, size_t size, ptrdiff_t stride,
            NAG_M01_FUN compare,
            Nag_SortOrder order, NagError NAG_HUGE *fail);
#else
extern void m01ctc();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL m01cuc(Pointer NAG_HUGE *base, ptrdiff_t offset,
             NAG_M01_FUN compare,
             Nag_SortOrder order, NagError NAG_HUGE *fail);
#else
extern void m01cuc();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL m01dsc(Pointer vec, size_t n, ptrdiff_t stride,
            NAG_M01_FUN compare,
            Nag_SortOrder order, size_t NAG_HUGE *rank, NagError NAG_HUGE *fail);
#else
extern void m01dsc();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL m01esc(Pointer vec, size_t n, size_t size, 
            ptrdiff_t stride, size_t NAG_HUGE *indices, NagError NAG_HUGE *fail);
#else
extern void m01esc();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP Boolean NAG_CALL m01fsc(Pointer key, Pointer vec, size_t n, ptrdiff_t stride,
               NAG_M01_FUN compare,
               Nag_SortOrder order, Nag_SearchMatch final, Pointer NAG_HUGE *match, NagError NAG_HUGE *fail);
#else
extern Boolean m01fsc();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL m01zac(size_t NAG_HUGE rank[], size_t n, NagError NAG_HUGE *fail);
#else
extern void m01zac();
#endif

#ifdef NAG_PROTO
extern void m01zbc(Integer *iperm, Integer m1, Integer m2, Integer *ifail);
#else
extern void m01zbc();
#endif

#ifdef __cplusplus
}
#endif
#endif /* not NAGM01 */
