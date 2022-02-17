#ifndef NAGX04
#define NAGX04
#ifdef __cplusplus
extern "C"
{
#endif


/* <nagx04.h>
 *
 * Copyright 1996 Numerical Algorithms Group
 *
 * Include file for NAG C Library x04 Chapter
 *
 * Mark 4, revised, 1996.
 * Mark 5 revised. IER-2169 (Feb 1998).
 */
#include <nag_stddef.h>
#include <stdio.h>


#ifdef NAG_PROTO
extern void x04aaz(char NAG_HUGE *str);
#else
extern void x04aaz();
#endif

#ifdef NAG_PROTO
extern FILE *x04bax(const char *input_file_name);
#else
extern FILE *x04bax();
#endif


#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL x04bay(const char *output_file_name, const char *string);
#else
extern void x04bay();
#endif

#ifdef NAG_PROTO
extern void x04baz(Nag_FileSt NAG_HUGE *stream, char NAG_HUGE *str, int code, Nag_Mesg NAG_HUGE *mesg);
#else
extern void x04baz();
#endif

#ifdef NAG_PROTO
extern void x04bbz(char NAG_HUGE str[]);
#else
extern void x04bbz();
#endif

#ifdef NAG_PROTO
extern void x04bcz(char NAG_HUGE str[]);
#else
extern void x04bcz();
#endif

#ifdef NAG_PROTO
extern char NAG_HUGE *x04caz(char NAG_HUGE *s1, char NAG_HUGE *s2, int n);
#else
extern char *x04caz();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP Pointer NAG_CALL x04bbc(size_t size);
#else
 extern Pointer x04bbc();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP Pointer NAG_CALL x04bcc(Pointer ptr, size_t newsize);
#else
 extern Pointer x04bcc();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL x04bdc(Pointer *ptr);
#else
 extern void x04bdc();
#endif

#ifdef NAG_PROTO
extern Nag_HashError x04fag(Nag_HashTable *htable, unsigned long approx_size);
#else
extern Nag_HashError x04fag();
#endif

#ifdef NAG_PROTO
extern void x04fah(Nag_HashTable *htable);
#else
extern void x04fah();
#endif

#ifdef NAG_PROTO
extern Pointer x04faj(Nag_HashTable *htable, char *name);
#else
extern Pointer x04faj();
#endif

#ifdef NAG_PROTO
extern void x04fak(Nag_HashTable *htable, Pointer data, char *name, 
		   Nag_HashError *error);
#else
extern void x04fak();
#endif

#ifdef NAG_PROTO
extern Pointer x04fal(Nag_HashTable *htable, Pointer data, char *name, 
		      Nag_HashError *error);
#else
extern Pointer x04fal();
#endif

#ifdef NAG_PROTO
extern void x04fam(Nag_HashTable *htable, Pointer old_base, Pointer new_base,
		   size_t offset);
#else
extern void x04fam();
#endif

#ifdef __cplusplus
}
#endif
#endif /* not NAGX04 */
