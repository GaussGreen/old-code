//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : Malloc.hpp
//
//   Description : Class for representing addins
//
//   Author      : Mark A Robson
//
//   Date        : 12 Feb 2001
//
//
//----------------------------------------------------------------------------

#ifndef EDG_MALLOC_H
#define EDG_MALLOC_H

#include <cstddef>


CORE_BEGIN_NAMESPACE
/** Wrapper class for malloc that throws an exception on failure.
    Main reasons for using malloc: need realloc.

    Note: Memory obtained through Malloc must not be freed with delete and
    vice versa */
class CORE_DLL Malloc{
public:
    /** malloc amount number of bytes. Throws exception on failure */
    static void* allocate(size_t amount);

    /** malloc amount number of bytes and set to 0.
        Throws exception on failure */
    static void* callocate(size_t number, size_t size);

    /** free memory */
    static void deallocate(void* ptr);

    /** realloc - resize ptr to amount number of bytes.
        Throws exception on failure */
    static void* reallocate(void* ptr, size_t amount);
};


// Macros are a bit essential here
#ifndef NEW
#define NEW(t)              (t *) Malloc::callocate(1, sizeof(t))
#endif
#ifndef NEW_ARRAY
#define NEW_ARRAY(t,n)      (t *) Malloc::callocate((n), sizeof(t))
#endif
#ifndef FREE
#define FREE(x) Malloc::deallocate((void *) (x))
#endif
#ifndef REALLOC
#define REALLOC(ptr, type, size) \
    (type *) Malloc::reallocate((ptr), (size)*sizeof(type))
#endif

CORE_END_NAMESPACE

#endif

