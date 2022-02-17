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
#include "edginc/coreConfig.hpp"
#include "edginc/Malloc.hpp"
#include "edginc/ModelException.hpp"
#include "edginc/Format.hpp"

CORE_BEGIN_NAMESPACE


/** Wrapper class for malloc that throws an exception on failure.
    Main reasons for using malloc: need realloc.

    Note: Memory obtained through Malloc must not be freed with delete and
    vice versa */

/** malloc amount number of bytes. Throws exception on failure */
void* Malloc::allocate(size_t amount){
    void*  ptr = malloc(amount);
    if (!ptr){
        throw ModelException("Malloc::allocate", "Couldn't allocate "+
                             Format::toString((int)amount)+" bytes");
    }
    return ptr;
}

/** malloc amount number of bytes and set to 0.
    Throws exception on failure */
void* Malloc::callocate(size_t number, size_t size){
    void*  ptr = calloc(number, size);
    if (!ptr){
        throw ModelException("Malloc::callocate", "Couldn't allocate "+
                             Format::toString((int)(number*size))+" bytes");
    }
    return ptr;
}


/** free memory */
void Malloc::deallocate(void* ptr){
    free(ptr);
}


/** realloc - resize ptr to amount number of bytes.
    Throws exception on failure */
void* Malloc::reallocate(void* ptr, size_t amount){
    void*  newPtr = realloc(ptr, amount);
    if (!newPtr){
        throw ModelException("Malloc::allocate", "Couldn't reallocate "+
                             Format::toString((int)(amount))+" bytes");
    }
    return newPtr;
}
CORE_END_NAMESPACE
