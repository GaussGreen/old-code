#include "edginc/config.hpp"
#include "edginc/Lock.hpp"

#ifdef WIN32
#include <windows.h>
#endif

DRLIB_BEGIN_NAMESPACE

#ifdef WIN32

Lock::Lock() 
{ 
    _mutex = CreateMutex(NULL, false, NULL); 
}

void Lock::lock() 
{ 
    WaitForSingleObject(_mutex, INFINITE); 
}

void Lock::unlock() 
{ 
    ReleaseMutex(_mutex); 
}

Lock::~Lock() 
{ 
    CloseHandle(_mutex); 
}

#else

Lock::Lock() 
{
    pthread_mutex_init(&_mutex, NULL); 
}

void Lock::lock() 
{ 
    pthread_mutex_lock(&_mutex); 
}

void Lock::unlock() 
{ 
    pthread_mutex_unlock(&_mutex); 
}

Lock::~Lock() 
{ 
    pthread_mutex_destroy(&_mutex); 
}

#endif

DRLIB_END_NAMESPACE


