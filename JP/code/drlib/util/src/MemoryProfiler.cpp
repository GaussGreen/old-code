//----------------------------------------------------------------------------
//
//   Group       : QRD
//
//   Filename    : MemoryProfiler.cpp
//
//   System dependent classes to gather current memory consumption
//
//   Author      : Vladimir Grebinskiy
//
//   Date        : 03/27/2006
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/MemoryProfiler.hpp"
#include <iostream>
#include <cstdio>

DRLIB_BEGIN_NAMESPACE

void MemoryProfiler::print(size_t /*kb*/)
{
}

MemoryProfiler::MemoryProfiler(bool _printInDtor)
    : last(0), printInDtor(_printInDtor)
{
}

size_t      MemoryProfiler::getCurrent() // get current memory consumption
{
    return last = getVMsize();
}

size_t      MemoryProfiler::getDelta()  //  get delta between getCurrent and the previous one
{
    size_t prev = last;
    return getCurrent() - prev;
}

// When we use  Glibc based malloc, some information about heap can be retrieved
// See mallinfo in <malloc.h> for more details. We use to get current heap consumption
class GlibcMallocInfo::Impl {
    public:
        size_t getVMsize(void);
};

// FIXME find a better test to check if mallinfo is available
#if defined (LINUX) || defined (LINUX64)
#include <malloc.h>
size_t GlibcMallocInfo::Impl::getVMsize()
{
    struct mallinfo vm = mallinfo();
    return vm.arena /1024  + vm.hblkhd / 1024; /* total allocated space in KB */
}
#else // for other platform I don't know of a good way to retreive it
size_t GlibcMallocInfo::Impl::getVMsize()
{
    return 0;
}
#endif


GlibcMallocInfo::GlibcMallocInfo(string s, bool  _printInDtor) :
    MemoryProfiler(_printInDtor),
    name(s),
    impl(new Impl())
{
    getDelta();
}

GlibcMallocInfo::~GlibcMallocInfo()
{
    if (printInDtor) {
        print(getDelta());
        print(getCurrent());
    }
    delete impl;
}

size_t GlibcMallocInfo::getVMsize()
{
    return impl->getVMsize();
}

void GlibcMallocInfo::print(size_t kb)
{
    clog << name << " : " << kb << endl;
}


class ProcessVMSize::Impl {
    public:    
        size_t getVMsize();
    private:
        int readIntField(string name); // works even on Windows, where it will always return -1
};
#if defined (LINUX) || defined (LINUX64)
int ProcessVMSize::Impl::readIntField(string name)
{

    const char * info_file = "/proc/self/status"; // Linux specific
    FILE * f = fopen(info_file, "r");

    string fmt = name + ": %d";
    const char * cfmt = fmt.c_str();
    int result = -1;
    if (f == NULL)
        return result;

    char * s = NULL;
    size_t n = 0;
    while(getline(&s, &n, f) != -1) // Glibc specific
    {
        if (sscanf(s, cfmt, & result) == 1) // Glibc specific
            break;
    }
    fclose(f);
    if (s)
        free(s);

    return result;
}
#else // defined (LINUX) || defined (LINUX64)
int ProcessVMSize::Impl::readIntField(string name)
{
    return 0;
}
#endif

size_t ProcessVMSize::Impl::getVMsize()
{
    return readIntField("VmSize");
}

size_t ProcessVMSize::getVMsize()  // returns amount of memory consumed by malloc : see mallinfo in <malloc.h>
{
    return impl->getVMsize();
}

void ProcessVMSize::print(size_t kb)
{
    clog << name << " : " << kb << endl;
}

ProcessVMSize::ProcessVMSize(string s , bool  _printInDtor) :
        MemoryProfiler(_printInDtor),
        impl(new Impl()),
        name(s)
{
    getDelta();
}

ProcessVMSize::~ProcessVMSize()
{
    if (printInDtor) {
        print(getDelta());
        print(getCurrent());
    }
    delete impl;
}


DRLIB_END_NAMESPACE
