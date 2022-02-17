//----------------------------------------------------------------------------
//
//   Group       : QRD
//
//   Filename    : MemoryProfiler.hpp
//
//   System dependedn class to gather current memory consumption
//
//   Author      : Vladimir A grebinskiy
//
//   Date        : 03/27/2006
//
//
//----------------------------------------------------------------------------

#ifndef QLIB_MEMORYPROFILER_HPP
#define QLIB_MEMORYPROFILER_HPP

#include "edginc/config.hpp"
#include <time.h>

DRLIB_BEGIN_NAMESPACE

class UTIL_DLL MemoryProfiler {
    private:
        size_t     last;
    protected:
        bool        printInDtor; // print in destructor
        virtual     size_t getVMsize() = 0; // read some info about memory consumption
        virtual     void print(size_t kb);  // does nothing by default
    public:
        MemoryProfiler(bool printInDtor); //
        virtual     ~MemoryProfiler() {}; // we cannot print here, as it is too late by now
        size_t      getCurrent(); // get current memory consumption
        size_t      getDelta();  //  get delta between getCurrent and the previous one
};

// When we use  Glibc based malloc, some information about heap can be retrieved
// See mallinfo in <malloc.h> for more details. We use to get current heap consumption
class UTIL_DLL GlibcMallocInfo : public MemoryProfiler {
    string name;
    class Impl; // platform specific
    Impl * impl;
protected:
    virtual     size_t getVMsize(); 
    virtual     void print(size_t kb);
public:
    GlibcMallocInfo(string objectName = "malloc", bool  printInDtor = true);
    ~GlibcMallocInfo();
};

// The Linux kernel allows to read current consumption of the virtual memory via reading /proc/self/

class UTIL_DLL ProcessVMSize : public MemoryProfiler {
    string name;
    class Impl; // platform specific
    Impl * impl;
protected:
    virtual     size_t getVMsize(); // returns amount of memory consumed by malloc : see mallinfo in <malloc.h>
    virtual     void print(size_t kb);
public:
    ProcessVMSize(string objName = "vmsize", bool  printInDtor = true);
    ~ProcessVMSize();
};

DRLIB_END_NAMESPACE

#endif