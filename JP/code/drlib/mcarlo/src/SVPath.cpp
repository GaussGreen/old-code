//----------------------------------------------------------------------------
//
//   Group       : QR&D
//
//   Filename    : SVPath.cpp
//
//   Description : A Legacy non-type-specific access mode to the Spot Variable 
//
//   Date        : March 2004
//
//
//----------------------------------------------------------------------------


#include "edginc/config.hpp"
#define QLIB_SVPATH_CPP
#include "edginc/SVPath.hpp"
//#include "edginc/MCPathConfigSRMGenSV.hpp"

USING_DRLIB_NAMESPACE

SVPath::~SVPath(){}
SVPath::SVPath(const double*      ppath, 
                   const vector<int>& offsets,
                   int                beginIdx, 
                   int                endIdx):
            ppath(ppath), 
            offsets(offsets), 
            beginIdx(beginIdx), 
            endIdx(endIdx), 
            currentIdx(beginIdx), 
            isTrivial(false)
{
                int i = 1;
}

SVPath::SVPath(const double* ppath, int beginIdx, int endIdx):
            ppath(ppath), 
            offsets(0), 
            beginIdx(beginIdx), 
            endIdx(endIdx), 
            currentIdx(beginIdx), 
            isTrivial(true) 
{}

SVPath::SVPath():
            ppath(0), 
            offsets(0), 
            beginIdx(0), 
            endIdx(0), 
            currentIdx(0), 
            isTrivial(true) 
{}

void SVPath::initialize(
                         const double*      _path, 
                         const vector<int>& _offsets,
                         int                _beginIdx, 
                         int                _endIdx)
{
    ppath = _path;
    offsets = _offsets; 
    beginIdx = _beginIdx; 
    endIdx = _endIdx; 
    currentIdx = _beginIdx;
    isTrivial = false;
}

void SVPath::initialize(
                         const double*      _path, 
                         int                _beginIdx, 
                         int                _endIdx)
{
    ppath = _path;
    offsets.clear(); 
    beginIdx = _beginIdx; 
    endIdx = _endIdx; 
    currentIdx = _beginIdx;
    isTrivial = true;
}

void SVPath::initialize(
                         int                _beginIdx, 
                         int                _endIdx)
{
    ppath = NULL;
    offsets.clear(); 
    beginIdx = _beginIdx; 
    endIdx = _endIdx; 
    currentIdx = _beginIdx;
    isTrivial = true;
}


double SVPath::element(int idx) const
{ 
    throw ModelException(
        "SVPath::element()", 
        "This function shall not be called in the base class. "
        "Most likely cause of this error is that Path has NULL pointer inside."); 
}


