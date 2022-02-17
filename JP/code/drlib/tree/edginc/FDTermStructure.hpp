//----------------------------------------------------------------------------
//
//   Group       : Convertibles Research
//
//   Filename    : FDTermStructure.hpp
//
//   Description : Defines a term structure used in the generic 1F FD Engine
//
//   Author      : André Segger
//
//   Date        : 10 October 2003
//
//----------------------------------------------------------------------------

#ifndef TERM_STRUCTURE_HPP
#define TERM_STRUCTURE_HPP

#include "edginc/DateTime.hpp"

DRLIB_BEGIN_NAMESPACE

/** need to make this an abstract base class with flat and term structure
    derived classes when proven to be suitable */
class TREE_DLL FDTermStructure : public CObject {
public:
    static CClassConstSP const TYPE;
    friend class FDTermStructureHelper;
    friend class FD1FE2C;

    virtual ~FDTermStructure () {}

    FDTermStructure (double flatTerm);
    FDTermStructure (DoubleArraySP termArray);

    void scale(const double& scalingFactor);

    bool            isFlat() { return (dimension==0); }

    double	    operator()(const int r) const;

private:
    FDTermStructure (): CObject(TYPE), dimension(0), term(0.0) {}

    int           dimension; // $unregistered
    double        term; // $unregistered
    DoubleArraySP termStructure; // $unregistered
};


typedef smartConstPtr<FDTermStructure >             FDTermStructureConstSP;
typedef smartPtr<FDTermStructure >                  FDTermStructureSP;

typedef array<FDTermStructureSP, FDTermStructure >  FDTermStructureArray;
typedef smartPtr<FDTermStructureArray>              FDTermStructureArraySP;
typedef smartConstPtr<FDTermStructureArray>         FDTermStructureArrayConstSP;

DRLIB_END_NAMESPACE
#endif
