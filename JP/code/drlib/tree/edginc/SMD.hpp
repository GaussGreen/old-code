//----------------------------------------------------------------------------
//
//   Group       : Global QR&D
//
//   Filename    : SMD.hpp
//
//----------------------------------------------------------------------------

#ifndef _SMD_HPP
#define _SMD_HPP

#include "edginc/SingleCcyTree.hpp"

DRLIB_BEGIN_NAMESPACE

class TREE_DLL SMD : public SingleCcyTree {
public:

    /*************** FDModel interface ****************/
protected:
    virtual void retrieveFactor();

    /********************** local methods  **************************/
protected:
    SMD(const CClassConstSP &type = TYPE) : SingleCcyTree(type) {}
private:
    static IObject* defaultConstructor(void) { return new SMD(); }
    static void load(CClassSP& clazz); /* declare what is being exposed */

    /***********************export variables ********************************/

    // SMD parameters - ??? make a separate structure and market family -
    double correlationSkew;
    double correlationCurvature;
    double correlationLevel;
    double correlationTermStructure;

public:
    static CClassConstSP const TYPE;

};

typedef smartPtr<SMD> SMDSP;

DRLIB_END_NAMESPACE

#endif

