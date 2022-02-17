//----------------------------------------------------------------------------
//
//   Group       : Global QR&D
//
//   Filename    : Fix3TD.hpp
//
//----------------------------------------------------------------------------

#ifndef _FIX3TD_HPP
#define _FIX3TD_HPP

#include "edginc/SingleCcyTree.hpp"

DRLIB_BEGIN_NAMESPACE

class TREE_DLL Fix3TD : public SingleCcyTree {
public:

    /*************** FDModel interface ****************/
protected:
    virtual void retrieveFactor();

    /********************** local methods  **************************/

    virtual void initTreeData();  // SingleCcyTree::initTreeData

    void Populate2QSmile(const string& smileKey, IRExoticParamTable& smileTable,
                         FIX3_TREE_DATA& treeData, MKTVOL_DATA& mktVolData);
    void PopulateEngineParams(const string& engineKey, IRModelConfigTable& engineTable,
                              FIX3_TREE_DATA& treeData, MKTVOL_DATA& mktVolData);
    void PopulateModelParams(const string& modelKey, IRExoticParamTable& modelTable,
                             FIX3_TREE_DATA& treeData, MKTVOL_DATA& mktVolData,
                             bool termStructure, int nbFactors);

protected:
    Fix3TD(const CClassConstSP &type = TYPE) : SingleCcyTree(type) {}
private:
    static IObject* defaultConstructor(void) { return new Fix3TD(); }
    static void load(CClassSP& clazz); /* declare what is being exposed */

    /***************************** variables ********************************/
public:
    static CClassConstSP const TYPE;

};

typedef smartPtr<Fix3TD> Fix3TDSP;

DRLIB_END_NAMESPACE

#endif

