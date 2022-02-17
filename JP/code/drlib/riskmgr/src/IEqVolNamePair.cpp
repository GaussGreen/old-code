/**
 * @file IEqVolNamePair.cpp
 */

#include "edginc/config.hpp"
#include "edginc/IEqVolNamePair.hpp"
#include "edginc/ObjectIteration.hpp"

DRLIB_BEGIN_NAMESPACE

IEqVolNamePair::~IEqVolNamePair() {}

class GetNamePair : virtual public ObjectIteration::IActionConst {
public:
    GetNamePair(map<string, OutputNameSP>& namePairs) : namePairs(namePairs){}
    virtual bool invoke(const ObjectIteration::State& state,IObjectConstSP obj){
        const IEqVolNamePair& eq = dynamic_cast<const IEqVolNamePair&>(*obj);
        string eqName, volName;
        bool recurse = eq.getNamePairs(eqName, volName);
        namePairs[eqName] = OutputNameSP(new OutputName(volName));
        return recurse;
    }
private:
    map<string, OutputNameSP>& namePairs;
};

map<string, OutputNameSP> IEqVolNamePair::namePairs(IObjectConstSP root) {
    map<string, OutputNameSP> them;
    GetNamePair action(them);
    ObjectIteration iter(IEqVolNamePair::TYPE);
    iter.recurse(action, root);
    return them;
}

CClassConstSP const IEqVolNamePair::TYPE = CClass::registerInterfaceLoadMethod(
    "IEqVolNamePair", typeid(IEqVolNamePair), 0);

DRLIB_END_NAMESPACE
