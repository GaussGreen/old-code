#ifndef SVGEN_IPROD_CREATOR_HPP
#define SVGEN_IPROD_CREATOR_HPP

#include "edginc/StateVariableClient.hpp"

DRLIB_BEGIN_NAMESPACE

class IRPRODUCTS_DLL SVIProdCreator : public virtual IStateVariable {
	 public:
		  virtual double elem(size_t index) = 0;
};
DECLARE_REF_COUNT(SVIProdCreator);

class  IRPRODUCTS_DLL SVGenIProdCreator : public virtual IStateVariableGen,
	                                 public virtual IStateVariableClient {
public:    
    virtual SVIProdCreatorSP getSVIProdCreator(IStateVariableGen::IStateGen* pathGen) const = 0;    
    virtual IStateVariableSP create(IStateVariableSP oldStateVar, IStateVariableGen::IStateGen* pathGen) const{
        return IStateVariableSP(getSVIProdCreator(pathGen).get());
    }
    virtual void collectStateVars(IStateVariableCollectorSP svCollector) const = 0;
};
DECLARE_REF_COUNT(SVGenIProdCreator);

DRLIB_END_NAMESPACE
#endif
