//
//   Filename    : IElemStateVariableGenVisitor.hpp
//
//   Description : Interface class for IElemStateVariableGen Generators
//
//   Date        : 6 Apr 2006
//
//----------------------------------------------------------------------------
#ifndef IELEMSTATEVARIABLEGENVISITOR_HPP
#define IELEMSTATEVARIABLEGENVISITOR_HPP

DRLIB_BEGIN_NAMESPACE

class IElemStateVariableGen;

class SVGenDiscFactor;
class SVGenExpectedDiscFactor;

class SVGenSurvivalDiscFactor;
class SVGenExpectedSurvivalDiscFactor;
class SVGenAggregatedSurvivalDiscFactor;

class SVGenSpot;
class SVGenExpectedSpot;
class SVGenExpectedEnergyFuture;
class SVGenExpectedBasisFwdSpread;
class SVGenPathWeight;
class SVGenDateOfDefault;

class IElemStateVariableGenVisitor {
    public:
        virtual void processSVGen(const IElemStateVariableGen* base) = 0;
        virtual void processSVGen(const SVGenDiscFactor* df) = 0;
        virtual void processSVGen(const SVGenExpectedDiscFactor* edf) = 0;
        virtual void processSVGen(const SVGenSurvivalDiscFactor* sdf) = 0;
        virtual void processSVGen(const SVGenExpectedSurvivalDiscFactor* esdf) = 0;
        virtual void processSVGen(const SVGenAggregatedSurvivalDiscFactor* asdf) = 0;
        virtual void processSVGen(const SVGenExpectedSpot* es) = 0;
        virtual void processSVGen(const SVGenSpot* s) = 0;
        virtual void processSVGen(const SVGenExpectedBasisFwdSpread* s) = 0;
        virtual void processSVGen(const SVGenDateOfDefault* dod) = 0;
        virtual void processSVGen(const SVGenExpectedEnergyFuture* en) = 0;
        virtual void processSVGen(const SVGenPathWeight* pw) = 0;
    virtual ~IElemStateVariableGenVisitor() {}
};

DRLIB_END_NAMESPACE
#endif // IELEMSTATEVARIABLEGENVISITOR_HPP
