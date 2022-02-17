//----------------------------------------------------------------------------
//
//   Group       : QR&D
//
//   Filename    : QMCStrata.hpp
//
//   Description : A definition of strata for stratified sampling
//
//   Date        : 22 Nov 2006
//
//----------------------------------------------------------------------------

#ifndef QMCStrata_HPP
#define QMCStrata_HPP
#include "edginc/VirtualDestructorBase.hpp"
#include "edginc/Object.hpp"
#include "edginc/FORWARD_DECLARE.hpp"
#include <map>

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE(QMCStrata);

/** An interface to a single strata */
class MCARLO_DLL QMCStrata : public CObject
{
public:
    static CClassConstSP const TYPE;
    virtual string getName() const { return name; }
    virtual bool notIntersecting(QMCStrataConstSP strata2) const;
    void checkTypes(QMCStrataConstSP strata2) const;

protected:
    QMCStrata(const CClassConstSP& clazz) : CObject(clazz) {}
    QMCStrata() : CObject(TYPE) {}
    static void load(CClassSP& clazz);
private:
    static IObject* defaultConstructor() { return new QMCStrata; }
    string name;
    
};
DECLARE(QMCStrata);


class MCARLO_DLL QMCStratificationRecord : public CObject
{
public:
    static CClassConstSP const TYPE;

    virtual string getName() const { return name; }
    virtual double getWeight() const { return weight; }
    QMCStrataConstSP getStrata() const { return strata; }

    QMCStratificationRecord(double _weight, QMCStrataConstSP _strata) : 
        CObject(TYPE), weight(_weight), strata(_strata) {}

protected:
    QMCStratificationRecord(const CClassConstSP& clazz) : CObject(clazz), weight(0.0) {}

private:
    static void load(CClassSP& clazz);
    static IObject* defaultConstructor() { return new QMCStratificationRecord; }
    QMCStratificationRecord(): CObject(TYPE), weight(0.0) {}

    string name;
    double weight;
    QMCStrataConstSP strata;
};
DECLARE(QMCStratificationRecord);

#ifndef QLIB_QMCSTRATA_CPP
EXTERN_TEMPLATE(class MCARLO_DLL_SP smartConstPtr<QMCStratificationRecord>);
EXTERN_TEMPLATE(class MCARLO_DLL_SP smartPtr<QMCStratificationRecord>);
#else
INSTANTIATE_TEMPLATE(class MCARLO_DLL smartConstPtr<QMCStratificationRecord>);
INSTANTIATE_TEMPLATE(class MCARLO_DLL smartPtr<QMCStratificationRecord>);
#endif



// specific stratas implemented for use:
// CID total number of jumps control
class MCARLO_DLL QMCStrataCIDJumps : public QMCStrata
{
public:
    static CClassConstSP const TYPE;
    //virtual void getMarket(const IModel* model, const MarketData* market) {}
    virtual pair<int, int> getNJumps() const { return make_pair(lowestNJumps, highestNJumps);}
    virtual bool notIntersecting(QMCStrataConstSP strata2) const;

    QMCStrataCIDJumps(int _lowestNJumps, int _highestNJumps) : 
        QMCStrata(TYPE), lowestNJumps(_lowestNJumps), highestNJumps(_highestNJumps) {}

protected:
    QMCStrataCIDJumps(const CClassConstSP& clazz) : QMCStrata(clazz), lowestNJumps(0), highestNJumps(-1) {}

private:
    static void load(CClassSP& clazz);
    static IObject* defaultConstructor() { return new QMCStrataCIDJumps; }
    QMCStrataCIDJumps(): QMCStrata(TYPE) {}

    int    lowestNJumps;
    int    highestNJumps;
};
DECLARE(QMCStrataCIDJumps);


DRLIB_END_NAMESPACE

#endif // QMCStrata
