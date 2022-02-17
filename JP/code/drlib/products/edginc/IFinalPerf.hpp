//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : CallableEquityKOSwap.hpp
//
//   Description   Callable Equity KO Swap 
//
//
//   $Log: CallableEquityKOSwap.hpp,v $
//----------------------------------------------------------------------------

#ifndef EDG_FNLP_HPP
#define EDG_FNLP_HPP

#include "edginc/Object.hpp"
#include "edginc/IDoubleArray.hpp"
#include "edginc/Average.hpp"
#include "edginc/Schedule.hpp"
#include "edginc/Model1F.hpp"
#include "edginc/TreeSlice.hpp"
#include "edginc/Tree1f.hpp"

DRLIB_BEGIN_NAMESPACE

//---------------------------------------------------------//
// Allow some flexibility in styles of FinalPerf by providing interface and 
// opportunity to extend choice of implementers. 

class IFinalPerf;
typedef refCountPtr<IFinalPerf> IFinalPerfSP;

class PRODUCTS_DLL IFinalPerf {
public:
    virtual double getValue(int iPrice, int simDateIdx, double s, double k) = 0;

    class getValue_oper : public SliceMarker< getValue_oper >
    {
    public:
        getValue_oper(const TreeSlice & spot,
                 //const TreeSlice & price,        
                 const IFinalPerfSP fp,
                 int iPrice, 
                 int simDateIdx, 
                 double k)
                 :
                 spot( spot ),
                 //price( price ),
                 fp( fp ),
                 iPrice( iPrice ),
                 simDateIdx( simDateIdx ),
                 k( k ){};

        // TreeSlice "expression template" primitives
        static const int sliceCount = 1;
        template< typename S >
        const S** listSlices(const S** list) const
        {
            //return price.listSlices( spot.listSlices( list ) );
            return spot.listSlices( list ) ;
        }
        inline double calc() const
        {
            return apply( spot.calc() );
        }
        void printDebug(char *s) const
        {
            strcat(s, "(PenultSmooth)");
        }

    private:
        double apply( double s ) const
        {
            return fp->getValue(iPrice, simDateIdx, s , k);
        }

        const TreeSlice & spot;
        //const TreeSlice & price;
        IFinalPerfSP fp;
        int iPrice;
        int simDateIdx;
        double k;
    };

    // return intrinsic value at terminal (assuming all sampling are done.)
    virtual double getIntrinsic(double s, double k) = 0;
    
    virtual void upDatePerf(CTree1f* tree1f, int step){
                // nothing to do
    };

    virtual void upDateValues(int simDateIdx, const double s, double & p1, const double p2){
        // nothing to do, usually.
    };
    
    class upDateValues_oper : public SliceMarker< upDateValues_oper >
    {
    public:
        upDateValues_oper(
                const IFinalPerfSP fp,
                int simDateIdx,
                const TreeSlice & spot, 
                const TreeSlice & price1, 
                const TreeSlice & price2)
                :
                fp( fp ),
                simDateIdx( simDateIdx ),
                spot( spot ),
                price1( price1 ),
                price2( price2 )
                {};

        // TreeSlice "expression template" primitives
        static const int sliceCount = 2;
        template< typename S >
        const S** listSlices(const S** list) const
        {
            return price2.listSlices( price1.listSlices( spot.listSlices( list ) ) );
        }
        inline double calc() const
        {
            return apply( spot.calc(), price1.calc(), price2.calc() );
        }
        void printDebug(char *s) const
        {
            strcat(s, "(PenultSmooth)");
        }

    private:
        double apply( const double s, double p1, const double p2 ) const
        {
            fp->upDateValues(simDateIdx, s , p1, p2);
            return p1;
        }

        const TreeSlice & spot;
        const TreeSlice & price1;
        const TreeSlice & price2;
        IFinalPerfSP fp;
        int iPrice;
        int simDateIdx;
        double k;
    };

    // return true if it should have insert node level.
    virtual bool insNodeLevel(int simDateIdx, double& insLvl){
        return false;
    }

    virtual ~IFinalPerf() {};
};


// Build IFinalPerf via a Maker class.
// It's a bit like building the thing we need in 2 stages.
class PRODUCTS_DLL IFinalPerfMaker : virtual public IObject {
public:
    static CClassConstSP const TYPE;
    
    virtual IFinalPerfSP getFinalPerf(const CAsset*        asset,
                                      const DateTime       valueDate,
                                      const double         refLevel,
                                      const YieldCurveConstSP   discount,
                                      const InstrumentSettlementSP instSettle,
                                      const InstrumentSettlementSP premiumSettle,
                                      string                       ccyTreatment,
                                      const DateTimeArray& simDates) = 0;

    virtual DateTimeArray getMatDates() = 0;

    virtual DateTime getMatStartDate() = 0;

    virtual bool sensShift(Theta* shift, CAsset* asset, DateTime valueDate) = 0;

    // how many insert node is required for this performer?
    virtual int addInsNode(){       
        return 0;
    }

    // how many number of price state needs.
    virtual int addNumPrices() {
        return 0;
    }

    virtual bool hasCritDates(DateTimeArray& critDates){
        return false;
    }

    // return the simulation is still necessary or not.
    // return also endDate & barLvel.  Mainly for getEvents for Barrier!
    virtual bool isEnd(DateTime valDate, double spot, DateTime* endDate, double &barLvl) = 0;

    // has the barrier???
    virtual bool hasBarrier(){
        return false;
    }
    
    // for Legal Term Shift
    virtual bool useEcoBar() = 0;

    // return barrier schedule & intraday if it has.
    virtual ScheduleSP getBarrier(bool isEconomic, bool &isContinuous) const{
        return ScheduleSP(0);
    }

    virtual ~IFinalPerfMaker() {};

private:
    static void load(CClassSP& clazz);
};

typedef smartPtr<IFinalPerfMaker> IFinalPerfMakerSP;

//---------------------------------------------------------//
//  Vanilla Out at final
class PRODUCTS_DLL VanillaPerfMaker : public CObject,
                        virtual public IFinalPerfMaker {
public:
    static CClassConstSP const TYPE;
    friend class VanillaPerf;
    friend class VanillaPerfMakerHelper;
    
    // constructor, can access internally
    VanillaPerfMaker(IDoubleArrayModifierMakerSP  matPerformance,
                     DateTime               matDateTime); 

    virtual IFinalPerfSP getFinalPerf(const CAsset*        asset,
                                      const DateTime       valueDate,
                                      const double         refLevel,
                                      const YieldCurveConstSP   discount,
                                      const InstrumentSettlementSP instSettle,
                                      const InstrumentSettlementSP premiumSettle,
                                      string                       ccyTreatment,
                                      const DateTimeArray& simDates);

    virtual DateTimeArray getMatDates();

    virtual DateTime getMatStartDate();

    virtual bool sensShift(Theta* shift, CAsset* asset, DateTime valueDate);

    // validation & and make aop for the calling from the public interface.
    void validatePop2Object();

    // for Legal Term Shift
    virtual bool useEcoBar(){
        return true;  // no tweak, but return true to ignore in General.
    };

    // not implemented, yet.  Maybe no need....
    virtual bool isEnd(DateTime valDate, double spot, DateTime* endDate, double &barLvl){
        return false;        
    };

private:
    VanillaPerfMaker(): CObject(TYPE) {} // for reflection
    VanillaPerfMaker(const VanillaPerfMaker& rhs); // not implemented
    VanillaPerfMaker& operator=(const VanillaPerfMaker& rhs); // not implemented

public:
    IDoubleArrayModifierMakerSP  matPerformance;    // option at maturity
    DateTime                     matDateTime;
};
typedef smartPtr<VanillaPerfMaker> VanillaPerfMakerSP;

//---------------------------------------------------------//
//  Average Out at final
class PRODUCTS_DLL AveragePerfMaker : public CObject,
                        virtual public IFinalPerfMaker {
public:
    static CClassConstSP const TYPE;
    friend class AveragePerf;
    friend class AveragePerfMakerHelper;
    
    AvgOutPerfSP getAvgOutPerf();

    AveragePerfMaker(AvgOutPerfMakerSP aopMaker);

    virtual IFinalPerfSP getFinalPerf(const CAsset*        asset,
                                      const DateTime       valueDate,
                                      const double         refLevel,
                                      const YieldCurveConstSP   discount,
                                      const InstrumentSettlementSP instSettle,
                                      const InstrumentSettlementSP premiumSettle,
                                      string                       ccyTreatment,
                                      const DateTimeArray& simDates);

    virtual DateTimeArray getMatDates();

    virtual DateTime getMatStartDate();

    virtual bool sensShift(Theta* shift, CAsset* asset, DateTime valueDate);

    // validation & and make aop for the calling from the public interface.
    void validatePop2Object();

    // for Legal Term Shift
    virtual bool useEcoBar(){
        return true;  // no tweak, but return true to ignore in General.
    };

    // not implemented, yet.  Maybe no need....
    virtual bool isEnd(DateTime valDate, double spot, DateTime* endDate, double &barLvl){
        return false;        
    };

private:
    AveragePerfMaker(): CObject(TYPE) {} // for reflection
    AveragePerfMaker(const AveragePerfMaker& rhs); // not implemented
    AveragePerfMaker& operator=(const AveragePerfMaker& rhs); // not implemented

    AvgOutPerfMakerSP aopMaker;    

    AvgOutPerfSP           aop;     // need to keep this for roll(). $unregistered
public:

    double                 refLevel; // $unregistered
};
typedef smartPtr<AveragePerfMaker> AveragePerfMakerSP;


//---------------------------------------------------------//
//  Knock-In option
class PRODUCTS_DLL KnockInPerfMaker : public CObject,
                        virtual public IFinalPerfMaker {
public:
    static CClassConstSP const TYPE;
    friend class KnockInPerf;
    friend class KnockInPerfMakerHelper;
    
//    KnockInPerfMaker();

    virtual IFinalPerfSP getFinalPerf(const CAsset*        asset,
                                      const DateTime       valueDate,
                                      const double         refLevel,
                                      const YieldCurveConstSP   discount,           //not used
                                      const InstrumentSettlementSP instSettle,      //not used
                                      const InstrumentSettlementSP premiumSettle,   //not used
                                      string                       ccyTreatment,    //not used
                                      const DateTimeArray& simDates);

    virtual DateTimeArray getMatDates();

    virtual DateTime getMatStartDate();

    virtual bool sensShift(Theta* shift, CAsset* asset, DateTime valueDate);

    // how many insert node is required for this performer?
    virtual int addInsNode();

    // how many number of price state needs.
    virtual int addNumPrices();

    virtual bool hasCritDates(DateTimeArray& critDates);

    // has the barrier???
    virtual bool hasBarrier(){
        return true;
    }

    virtual bool useEcoBar();

    virtual ScheduleSP getBarrier(bool isEconomic, bool &isContinuous) const;

    // return the simulation is still necessary or not.
    // return also endDate & barLvel.  Mainly for getEvents for Barrier!
    virtual bool isEnd(DateTime valDate, double spot, DateTime* endDate, double &barLvl);

private:
    KnockInPerfMaker(): CObject(TYPE) {} // for reflection
    KnockInPerfMaker(const KnockInPerfMaker& rhs); // not implemented
    KnockInPerfMaker& operator=(const KnockInPerfMaker& rhs); // not implemented
  
    ScheduleSP barrier;
    ScheduleSP ecoBarrier;
    bool       intraDayMonitor;
    DateTime   matDateTime;
    double     refLevel; // $unregistered
    bool       isBreached;
    DateTime   breachDate;
        
    IDoubleArrayModifierMakerSP  matPerformance;    // option at maturity

};
typedef smartPtr<KnockInPerfMaker> KnockInPerfMakerSP;

DRLIB_END_NAMESPACE
#endif
