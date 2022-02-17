//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : VolMertonLV.cpp
//
//   Description : VolMertonLV
//
//   Author      : Francois Lu
//
//   Date        : 2 Feb 2004
//
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#include "edginc/VolMertonLVProcessed.hpp"
#include "edginc/EquityBase.hpp"
#include "edginc/VolRequestTime.hpp"

DRLIB_BEGIN_NAMESPACE

/*
void VolMertonLV::validatePop2Object(){
}
*/

/** Returns name of vol */
string VolMertonLV::getName() const{
    return name;
}

void VolMertonLV::acceptValueDateCollector(VolMertonLV* vol, 
                                    CValueDateCollector* collector)
{
    collector->valueDateValidate(vol->baseDate, 
                                 vol->getName());
}

/*** Creates a "Matrix Vol Surface" on provided Strikes  ***/
void VolMertonLV::createMertonVolSurface(VolSurfaceSP& volSurface){
    static const string routine("VolMertonLV::createMertonVolSurface");
    try{
        const DateTimeArray& dates = volSurface->getDates();
		const CDoubleArray& strikes = volSurface->getStrikes();
        CDoubleMatrix matrix(strikes.size(),
                             dates.size());
        for (int iStrike = 0; iStrike < strikes.size(); iStrike++) {
            for (int iMat = 0; iMat < dates.size(); iMat++) {
                matrix[iStrike][iMat] = impliedMertonVols[iStrike][iMat];
            }
        }
        /** for performance need constructor that takes in
            cached values (to do) */
        mertonSurface = VolSurfaceSP(
				            new VolSurface(getName(),
                           timeMetric.get(),
                           strikes,
                           matrix,
                           volSurface->getExpiries().get(),
                           baseDate));

        mertonSurfaceSpline = VolSplineSP(new VolSpline(*mertonSurface));

    } catch (exception& e){
        throw ModelException(e, routine);
    }
} 

/** Combines market and instrument data together to give a Processed Vol */
CVolProcessed* VolMertonLV::getProcessedVol(
    const CVolRequest* request,
    const CAsset* myAsset) const{
    static const string  routine("VolMerton::getProcessedVol");

    const static string method = "VolMertonLV::getProcessedVol";
    const double MAX_JUMP_WIDTH = 0.7, MIN_JUMP_WIDTH = 0.0;
    const double MAX_JUMP_MEAN = 2.0, MIN_JUMP_MEAN = -2.0;
    const double MAX_JUMP_RATE = 20.0, MIN_JUMP_RATE = 0.0;

    if (VolRequestTime::TYPE->isInstance(request)){
        return VolRequestTime::createVolProcessed(getName(), timeMetric);
    }
    
    // validate jump parameters
    if (JumpRate > MAX_JUMP_RATE || JumpRate < MIN_JUMP_RATE)
    {
        throw ModelException(method, "jump rate out of range : " +
                             Format::toString(MIN_JUMP_RATE) +
                             " to " + Format::toString(MAX_JUMP_RATE));
    }
    //JumpMean
    if (JumpMean > MAX_JUMP_MEAN || JumpMean < MIN_JUMP_MEAN)
    {
        throw ModelException(method, "jump mean out of range : " + 
                             Format::toString(MIN_JUMP_MEAN) +
                             " to " + Format::toString(MAX_JUMP_MEAN));
    }
    //JumpWidth
    if (JumpWidth > MAX_JUMP_WIDTH ||  JumpWidth < MIN_JUMP_WIDTH)
    {
        throw ModelException(method, "jump width out of range : " +
                             Format::toString(MIN_JUMP_WIDTH) +
                             " to " + Format::toString(MAX_JUMP_WIDTH));
    }


	const LocVolRequest* lv = dynamic_cast<const LocVolRequest*>(request);

	// keep data
    CVolParamConstSP volParam =CVolParamConstSP(mertonSurfaceSpline->createVolParam());

	VolMertonLVProcessed* ptr= new VolMertonLVProcessed(mertonSurfaceSpline.get(),
                                volParam,
                                lv, 
                                myAsset,
                                mertonSurface.get(),
								VolMertonLVSP(copy(this)));

    ptr->startDate=mertonSurface->getBaseDate();
    ptr->asset = AssetConstSP(copy(myAsset));

	const EquityBase* eq = dynamic_cast<const EquityBase*>(myAsset);
	if (!eq)
        throw ModelException(method, "only equity asset supported");

	ptr->equityYC = eq->getEquity()->getYC().getSP();
    return ptr;
}

// version for quanto and struck
CVolProcessed* VolMertonLV::getProcessedVol(
    const CVolRequest*,
    const CAsset*,
    const FXAsset*,
    const Correlation*) const
{
    throw ModelException("VolMertonLV::getProcessedVol",
						"quanto and ccy struck version not implemented");
}

IObject* VolMertonLV::defaultCtor(){
    return new VolMertonLV();
}

/** Invoked when Class is 'loaded' */
void VolMertonLV::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(VolMertonLV, clazz);
    SUPERCLASS(CVolBase);
    EMPTY_SHELL_METHOD(defaultCtor);
    FIELD(name, "Vol identifier");
    FIELD(ATMVol, "");
    FIELD(JumpRate, "");
    FIELD(JumpMean, "");
    FIELD(JumpWidth, "");
    FIELD(impliedMertonVols,  "impliedMertonVols");
    FIELD_MAKE_OPTIONAL(impliedMertonVols);

    // transient        
    FIELD(timeMetric, "");
    FIELD_MAKE_TRANSIENT(timeMetric);
    FIELD(baseDate, "");
    FIELD_MAKE_TRANSIENT(baseDate);
    FIELD(mertonSurface, "");
    FIELD_MAKE_TRANSIENT(mertonSurface);
    FIELD(mertonSurfaceSpline, "");
    FIELD_MAKE_TRANSIENT(mertonSurfaceSpline);
}

VolMertonLV::VolMertonLV(): CVolBase(TYPE),
                            ATMVol(0.2), JumpRate(0.1),
                            JumpMean(-0.8), JumpWidth(0.4){}

void VolMertonLV::getMarket(const IModel*     model, 
							const MarketData* market){

    /** returns a constant reference to surface to be used for the backbone */
    VolSurfaceSP volSurface = VolSurfaceSP(VolSurfaceSP::dynamicCast(
					market->GetData(getName(),VolSurface::TYPE)));
    volSurface->getMarket(model,market);  // put holiday and so on into volsurface.
    baseDate = volSurface->getBaseDate();
    timeMetric = volSurface->getTimeMetric();

    createMertonVolSurface(volSurface);
}

CClassConstSP const VolMertonLV::TYPE =
CClass::registerClassLoadMethod("VolMertonLV", typeid(VolMertonLV), load);

bool VolMertonLVLinkIn(){
    return (VolMertonLV::TYPE && true);
}

DRLIB_END_NAMESPACE
