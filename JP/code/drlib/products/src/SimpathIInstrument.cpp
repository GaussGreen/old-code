//----------------------------------------------------------------------------
//
//   Group       : QR&D
//
//   Filename    : SimPathIInstrument.cpp
//
//   Description : Pseudo-instrument for RM applications
//
//   Author      : Afshin Bayrooti
//
//   Date        : October 2005
//
//
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#include "edginc/SimpleEquity.hpp"
#include "edginc/SimpathIInstrumentMC.hpp"

#include "edginc/FXAsset.hpp"
#include "edginc/ICDSParSpreads.hpp"
#include "edginc/SVGenDiscFactor.hpp"
#include "edginc/SVGenExpectedDiscFactor.hpp"
#include "edginc/SVGenSurvivalDiscFactor.hpp"
#include "edginc/SVGenExpectedSurvivalDiscFactor.hpp"
#include "edginc/SVGenAggregatedSurvivalDiscFactor.hpp"
#include "edginc/SVGenExpectedEnergyFuture.hpp"
#include "edginc/SVGenExpectedBasisForward.hpp"

#include "edginc/SimpathIInstrument.hpp"

DRLIB_BEGIN_NAMESPACE

using namespace RM_Assets;

SimpathICurveInfo::SimpathICurveInfo() : CObject(TYPE){};

IObject* SimpathICurveInfo::defaultSimpathICurveInfo()
{
    return new SimpathICurveInfo();
}

// Invoked when Class is 'loaded'
void SimpathICurveInfo::load(CClassSP& clazz)
{
    REGISTER(SimpathICurveInfo, clazz);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(defaultSimpathICurveInfo);
    FIELD(name, "Curve name");
    FIELD(rngSlot, "Slot in random number generator");
    FIELD(maxDiffDate, "Last diffusion date");
    FIELD(maxCurveMaturity, "Curves will not be available beyond this date");
    clazz->setPublic(); // make visible to EAS/spreadsheet
    Addin::registerConstructor(
        "SimpathI_Curve_Info",
        Addin::UTILITIES,
        "Create SimpathI Curve Info object",
        TYPE);
}

CClassConstSP const SimpathICurveInfo::TYPE = CClass::registerClassLoadMethod(
    "SimpathICurveInfo", typeid(SimpathICurveInfo),
    SimpathICurveInfo::load);

DEFINE_TEMPLATE_TYPE(SimpathICurveInfoArray);

SimpathISpotInfo::SimpathISpotInfo() : CObject(TYPE) {};

IObject* SimpathISpotInfo::defaultSimpathISpotInfo()
{
    return new SimpathISpotInfo();
}

// Invoked when Class is 'loaded'
void SimpathISpotInfo::load(CClassSP& clazz)
{
    REGISTER(SimpathISpotInfo, clazz);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(defaultSimpathISpotInfo);
    FIELD(name, "Spot asset name");
    FIELD(rngSlot, "Slot in random number generator");
    FIELD(maxDiffDate, "Last diffusion date");
    clazz->setPublic(); // make visible to EAS/spreadsheet
    Addin::registerConstructor(
        "SimpathI_Spot_Info",
        Addin::UTILITIES,
        "Create SimpathISpotInfo object",
        TYPE);
}

CClassConstSP const SimpathISpotInfo::TYPE = CClass::registerClassLoadMethod(
    "SimpathISpotInfo", typeid(SimpathISpotInfo),
    SimpathISpotInfo::load);

DEFINE_TEMPLATE_TYPE(SimpathISpotInfoArray);

SimpathIInstrument::SimpathIInstrument() : GenericNFactor(TYPE) {};

void SimpathIInstrument::initializeMaps() const
{
    // Store info structures in nice form for lookups
    for ( int i = 0; i < swapCurveInfo->size(); ++i ) {
        string key = (*swapCurveInfo)[i]->name;
        swapCurveInfoMap[key] = (*swapCurveInfo)[i];
        asset2rng[key] = swapCurveInfoMap[key]->getRngSlot();
    }
    for ( int j = 0; j < cdsCurveInfo->size(); ++j ) {
        string key = (*cdsCurveInfo)[j]->name;
        cdsCurveInfoMap[key] = (*cdsCurveInfo)[j];
        asset2rng[key] = cdsCurveInfoMap[key]->getRngSlot();
    }
    for ( int k = 0; k < fxInfo->size(); ++k ) {
        string key = (*fxInfo)[k]->name;
        fxInfoMap[key] = (*fxInfo)[k];
        asset2rng[key] = fxInfoMap[key]->getRngSlot();
    }
    for ( int m = 0; m < eqInfo->size(); ++m ) {
        string key = (*eqInfo)[m]->name;
        eqInfoMap[key] = (*eqInfo)[m];
        asset2rng[key] = eqInfoMap[key]->getRngSlot();
    }
    if (enrgCurveInfo.get()) { // energy curve info is optional
	    for ( int o = 0; o < enrgCurveInfo->size(); ++o ) {
		    string key = (*enrgCurveInfo)[o]->name;
            enrgCurveInfoMap[key] = (*enrgCurveInfo)[o];
            asset2rng[key] = enrgCurveInfoMap[key] ->getRngSlot();
	    }
    }
    if (basisCurveInfo.get()) { // basis curve info is optional
        for ( int p = 0; p < basisCurveInfo->size(); ++p ) {
            string key = (*basisCurveInfo)[p]->name;
            basisCurveInfoMap[key] = (*basisCurveInfo)[p];
            asset2rng[key] = basisCurveInfoMap[key]->getRngSlot();
        }
    }
}

void SimpathIInstrument::classifyAssets(
    IMultiMarketFactorsSP   assets, // (I)
    vector<YieldCurveConstSP>& ycAsset,
    vector<FXAssetConstSP>& fxAsset,
    vector<ICDSParSpreadsConstSP>& cdsAsset,
	vector<SimpleEquityConstSP>& eqAsset,
	vector<EnergyFuturesCurveConstSP>& enrgAsset,
    vector<IBasisIndexCurveConstSP>& basisAsset,
    map<string,int>& fxAssetIDMap)
{
    // Classify the assets
    static CClassConstSP iYCType = CClass::forName("IYieldCurve");
    static CClassConstSP iFXAsset = CClass::forName("FXAsset");
    static CClassConstSP iCDS = CClass::forName("ICDSParSpreads");
	static CClassConstSP iEQAsset = CClass::forName("SimpleEquity");
	static CClassConstSP iENGAsset = CClass::forName("EnergyFuturesCurve");
    static CClassConstSP iBSSAsset = CClass::forName("IBasisIndexCurve");

    for ( int i = 0; i < assets->NbAssets(); ++i ) {
        IMarketFactorConstSP factor = assets->getFactor(i);
        CClassConstSP type = factor->getClass();
        if ( iYCType->isAssignableFrom( type ) ) {
            YieldCurveConstSP yc = YieldCurveConstSP::dynamicCast( factor );
            ycAsset.push_back( yc );
        } else if ( iFXAsset->isAssignableFrom( type ) ) {
            FXAssetConstSP fx = FXAssetConstSP::dynamicCast( factor );
            fxAsset.push_back( fx );
            ycAsset.push_back( fx->getRiskCcy() );
            assert(fxAssetIDMap.find(fx->getRiskCcy()->getCcy()) == fxAssetIDMap.end()); // assert we don't override
            fxAssetIDMap.insert(make_pair(fx->getRiskCcy()->getCcy(), fxAssetIDMap.size()));
        } else if ( iCDS->isAssignableFrom( type ) ) {
            ICDSParSpreadsConstSP parCDS = ICDSParSpreadsConstSP::dynamicCast( factor );
            cdsAsset.push_back( parCDS );
		} else if ( iEQAsset->isAssignableFrom( type ) ) {
			SimpleEquityConstSP eq = SimpleEquityConstSP::dynamicCast( factor );
			eqAsset.push_back( eq );
			// TO DO:  look this over
			fxAssetIDMap.insert(make_pair(eq->getName(), fxAssetIDMap.size()));
		} else if ( iENGAsset->isAssignableFrom( type ) ) {
			EnergyFuturesCurveConstSP enrgFutureCurve = EnergyFuturesCurveConstSP::dynamicCast( factor );
			enrgAsset.push_back( enrgFutureCurve );
        } else if ( iBSSAsset->isAssignableFrom( type ) ) {
            IBasisIndexCurveConstSP basisCurve = IBasisIndexCurveConstSP::dynamicCast( factor );
            basisAsset.push_back( basisCurve );
		} else {
            throw ModelException("SimpathIInstrument::createProduct",
                "Asset type " + type->getName() +
                " is not supported");
        }
    }
}

void SimpathIInstrument::insertDFSVGen(
    const DateTimeArray& timeLine,
    const vector<YieldCurveConstSP>& ycAsset,
    vector<IStateVariableGen*>& gens ) const
{
    // For each yield curve, request discount factors from t0 to each
    // required date in time line
    for ( size_t i = 0; i < ycAsset.size(); ++i ) {
        SimpathICurveInfoMap::const_iterator I =
            swapCurveInfoMap.find( ycAsset[i]->getName() );
        if ( I == swapCurveInfoMap.end() ) {
            throw ModelException( "No swap curve info available for "
                + ycAsset[i]->getName() );
        }
        const SimpathICurveInfoSP& info = I->second;
        // Only diffuse as far as needed
        DateTimeArray dates( info->maxDiffDate.getPastDates( timeLine ) );
        gens.push_back(
            new SVGenDiscFactor(
                valueDate,
                ycAsset[i],
                dates ) );
    }
}


void SimpathIInstrument::insertExpectedDFSVGen(
    const DateTime& currentDate,
    vector<YieldCurveConstSP>& ycAsset,
    vector<IStateVariableGen*>& gens ) const
{
    DateTimeArray ycPillarDate( yieldCurveOffsets->size() );
    for ( int j = 0; j < yieldCurveOffsets->size(); ++j ) {
        MaturityPeriod period( (*yieldCurveOffsets)[j] );
        ycPillarDate[j] = period.toDate( currentDate );
    }

    for ( size_t n = 0; n < ycAsset.size(); ++n ) {
        SimpathICurveInfoMap::const_iterator I =
            swapCurveInfoMap.find( ycAsset[n]->getName() );
        if ( I == swapCurveInfoMap.end() ) {
            throw ModelException( "No swap curve info available for "
                + ycAsset[n]->getName() );
        }
        const SimpathICurveInfoSP& info = I->second;
        if ( currentDate <= info->maxDiffDate ) { // Check if need to diffuse further
            // Only diffuse as far as needed
            DateTimeArray pillarDates( info->maxCurveMaturity.getPastDates( ycPillarDate ) );
            gens.push_back(
                new SVGenExpectedDiscFactor(
                currentDate,
                currentDate,
                ycAsset[n],
                pillarDates,
                true) ); // computeLog
        }
    }
}

void SimpathIInstrument::insertSurvivalDFSVGen(
    const DateTimeArray& timeLine,
    const vector<ICDSParSpreadsConstSP>& cdsAsset,
    vector<IStateVariableGen*>& gens ) const
{
    // For each CDS curve, request survival discount factors from t0 to each
    // date in time line
    for ( size_t i = 0; i < cdsAsset.size(); ++i ) {
        SimpathICurveInfoMap::const_iterator I =
            cdsCurveInfoMap.find( cdsAsset[i]->getName() );
        if ( I == cdsCurveInfoMap.end() ) {
            throw ModelException( "No cds curve info available for "
                + cdsAsset[i]->getName() );
        }
        const SimpathICurveInfoSP& info = I->second;
        // Only diffuse as far as needed
        DateTimeArray dates( info->maxDiffDate.getPastDates( timeLine ) );
        gens.push_back(
            new SVGenSurvivalDiscFactor(
                valueDate,
                cdsAsset[i],
                dates ) );
    }
}

void SimpathIInstrument::insertExpectedSurvivalDFSVGen(
    const DateTime& currentDate,
    const vector<ICDSParSpreadsConstSP>& cdsAsset,
    vector<IStateVariableGen*>& gens ) const
{
    DateTimeArray cdsPillarDate( cdsOffsets->size() );
    for ( int k = 0; k < cdsOffsets->size(); ++k ) {
        MaturityPeriod period( (*cdsOffsets)[k] );
        cdsPillarDate[k] = period.toDate( currentDate );
    }
    // For each CDS curve, request expected survival discount factors
    // between t and each cds pillar date
    for ( size_t m = 0; m < cdsAsset.size(); ++m ) {
        SimpathICurveInfoMap::const_iterator I =
            cdsCurveInfoMap.find( cdsAsset[m]->getName() );
        if ( I == cdsCurveInfoMap.end() ) {
            throw ModelException( "No cds curve info available for "
                + cdsAsset[m]->getName() );
        }
        const SimpathICurveInfoSP& info = I->second;
        if ( currentDate <= info->maxDiffDate ) { // Check if need to diffuse further
            // Only diffuse as far as needed
            DateTimeArray pillarDates( info->maxCurveMaturity.getPastDates( cdsPillarDate ) );
            gens.push_back(
                new SVGenExpectedSurvivalDiscFactor(
                currentDate,
                currentDate,
                cdsAsset[m],
                pillarDates,
                true) ); // computeLog
        }
    }
}

void SimpathIInstrument::insertExpectedEnergyFutureSVGen(
    const DateTime& currentDate,
    const vector<EnergyFuturesCurveConstSP>& enrgAsset,
    vector<IStateVariableGen*>& gens ) const
{
	// For each energy future curve, request expected future prices between t and each pillar date
    for ( size_t m = 0; m < enrgAsset.size(); ++m )
	{
		EnergyFuturesCurveConstSP thisEnergy = enrgAsset[m];
		EnergyUnderlyerConstSP thisEnergyUnderlyer = thisEnergy->getEnergyUnderlyer();

        SimpathICurveInfoMap::const_iterator I = enrgCurveInfoMap.find( enrgAsset[m]->getName() );
        if ( I == enrgCurveInfoMap.end() ) {
            throw ModelException( "No energy future curve info available for " + enrgAsset[m]->getName() );
        }

        const SimpathICurveInfoSP& info = I->second;
        if ( currentDate <= info->maxDiffDate ) { // Check if need to diffuse further

			// prepare the energy future pillar dates
			//DateTimeArray enrgPillarDates( thisEnergy->getSize() );
			//for (int i = 0; i < thisEnergy->getSize(); ++i) {
			//	EnergyContractLabel label(thisEnergy->getExpiryLabel(i));
			//	enrgPillarDates[i] = thisEnergyUnderlyer->expiryDate(label);
			//}
            const DateTimeArray& enrgPillarDates = thisEnergy->getFutureMaturityDates();

            // Only diffuse as far as needed
			DateTimeArray pillarDates = DateTime::getInclusiveDates(
				currentDate,
				info->maxCurveMaturity,
				enrgPillarDates);
            gens.push_back(
                new SVGenExpectedEnergyFuture(
					currentDate,
					enrgAsset[m],
					pillarDates,
					true)
				); // computeLog
        }
    }
}

//takes current date; adds a list of offsets, but stops at maxDate
static DateTimeArray addOffsets(  const DateTime& currentDate,
                                  const StringArray&  offsets)
{
    DateTimeArray dates;
    for (int t = 0; t < offsets.size(); ++t )
    {
        MaturityPeriod period( offsets[t] );
        dates.push_back(period.toDate( currentDate ));
    }
    return dates;
}

static DateTime shiftDate(const DateTime& date, long offset)
{
    return DateTime(date.getDate() + offset, date.getTime());
}
//takes current date; adds a list of offsets, but stops at maxDate
static DateTimeArray addOffsets(  const DateTime& currentDate,
                                  const vector<long>&  offsets)
{
    DateTimeArray dates;
    for (size_t t = 0; t < offsets.size(); ++t )
        dates.push_back(shiftDate(currentDate, offsets[t]));
    return dates;
}

// same function but works on DateTimeArray
static DateTimeArray addOffsets(  const DateTimeArray& dates,
                                  const vector<long>&  offsets)
{
    vector<DateTimeArray> storage(dates.size());
    for(size_t i=0; i != dates.size(); ++i)
        storage[i]=addOffsets(dates[i], offsets);

    return DateTime::merge(storage);
}

//takes current date; adds a list of offsets, but stops at maxDate
DateTimeArray SimpathIInstrument::calcForwardDates(  const DateTime& currentDate,
                                        const StringArray&  offsets,
                                        const DateTime& maxDate)
{
    DateTimeArray result;
    for(int i=0; i < offsets.size(); ++i)
    {
        MaturityPeriod period( offsets[i] );
        DateTime date(period.toDate( currentDate ));

        if (date <= maxDate)
            result.push_back(date);
        else
            break;
    }
    return result;
}

//takes current date; adds a list of offsets, but stops at maxDate
DateTimeArray SimpathIInstrument::calcForwardDates(  const DateTime& currentDate,
                                        const vector<long>&  offsets,
                                        const DateTime& maxDate)
{
    DateTimeArray result;
    for(size_t i=0; i < offsets.size(); ++i)
    {
        DateTime date = shiftDate(currentDate, offsets[i]);
        if (date <= maxDate)
            result.push_back(date);
        else
            break;
    }
    return result;
}

static unsigned long parseForMode(const char* spec)
{
    unsigned long mode = 0;

    if (*spec == ':')
    {
        string s(++spec);
        string::size_type head=0, tail;
        while (head < s.length())
        {
            tail = s.find('+', head);
            if (tail == string::npos)
                tail = s.length();

            if (strncmp(s.c_str()+head, "PATH", tail-head) == 0)
                mode |= IMCWriter::PathMode;
            else if (strncmp(s.c_str()+head, "SLICE", tail-head) == 0)
                mode |= IMCWriter::SliceMode;
            head = tail + 1;
        }
    }

    return mode ? mode : IMCWriter::PathMode;
}

IMCProduct* SimpathIInstrument::createProduct(const MonteCarlo* ) const
{
//     int i;
//     size_t j;

    vector<YieldCurveConstSP> ycAsset;
    vector<FXAssetConstSP> fxAsset;
    vector<ICDSParSpreadsConstSP> cdsAsset;
	vector<SimpleEquityConstSP> eqAsset;
	vector<EnergyFuturesCurveConstSP> enrgAsset;
    vector<IBasisIndexCurveConstSP> basisAsset;
    map<string,int> fxAssetIDMap;
    // go over all assets in this->assets and classify them
    // creates a map to lookup FX by its name in fxAsset[]
    classifyAssets( this->assets, ycAsset, fxAsset, cdsAsset, eqAsset, enrgAsset, basisAsset, fxAssetIDMap);

    // Hold data in more convenient form
    // initialize several maps: swapCurveInfoMap, cdsCurveInfoMap, enrgCurveInfoMap, basisCurveInfoMap, fxInfoMap, eqInfoMap to lookup info by name
    initializeMaps();
    DateTime zeroDate(0, 0);
    // Create the timeline
    DateTimeArray timeLine( timeLineOffsets->size() );
    for (int i = 0; i < timeLineOffsets->size(); ++i ) {
        MaturityPeriod tpPeriod( (*timeLineOffsets)[i]);
        timeLine[i] = tpPeriod.toDate( valueDate );
    }
    
    ASSERT(timeLine.size() > 1); // For ccy/basis fix we need [1] element (as [0] maybe today)
    
//  TODO: this does not work properly yet...
//    if (timeLine[0] != valueDate)
//        timeLine.insert(timeLine.begin(), valueDate);

    vector<long> ycdayoffsets( yieldCurveOffsets->size() );
    for (int i = 0; i < yieldCurveOffsets->size(); ++i ) {
        MaturityPeriod period( (*yieldCurveOffsets)[i] );
        ycdayoffsets[i] = period.toDate( zeroDate ).daysDiff( zeroDate );
    }

    vector<long> cdsdayoffsets( cdsOffsets->size() );
    for (int i = 0; i < cdsOffsets->size(); ++i ) {
        MaturityPeriod period( (*cdsOffsets)[i] );
        cdsdayoffsets[i] = period.toDate( zeroDate ).daysDiff( zeroDate );
    }

    vector<long> basisdayoffsets(basisCurveOffsets.get() ? basisCurveOffsets->size() : 0);
    for (size_t i = 0; i < basisdayoffsets.size(); ++i ) {
        MaturityPeriod period( (*basisCurveOffsets)[i] );
        basisdayoffsets[i] = period.toDate( zeroDate ).daysDiff( zeroDate );
    }

    SimPathICollectionMap collectionMap;
    
    // Request fx rate for each date in time line.
    // CURRENTLY, WILL NEED A SEPARATE SVGenSpot FOR EACH FX RATE IF WE WANT TO SPECIFY
    // HOW FAR EACH ASSET SHOULD BE DIFFUSED.
    SVGenSpotSP fxGenGlobal (new SVGenSpot(
                                fxAsset.size() + eqAsset.size(),
                                timeLine ));

    // Loop thru all energy future curves:
	for (size_t j = 0; j < enrgAsset.size(); ++j )
	{
		EnergyFuturesCurveConstSP thisEnergy = enrgAsset[j];
		EnergyUnderlyerConstSP thisEnergyUnderlyer = thisEnergy->getEnergyUnderlyer();
        const string & thisEnergyName = thisEnergy->getName();

        SimpathICurveInfoMap::const_iterator I = enrgCurveInfoMap.find( thisEnergyName );
        if ( I == enrgCurveInfoMap.end() )
            throw ModelException( "No energy future curve info available for " + thisEnergyName );

        const SimpathICurveInfoSP& info = I->second;
#if 0        
        cerr << "Adding energy: " << thisEnergyName << " maxDiff= " << info->maxDiffDate.toString()  << " (" << info->maxDiffDate.getDate() << ")  " \
			 << " maxCurveMat= " << info->maxCurveMaturity.toString() << " (" << info->maxCurveMaturity.getDate() << ")" << endl;
#endif		
        // Only diffuse as far as needed
        DateTimeArray dates( info->maxDiffDate.getPastDates( timeLine ) );

		// prepare the energy future pillar date vectors
        //const DateTimeArray& enrgPillarDates = thisEnergy->getFutureMaturityDates();
        DateTimeArray pillarDatesAll(thisEnergy->getFutureMaturityDates());
        DateTimeArray enrgPillarDates = info->maxCurveMaturity.getPastDates(pillarDatesAll); // chop by max curve maturity
        vector<TDate> enrgPillarDatesAsTDateVec( enrgPillarDates.getLength() );
        for (int i = 0; i < enrgPillarDates.getLength(); ++i)
            enrgPillarDatesAsTDateVec[i] = (TDate)enrgPillarDates[i].getDate();

        // no records created yet, so create one and alias it:
        SimpathIEnergyCollectionSP pCollection (new SimpathIEnergyCollection(*info, enrgPillarDatesAsTDateVec, enrgAsset[j].get()));

        // For each energy future curve, request future price between t and each pillar date
        for (int i = 0; i < dates.size() ; ++i )
        {
			// Only diffuse as far as needed
            DateTime currentDate = dates[i];
			DateTimeArray pillarDates = DateTime::getInclusiveDates(
				currentDate,
				info->maxCurveMaturity,
				enrgPillarDates);

            // check if current date pass maxCurveMaturity or not:
            if (!pillarDates.empty())
            {
                // create SV generator
                SVGenExpectedEnergyFutureSP futureGenSP(
					    new SVGenExpectedEnergyFuture(
						    currentDate,
						    enrgAsset[j],
						    pillarDates,
						    true) // compute log
					    );

    			// add the SV generator to the energy collection object
                pCollection->vExpFpGen.push_back(futureGenSP);
            }
        }

        if(collectionMap.count(thisEnergyName))
            throw ModelException("Energy future name " + thisEnergyName + " is repeated.");

        collectionMap[thisEnergyName] = pCollection;
	}


    for(size_t j = 0; j<ycAsset.size(); ++j)
    {
        SimpathICurveInfoMap::const_iterator I =
            swapCurveInfoMap.find( ycAsset[j]->getName() );
        if ( I == swapCurveInfoMap.end() ) {
            throw ModelException( "No swap curve info available for "
                + ycAsset[j]->getName() );
        }
        const SimpathICurveInfoSP& info = I->second;

        vector<FXAssetConstSP>::const_iterator found = fxAsset.begin();
        while ( found != fxAsset.end() && (*found)->getRiskCcyIsoCode() != ycAsset[j]->getCcy())
            ++found;

//         cout << "FX " << j << " maxDiffDate= " << info->maxDiffDate.toString() << endl;
//         cout << "FX " << j << " maxCurve= " << info->maxCurveMaturity.toString() << endl;

        // Only diffuse as far as needed
        DateTimeArray dates( info->maxDiffDate.getPastDates( timeLine ) );

        // It is possible that maxMaturity was set to 0 and we don't have any elements here. As we will need this asset, we create it by asking for one fake date (well, the date is not totally fake, we need something on the timeline, but not "today"; we don't use .back() as that would incur diffusing more than needed in the worst case.
        if (dates.empty()) {
            std::cerr << ( "All dates are past maxDiffDate (" + info->maxDiffDate.toString() +") [" + Format::toString(info->maxDiffDate.getDate()) + "] for ycAsset= " + ycAsset[j]->getName() );
            std::cerr << " Fix activated." << std::endl;
            dates.push_back(timeLine[1]); 
        }
        // no records created yet, so create one and alias it
        SimpathICurrencyCollectionSP pCollection (new SimpathICurrencyCollection(*info, ycdayoffsets, ycAsset[j].get(), found == fxAsset.end() ? 0 : found->get()));

        pCollection->dfGen = SVGenDiscFactorSP(new SVGenDiscFactor(
                                    valueDate,
                                    ycAsset[j],
                                    dates )) ;

        pCollection->setFxEqGen(fxGenGlobal);

        if (fxAssetIDMap.count(ycAsset[j]->getName()))
            pCollection->fxPosition = fxAssetIDMap[ycAsset[j]->getName()];
        else
            pCollection->fxPosition = -1;


        // For each yield curve, request expected discount factors between t and
        // each pillar date
        for (int i = 0; i < dates.size() ; ++i )
        {
            DateTime currentDate = dates[i];
            DateTimeArray pillarDates = calcForwardDates(currentDate, ycdayoffsets, info->maxCurveMaturity);
            if (!pillarDates.empty())
            {
                pCollection->vExpDFGen.push_back(SVGenExpectedDiscFactorSP(new SVGenExpectedDiscFactor(
                                                            currentDate,
                                                            currentDate,
                                                            ycAsset[j],
                                                            pillarDates,
                                                            true) ) ); // computeLog
            }
        }

        collectionMap[ycAsset[j]->getName()] = pCollection;
    }

    // Loop thru all basis curves:
    for(size_t j = 0; j < basisAsset.size(); ++j)
    {
        SimpathICurveInfoMap::const_iterator I =
            basisCurveInfoMap.find( basisAsset[j]->getName() );
        if ( I == basisCurveInfoMap.end() ) {
            throw ModelException( "No basis curve info available for "
                + basisAsset[j]->getName() );
        }
        const SimpathICurveInfoSP& info = I->second;

        // Only diffuse as far as needed
        DateTimeArray dates( info->maxDiffDate.getPastDates( timeLine ) );
        
        if (dates.empty()) {
            std::cerr << ( "All dates are past maxDiffDate (" + info->maxDiffDate.toString() +") [" + Format::toString(info->maxDiffDate.getDate()) + "] for basisAsset= " + basisAsset[j]->getName() );
            std::cerr << " Fix activated." << std::endl;
            dates.push_back(timeLine[1]); 
        }

        // no records created yet, so create one and alias it
        SimpathIBasisCollectionSP pCollection (new SimpathIBasisCollection(*info, basisdayoffsets, basisAsset[j].get()));

        // For each basis curve, request expected basis fwd sprds between t
        // and each pillar date
        for (int i = 0; i < dates.size() ; ++i )
        {
            DateTime currentDate = dates[i];
            DateTimeArray pillarDates = calcForwardDates(currentDate, basisdayoffsets, info->maxCurveMaturity);
            if (!pillarDates.empty())
            {
                pCollection->vExpBasisFwdSprdGen.push_back(SVGenExpectedBasisFwdSpreadSP(
                    new SVGenExpectedBasisFwdSpread(
                    basisAsset[j],
                    currentDate,
                    pillarDates))); 
            }
        }

        collectionMap[basisAsset[j]->getName()] = pCollection;
    }

    // inserting equities
    DateTimeArraySP timeLineSP(new DateTimeArray(timeLine));
    DateTimeArraySP forwardDatesSP(new DateTimeArray(addOffsets(timeLine, cdsdayoffsets)));
    for (size_t j = 0; j < eqAsset.size(); ++j ) {
        SimpathISpotInfoMap::const_iterator I =
            eqInfoMap.find( eqAsset[j]->getName() );
        if ( I == eqInfoMap.end() ) {
            throw ModelException( "No equity info available for "
                + eqAsset[j]->getName() );
        }
        const SimpathISpotInfoSP& info = I->second;
        SimpathIEquityCollectionSP pCollection(
            new SimpathIEquityCollection(*info, eqAsset[j].get/*->YCIsoCode*/()));

        pCollection->setFxEqGen(fxGenGlobal);

        if (fxAssetIDMap.count(eqAsset[j]->getName()))
            pCollection->eqPosition = fxAssetIDMap[eqAsset[j]->getName()];
        else {
            throw ModelException("No equity info in internal map for "
                + eqAsset[j]->getName() );
        }
        collectionMap[eqAsset[j]->getName()] = pCollection;
    }

    SpotIdxArraySP    sdfIdxSP(new SpotIdxArray);
    FwdIdxArraySP     esdfReqIdxSP(new FwdIdxArray);
    FwdIdxArraySP     esdfForIdxSP(new FwdIdxArray);

    for(size_t j = 0; j<cdsAsset.size(); ++j)
    {
        SimpathICurveInfoMap::const_iterator I =
            cdsCurveInfoMap.find( cdsAsset[j]->getName() );
        if ( I == cdsCurveInfoMap.end() ) {
            throw ModelException( "No cds curve info available for "
                + cdsAsset[j]->getName() );
        }
        const SimpathICurveInfoSP& info = I->second;
/*        cout << "FX " << j << " maxDiffDate= " << info->maxDiffDate.toString() << endl;
        cout << "FX " << j << " maxCurve= " << info->maxCurveMaturity.toString() << endl;*/

        // no records created yet, so create one and alias it
        SimpathICreditCollectionSP pCollection (new SimpathICreditCollection(*info, cdsdayoffsets, cdsAsset[j].get()));
        pCollection->aggSDFGen = SVGenAggregatedSurvivalDiscFactorSP(new SVGenAggregatedSurvivalDiscFactor(
                            timeLineSP,//common for all CR
                            sdfIdxSP,
                            timeLineSP,//common for all CR
                            esdfReqIdxSP,
                            forwardDatesSP,
                            esdfForIdxSP,
                            info->maxDiffDate, // specific for assets
                            info->maxCurveMaturity,// specific for assets
                            cdsAsset[j], // specific for assets
                            true));

        if(collectionMap.count(cdsAsset[j]->getName()))
            throw ModelException("Asset Name " + cdsAsset[j]->getName() + " is repeated.");
        collectionMap[cdsAsset[j]->getName()] = pCollection;
    }
    // setup terminal writer

    MCWriterDataHolder writableData;

        for (SimPathICollectionMapIterator it = collectionMap.begin();
             it != collectionMap.end();
             ++it)
            {
                ISimPathICollectionSP & collection = it->second;
                collection->initAssetInfo(writableData);
            }

        std::vector<TDate> dates;
        dates.reserve(1 + timeLine.size());
        dates.push_back(valueDate.getDate()); // also writes out initial data
        for(int i=0; i<timeLine.size(); ++i)
            dates.push_back(timeLine[i].getDate());
#if 0
        // The MCWriter follows Decorator pattern
        // that means that you can chain writers to produce different combinations of output
        // Produce XML output to adhere rules of models.exe
        // But it is possible to to chain it with other filters.
        // In any case, the IMCWriterSP(new MCWriter(writableData, dates)) should be the last in the chain of filters
        stateIMCWriter = IMCWriterSP(new MCWriter(writableData, dates));

        string filename = (filenamePrefix == "") ? "SimpathIDiffusion" : filenamePrefix;

        string::size_type head=0, tail;
        while (head < format.length())
        {
            tail = format.find('|', head);
            if (tail == string::npos)
                tail = format.length();

            string spec = format.substr(head, tail - head);
            if (spec.find("XML", 0) == 0)
            {
                stateIMCWriter = IMCWriterSP(new MCXMLWriter(filename, stateIMCWriter));
            }
            else if (spec.find("DEBUG") == 0)
            {
                unsigned  long mode = parseForMode(spec.c_str() + strlen("DEBUG"));
                stateIMCWriter = IMCWriterSP(new MCDebugWriter(filename, stateIMCWriter, mode));
            }
            else if (spec.find("SANE") == 0)
            {
                unsigned  long mode = parseForMode(spec.c_str() + strlen("SANE"));
                stateIMCWriter = IMCWriterSP(new MCDebugWriter(filename, stateIMCWriter, mode, true));
            }
            else if (spec.find("BINARY") == 0)
            {
                unsigned  long mode = parseForMode(spec.c_str() + strlen("BINARY"));
                stateIMCWriter = IMCWriterSP(new MCBinaryWriter(filename, stateIMCWriter, mode));
            }
            
            head = tail + 1;
        }
#endif

#if 0
    std::vector<TDate> dates;
    for(int i=0; i<timeLine.size(); ++i)
        dates.push_back(timeLine[i].getDate());
#endif
    
    std::pair<IMCWriterSP, bool> stateIMCWriter = createWriter(writableData, dates);

    // create asset data entries (one per equity)
    // TODO: Gens for equity
    // insertEquityAD ( timeLine, ycAsset, mcAssetData );

    return new SimpathIInstrumentMC(
        assets.get(),
        valueDate,
        discount.get(),
        collectionMap,
        timeLine,
        stateIMCWriter.first,   // IMCWriterSP
        asset2rng,
		stateIMCWriter.second); // doPrice
}

// The MCWriter follows Decorator pattern
// that means that you can chain writers to produce different combinations of output
// Produce XML output to adhere rules of models.exe
// But it is possible to to chain it with other filters.
// In any case, the IMCWriterSP(new MCWriter(writableData, dates)) should be the last in the chain of filters

std::pair<IMCWriterSP, bool> SimpathIInstrument::createWriter(const MCWriterDataHolder& writableData, const std::vector<TDate>& dates) const
{
    IMCWriterSP stateIMCWriter = IMCWriterSP(new MCWriter(writableData, dates));
	bool doPrices = false; // true => XML writer was requested; 
    string filename = (filenamePrefix == "") ? "SimpathIDiffusion" : filenamePrefix;

    string::size_type head=0, tail;
    while (head < format.length())
    {
        tail = format.find('|', head);
        if (tail == string::npos)
            tail = format.length();

        string spec = format.substr(head, tail - head);
        if (spec.find("XML", 0) == 0)
        {
            // stateIMCWriter = IMCWriterSP(new MCXMLWriter(filename, stateIMCWriter));
            // stateIMCWriter = IMCWriterSP(new MCRegTestWriter(stateIMCWriter));
			doPrices = true;
        }
        else if (spec.find("DEBUG") == 0)
        {
            unsigned  long mode = parseForMode(spec.c_str() + strlen("DEBUG"));
            stateIMCWriter = IMCWriterSP(new MCDebugWriter(filename, stateIMCWriter, mode));
        }
	else if (spec.find("SANE") == 0)
	{
	    unsigned  long mode = parseForMode(spec.c_str() + strlen("SANE"));
	    stateIMCWriter = IMCWriterSP(new MCDebugWriter(filename, stateIMCWriter, mode, true));
	}
        else if (spec.find("BINARY") == 0)
        {
            unsigned  long mode = parseForMode(spec.c_str() + strlen("BINARY"));
            stateIMCWriter = IMCWriterSP(new MCBinaryWriter(filename, stateIMCWriter, mode));
        }
        
        head = tail + 1;
    }

    return make_pair(stateIMCWriter, doPrices);
}
IObject* SimpathIInstrument::defaultSimpathIInstrument()
{
    return new SimpathIInstrument();
}

// Invoked when Class is 'loaded'
void SimpathIInstrument::load(CClassSP& clazz)
{
    REGISTER(SimpathIInstrument, clazz);
    SUPERCLASS(GenericNFactor);
    EMPTY_SHELL_METHOD(defaultSimpathIInstrument);
    FIELD(timeLineOffsets, "rigid timeline specified as offsets to today: 1D, 2W, etc.");
    FIELD(yieldCurveOffsets, "yield curve offsets relative to each point in time line: 1M, 5Y, etc.");
    FIELD(basisCurveOffsets, "basis curve offsets relative to each point in time line: 1M, 5Y, etc.");
    FIELD(cdsOffsets, "cds offsets relative to each point in time line: 1M, 5Y, etc.");
    FIELD(swapCurveInfo, "Swap curve information");
    FIELD(basisCurveInfo, "Basis spread curve information");
    FIELD(cdsCurveInfo, "CDS curve information");
	FIELD(enrgCurveInfo, "Energy future curve information");
    FIELD(fxInfo, "fx information");
    FIELD(eqInfo, "eq information");
    FIELD(filenamePrefix, "filename to dump results of diffusion per path and datetime");
    FIELD(format, "Valid choices are DEBUG/XML/BINARY");

    FIELD_MAKE_OPTIONAL(enrgCurveInfo);
    FIELD_MAKE_OPTIONAL(basisCurveInfo);
    FIELD_MAKE_OPTIONAL(basisCurveOffsets);
    FIELD_MAKE_OPTIONAL(filenamePrefix);
    FIELD_MAKE_OPTIONAL(format);
    clazz->setPublic(); // make visible to EAS/spreadsheet
}


CClassConstSP const SimpathIInstrument::TYPE = CClass::registerClassLoadMethod(
    "SimpathIInstrument", typeid(SimpathIInstrument),
    SimpathIInstrument::load);

// for class loading (avoid having header file)
bool SimpathIInstrumentLoad() {
    return (SimpathIInstrument::TYPE != 0);
}

DRLIB_END_NAMESPACE

