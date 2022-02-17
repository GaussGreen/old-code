//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : CreditEngine.cpp
//
//   Description : Credit engine
//
//   Author      : Ning Shen
//
//   Date        : 10 June 2002
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/Nrfns.hpp"
#include "edginc/Maths.hpp"
#include "edginc/CreditEngine.hpp"
#include "edginc/ATMVolRequest.hpp"
#include <time.h>
#include <algorithm>

DRLIB_BEGIN_NAMESPACE

CClassConstSP const CreditEngine::TYPE = CClass::registerClassLoadMethod(
							"CreditEngine", typeid(CreditEngine), load);

CreditEngine::CreditEngine() : CObject(TYPE) 
{
    useThetaFwdRate = false;
}

void CreditEngine::load(CClassSP& clazz)
{
        REGISTER(CreditEngine, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultCreditEngine);
        clazz->setPublic(); // make visible to EAS/spreadsheet
        FIELD(creditData, "credit data");
    	FIELD(market, "market cache"); 
        FIELD(marginAcct,"margin collateral account");
        FIELD(percentile, "percentile used for peak computation");
        FIELD(creditCcyName, "name of currency we wish to compute exposure in");
		FIELD(creditSpreads, "optional field: required if cvr is to be computed.");
		FIELD_MAKE_OPTIONAL(creditSpreads);
		FIELD(outputGrids, "if true, returns all values on all paths");
		FIELD(useThetaFwdRate, "optional field: true = use fwd rate when shift value date; false(default)=shift whole yield curve.");
		FIELD_MAKE_OPTIONAL(useThetaFwdRate);
        FIELD(model, "used for getting yield curve from market");
		FIELD_MAKE_OPTIONAL(model);

        Addin::registerClassObjectMethod("CREDIT_ENGINE",
                                         Addin::RISK, // CREDIT ??
                                         "simulate path value grids, netting performed for array of instruments",
                                         TYPE,
                                         false,
                                         Addin::returnHandle,
                                         (Addin::ObjMethod*)computePathValuesAddin);
}

// EdrAction version of addin
IObjectSP CreditEngine::run()
{
	return runTest();
}

/** for regression run */
IObjectSP CreditEngine::runTest() const
{
    return computePathValuesAddin( this );
}

IObjectSP CreditEngine::computePathValuesAddin(const CreditEngine* obj)
{
	clock_t startTime = clock();
    
    try{

    CreditPathValuesOutSP out;
    out = obj->computePathValues();
	if (!!out)
	{
		out->calcPeakAverageStdAndDRE(obj->market, obj->marginAcct, obj->percentile, 
			obj->creditCcyName,	obj->creditSpreads, obj->model);
	}
	else
	{
		throw ModelException("CreditEngine::computePathValuesAddin",
			"Unable to compute credit for one or more instruments");
	}
    if (!obj->outputGrids)
    {
        out->flushGrids();
		out->creditDebug->calcTime = 0; // hard code to avoid regtest difference
    }
	else
	{
		clock_t endTime = clock();

		// set calcTime field in out object
		out->creditDebug->calcTime = (double)(endTime-startTime)/CLOCKS_PER_SEC;
    } 

 

    
    return (IObjectSP) out;

    }     catch(ModelException& e){
        throw ModelException(e, "CreditEngine::computePathValuesAddin");
    }

}

/************************* local helpers **********/
class CreditDateSort
{
public:
    class DataPair{
    public:
        DateTime date;
        int      type;
        bool operator < (const DataPair& v) const{
            return (date < v.date); // time compared !!
        }
    };

    CreditDateSort(DateTimeArraySP pathDates, CIntArraySP pathDateTypes)
    {
        if (pathDates->size() == 0)
            return;
        if (pathDates->size() != pathDateTypes->size())
            throw ModelException("CreditDateSort", "path date array and type array are different in size");

        int i;
        vector<DataPair> arr(pathDates->size());
        for (i=0; i<pathDates->size();i++)
        {
            arr[i].date = (*pathDates)[i];
            arr[i].type = (*pathDateTypes)[i];
        }
        // sort the array
        sort(arr.begin(), arr.end());

        // remove duplicate
        vector<DataPair>::iterator iter = arr.begin();
        if (arr.size() > 0)
        {
            for (iter = arr.begin()+1; iter != arr.end(); iter++)
            {
                if (iter->date == (iter-1)->date) // time compared !!!
                {
                    if (iter->type != (iter-1)->type)
                        (iter-1)->type = 2; // this must be both mtm and exposure date !!

                    arr.erase(iter);
                    iter --;
                }
            }
        }
        // return results
        pathDates->resize(arr.size());
        pathDateTypes->resize(arr.size());
        for (i=0; i<pathDates->size();i++)
        {
            (*pathDates)[i] = arr[i].date;
            (*pathDateTypes)[i] = arr[i].type;
        }
    }
};

/** generate n normal deviates */
void randGen(int seed, int n, vector<float>::iterator out)
{
    long idum = -seed;

    for (int i=0; i<n; i++)
        out[i] = gasdev(&idum);
}

/************************* end local helpers **********/

/** create paths dates that require sample generation */ 
void CreditEngine::createPathDates(const DateTime & valueDate, DateTimeArraySP& pathDates, CIntArraySP& pathDateTypes) const
{

    // to do: use minCreditPPY; leave this for spreadsheet ?

    // create new path dates from exposure dates
    pathDates = DateTimeArraySP(new DateTimeArray(creditData->exposureDates));
    pathDateTypes = CIntArraySP(new CIntArray(pathDates->size(), 0));
    // mtm date type = 1
    CIntArray mtmTypes(creditData->mtmDatesInUse.size(), 1);
    // append from mtm dates
    pathDates->insert(pathDates->end(), creditData->mtmDatesInUse.begin(), creditData->mtmDatesInUse.end());
    pathDateTypes->insert(pathDateTypes->end(), mtmTypes.begin(), mtmTypes.end());
    // sort dates
    CreditDateSort(pathDates, pathDateTypes);

    // set all pathDates times to valueDate time
    for( int i=0; i < pathDates->size(); ++i )
        (*pathDates)[i] = DateTime( (*pathDates)[i].getDate(), valueDate.getTime() );
}

/** calculate fx fwds  */
CreditSupportSP CreditEngine::getCreditSupport(CInstrument* inst) const
{
    CreditSupport::Interface* creditInterface = 
        dynamic_cast<CreditSupport::Interface*>(inst);
    if (!creditInterface)
    {
        string msg = "instrument type :"
                    + inst->TYPE->getName()
                    + " does not support Credit Engine.";
        throw ModelException("CreditEngine::getCreditSupport", msg);
    }

    CreditSupportSP creditSupport = creditInterface->createCreditSupport(market);
    creditSupport->setValueDateShift(useThetaFwdRate);

    return creditSupport;
}

/** generate spot paths */
void CreditEngine::genPaths(CreditUndSP und,
                            const DateTimeArray& pathDates,
                            int numOfPaths,
                            const DateTime& valueDate, 
                            int seed,
                            int marketIndex,
                            double corr,
                            DoubleArray& spot1,
                            DoubleArray& spot2,
                            DoubleArray& atmFwd,
							DoubleArray& var_dt) const
{
    int i,j, k, m;

    // get vol
    ATMVolRequest atmReq;
    CVolProcessedBSSP atmVol = CVolProcessedBSSP::dynamicCast(
            CVolProcessedSP(und->getProcessedVol(&atmReq)));

    double spot = und->getSpot();

    ASSERT(2*(numOfPaths/2) == numOfPaths); // even path assumed for antitheptic
    int numOfDates = pathDates.size();
    int numRandom = numOfPaths*numOfDates/2; // the other half comes from anti path

    // get all random numbers in one go
    vector<float> rand(numRandom); 
    if (!Maths::isZero(corr-1)) // just use market index if corr is one, no need to simulate separately
        randGen(seed, numRandom, rand.begin());

    // process correlation if needed
    if(!Maths::isZero(corr))
    {
        if( (int)randMarket.size() <= marketIndex )
            randMarket.resize( marketIndex + 1 );
        
        // sim market index if needed
        if ((int)randMarket[ marketIndex ].size() != numRandom)
        {
            randMarket[ marketIndex ].resize(numRandom);
            randGen(creditData->marketSeeds[ marketIndex ], numRandom, randMarket[ marketIndex ].begin());
        }
        if (Maths::isZero(corr-1))
        {
            for (i=0; i<numRandom; i++)
                rand[i] = randMarket[ marketIndex ][i]*(float)corr;
        }
        else
        {
            float sqr_corr = (float) sqrt(1-corr*corr);
            for (i=0; i<numRandom; i++)
                rand[i] = rand[i]*sqr_corr + randMarket[ marketIndex ][i]*(float)corr;
        }
    }

   // calc drift and var for each step
    DoubleArray v_sqt(numOfDates);
    DoubleArray drift_dt(numOfDates);
    und->fwdValue(pathDates, atmFwd);
    atmVol->CalcVar(valueDate, pathDates, CVolProcessedBS::forward, var_dt);

	drift_dt = atmFwd;
    for (i=numOfDates-1; i>0; i--)
    {
		drift_dt[i] *= exp(-var_dt[i]/2.0)/drift_dt[i-1];
        v_sqt[i] = sqrt(var_dt[i]);
	}
    drift_dt[0] *= exp(-var_dt[0]/2.0)/spot;
    v_sqt[0] = sqrt(var_dt[0]);

    // simulate spot, for path and anti-path
    double deviate;

    for (i=0; i<numOfPaths/2; i++)
    {
        k = i*numOfDates;
		deviate = exp(rand[k]*v_sqt[0]);
        spot1[k] = spot*drift_dt[0]*deviate;
        spot2[k] = spot*drift_dt[0]/deviate;
        for (j=1; j<numOfDates; j++)
        {
            m = k + j;
			deviate = exp(rand[m]*v_sqt[j]);
            spot1[m] = spot1[m-1]*drift_dt[j] * deviate;
            spot2[m] = spot2[m-1]*drift_dt[j] / deviate;
        }
    }
}

/** local helper. check for consistent seed and underlier, correlation */
inline void checkOneUnd(const string& undName, const string& otherUndName, double corr1, double corr2)
{
    if (undName != otherUndName )
        throw ModelException("underlier '" + otherUndName +
                        "' and '" + undName + "' have same seed");

    if(!undName.empty() && fabs(corr1 - corr2) > 1e-6) // good enough and Format is happy
        throw ModelException("underlier '" + undName + "' has two correlation values to market " +
                              Format::toString(corr1) +" and " + Format::toString(corr2));
}

/** returns the maximum of all instrument's last exposure date. */
DateTime CreditEngine::getLastExposureDateFromImnts() const
{
    CInstrumentArray& instArr = *creditData->instArr;
    DateTime lastExpDate = DateTime(0,0);
    for (int i = 0; i < instArr.size(); i++)
    {
        DateTime thisLastExpDate = getCreditSupport(instArr[i].get())->getInstLastExposureDate();
        if (thisLastExpDate > lastExpDate)
        {
            lastExpDate = thisLastExpDate;
        }
    }

    if (lastExpDate.getDate() == 0)
    {
        throw ModelException("CreditEngine::getLastExposureDates", "All supplied instruments have uninitialized "
            "lastExposureDates.");
    }

    return lastExpDate;
}

/** generate paths for one assset and call price() on all instruments with the same underlier */
void CreditEngine::calcOneUndPath(CreditPathValuesOutSP target, 
                                  map<string, int> & usedUnd,
                                  const CIntArray& instIndex,
                                  int seed, int market, double corr,
								  const DoubleArray& positions) const
{
    const string method = "CreditEngine::calcOneUndPath";

    try{

    if (instIndex.size() == 0)
        return; // nothing to do

    int i,j, k;
	
    CreditSupportSP creditSupport = getCreditSupport((*creditData->instArr)[instIndex[0]].get());

    DateTimeArrayConstSP pathDates = target->getPathDates();
    int numOfPaths = (*target->getPathValues()).size();
    int numOfDates = pathDates->size();
    int numRandom = numOfPaths*numOfDates/2; // the other half comes from anti path
    
    DateTime valueDate = (*creditData->instArr)[instIndex[0]]->getValueDate();
    double spot = 0.0;

    DoubleArray atmFwd(numOfDates);
    DoubleArray var_dt(numOfDates);
    DoubleArray spot1(numRandom);
    DoubleArray spot2(numRandom);

    CreditUndSP und = creditSupport->getUnderlier();
    string undName;

    // generate spot paths
    if( und.get() )
    {
        undName = und->getName();

    	// check that underlier has diff seed from market 
	    if( seed == creditData->marketSeeds[ market ] )
		    throw ModelException("underlier " + undName + " has same seed as market");

        spot = und->getSpot();
        genPaths(und,
                 *pathDates,
                 numOfPaths,
                 valueDate, 
                 seed,
                 market,
                 corr,
                 spot1, // output
                 spot2, // output
			     atmFwd, // output
                 var_dt); // output

        // check if this underlier is not processed before (different seed must correspond to different underlier)
        if (usedUnd.find(undName) != usedUnd.end())
            throw ModelException(method, "underlier " + undName + " has two different seeds " +
                                 Format::toString(usedUnd[undName]) +" and "+ Format::toString(seed));
        else
            usedUnd[undName] = seed;
    }

    // prepare local result
    CreditPathValuesOutSP oneResult;
    if (instIndex.size()>1)
    {
        oneResult = CreditPathValuesOutSP::dynamicCast((IObjectSP)target->clone());

        // clean up target content to prepare for add() later
        for (j=0; j<numOfPaths; j++)
        {
            for (k=0; k<numOfDates; k++)
            {
                (*target->getPathValues())[j][k] =0;
            }
        }
    }
    else
        oneResult = target;

    DoubleArray fx(numOfDates, 1.0);
    string creditCcy = creditData->exposureCcy;
    string prevInstCcy = "";
    string instCcy;
    // loop through instruments for all paths
    for (i=0; i<instIndex.size(); i++)
    {
        if (i>0)
        {
            creditSupport = getCreditSupport( (*creditData->instArr)[instIndex[i]].get() );

            // check if underlier is really the same
            CreditUndSP otherUnd = creditSupport->getUnderlier();
            string otherUndName;
            if( otherUnd.get() )
                otherUndName = otherUnd->getName();

            checkOneUnd(undName, otherUndName, corr, creditData->correlationArr[instIndex[i]]);
        }
        // perform any pre calc ( caching)
        creditSupport->preProcess(*oneResult->getPathDates(), atmFwd, var_dt);

        // see if fx needs taken care of
        instCcy = creditSupport->getInstCcyCode();
        if(creditCcy != instCcy)
        {
            if (instCcy == "")
                throw ModelException(method, "empty name for instrument ccy.");
            
            if (prevInstCcy != instCcy)
            {
                calcFXFwd(creditCcy, instCcy, *oneResult->getPathDates(), fx);
                prevInstCcy = instCcy;
            }
        }

		// for debug: compute DoubleArrayArray of spots(path, date)
		if (outputGrids && und.get())
		{
			DoubleArrayArrayArraySP& pathSpots = target->creditDebug->pathSpots;
			int n = pathSpots->size();
			pathSpots->resize(n + 1);
			(*pathSpots)[n].resize(numOfPaths);
			for (int j = 0; j < numOfPaths/2; j++)
			{
				(*pathSpots)[n][j*2].resize(numOfDates);
				(*pathSpots)[n][j*2+1].resize(numOfDates);
				for (int k = 0; k< numOfDates; k++)
				{
					(*pathSpots)[n][j*2][k] = spot1[j*numOfDates + k];
					(*pathSpots)[n][j*2+1][k] = spot2[j*numOfDates + k];
				}
			}
		}

        // calc prices for each path
        for (j=0; j<numOfPaths/2; j++)
        {                
            creditSupport->calcPathValues((*oneResult->getPathValues())[j*2], 
                                          *oneResult->getPathDates(), 
                                          &spot1[j*numOfDates],
                                          spot);
            creditSupport->calcPathValues((*oneResult->getPathValues())[j*2+1], 
                                          *oneResult->getPathDates(), 
                                          &spot2[j*numOfDates],
                                          spot);
        }

		// multiply by position and fx
		double currPosition = positions[instIndex[i]];
		for (j=0; j<numOfPaths; j++)
		{
			for (k=0; k<numOfDates; k++)
			{
				(*oneResult->getPathValues())[j][k] *= currPosition*fx[k];
			}
		}

		// for debug: compute DoubleArrayArray of values(path, date)
		if (outputGrids && und.get())
		{
			smartPtr< DateTimeArrayArray > & pathDates = target->creditDebug->pathDates;
			int nDates = pathDates->size();
            pathDates->resize(nDates + 1);
			(*pathDates)[nDates].resize(numOfDates);
			for (int k = 0; k< numOfDates; k++)
				(*pathDates)[nDates][k] = (*oneResult->getPathDates())[k];

			DoubleArrayArrayArraySP & pathValues = target->creditDebug->pathValues;
			int nValues = pathValues->size();
            pathValues->resize(nValues + 1);
			(*pathValues)[nValues].resize(numOfPaths);
			for (int j = 0; j < numOfPaths; j++)
			{
				(*pathValues)[nValues][j].resize(numOfDates);
				for (int k = 0; k< numOfDates; k++)
					(*pathValues)[nValues][j][k] = (*oneResult->getPathValues())[j][k];
			}
		}

        // net result if needed
        if(instIndex.size() >1)
            target->add(oneResult.get());
    }

    } catch(ModelException& e){
        throw ModelException(e, "CreditEngine::calcOneUndPath");
    }


}

/** takes credit data, group instruments according to underliers, 
    call calcOneUndPath() with an array of indicators that signal the instruments using the same underlier. */
CreditPathValuesOutSP CreditEngine::computePathValues() const
{
    try{
        int i, j;
        int num = creditData->instArr->size();
        map<string, int> usedUnd; // keep track of processed underlier and seed

        // set mtmDatesInUse
        creditData->createMTMDates( marginAcct->marginCallHol, marginAcct->marginCallDelay,
				market->GetReferenceDate());

        // underlier group index
        CIntArray   undIndex;
        CBoolArray  undDone(num);

        CreditPathValuesOutSP results(new CreditPathValuesOut());
        // copy input data 
        results->exposureCcy = creditData->exposureCcy;
        results->gemBandDates = DateTimeArray(creditData->gemBandDates);
        
        CreditPathValuesOutSP netResults;

        // simulation dates and types
        createPathDates(market->GetReferenceDate(), results->getPathDates(), results->getPathDateTypes());
        // set the size of value arrays
        results->getPathValues() = DoubleArrayArraySP(new DoubleArrayArray(creditData->nPaths));
		results->creditDebug->pathDates = smartPtr< DateTimeArrayArray >(new DateTimeArrayArray(0));
		results->creditDebug->pathSpots = DoubleArrayArrayArraySP(new DoubleArrayArrayArray(0));
		results->creditDebug->pathValues = DoubleArrayArrayArraySP(new DoubleArrayArrayArray(0));
        results->lastExposureDate = getLastExposureDateFromImnts();

        for (i=0; i<creditData->nPaths; i++)
        {
            (*results->getPathValues())[i].resize(results->getPathDates()->size());
        }

        // set all index to be processed
        for (i=0; i<num; i++)
            undDone[i] = false;

        // loop underliers
        for (i=0; i<num; i++)
        {
            if (!undDone[i])
            {
                undIndex.clear();
                undIndex.push_back(i);
			    undDone[i] = true;
               // find instruments with same underlying
                for (j=i+1; j<num; j++)
                {
                    if (!undDone[j] && creditData->seedArr[j] == creditData->seedArr[i])
                    {
                        undIndex.push_back(j);
                        undDone[j] = true;
                    }
                }
                // process all instrument of the same underlying.
                calcOneUndPath(results, usedUnd, undIndex,
                               creditData->seedArr[i], creditData->marketIndexes[i],
                               creditData->correlationArr[i], creditData->positions);
        
			    if (!netResults)
			    {
				    if (num >1) // copy
					    netResults = CreditPathValuesOutSP::dynamicCast((IObjectSP)results->clone());
				    else
					    netResults = results;
			    }
			    else // need netting
				    netResults->add(results.get());
		    }
        }

        return netResults;
    }
    catch(ModelException& e){
        throw ModelException(e, "CreditEngine::computePathValues");
    }
}

/** convert values to exposure ccy by simply multiplying fwd fx */
void CreditEngine::calcFXFwd(const string& creditCcy, const string& instCcy,
                             const DateTimeArray& dates, DoubleArray& fx) const
{
	// fetch fxAsset if it's not set or if it's different from what is requested
    if (!fxAsset || 
		fxAsset->getRiskCcyIsoCode() != instCcy || 
		fxAsset->getBaseCcyIsoCode() != creditCcy )
    { // get fx 
        if(!market->hasFXData(instCcy, creditCcy))
            throw ModelException("CreditEngine::calcFXFwd", "cannot find fx data in market for "
                            + instCcy + " vs " + creditData->exposureCcy);

        string name = market->getFXName(instCcy, creditCcy);
        fxAsset = FXAssetSP::dynamicCast((IObjectSP)market->GetData(name, FXAsset::TYPE));
    }
    fxAsset->fwdValue(dates, fx);
}

DRLIB_END_NAMESPACE

