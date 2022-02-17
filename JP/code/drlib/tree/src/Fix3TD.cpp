//----------------------------------------------------------------------------
//
//   Group       : Global QR&D
//
//   Filename    : Fix3TD.cpp
//
//   Description : temporary version of fix3 to work on the claim/bank and
//                 tree starting today without cluttering up the existing Fix3
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/Fix3TD.hpp"
#include "edginc/IrConverter.hpp"
#include "edginc/IRSmile2Q.hpp"
#include "edginc/IREngineTree.hpp"
#include "edginc/IRModelVNFM.hpp"

DRLIB_BEGIN_NAMESPACE

/************************************** Fix3TD *************************************/

/**  collect model initialisation data, set up timeline  */

void Fix3TD::retrieveFactor() {
    try{
        SingleCcyTree::populateTreeCurves();

        IRVolSelector volSelector(getIRVolRaw(IRParams), treeCurves[0],
                                  IRParams->volCalibIndex,
                                  cetSmoothing, IRParams);

        volSelector.getMktVolData(mktVolData);

        mktVolData.ModelChoice = FIX3_TIMEDEP;
        mktVolData.IsNmrModel = false;
    }
    catch (exception& e){
        throw ModelException(e, __FUNCTION__);
    }
}

// expects time dependent VNFM parameters, and standard 2Q smile parameters
void Fix3TD::initTreeData() {
    try {

        treeData.NbFactor = nbFactors;
        mktVolData.NbFactor = nbFactors;

        if (!IRParams->irCalib.isEmpty())
            throw ModelException("Model does not support legacy IRCalib structure "
                "(field name irCalib) to define model and smile parameters - must "
                "use new market structures");

        // retrieve and populate 2Q smile parameters
        if (!IRParams->smileTable.get())
            throw ModelException("Model requires ModelTable structure to be provided "
                "(field name modelTable) for IR definition");

        Populate2QSmile(IRParams->smileSet, *(IRParams->smileTable.get()),
                        treeData, mktVolData);

        // retrieve and populate VNFM IR model parameters
        if (!IRParams->modelTable.get())
            throw ModelException("Model requires ModelTable structure to be provided "
                "(field name modelTable) for IR definition");

        PopulateModelParams(IRParams->modelSet, *(IRParams->modelTable.get()), treeData,
                            mktVolData, true, nbFactors);
        
        if (!engineTable.get())
            throw ModelException("Model requires EngineTable structure to be provided ");

        PopulateEngineParams(engineSet, *(engineTable.get()), treeData, mktVolData);

        // this follows the convention of the allocation of curves in retrieveFactors
        treeData.CvDiff = 0;
        treeData.CvIdx1 = 1;
        treeData.CvIdx2 = 2;
        treeData.CvDisc = 1;


          /* check validity of input */
        if (Fix3_Param_Check_TimeDep(treeData.NbFactor, &mktVolData, &treeData) != SUCCESS)
            throw ModelException("Fix3_Param_Input falied: "+IrConverter::slog.pop());
    }
    catch (exception& e) {
        throw ModelException(e, __FUNCTION__);
    }
}


void Fix3TD::Populate2QSmile(const string& smileKey, IRExoticParamTable& smileTable,
                             FIX3_TREE_DATA& treeData, MKTVOL_DATA& mktVolData) {
    try {

        if (smileKey.empty())
            throw ModelException("supplied smileKey is empty string");

        MarketObjectSP ptrSmile = smileTable.getMarketObject(smileKey);
        MarketObject*  pSmile = ptrSmile.get();
        if (!pSmile)
            throw ModelException("Smile object " + smileKey + 
                " not found in Smile Table " + smileTable.getName() + "!");
        const IRSmile2Q* pIRSmile2Q = dynamic_cast<const IRSmile2Q*>(pSmile);
        if (!pIRSmile2Q)
            throw ModelException("Object " + smileKey + " is not of type IRSmile2Q!");
            
        // Populate smile parameters
        mktVolData.QLeft = pIRSmile2Q->qLeft();
        mktVolData.QRight = pIRSmile2Q->qRight();
        mktVolData.FwdShift = pIRSmile2Q->fwdShift();

        // At input level q=0 means log-normal whereas
        // internally q=0 means normal: we switch here
        mktVolData.QLeft  = 1. - mktVolData.QLeft;
        mktVolData.QRight = 1. - mktVolData.QRight;

        // populate time dependent smile params
        mktVolData.NbSmileDates = 1;
        mktVolData.SmileDate[0] = mktVolData.SwapSt[0];
        mktVolData.QLeftTD[0] = mktVolData.QLeft;
        mktVolData.QRightTD[0] = mktVolData.QRight;
        mktVolData.FwdShiftTD[0] = mktVolData.FwdShift;
    }
    catch (exception& e) {
        throw ModelException(e, __FUNCTION__);
    }
}

void Fix3TD::PopulateEngineParams(const string& engineKey, IRModelConfigTable& engineTable,
                                  FIX3_TREE_DATA& treeData, MKTVOL_DATA& mktVolData) {
   try {

       if (engineKey.empty())
            throw ModelException("supplied engineKey is empty string");

       MarketObjectSP ptrEngine = engineTable.getMarketObject(engineKey);
        MarketObject*  pEngine = ptrEngine.get();
        if (!pEngine)
            throw ModelException("Engine object " + engineKey + 
                " not found in Engine Table " + engineTable.getName() + "!");
        const IREngineTree* pIREngineTree = dynamic_cast<const IREngineTree*>(pEngine);
        if (!pIREngineTree)
            throw ModelException("Engine object " + engineKey + 
                " is not of type IREngineTree!");

        // Populate engine parameters
        mktVolData.CetNbIter = pIREngineTree->nbCET;
        treeData.Ppy = pIREngineTree->nbPPY;
        treeData.NbSigmaMax = pIREngineTree->nbStdDevs;
        treeData.NbStdDevStates = pIREngineTree->nbStateVarStdDevs;
        treeData.NbStates = pIREngineTree->nbStateVars;
    }
    catch (exception& e) {
        throw ModelException(e, __FUNCTION__);
    }
}

// this only populates the structures relevant for TIMEDEP version of fix3 and for
// now should only be called for this model - not other fix3 versions                                  
void Fix3TD::PopulateModelParams(const string& modelKey, IRExoticParamTable& modelTable,
                                 FIX3_TREE_DATA& treeData, MKTVOL_DATA& mktVolData,
                                 bool termStructure, int nbFactors) {
    try {

        if (modelKey.empty())
            throw ModelException("supplied modelKey is empty string");

        MarketObjectSP ptrModel = modelTable.getMarketObject(modelKey);
        MarketObject*  pModel = ptrModel.get();
        if (!pModel)
            throw ModelException("Model object key name " + modelKey + 
                " not found in Model Table " + modelTable.getName());
        const IRModelVNFM* pIRModelVNFM = dynamic_cast<const IRModelVNFM*>(pModel);
        if (!pIRModelVNFM)
            throw ModelException("Object " + modelKey + " is not of type IRModelVNFM!");
        
        // Populate model parameters
        if (pIRModelVNFM->numFactors() != nbFactors)
            throw ModelException("Number of factors defined in the model " +
            Format::toString(nbFactors) + " does not equal number of factors defined in "
            "modelTable (" + modelTable.getName() + "), modelKey " + 
            modelKey + ", which defines " + Format::toString(pIRModelVNFM->numFactors()) + 
            " factors");

        const DoubleMatrix& mMR = pIRModelVNFM->meanReversion();
        const DoubleMatrix& mWeight = pIRModelVNFM->weight();
        const DoubleMatrix& mCorr = pIRModelVNFM->correlation();
        
        int i, t;
        if (termStructure)
        {
            // Get vnfm dates
            const ExpiryArray& dates = pIRModelVNFM->dates();
            if (dates.empty())
                throw ModelException("Internal consistency error - termStructure VNFM model "
                    "expected but list of expiryDates is empty in modelTable (" + 
                    modelTable.getName() + "), modelKey " + modelKey);

            int nbDates = dates.size();
            mktVolData.NbTDInp = nbDates;
            DateTime marketBaseDate = getToday();

            for (t = 0; t < nbDates; t++)
            {
                mktVolData.TDInpDate[t] = dates[t]->toDate(marketBaseDate).toIrDate();
                // ??? Ensure dates are in order when supplied to VNFMParams structure

                for (i = 0; i < nbFactors; ++i)
                {
                    mktVolData.AlphaTD[i][t] = mWeight[i][t];
                    mktVolData.BetaTD[i][t] = mMR[i][t];
                }
            
                for (i = nbFactors + 1; i < 3; ++i)
                {  // to -999 for unused bits
                    mktVolData.AlphaTD[i][t] = -999.;
                    mktVolData.BetaTD[i][t] = -999.;
                }

                for (i = 0; i < mCorr.numCols(); ++i)
                    mktVolData.RhoTD[i][t] = mCorr[i][t];
                for (i = mCorr.numCols() + 1; i < 3; ++i)
                    mktVolData.RhoTD[i][t] = -999.;
            }
        }
        else {
            if (mMR.numRows() != 1)
                throw ModelException("Only expecting 1 meanReversion row for non "
                    "time dependent VNFM model mode - modelKey " + modelKey + 
                    " from modelTable " + modelTable.getName() + " contains " +
                    Format::toString(mMR.numRows()) + " entries");

            if (mWeight.numRows() != 1)
                throw ModelException("Only expecting 1 factor weight row for non "
                    "time dependent VNFM model mode - modelKey " + modelKey + 
                    " from modelTable " + modelTable.getName() + " contains " +
                    Format::toString(mWeight.numRows()) + " entries");

            if (mCorr.numRows() != 1)
                throw ModelException("Only expecting 1 factor correlation row for non "
                    "time dependent VNFM model mode - modelKey " + modelKey + 
                    " from modelTable " + modelTable.getName() + " contains " +
                    Format::toString(mCorr.numRows()) + " entries");

            mktVolData.NbTDInp = 1;

            // if non time dependent mode, ensure no dates are supplied
            const ExpiryArray& dates = pIRModelVNFM->dates();
            if (!dates.empty())
                throw ModelException("To avoid confusion, if running in standard mode "
                    "ie. non TimeDependent, model parameters in modelTable (" + 
                    modelTable.getName() + "), modelKey " + modelKey +
                    " should not supply zero dates - " + Format::toString(dates.size()) +
                    " received");

            for (i = 0; i < nbFactors; ++i) {
                mktVolData.AlphaTD[i][0] = mWeight[i][0];
                mktVolData.BetaTD[i][0] = mMR[i][0];
            }
            for (i = nbFactors + 1; i < 3; ++i) { // to -999 for unused bits
                mktVolData.AlphaTD[i][0] = -999.;
                mktVolData.BetaTD[i][0] = -999.;
            }
            for (i = 0; i < mCorr.numCols(); ++i)
                mktVolData.RhoTD[i][0] = mCorr[i][0];
            for (i = mCorr.numCols() + 1; i < 3; ++i)
                mktVolData.RhoTD[i][0] = -999.;
        }

        mktVolData.Bbq = pIRModelVNFM->backbone();
        mktVolData.Bbq = 1. - mktVolData.Bbq;

    }   
    catch (exception& e) {
        throw ModelException(e, __FUNCTION__);
    }
}


//------------------------------------------------------------------------------

void Fix3TD::load(CClassSP& clazz) {
    clazz->setPublic();
    REGISTER(Fix3TD, clazz);
    SUPERCLASS(SingleCcyTree);
    EMPTY_SHELL_METHOD(defaultConstructor);
}

// interface reflection/class loading mechanism
CClassConstSP const Fix3TD::TYPE = CClass::registerClassLoadMethod(
    "Fix3TD", typeid(Fix3TD), Fix3TD::load);

bool Fix3TDLoad(void) {
    return (Fix3TD::TYPE != 0);
}

DRLIB_END_NAMESPACE
