//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : CalibratorAddin.cpp
//
//   Description : CalibratorAddin
//
//   Author      : regis Guichard
//
//   Date        : 21 May 2002
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/Class.hpp"
#include "edginc/Calibrator.hpp"
#include "edginc/SensMgr.hpp"
#include "edginc/IInstrumentCollection.hpp"
#include "edginc/Format.hpp"
#include "edginc/Results.hpp"
#include "edginc/ObjectIteration.hpp"
#include "edginc/Addin.hpp"
#include "edginc/Maths.hpp"
#include "edginc/NotApplicable.hpp"
#include "edginc/NonPricingModel.hpp"
#include "edginc/IBootstrapper.hpp"
#include <set>
#include ext_hash_map
#include "edginc/ClientRunnable.hpp"
#include "edginc/SpreadSheetMode.hpp"
#include "edginc/Timer.hpp"
#include "edginc/SRMInitialGuess.hpp"

DRLIB_BEGIN_NAMESPACE
class CalibratorAddin2: public CObject,
                        public ClientRunnable{
    OptimizerNDSP               optimizer;
    Calibrator::ObjFuncSP       objFunc;
    Calibrator::InstanceIDArray ids;
    CMarketDataSP               market;
    bool                        doBootstrapping;
	string						initialGuess;
	string						initGuessVolType;

public:
    static CClassConstSP const TYPE;

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(CalibratorAddin2, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(ClientRunnable);
        EMPTY_SHELL_METHOD(defaultCtor);
        FIELD(optimizer, "N-dimensional optimizer");
        FIELD(objFunc, "Objective function");
        FIELD(ids, "Identifies what to calibrate");
        FIELD(market, "Market");
        FIELD_MAKE_OPTIONAL(market);
        FIELD(doBootstrapping, "Indicates if the bootstrapping routine is used");
        FIELD_MAKE_OPTIONAL(doBootstrapping);
		FIELD(initialGuess, "Indicates if initial guess should be used (SRM specific)");
        FIELD_MAKE_OPTIONAL(initialGuess);
		FIELD(initGuessVolType, "Indicates the type of vol used by the guess (SRM specific)");
        FIELD_MAKE_OPTIONAL(initGuessVolType);
                
        // registration for addin function
        Addin::registerClassObjectMethod(
            "CALIBRATOR2",
            Addin::RISK,
            "Calibrates parameters against an objective function",
            CalibratorAddin2::TYPE,
            true,
            Addin::returnHandle,
            (Addin::ObjMethod*)go);        
    }

    
    

//private: // shuts the compiler up
    CalibratorAddin2(const CalibratorAddin2& rhs);
    CalibratorAddin2& operator=(const CalibratorAddin2& rhs);

    static IObjectSP go(CalibratorAddin2* params){
        return params->run();
    }


    /************************************************************************/
    /* Main method to run calibration                                       */
    /************************************************************************/
    virtual IObjectSP run()
    {
        static const string method("CalibratorAddin2::run");

        // if there's a market, give objFunc a chance to get it
        // otherwise, trust (!) that objFunc already has all the market
        // data it needs 
        if (market.get()){
            objFunc->getMarket(market.get());
        }

        
        // adjust optimizer if needed ////////////////////////////////////
        OptimizerNDSP tempOptimizer;
        double MEmyCalcTime;
        // start clock
        Timer timerAdjustment;
        if (!initialGuess.empty())
         {
            // SRM specific special case to adjust optimizer    
            tempOptimizer = getSRMoptimizer();
            MEmyCalcTime = timerAdjustment.calcTime();
        }
        else {
            // further validation irrespective of getMarket having been
            // called or not
            objFunc->validate();
            tempOptimizer = optimizer;
        }
        
        // set up calibrator ////////////////////////////////////////
        Calibrator calibrator(tempOptimizer);
        
    
        // calibrate ////////////////////////////////////////////
        if (!doBootstrapping)
        {
            // Global calibration -----------------------------------------
            // run the calibrator
            return calibrator.run(*objFunc, ids);
        } 
        else if (Calibrator::ObjFunc::IGenericBootstrappable::TYPE->isInstance(objFunc.get())) 
        {
			// "Generic Bootstrapping" --------------------------------------
            //(this does NOT include equity vols and SRM bootstrapping for the moment...)
            return genericBoostrapping(
                tempOptimizer,
                Calibrator::ObjFunc::IGenericBootstrappableSP::dynamicCast(objFunc));

        } 
        else
        {
            /// SRM and equity bootstrapping ------------------
            // Ideally bespoke calibration such as the SRM should fall into the generic framework above.
            // Here we should just have a basic default bootstrapping algorithm
            return equitySrmBootstrapping(calibrator, MEmyCalcTime);
        }
    }

    /************************************************************************/
    /* Method to perform 'generic bootstrapping                             */
    /* This can be completely flexible user defined calibration algorithm   */
    /* To use need to make objective function derive off                    */
    /* 1)Calibrator::objFunc::IGenericBootstrappable                        */
    /* 2) Create getBoostrapper(..) method in objFunc                       */
    /* 3) create instance of IBootrapper which contains the calibration     */
    /*  algorithm                                                           */
    /************************************************************************/
    CResultsSP genericBoostrapping(
        OptimizerNDSP tempOptimizer,
        Calibrator::ObjFunc::IGenericBootstrappableSP bootstrappable)
    {
        static const string method = "CalibratorAddin2::genericBoostrapping()";
        try
        {
            // start clock
            Timer timer;

            // intialise the ids
            IObjectSP adjGroup(objFunc->getAdjustableGroup());
            int i;
            for(i=0;i< ids.size(); i++)
            {
                ids[i]->initialise(adjGroup);
            }
            // Get the "bootstrapper" 
            IBootstrapperSP bootstrapper = bootstrappable->getBootstrapper(ids);


            // ids for the calibrated values
            Calibrator::InstanceIDArray aggregIds;

            // DoubleArray used to aggregate the calibrated values
            DoubleArray aggregCalibValues;

            // number of variables calibrated
            int nbVars;

            // array used to aggregate objective function values
            DoubleArraySP objFuncValues(new DoubleArray());

            // do bootstrapping
            bootstrapper->bootstrap(
                tempOptimizer,                     // (I) calibrator to be used for calibration
                nbVars,                         // (O) number of variables calibrated
                aggregIds,                      // (0)
                aggregCalibValues,                // (O) calibrated values
                *objFuncValues                    // (O) objective function values
                );


            // record time
            double calcTime = timer.calcTime();

            // store results ---------------------------------------
            CResultsSP res(new Results());

            // store the nber of vars
            res->storeScalarGreek(
                nbVars, 
                "Calibrator",
                OutputNameSP(new OutputName("nbVars")));

            // store the calibrated values
            res->storeGreek(
                Calibrator::InstanceID::writeResults(aggregIds, aggregCalibValues), 
                "Calibrator", 
                OutputNameSP(new OutputName("calibratedVars")));

            // store the optimizer statistics (skipped for the moment)
            res->storeGreek(
                IObjectSP(CBool::create(false)), 
                "Calibrator",
                OutputNameSP(new OutputName("leastSquareStats", "status")));

            // export objective function values corresponding
            // to each calibrated step
            res->storeGreek(
                objFuncValues,
                "Calibrator",
                OutputNameSP(new OutputName("calibratedObjFuncValue")));

            // export calc time
            res->storeScalarGreek(
                calcTime, 
                "Calibrator",
                OutputNameSP(new OutputName("CALC_TIME")));   // CALC_TIME will be ignored by testcmp.pl

            
            return res;

        }
        catch(exception & e)
        {
            throw ModelException(e, method);
        }
        
    }

    /************************************************************************/
    /* SRM and equity bootstrapping                                         */
    /* This should really live in a separate file since very asset class    */
    /* specific                                                             */
    /************************************************************************/
    CResultsSP equitySrmBootstrapping(Calibrator & calibrator, double MEmyCalcTime)
    {
        static const string method = "CalibratorAddin2::equitySrmBootstrapping";

        // for SRM
        
        string normalizeWeights;
        string normalizeObjFunc;

        // add to the  output file 
        ExpiryArraySP expiries;

        //export objFunc values for each maturity 
        DoubleArraySP objFuncValues;
        bool exportObjFuncValues;
        bool exportAdjustmentStatus;
        bool isMyAdjustmentSucceeded;
        bool MEmyAdjustmentFailedOtherOption;

        DoubleArraySP objFuncVals;
        CDoubleMatrixSP calibVarsGuess;
        double abstolGuess;
        int count1 = 0;
        int count2 = 0;
        int count3 = 0;
        


        ///////////////////////////////////////////////////////////////////////////////////////////////
        // use initial guess ? 

        // author : Michel Elmaleh
        //          several changes : 
        //          1) the weightMatrix is normalized, i.e the weights sum to one for each maturity.
        //          2) I make some adjustments on the initialGuess before running the full calibration, i.e I 
        //             apply to the objective function a series of one dimensional gradient descent minimizations
        //             in order to get closer to "the minimum".
        //          3) Define a new stopping criterion for the simplex.
        ///////////////////////////////////////////////////////////////////////////////////////////////
        if (!initialGuess.empty()){

            

            // 1) do not normalize the objective function 
            Calibrator::ObjFuncLeastSquareSP objFuncLS = Calibrator::ObjFuncLeastSquareSP::dynamicCast(objFunc);
            objFuncLS->isObjFuncNormalized("do_not_normalize");
            normalizeObjFunc = objFuncLS->getNormalizedWeights();

            // 2) normalize the weights, i.e they sum to 1 
            VanillaGrid::LeastSquareSimpleSP objFuncVGLS;
            try{
                objFuncVGLS = VanillaGrid::LeastSquareSimpleSP::dynamicCast(objFunc);
            } catch (exception&) {
                Calibrator::ObjFuncLeastSquareComboSP objFuncLSCombo = Calibrator::ObjFuncLeastSquareComboSP::dynamicCast(objFunc);
                objFuncVGLS = VanillaGrid::LeastSquareSimpleSP::dynamicCast(objFuncLSCombo->getObjFuncArray()[0]);
            }
            objFuncVGLS->useNormalizedWeights("normalize");
            normalizeWeights = objFuncVGLS->getNormalizedWeights();

            //get number of mc simulations
            IModelSP model = objFuncVGLS->getModel();
            MonteCarlo *mc = dynamic_cast<MonteCarlo *>(model.get()); 
            int nbIter = mc->getNbIters();
            int minNbIterAllowed = 5000;


            // make validations on the objFunc and on the weightMatrix
            static const string method = "invalid normalization for the weights or for the objective function";

            if (!CString::equalsIgnoreCase(normalizeObjFunc,"do_not_normalize") && 
                CString::equalsIgnoreCase(normalizeWeights,"normalize")){
                    throw ModelException(method, "can't normalize both the weights and the objective function !");
                }
                if (CString::equalsIgnoreCase(normalizeObjFunc,"do_not_normalize") && 
                    !CString::equalsIgnoreCase(normalizeWeights,"normalize")){
                        throw ModelException(method, "the weights are not normalized : objFunc should be normalized !");
                    }

                    // create new weightMatrix
                    objFunc->validate();

                    // 3) compute the initialGuess 
                    SRMInitialGuess guess(market,objFunc,ids,initialGuess,initGuessVolType);
                    guess.saveGuess(objFunc);

                    // add the expiries to the output file 
                    expiries = guess.getExpiries();
                    int nbMat = expiries->size();
                    int nbStrikes = guess.getNbStrikes();

                    //export objFuncValues for each maturity
                    objFuncValues = DoubleArraySP(new DoubleArray(nbMat));
                    exportObjFuncValues = false;
                    exportAdjustmentStatus = false;

                    int center;
                    CDoubleMatrix valGuess(nbMat,nbStrikes);
                    DoubleArraySP guessAtmVal = DoubleArraySP(new DoubleArray(nbMat)),
                        guessSmileA1 = DoubleArraySP(new DoubleArray(nbMat)),
                        guessSmileA2 = DoubleArraySP(new DoubleArray(nbMat)),
                        guessSmileA3 = DoubleArraySP(new DoubleArray(nbMat));

                    // make some adjustments on the initialGuess 
                    try{
                        center = guess.getCenterIndex(objFuncVGLS);

                        //make some adjustments on the initialGuess spotVols
                        guessAtmVal  = guess.getAtmVol(),
                            guessSmileA1 = guess.getSmileA1(),
                            guessSmileA2 = guess.getSmileA2(),
                            guessSmileA3 = guess.getSmileA3();


                        DoubleArray shiftSpotVols(nbMat),
                            shiftSmileA1(nbMat),
                            shiftSmileA2(nbMat);
                        DoubleArraySP tweakSpotVols = DoubleArraySP(new DoubleArray(nbMat)),
                            tweakSmileA1 = DoubleArraySP(new DoubleArray(nbMat)),
                            tweakSmileA2 = DoubleArraySP(new DoubleArray(nbMat));
                        CDoubleMatrix sensitivitiesSpotVols(nbMat,nbStrikes),
                            sensitivitiesSmileA1(nbMat,nbStrikes),
                            sensitivitiesSmileA2(nbMat,nbStrikes);
                        objFuncVals = DoubleArraySP(new DoubleArray(nbMat));
                        double _shift = 0.01;
                        double tol = 4.0E-06;
                        int count;

                        //make  adjustment on spotVols
                        try{
                            if(nbIter < minNbIterAllowed){
                                throw ModelException("the number of paths should be at least " 
                                    + Format::toString(minNbIterAllowed) + " ");
                            }
                            objFuncVGLS->calcInitialVals(valGuess,center);
                            guess.nbrOfExpiriesToImprove(valGuess,tol,count,objFuncVals);
                            if(count > 0){
                                guess.tweakOneParameter(_shift,guessAtmVal,shiftSpotVols,tweakSpotVols);
                                guess.adjustAtmVol(tweakSpotVols);
                                guess.saveGuess(objFunc);
                                objFuncVGLS->getSensitivities(sensitivitiesSpotVols,valGuess,shiftSpotVols,center);
                                guess.getAdjustedGuess(sensitivitiesSpotVols,valGuess,guessAtmVal);
                                guess.adjustAtmVol(guessAtmVal);

                                guess.saveGuess(objFunc);
                                objFuncVGLS->calcInitialVals(valGuess,center);
                                guess.nbrOfExpiriesToImprove(valGuess,tol,count,objFuncVals);
                            }
                        }catch(exception& e ){
                            throw ModelException(e, "my adjustments failed!");
                        }


                        //make  adjustment on smileA1
                        try{
                            if(count > 0){
                                guess.tweakOneParameter(_shift,guessSmileA1,shiftSmileA1,tweakSmileA1);
                                guess.adjustSmileA1(tweakSmileA1);
                                guess.saveGuess(objFunc);
                                objFuncVGLS->getSensitivities(sensitivitiesSmileA1,valGuess,shiftSmileA1,center);
                                guess.getAdjustedGuess(sensitivitiesSmileA1,valGuess,guessSmileA1);
                                guess.adjustSmileA1(guessSmileA1);

                                guess.saveGuess(objFunc);
                                objFuncVGLS->calcInitialVals(valGuess,center);
                                guess.nbrOfExpiriesToImprove(valGuess,tol,count,objFuncVals);
                            }
                        }catch(exception& e){
                            throw ModelException(e, "my adjustments failed !");
                        }

                        //make new adjustment on spotVols
                        try{
                            if(count > 0){
                                guess.tweakOneParameter(_shift,guessAtmVal,shiftSpotVols,tweakSpotVols);
                                guess.adjustAtmVol(tweakSpotVols);
                                guess.saveGuess(objFunc);
                                objFuncVGLS->getSensitivities(sensitivitiesSpotVols,valGuess,shiftSpotVols,center);
                                guess.getAdjustedGuess(sensitivitiesSpotVols,valGuess,guessAtmVal);
                                guess.adjustAtmVol(guessAtmVal);

                                guess.saveGuess(objFunc);
                                objFuncVGLS->calcInitialVals(valGuess,center);
                                guess.nbrOfExpiriesToImprove(valGuess,tol,count,objFuncVals);
                            }
                        }catch(exception& e){
                            throw ModelException(e, "my adjustments failed !");
                        }


                        //make new adjustment on smileA1
                        try{
                            if(count > 0){
                                guess.tweakOneParameter(_shift,guessSmileA1,shiftSmileA1,tweakSmileA1);
                                guess.adjustSmileA1(tweakSmileA1);
                                guess.saveGuess(objFunc);
                                objFuncVGLS->getSensitivities(sensitivitiesSmileA1,valGuess,shiftSmileA1,center);
                                guess.getAdjustedGuess(sensitivitiesSmileA1,valGuess,guessSmileA1);
                                guess.adjustSmileA1(guessSmileA1);

                                guess.saveGuess(objFunc);
                                objFuncVGLS->calcInitialVals(valGuess,center);
                                guess.nbrOfExpiriesToImprove(valGuess,tol,count,objFuncVals);
                            }
                        }catch(exception& e){
                            throw ModelException(e, "my adjustments failed !");
                        }


                        //make new adjustment on spotVols 
                        try{
                            if(count > 0){
                                guess.tweakOneParameter(_shift,guessAtmVal,shiftSpotVols,tweakSpotVols);
                                guess.adjustAtmVol(tweakSpotVols);
                                guess.saveGuess(objFunc);
                                objFuncVGLS->getSensitivities(sensitivitiesSpotVols,valGuess,shiftSpotVols,center);
                                guess.getAdjustedGuess(sensitivitiesSpotVols,valGuess,guessAtmVal);
                                guess.adjustAtmVol(guessAtmVal);

                                guess.saveGuess(objFunc);
                                objFuncVGLS->calcInitialVals(valGuess,center);
                                guess.nbrOfExpiriesToImprove(valGuess,tol,count,objFuncVals);
                            }
                        }catch(exception& e){
                            throw ModelException(e, "my adjustments failed !");
                        }



                        //make new adjustment on spotVols 
                        try{
                            if(count > 0){
                                guess.tweakOneParameter(_shift,guessAtmVal,shiftSpotVols,tweakSpotVols);
                                guess.adjustAtmVol(tweakSpotVols);
                                guess.saveGuess(objFunc);
                                objFuncVGLS->getSensitivities(sensitivitiesSpotVols,valGuess,shiftSpotVols,center);
                                guess.getAdjustedGuess(sensitivitiesSpotVols,valGuess,guessAtmVal);
                                guess.adjustAtmVol(guessAtmVal);

                                guess.saveGuess(objFunc);
                                objFuncVGLS->calcInitialVals(valGuess,center);
                                guess.nbrOfExpiriesToImprove(valGuess,tol,count,objFuncVals);
                            }
                        }catch(exception& e){
                            throw ModelException(e, "my adjustments failed !");
                        }


                        //make adjustment on smileA2
                        try{
                            if(count > 0){
                                guess.tweakOneParameter(_shift,guessSmileA2,shiftSmileA2,tweakSmileA2);
                                guess.adjustSmileA2(tweakSmileA2);
                                guess.saveGuess(objFunc);
                                objFuncVGLS->getSensitivities(sensitivitiesSmileA2,valGuess,shiftSmileA2,center);
                                guess.getAdjustedGuess(sensitivitiesSmileA2,valGuess,guessSmileA2);
                                guess.adjustSmileA2(guessSmileA2);

                                guess.saveGuess(objFunc);
                                objFuncVGLS->calcInitialVals(valGuess,center);
                                guess.nbrOfExpiriesToImprove(valGuess,tol,count,objFuncVals);
                            }
                        }catch(exception& e){
                            throw ModelException(e, "my adjustments failed !");
                        }



                        //make new adjustment on spotVols 
                        try{
                            if(count > 0){
                                guess.tweakOneParameter(_shift,guessAtmVal,shiftSpotVols,tweakSpotVols);
                                guess.adjustAtmVol(tweakSpotVols);
                                guess.saveGuess(objFunc);
                                objFuncVGLS->getSensitivities(sensitivitiesSpotVols,valGuess,shiftSpotVols,center);
                                guess.getAdjustedGuess(sensitivitiesSpotVols,valGuess,guessAtmVal);
                                guess.adjustAtmVol(guessAtmVal);
                                guess.saveGuess(objFunc);
                            }
                        }catch(exception& e){
                            throw ModelException(e, "my adjustments failed !");
                        }



                        ////////////////////////////////
                        try{
                            objFuncVGLS->calcInitialVals(valGuess,center);
                            tol = 5.0E-06;
                        }catch(exception& e){
                            throw ModelException(e, "my adjustments failed !");
                        }


                        //make new adjustment on smileA2
                        try{
                            guess.nbrOfExpiriesToImprove(valGuess,tol,count,objFuncVals);
                            if(count > 0){
                                guess.tweakOneParameter(_shift,guessSmileA2,shiftSmileA2,tweakSmileA2);
                                guess.adjustSmileA2(tweakSmileA2);
                                guess.saveGuess(objFunc);
                                objFuncVGLS->getSensitivities(sensitivitiesSmileA2,valGuess,shiftSmileA2,center);
                                guess.getAdjustedGuess(sensitivitiesSmileA2,valGuess,guessSmileA2);
                                guess.adjustSmileA2(guessSmileA2);

                                guess.saveGuess(objFunc);
                                objFuncVGLS->calcInitialVals(valGuess,center);
                                guess.nbrOfExpiriesToImprove(valGuess,tol,count,objFuncVals);

                            }
                        }catch(exception& e){
                            throw ModelException(e, "my adjustments failed !");
                        }



                        //make new adjustment on smileA2
                        try{
                            if(count > 0){
                                guess.tweakOneParameter(_shift,guessSmileA2,shiftSmileA2,tweakSmileA2);
                                guess.adjustSmileA2(tweakSmileA2);
                                guess.saveGuess(objFunc);
                                objFuncVGLS->getSensitivities(sensitivitiesSmileA2,valGuess,shiftSmileA2,center);
                                guess.getAdjustedGuess(sensitivitiesSmileA2,valGuess,guessSmileA2);
                                guess.adjustSmileA2(guessSmileA2);

                                guess.saveGuess(objFunc);
                                objFuncVGLS->calcInitialVals(valGuess,center);
                                guess.nbrOfExpiriesToImprove(valGuess,tol,count,objFuncVals);
                            }
                        }catch(exception& e){
                            throw ModelException(e, "my adjustments failed !");
                        }



                        //make new adjustment on spotVols
                        try{
                            if(count > 0){
                                guess.tweakOneParameter(_shift,guessAtmVal,shiftSpotVols,tweakSpotVols);
                                guess.adjustAtmVol(tweakSpotVols);
                                guess.saveGuess(objFunc);
                                objFuncVGLS->getSensitivities(sensitivitiesSpotVols,valGuess,shiftSpotVols,center);
                                guess.getAdjustedGuess(sensitivitiesSpotVols,valGuess,guessAtmVal);
                                guess.adjustAtmVol(guessAtmVal);

                            }
                        }catch(exception& e){
                            throw ModelException(e, "my adjustments failed !");
                        }



                        // make sure that the adjusted parameters are in the open range defined for the calibration
                        calibVarsGuess = CDoubleMatrixSP(new CDoubleMatrix(nbMat,4));
                        for(int t=0;t<nbMat;++t){
                            (*calibVarsGuess)[t][0] = (*guessSmileA1)[t];
                            (*calibVarsGuess)[t][1] = (*guessSmileA2)[t];
                            (*calibVarsGuess)[t][2] = (*guessSmileA3)[t];
                            (*calibVarsGuess)[t][3] = (*guessAtmVal)[t];
                        }

                        guess.validate(calibVarsGuess);


                        //save the guess
                        guess.saveGuess(objFunc);

                        //compute the objFunc for each maturity
                        objFuncVals = DoubleArraySP(new DoubleArray(nbMat));
                        try{
                            objFuncVGLS->calcInitialVals(valGuess,center);
                            guess.getObjFuncArray(valGuess,objFuncVals,count1,count2,count3);
                        }catch(exception& e){
                            throw ModelException(e, "my adjustments failed !");
                        }

                        isMyAdjustmentSucceeded = true;

                    }catch(exception&){
                        //// save the initialGuess : no adjustment on the initialGuess spotVols
                        isMyAdjustmentSucceeded = false;
                        objFuncVals = DoubleArraySP(new DoubleArray(nbMat));
                        SRMInitialGuess guess(market,objFunc,ids,initialGuess,initGuessVolType);
                        guess.saveGuess(objFunc);

                        try{
                            if(nbIter < minNbIterAllowed){
                                throw ModelException("the number of paths should be at least " 
                                    + Format::toString(minNbIterAllowed) + " ");
                            }
                            MEmyAdjustmentFailedOtherOption = true;
                            objFuncVals = DoubleArraySP(new DoubleArray(nbMat));
                            objFuncVGLS->calcInitialVals(valGuess,center);
                            guess.getObjFuncArray(valGuess,objFuncVals,count1,count2,count3);

                            calibVarsGuess = CDoubleMatrixSP(new CDoubleMatrix(nbMat,4));
                            for(int t=0;t<nbMat;++t){
                                (*calibVarsGuess)[t][0] = (*guess.getSmileA1())[t];
                                (*calibVarsGuess)[t][1] = (*guess.getSmileA2())[t];
                                (*calibVarsGuess)[t][2] = (*guess.getSmileA3())[t];
                                (*calibVarsGuess)[t][3] = (*guess.getAtmVol())[t];
                            }

                        }catch(exception&){
                            MEmyAdjustmentFailedOtherOption = false;
                        }

                    }
        }

        // start clock
        Timer timer;

        // validate the ids 
        InstanceIDValidate();

        // validate obj func type
        if (!Calibrator::ObjFunc::IBootstrappable::TYPE->isInstance(objFunc.get())){
            throw ModelException(method, "The objective function has to be bootstrappable");
        }

        // no InstanceIDs to calibrate
        int nbExpiries;
        if (!ids.size() 
            || !(nbExpiries = ids[0]->numVariables())){
                // run the calibrator (just calculate the value of the objective function)
                return calibrator.run(*objFunc, ids);
            }
        else {
            // create a CResultsArray object to store the results
            // after each step in the bootstrapping routine
            CResultsArray resArray(nbExpiries);

            Calibrator::InstanceID::IBootstrappableArray bootstrapIds(ids.size());
            int iIds;
            for (iIds = 0 ; iIds < ids.size() ; iIds++){
                // dynamic ids to idsBootstrap (bootstrappable)
                bootstrapIds[iIds] = Calibrator::InstanceID::IBootstrappableSP::dynamicCast(ids[iIds]);
                if(bootstrapIds[iIds].get() == 0){
                    throw ModelException(method, "The array of variables have to be bootstrappable");
                }
            }

            // make some additional checks on the objective function
            Calibrator::ObjFunc::IBootstrappable& objFuncBootstrap =
                dynamic_cast<Calibrator::ObjFunc::IBootstrappable&>(*objFunc);
            objFuncBootstrap.validate(bootstrapIds);

            // storage of the calibrated values and optimizer statistics

            // DoubleArray used to aggregate the calibrated values
            Calibrator::InstanceIDArray aggregIds;
            DoubleArray aggregCalibValues;
            // double used to aggregate the variance
            double variance = 0.0;
            // DoubleArray used to aggregate the upper limits and lower limits
            DoubleArray upperLimits; // upper limits
            DoubleArray lowerLimits; // lower limits
            DoubleArray standardErrors; // standard errors
            // DoubleMatrix used to aggregate the correlation and covariance matrices
            int dimMat = nbExpiries * ids.size(); // dimension of the correlation and covariance matrices
            CDoubleMatrix correlationMatrix(dimMat, dimMat); // correlation matrix (initialized to 0.0)
            CDoubleMatrix covarianceMatrix(dimMat, dimMat); // covariance matrix (initialized to 0.0)

            bool leastSquareStatus = false; // are statistics requested?
            int nbVars = 0;

            if(!initialGuess.empty()){
                if(count1 > 0){ abstolGuess = 4.5E-06;}
                if(count2 > 0){ abstolGuess = 5.0E-06;}
                if(count3 == 1){ abstolGuess = 5.5E-06;}
                else if(count1 == 0 && count2 == 0 && count3 != 1){ abstolGuess = 4.0E-06;}
            }

            int iExpiry;
            double tempVal;
            for (iExpiry = 0 ; iExpiry < nbExpiries ; iExpiry++){
                try{
                    // if bootstrappable, the objective function might need some updates
                    // before each run of the calibrator
                    objFuncBootstrap.update(iExpiry);

                    // if bootstrappable, get the set of scalar InstanceID to calibrate
                    Calibrator::InstanceIDArray theseIds(ids.size());
                    for (iIds = 0 ; iIds < ids.size() ; iIds++){
                        theseIds[iIds] = bootstrapIds[iIds]->getInstanceID(iExpiry);
                        aggregIds.push_back(theseIds[iIds]);
                        nbVars += theseIds[iIds]->numVariables();
                    }


                    ///////////////////////////////////////////////////////////////////////////////////////////////
                    // for SRM
                    // author : Michel Elmaleh
                    ///////////////////////////////////////////////////////////////////////////////////////////////
                    if(!initialGuess.empty()){
                        double objFuncVal = (*objFuncVals)[iExpiry];
                        DoubleArray calibGuess(ids.size());
                        if(MEmyAdjustmentFailedOtherOption){
                            for(int i=0; i < ids.size(); ++i){
                                calibGuess[i] = (*calibVarsGuess)[iExpiry][i];
                            }
                            if(objFuncVal < abstolGuess){
                                resArray[iExpiry] = calibrator.run(*objFunc,theseIds,calibGuess,objFuncVal);
                                tempVal = sqrt(objFuncVal) * 10000.;
                            }
                            else{
                                // run the calibrator
                                resArray[iExpiry] = calibrator.run(*objFunc, theseIds);
                                tempVal = 10000. *  sqrt(resArray[iExpiry]->retrieveScalarGreek("Calibrator",
                                    OutputNameSP(new OutputName("calibratedObjFuncValue")))); 
                            }
                        }
                        else{
                            // run the calibrator
                            resArray[iExpiry] = calibrator.run(*objFunc, theseIds);
                            tempVal = 10000. *  sqrt(resArray[iExpiry]->retrieveScalarGreek("Calibrator",
                                OutputNameSP(new OutputName("calibratedObjFuncValue")))); 
                        }
                        (*objFuncValues)[iExpiry] = floor(tempVal) + floor((tempVal - floor(tempVal)) * 10.)/10.;
                    }
                    else{
                        // run the calibrator
                        resArray[iExpiry] = calibrator.run(*objFunc, theseIds);
                    }
                    ///////////////////////////////////////////////////////////////////////////////////////////////
                    ///////////////////////////////////////////////////////////////////////////////////////////////


                    // extract calibrated values (XXX this is a complete shit show !!!)
                    IObjectConstSP calibVars0 = 
                        resArray[iExpiry]->retrieveGreek("Calibrator",
                        OutputNameSP(new OutputName("calibratedVars")));
                    IObjectSP calibVars1 = IObjectSP::constCast(calibVars0);
                    if (!calibVars1 || !ObjectArray::TYPE->isInstance(calibVars1)){
                        throw ModelException(method, "format error for calibrated vars");
                    }
                    ObjectArraySP calibVars2 = ObjectArraySP::dynamicCast(calibVars1);
                    if (calibVars2->size() != Calibrator::InstanceID::NB_RES 
                        || !ObjectArray::TYPE->isInstance((*calibVars2)[Calibrator::InstanceID::VALUE])){
                            throw ModelException(method, "format error for calibrated values");
                        }
                        ObjectArraySP vals0 = ObjectArraySP::dynamicCast((*calibVars2)[Calibrator::InstanceID::VALUE]);
                        int nbVals = vals0->size();
                        DoubleArray vals1;
                        int iVal;
                        for (iVal = 0; iVal < nbVals; ++iVal){
                            if (!(*vals0)[iVal]){
                                throw ModelException(method, "format error for calibrated vars");
                            }
                            else if (CDouble::TYPE->isInstance((*vals0)[iVal])){
                                CDoubleSP vals2 = CDoubleSP::dynamicCast((*vals0)[iVal]);
                                vals1.push_back(vals2->doubleValue());
                                // aggregate the ids and the calibrated values
                                aggregCalibValues.push_back(vals2->doubleValue());
                            }
                            else{
                                // XXX should never get there (for now!)
                                throw ModelException(method, "format error for calibrated vars");
                            }
                        }
                        // apply calibrated values to obj func
                        Calibrator::InstanceID::applyAdjustment(theseIds,
                            objFunc->getAdjustableGroup(),
                            vals1);

                        // aggregate the optimizer statistics
                        // calibration status for the current run of the calibrator
                        IObjectConstSP calibrationStatus =
                            resArray[iExpiry]->retrieveGreek("Calibrator",
                            OutputNameSP(new OutputName("leastSquareStats", "status")));
                        if (!calibrationStatus || !CBool::TYPE->isInstance(calibrationStatus))
                            throw ModelException(method, "format error for status");
                        CBoolConstSP calibStatus = CBoolConstSP::dynamicCast(calibrationStatus);
                        bool success = calibStatus->boolValue();

                        if (success){
                            // succeeded at least once
                            leastSquareStatus = true;
                            // variance
                            IObjectConstSP calibrationVariance = 
                                resArray[iExpiry]->retrieveGreek("Calibrator",
                                OutputNameSP(new OutputName("leastSquareStats", "variance")));
                            if (!calibrationVariance || !CDouble::TYPE->isInstance(calibrationVariance))
                                throw ModelException(method, "format error for variance");
                            CDoubleConstSP calibVariance = CDoubleConstSP::dynamicCast(calibrationVariance);
                            variance += calibVariance->doubleValue();
                            // upper limits
                            IObjectConstSP upperLimits1 = 
                                resArray[iExpiry]->retrieveGreek("Calibrator",
                                OutputNameSP(new OutputName("leastSquareStats", "upperLimits")));
                            if (!upperLimits1 || !DoubleArray::TYPE->isInstance(upperLimits1))
                                throw ModelException(method, "format error for upper limits");                          
                            DoubleArrayConstSP upperLimits2 = DoubleArrayConstSP::dynamicCast(upperLimits1);
                            // lower limits
                            IObjectConstSP lowerLimits1 = 
                                resArray[iExpiry]->retrieveGreek("Calibrator",
                                OutputNameSP(new OutputName("leastSquareStats", "lowerLimits")));
                            if (!lowerLimits1 || !DoubleArray::TYPE->isInstance(lowerLimits1))
                                throw ModelException(method, "format error for lower limits");                      
                            DoubleArrayConstSP lowerLimits2 = DoubleArrayConstSP::dynamicCast(lowerLimits1);
                            // standard errors
                            IObjectConstSP stdErrs1 = 
                                resArray[iExpiry]->retrieveGreek("Calibrator",
                                OutputNameSP(new OutputName("leastSquareStats", "standardErrors")));
                            if (!stdErrs1 || !DoubleArray::TYPE->isInstance(stdErrs1))
                                throw ModelException(method, "format error for standard errors");                           
                            DoubleArrayConstSP stdErrs2 = DoubleArrayConstSP::dynamicCast(stdErrs1);                    
                            // correlation matrix
                            IObjectConstSP correls1 = 
                                resArray[iExpiry]->retrieveGreek("Calibrator",
                                OutputNameSP(new OutputName("leastSquareStats", "correlationMatrix")));
                            if (!correls1 || !CDoubleMatrix::TYPE->isInstance(correls1))
                                throw ModelException(method, "format error for correlation matrix");                    
                            CDoubleMatrixConstSP correls2 = CDoubleMatrixConstSP::dynamicCast(correls1);
                            // covariance matrix
                            IObjectConstSP covars1 = 
                                resArray[iExpiry]->retrieveGreek("Calibrator",
                                OutputNameSP(new OutputName("leastSquareStats", "covarianceMatrix")));
                            if (!covars1 || !CDoubleMatrix::TYPE->isInstance(covars1))
                                throw ModelException(method, "format error for covariance matrix");                         
                            CDoubleMatrixConstSP covars2 = CDoubleMatrixConstSP::dynamicCast(covars1);
                            // put in contiguous arrays / matrices
                            for (iVal = 0 ; iVal < nbVals ; ++iVal){   
                                upperLimits.push_back((*upperLimits2)[iVal]); // upper limits
                                lowerLimits.push_back((*lowerLimits2)[iVal]); // lower limits
                                standardErrors.push_back((*stdErrs2)[iVal]); // standard errors
                                int jVal;
                                for (jVal = 0 ; jVal < nbVals ; ++jVal){
                                    correlationMatrix[iExpiry * ids.size() + iVal][iExpiry * ids.size() + jVal] = 
                                        (*correls2)[iVal][jVal]; // correlations
                                    covarianceMatrix[iExpiry * ids.size() + iVal][iExpiry * ids.size() + jVal] = 
                                        (*covars2)[iVal][jVal]; // covariances
                                }
                            }
                        }
                        else{
                            variance += 1.0e+30;
                            for (iVal = 0 ; iVal < nbVals ; ++iVal){
                                upperLimits.push_back(1.0e+30);
                                lowerLimits.push_back(-1.0e+30);
                                standardErrors.push_back(1.0e+30);
                                int jVal;
                                for (jVal = 0 ; jVal < nbVals ; ++jVal){
                                    correlationMatrix[iExpiry * ids.size() + iVal][iExpiry * ids.size() + jVal] = 1.0e+30;
                                    covarianceMatrix[iExpiry * ids.size() + iVal][iExpiry * ids.size() + jVal] = 1.0e+30;
                                }
                            }
                        }
                        objFuncBootstrap.reset();
                }
                catch(exception& e){
                    string message = "CalibratorAddin2::run: Bootstrapping failed at "
                        + Format::toString(iExpiry + 1)
                        + "-th expiry";
                    throw ModelException::addTextToException(e, message);
                }
            }

            CResultsSP res(new Results());

            // store the nber of vars
            res->storeScalarGreek(nbVars, 
                "Calibrator",
                OutputNameSP(new OutputName("nbVars")));

            // store the calibrated values
            res->storeGreek(Calibrator::InstanceID::writeResults(aggregIds, aggregCalibValues), 
                "Calibrator", 
                OutputNameSP(new OutputName("calibratedVars")));

            // store the optimizer statistics
            if (leastSquareStatus){
                // variance
                res->storeScalarGreek(variance, 
                    "Calibrator",
                    OutputNameSP(new OutputName("leastSquareStats", "variance")));
                // upper limits
                res->storeGreek(IObjectSP(upperLimits.clone()),
                    "Calibrator",
                    OutputNameSP(new OutputName("leastSquareStats", "upperLimits")));
                // lower limits
                res->storeGreek(IObjectSP(lowerLimits.clone()),
                    "Calibrator",
                    OutputNameSP(new OutputName("leastSquareStats", "lowerLimits")));
                // standard errors
                res->storeGreek(IObjectSP(standardErrors.clone()),
                    "Calibrator",
                    OutputNameSP(new OutputName("leastSquareStats", "standardErrors")));
                // covariance matrix
                res->storeGreek(IObjectSP(covarianceMatrix.clone()),
                    "Calibrator",
                    OutputNameSP(new OutputName("leastSquareStats", "covarianceMatrix")));
                // correlation matrix
                res->storeGreek(IObjectSP(correlationMatrix.clone()),
                    "Calibrator",
                    OutputNameSP(new OutputName("leastSquareStats", "correlationMatrix")));
            }
            res->storeGreek(IObjectSP(CBool::create(leastSquareStatus)), 
                "Calibrator",
                OutputNameSP(new OutputName("leastSquareStats", "status")));

            // calculate obj func and export obj func value 
            if (!initialGuess.empty()){
                int nbMat = expiries->size();
                double f = objFunc->calcValue();
                if(CString::equalsIgnoreCase(normalizeObjFunc,"do_not_normalize")  && 
                    CString::equalsIgnoreCase(normalizeWeights,"normalize")){
                        f /= nbMat;
                    }
                    res->storeScalarGreek(f,
                        "Calibrator",
                        OutputNameSP(new OutputName("calibratedObjFuncValue")));
            }
            else{
                double f = objFunc->calcValue();
                res->storeScalarGreek(f,
                    "Calibrator",
                    OutputNameSP(new OutputName("calibratedObjFuncValue")));
            }

            //export objFunc values for each maturity if needed
            if (!initialGuess.empty()){
                if(exportObjFuncValues){
                    res->storeGreek(objFuncValues,
                        "Calibrator",
                        OutputNameSP(new OutputName("basis points error on volatilities")));
                }
            }

            // export calc time
            double calcTime = timer.calcTime();
            if(!initialGuess.empty()){
                calcTime += MEmyCalcTime;
            }
            res->storeScalarGreek(calcTime, 
                "Calibrator",
                OutputNameSP(new OutputName("CALC_TIME")));   // CALC_TIME will be ignored by testcmp.pl

            if (!initialGuess.empty()){
                // export maturities 
                res->storeGreek(
                    expiries, 
                    "Calibrator",
                    OutputNameSP(new OutputName("ExpiryArray")));
            }

            if(!initialGuess.empty()){
                //export adjustments status
                if(exportAdjustmentStatus){
                    res->storeGreek(
                        IObjectSP(CBool::create(isMyAdjustmentSucceeded)),
                        "Calibrator",
                        OutputNameSP(new OutputName("Adjustments", "Status")));
                }
            }

            return res;
        }

    }

    /************************************************************************/
    /* Method to define bespoke optimizer used in SRM                       */
    /* Ideally this would be part of a generic bootstrapper                 */
    /************************************************************************/
    OptimizerNDSP getSRMoptimizer()
    {
        static const string method = "CalibratorAddin2::getSRMoptimizer";

        // 4) create new simplex for SRM and define the optimizer
        try{
            SimplexSP newOptimizer = SimplexSP::dynamicCast(optimizer);
            // absolute tolerance : stopping criterion for the simplex 
            double absTol = 2.0E-06;
            static const string method = "set the absolute tolerance used by the composite simplex";
            if(absTol > 4.0E-06){
                throw ModelException(method, "the absolute tolerance is too high : decrease it !");
            }
            double myLengthScale = 0.1;

            newOptimizer->useMyMapping("relative");
            if(newOptimizer->getLengthScale() > 0.){
                newOptimizer->setLengthScale(myLengthScale);
            }
            newOptimizer->useSRMSimplex("use_amoeba2");
            newOptimizer->setStoppingCriterionType("composite");
            newOptimizer->setAbsTol(absTol);
            return newOptimizer;

        }catch (exception&) {
            return optimizer;
        }

    }

    /************************************************************************/
    /*  Makes some checks on the InstanceIDs before running the calibration */
    /*  when bootstrapping routine is used                                  */
    /************************************************************************/
    void InstanceIDValidate() const{
        static const string method("CalibratorAddin2::InstanceIDValidate");
        try{
            int nberIds = ids.size(); // number of InstanceIds 

            // Must call initialise on instanceIDs before calling any other
            // methods. So first:
            /** get hold of market data names to calibrate */
            IObjectSP adjGroup(objFunc->getAdjustableGroup());

            // checks that the InstanceIds are of type vector
            // and initialise them
            int iIds;
            for (iIds = 0 ; iIds < nberIds ; iIds++){
                if (!Calibrator::InstanceID::IBootstrappable::TYPE->isInstance(ids[iIds].get())){
                    throw ModelException(method, 
                        "the InstanceIds have to be of type Calibrator::InstanceID::IBootstrappable "
                        "if the bootstrapping routine is used.");
                }
                // initialise instanceIDs before we invoke any methods on them
                ids[iIds]->initialise(adjGroup);
            }

            // checks that the InstanceIds have the same size
            if (ids.size() > 0) {
                int sizeIds = ids[0]->numVariables();

                for (iIds = 1 ; iIds < nberIds ; iIds++)
                {
                    if (ids[iIds]->numVariables() != sizeIds)
                        throw ModelException(method, 
                        "The "
                        + Format::toString(iIds+1)
                        + "-th InstanceId is not of the same size ("
                        + Format::toString(ids[iIds]->numVariables())
                        + ") as the previous ones ("
                        + Format::toString(sizeIds)
                        + ")");
                }
            }
        }
        catch(exception& e){
            throw ModelException(e, method);
        }
    }


    CalibratorAddin2(): 
        CObject(TYPE),
        doBootstrapping(false) {}
                                                
    static IObject* defaultCtor(){
        return new CalibratorAddin2();
    }

};



CClassConstSP const CalibratorAddin2::TYPE = CClass::registerClassLoadMethod(
    "CalibratorAddin2", typeid(CalibratorAddin2), CalibratorAddin2::load);

bool CalibratorAddinLoad() {
	return CalibratorAddin2::TYPE != NULL;
}

/* Addin class */
DRLIB_END_NAMESPACE
