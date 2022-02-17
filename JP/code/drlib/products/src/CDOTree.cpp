//----------------------------------------------------------------------------
//
//   Group       : Credit QR
//
//   Filename    : CDOTree.cpp
//
//   Description : SpreadLossTree product for GeneralisedCDO
//
//   Author      : Matthias Arnsdorf
//
//   Date        : August, 2006
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#include "edginc/Maths.hpp"
#include "edginc/Addin.hpp"
#include "edginc/EventResults.hpp"
#include "edginc/CreditTrancheLossConfig.hpp"
#include "edginc/CDOTree.hpp"


DRLIB_BEGIN_NAMESPACE


/////////////////////////////////////////////////////////////////////////////////////
// CDOTree class
////////////////////////////////////////////////////////////////////////////////////

/** constructor */
CDOTree::CDOTree(GeneralisedCDOConstSP inst, FDModel* model) :
FDProduct(model), inst(inst) , protLegProd(), feeLegProd(), calcCondFwds(false)
{

	const static string method = "CDOTree::CDOTree";

	crTree = dynamic_cast<SpreadLossTree*>(model);
	if (!crTree) {
		throw ModelException(method,
			"Only SpreadLossTree models are supported!");
	}
	

    valueDate = inst->getValueDate();

	const string& discYCName = inst->discountYieldCurveName();
	if (discYCName.empty())
		throw ModelException(method, 
		"Discount YieldCurve name of component must be set.");
	

	/** creditLossConfig in GeneralisedCDO */
	CreditTrancheLossConfig* portfolio = dynamic_cast< CreditTrancheLossConfig* > (inst->portfolio.get());
	if(!portfolio)
	{
		throw ModelException("portfolio needs to be of type CreditTrancheLossConfig", method);
	}


	/** low strike as absolute cash amount */
	double lowStrike;
	/** high Strike as cash amount*/
	double highStrike;

	portfolio->getTrancheStrikes(lowStrike, highStrike);

	/** total portfolio notional */
	double notional = portfolio->portfolioNotional(); 

	/** primitive spreadLossTree contingent leg instrument */
	TrancheContingentLegSP contingentLeg(0);
	/** primitive spreadLossTree fee leg instrument */
	TrancheFeeLegSP feeLeg(0);

	if(!!inst->cLeg)
	{
		contingentLeg = TrancheContingentLegSP(new TrancheContingentLeg(
			inst->cLeg,
			lowStrike / notional,
			highStrike / notional,
			inst							/// bad day adjuster
			));
	}

	if(!!inst->fLeg)
	{
		feeLeg = TrancheFeeLegSP(new TrancheFeeLeg(
			inst->fLeg,
			lowStrike / notional,
			highStrike / notional,
			inst->getRecoverNotional(),
			inst								// bad day adjuster
			));
	}

	/// define underlying instruments
	if(!!contingentLeg)
	{
		protLegProd = crTree->createProduct(contingentLeg);
	}
	if(!!feeLeg)
	{
		feeLegProd = crTree->createProduct(feeLeg);
	}

}

void CDOTree::init(Control* ctrl) const
{
	// check if conditional fwds should be calculated
    // i.e if any condional quantity is needed
    OutputRequest* request = ctrl->requestsOutput(OutputRequest::CONDITIONAL_FWD);
    if (request)
    {
        // get date at which conditional fwds should be calculated
        // this is the minimum of the first coupon date or protection start date
        if(!!protLegProd)
        {
            condFwdDate = protLegProd->getStartDate();
        }
        if(!!feeLegProd)
        {
            if(feeLegProd->getStartDate().isLess(condFwdDate))
                condFwdDate = feeLegProd->getStartDate();
        }

        // check if condFwdDate is greater than valueDate in which case we turn on cond fwd calculation
        if(condFwdDate.isGreater(valueDate))
            calcCondFwds = true;
        
    }

    //add fwd date to crit dates if needed
    if(calcCondFwds)
        crTree->addCritDate(condFwdDate);
}

void CDOTree::initProd(void)
{

	string discCurve = inst->discountYieldCurveName();


	getValueSlice = crTree->createSlice();
	getValueSlice->name = inst->outputName + "_getValueSlice";

	instOutputName = inst->outputName;  // to help debugging, store local copy of instrument name


}

/************************************************************************/
/*   Extract value of instrument at given time step                     */
/************************************************************************/
const TreeSlice & CDOTree::getValue(int step) const 
{
	static const string method = "CDOTree::getValue";
	try
	{

		*getValueSlice = 0;
		if(!!protLegProd)
		{
			*getValueSlice += protLegProd->getValue(step);
		}
		if(!!feeLegProd)
		{
			*getValueSlice -= feeLegProd->getValue(step);
		}

        if(!inst->isLong)
        {
            //change sign of MTM
            *getValueSlice *= -1;
        }

		return *getValueSlice;
	}
	catch (exception & e)
	{
		throw ModelException(e, method);
	}

}



/************************************************************************/
/*  update                                                              */
/************************************************************************/
void CDOTree::update(int & step, FDProduct::UpdateType update) 
{
	try 
	{
		// get conditional fwd quantities if needed 
        if(calcCondFwds)
        {
            if(crTree->getDate(step).equals(condFwdDate))
            {
                condFwdStep = step;

                if(!!protLegProd)
                {
                    condContLegSlice = crTree->createSlice();
                    condContLegSlice->name = inst->outputName + "_condContLegSlice";
                    *condContLegSlice = protLegProd->getValue(step);
                }
                if(!!feeLegProd)
                {
                    condFeeLegSlice = crTree->createSlice();
                    condFeeLegSlice->name = inst->outputName + "_condFeeLegSlice";
                    *condFeeLegSlice = feeLegProd->getValue(step);
                }
            }
        }
	}
	catch (exception& e) {
		throw ModelException(e, "CDOTree::update");
	}
}

/************************************************************************/
/* getOutputName                                                        */
/************************************************************************/
string CDOTree::getOutputName() const
{
	return inst->outputName;
}

/************************************************************************/
/* record output                                                        */
/************************************************************************/
void CDOTree::recordOutput(Control* ctrl, YieldCurveConstSP, Results* results) 
{
    string ccy;
    if (inst->discount.get()) ccy = inst->discount->getCcy();

    // CDO value
    recordSliceToOutputName(
        ctrl, results, model, 
        inst->isTopComponent(), 
        inst->outputName, 
        ccy, getValue(0));


    string requestName;
    // TRANCHE_CONTINGENT_LEG_PRICE
    if(!!protLegProd)
    {
        requestName = OutputRequest::TRANCHE_CONTINGENT_LEG_PRICE;
        
        OutputRequest* request = ctrl->requestsOutput(requestName);
        if (request)
        {
            recordSliceToOutputName(
                ctrl, results, model, 
                false, // is main price 
                requestName, 
                ccy, 
                protLegProd->getValue(0));

        }
    }
    
    
    // TRANCHE_FEE_LEG_PRICE
    if(!!feeLegProd)
    {
        requestName = OutputRequest::TRANCHE_FEE_LEG_PRICE;

        OutputRequest*  request = ctrl->requestsOutput(requestName);
        if (request)
        {
            recordSliceToOutputName(
                ctrl, results, model, 
                false, // is main price 
                requestName, 
                ccy, 
                feeLegProd->getValue(0));

        }
    }

    // CONDITIONAL FWDS ----------------------------------------------------------
    if(calcCondFwds)
    {
        requestName = OutputRequest::CONDITIONAL_FWD;
        OutputRequest*  request = ctrl->requestsOutput(requestName);
        // sanity check
        if(!request) throw ModelException("invalid conditional fwd request","CDOTree::recordOutput");

        // CONDITIONAL FWD MTM added at end
        DoubleArraySP condFwdMTM;

        // default levels
        IntArraySP defLevels= crTree->getDefaultLevels();
        results->storeGreek(defLevels, Results::INSTRUMENT_PACKET, 
            OutputNameSP(new OutputName("DEFAULT_LEVELS")));

        // Marginal probabilities
        DoubleArraySP marginals = crTree->getMarginalDefaultProbs(condFwdStep);
        results->storeGreek(marginals, Results::INSTRUMENT_PACKET, 
            OutputNameSP(new OutputName("MARGINAL_DEFAULT_PROBS")));


        // Contingent leg
        if(!!protLegProd)
        {
            // get conditional values
            DoubleArraySP condFwdContLeg = crTree->getValConditionalOnLoss(*condContLegSlice, condFwdStep);
            // record output  as 'greek' in instrument packet 
            results->storeGreek(condFwdContLeg, Results::INSTRUMENT_PACKET, 
                OutputNameSP(new OutputName("CONDITIONAL_CONT_LEG")));
                
            // copy cont leg to MTM
            condFwdMTM = DoubleArraySP(new DoubleArray(*condFwdContLeg));
        }
        // Fee leg
        if(!!feeLegProd)
        {
            // get conditional values
            DoubleArraySP condFwdFeeLeg = crTree->getValConditionalOnLoss(*condFeeLegSlice, condFwdStep);
            // record output    
            results->storeGreek(condFwdFeeLeg, Results::INSTRUMENT_PACKET, 
                OutputNameSP(new OutputName("CONDITIONAL_FEE_LEG")));

            // subtract fee leg from MTM
            QLIB_VERIFY(condFwdFeeLeg->size() == condFwdMTM->size(), 
                "condional cont leg and fee leg have different size");

            int i;
            for(i = 0; i < condFwdMTM->size(); i++)
            {
                (*condFwdMTM)[i] -= (*condFwdFeeLeg)[i];
                if(!inst->isLong)
                {
                    (*condFwdMTM)[i] *= -1;
                }
                
            }

        }

        //store MTM
        results->storeRequestResult(request, condFwdMTM);
        

    } // end if calcCondFwds
    
    inst->recordExtraOutput(ctrl, results);
}




/************************************************************************/
/* printInfo                                                            */
/************************************************************************/
void CDOTree::printInfo(ostream& outputStream) const {
	string prodType = "CDOTree";
	string name = getOutputName();

	outputStream << "Product type: " << prodType << " name: " << name << endl << endl;

	/* TODO WRITE SECTION PRINTING OUT PRODUCT DETAILS */

	// internal slices and slice dimensions info - change slice type, though!
	outputStream << endl;
	outputStream << "Internal product slice details:" << endl;
	outputStream << "instanceName   dimension   devCurveName   descriptiveName" << endl;
	//outputStream << "getValueSlice   " << DYNAMIC_POINTER_CAST<TreeSliceRates>(getValueSlice)->getDim() << "   " <<
	//	getValueSlice->getCurveToDEV() << endl; // "   " << getValueSlice->getSliceName() << endl;

}



DRLIB_END_NAMESPACE
