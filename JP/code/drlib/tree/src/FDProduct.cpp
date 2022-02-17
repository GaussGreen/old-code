//----------------------------------------------------------------------------
//
//   Group       : Derivatives Research
//
//   Filename    : FDProduct.cpp
//
//   Description : FD/tree product base class.
//
//   Author      : Ning Shen
//
//   Date        : February 8, 2005
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"

#include "edginc/Control.hpp"
#include "edginc/Results.hpp"
#include "edginc/FDModel.hpp"
#include "edginc/FDProduct.hpp"

DRLIB_BEGIN_NAMESPACE

const vector< TreeSliceSP > & FDProduct::getSlicesToDEV() const
{
    return slicesToDEV;
}

bool FDProduct::startDEV( const TreeSliceSP & slice )
{
    vector< TreeSliceSP >::iterator iter = find( slicesToDEV.begin(), slicesToDEV.end(), slice );
    if( iter != slicesToDEV.end() )
        return false;
    slicesToDEV.push_back( slice );
    return true;
}

bool FDProduct::stopDEV( const TreeSliceSP & slice )
{
    vector< TreeSliceSP >::iterator iter = find( slicesToDEV.begin(), slicesToDEV.end(), slice );
    if( iter == slicesToDEV.end() )
        return false;
    slicesToDEV.erase( iter );
    return true;
}

void FDProduct::recordSliceToOutputName(Control* ctrl, Results* results, FDModel* model,
    bool isMainPrice, const string &outputName, string ccy, const TreeSlice & price )
{
    if (isMainPrice) 
        results->storePrice(model->getPrice0( price ), ccy);

    if (ctrl && ctrl->isPricing()) {
        double price0 = model->getPrice0( price );

        if (!outputName.empty()) {
            OutputRequest* request = ctrl->requestsOutput(outputName);
            if (request) {
                results->storeRequestResult(request, price0);
            }
            else if (ctrl->requestsOutput(OutputRequest::DBG)) {
                // store it as a DEBUG_PACKET
                results->storeScalarGreek(price0, Results::DEBUG_PACKET, 
                    OutputNameSP(new OutputName(outputName)));
            } // else ignored
        }
    }
}

const TreeSlice & FDProduct::getValue( int step ) const
{
    return getValue( step, model->getDate( step ) );
}
const TreeSlice & FDProduct::getValue( int step, DateTime eventDate ) const
{
    if( eventDate != model->getDate( step ) )
        throw ModelException( __FUNCTION__, "Cannot be valued at a date != step date" );

    return getValue( step );
}

/* default implementation - return product type and default message */
void FDProduct::printInfo(ostream& outputStream) const
{
    string prodType = typeid((*this)).name();
    string name = this->getOutputName();

    outputStream << "Product type: " << prodType << ", name: " << name << 
                 ", no product specific information" << endl;
}

const TreeSlice & FDProduct::getCashFlow( int step ) const {
    throw ModelException("FDProduct::getCashFlow", 
        "getCashFlow() called on a component that does not generate cashflows");
}

string FDProduct::getOutputName() const { 
    return "(" + string(typeid(*this).name()) + ")"; 
}


ModelException FDProduct::makeException(exception &e, const string &functionName) const {
    return ModelException(e, functionName+" "+getOutputName());
}

bool FDProduct::historicalValueAvailable(IProdCreator const &inst, TreeSlice &slice, DateTime eventDate) const {
    try {
        DateTime today = model->getToday();
        if (eventDate <= today) {
            CashflowInfo cfi;
            slice = inst.getValue(eventDate, cfi);
            if (cfi.amountType==CashflowInfo::AmountType::KNOWN)
                return true;
            if (eventDate < today)
                throw ModelException("Could not get value on "+eventDate.toString());
        }
        return false;
    }
    catch (exception& e) {
        throw makeException(e, __FUNCTION__);     
    }
}

void FDProduct::addModelResetDate(DateTime modelResetDate, DateTime lateResetDate) {
    DateTimeArray reset(1, modelResetDate);
    DateTimeArray lateReset(1, lateResetDate);
    addModelResetDates(reset, lateReset);
}

DRLIB_END_NAMESPACE
