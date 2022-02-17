//----------------------------------------------------------------------------
//
//   Group       : Credit Hybrids QR
//
//   Filename    : Instrument2TestMatrixesOfSPsRRs.cpp
//
//   Description : Instrument we use to test matrix Of SPs 
//
//   Date        : May 2006
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/Instrument2TestMatrixesOfSPsRRs.hpp"
#include "edginc/LessEqualGreaterEps.hpp"
#include "edginc/CIDClosedFormProduct2TestMatrixesOfSPsRRs.hpp"
#include "edginc/TemplateFunctions4Containers.hpp"

DRLIB_BEGIN_NAMESPACE

///////////////////////////////////////////////////////////////////////////////
// const string Instrument2TestMatrixesOfSPsRRs::COND_MC
//                                                = "Fast Monte-Carlo";
//const string Instrument2TestMatrixesOfSPsRRs::FULL_BLOWN_MC
//                                                = "Full-Blown Monte-Carlo";
///////////////////////////////////////////////////////////////////////////////
//MCStateVarSPsRRsSP Instrument2TestMatrixesOfSPsRRs::getStateVar() const
//{
//    static const string routine = "Instrument2TestMatrixOfSPs::getStateVar";
//    try
//    {
//        if (m_typeOfComputation == COND_MC)
//        {
//            return MCStateVarSPsRRsSP(new MCSVBuSPsRRsCondCID::StateVar(
//                                                                     m_dates, 
//                                                                     m_names));
//        }
//        if (m_typeOfComputation == FULL_BLOWN_MC)
//        {
//            return MCStateVarSPsRRsSP(new MCSVBuSPsRRsFullCID::StateVar(
//                                                                     m_dates, 
//                                                                     m_names));
//        }
//        // there are only 2 CID Monte-Carlos: full-blown MC and conditional one
//        throw ModelException("Unrecognized Monte-Carlo type", routine);
//    }
//    catch (exception& e)
//    {
//        throw ModelException(e, routine);
//    }
//}
///////////////////////////////////////////////////////////////////////////////
Instrument2TestMatrixesOfSPsRRs::Instrument2TestMatrixesOfSPsRRs()
    :   CInstrument(Instrument2TestMatrixesOfSPsRRs::TYPE)
    ,   names(0)
    ,   dates(0)
    ,   valueDate(0,0)
{}
///////////////////////////////////////////////////////////////////////////////
CClassConstSP const Instrument2TestMatrixesOfSPsRRs::TYPE 
        = CClass::registerClassLoadMethod(
                                "Instrument2TestMatrixesOfSPsRRs", 
                                typeid(Instrument2TestMatrixesOfSPsRRs), load);
///////////////////////////////////////////////////////////////////////////////
void Instrument2TestMatrixesOfSPsRRs::load(CClassSP& clazz)
{
    clazz->setPublic();
    REGISTER(Instrument2TestMatrixesOfSPsRRs, clazz);
    SUPERCLASS(CInstrument);
    EMPTY_SHELL_METHOD(defaultConstructor);
    IMPLEMENTS(CIDClosedFormModel::IIntoProduct);
    FIELD(names,     "Names");
    FIELD(dates,     "Dates");
    FIELD(valueDate, "Value date");
}
///////////////////////////////////////////////////////////////////////////////
IObject* Instrument2TestMatrixesOfSPsRRs::defaultConstructor()
{
	static const string routine 
                       = "Instrument2TestMatrixesOfSPsRRs::defaultConstructor";
	try
    {
        return new Instrument2TestMatrixesOfSPsRRs();
	}
	catch (exception& e)
	{
		throw ModelException(e, routine);
	}
}
///////////////////////////////////////////////////////////////////////////////
StringArrayConstSP   Instrument2TestMatrixesOfSPsRRs::getNames()    const
{
    static const string routine = "Instrument2TestMatrixesOfSPsRRs::getNames";
    try
    {
        return names;
    }
    catch (exception& e)
    {
        throw ModelException(e, routine);
    }
}
///////////////////////////////////////////////////////////////////////////////
DateTimeArrayConstSP Instrument2TestMatrixesOfSPsRRs::getDates()    const
{
    static const string routine = "Instrument2TestMatrixesOfSPsRRs::getDates";
    try
    {
        return dates;
    }
    catch (exception& e)
    {
        throw ModelException(e, routine);
    }
}
///////////////////////////////////////////////////////////////////////////////
void Instrument2TestMatrixesOfSPsRRs::validatePop2Object()
{
    static const string method 
                       = "Instrument2TestMatrixesOfSPsRRs::validatePop2Object";
    try
    {   
        // NULLiness and emptiness checkings ----------------------------------
        QLIB_VERIFY(!((!names) || (!dates)), "Inputs are NULL");
        QLIB_VERIFY(!((names->empty()) || (dates->empty())),
                                                 "Names or dates are empty");
        // CHECK NAMES ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        // nothing to check
        // CHECK DATES ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        // dates are after today's date ---------------------------------------
        QLIB_VERIFY(allElemsOf1AreGreaterThen2(*dates, valueDate),
                    "Attempt to compute SPs for the dates in the past"); 
    }
    catch(exception& e)
    {
        throw ModelException(e, method);
    }
}
///////////////////////////////////////////////////////////////////////////////
CIDClosedFormModel::IProductSP Instrument2TestMatrixesOfSPsRRs::createProduct(
                                        const CIDClosedFormModel* gModel) const
{                                                
	static const string routine 
                            = "Instrument2TestMatrixesOfSPsRRs::createProduct";
	try
    {
        return CIDClosedFormModel::IProductSP(
                 new CIDClosedFormProduct2TestMatrixesOfSPsRRs(this));
	}
	catch (exception& e)
	{
		throw ModelException(e, routine);
	}
}
///////////////////////////////////////////////////////////////////////////////
DateTime    Instrument2TestMatrixesOfSPsRRs::getValueDate() const
{
	static const string routine 
                             = "Instrument2TestMatrixesOfSPsRRs::getValueDate";
	try
    {
        return valueDate;
	}
	catch (exception& e)
	{
		throw ModelException(e, routine);
	}
}
///////////////////////////////////////////////////////////////////////////////
/* for class loading (avoid having header file) */
bool Instrument2TestMatrixesOfSPsRRsLinkIn()
{
    return (Instrument2TestMatrixesOfSPsRRs::TYPE != 0);
}
///////////////////////////////////////////////////////////////////////////////
DRLIB_END_NAMESPACE
