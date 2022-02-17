//----------------------------------------------------------------------------
//
//   Group       : Credit Hybrids QR
//
//   Filename    : CIDClosedFormProduct2TestMatrixesOfSPsRRs.cpp
//
//   Description : Instrument we use to test matrices Of SPs and RRs
//                                                     + CID Closed Form Model
//
//   Date        : Aug 2006
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/CIDClosedFormProduct2TestMatrixesOfSPsRRs.hpp"
#include "edginc/OutputRequest.hpp"
#include "edginc/Results.hpp"

DRLIB_BEGIN_NAMESPACE
///////////////////////////////////////////////////////////////////////////////
CIDClosedFormProduct2TestMatrixesOfSPsRRs
                  ::CIDClosedFormProduct2TestMatrixesOfSPsRRs(
                           const Instrument2TestMatrixesOfSPsRRs * gInstrument)
    : instrument(gInstrument)
{;};
///////////////////////////////////////////////////////////////////////////////
void CIDClosedFormProduct2TestMatrixesOfSPsRRs::price( 
                                        CIDClosedFormModel * gModel,
                                        CControl           * gControl, 
                                        CResults           * gResults) const
{
    static const string routine 
                          = "CIDClosedFormProduct2TestMatrixesOfSPsRRs::price";
    try
    {
        OutputRequest* _outputRequest = 0; 
        if (gControl->requestsOutput(OutputRequest::DBG, _outputRequest))
        {
            DateTimeArrayConstSP _dates = instrument->getDates();
            const StringArray  & _names(*instrument->getNames());
            const int _qtyOfDates = _dates->size();
            const int _qtyOfNames = _names.size();
            CDoubleMatrixSP _SPs(new DoubleMatrix(_qtyOfNames, _qtyOfDates));
            const CIDParameters & _CIDparams(gModel->getCIDParameters());
            int _n; for(_n = 0; _n < _qtyOfNames; _n++)
            {
                CIDNameRecordConstSP pNameRecord = _CIDparams.getSingleNameRecord(_names[_n]);
                int _d; for(_d = 0; _d < _qtyOfDates; _d++)
                {
                      (*_SPs)[_n][_d]
                                  =pNameRecord->calcProperSurvProbability(instrument->getValueDate(),(*_dates)[_d]);
                }
            }
            gResults->storeRequestResult(_outputRequest, _SPs);
        }
    }
    catch (exception& e)
    {
        throw ModelException(e, routine);
    }
}
///////////////////////////////////////////////////////////////////////////////
DRLIB_END_NAMESPACE
