/**
 * @file Instrument_forward.hpp
 */

#ifndef DRLIB_Instrument_forward_H
#define DRLIB_Instrument_forward_H

#include "edginc/Array.hpp"

DRLIB_BEGIN_NAMESPACE

class CInstrument;
typedef CInstrument Instrument;
typedef smartConstPtr<CInstrument> CInstrumentConstSP;
typedef smartPtr<CInstrument> CInstrumentSP;
#ifndef QLIB_INSTRUMENT_CPP
EXTERN_TEMPLATE(class RISKMGR_DLL_SP smartConstPtr<CInstrument>);
EXTERN_TEMPLATE(class RISKMGR_DLL_SP smartPtr<CInstrument>);
#else
INSTANTIATE_TEMPLATE(class RISKMGR_DLL smartConstPtr<CInstrument>);
INSTANTIATE_TEMPLATE(class RISKMGR_DLL smartPtr<CInstrument>);
#endif
typedef CInstrumentConstSP InstrumentConstSP;
typedef CInstrumentSP InstrumentSP;
typedef array<CInstrumentSP, CInstrument> CInstrumentArray;
typedef smartPtr<CInstrumentArray> CInstrumentArraySP;
typedef smartConstPtr<CInstrumentArray> CInstrumentArrayConstSP;
#ifndef QLIB_INSTRUMENT_CPP
EXTERN_TEMPLATE(class RISKMGR_DLL array<CInstrumentSP _COMMA_ CInstrument>);
EXTERN_TEMPLATE(class RISKMGR_DLL_SP smartConstPtr<CInstrumentArray>);
EXTERN_TEMPLATE(class RISKMGR_DLL_SP smartPtr<CInstrumentArray>);
#else
INSTANTIATE_TEMPLATE(class RISKMGR_DLL array<CInstrumentSP _COMMA_ CInstrument>);
INSTANTIATE_TEMPLATE(class RISKMGR_DLL smartConstPtr<CInstrumentArray>);
INSTANTIATE_TEMPLATE(class RISKMGR_DLL smartPtr<CInstrumentArray>);
#endif

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***

#endif
