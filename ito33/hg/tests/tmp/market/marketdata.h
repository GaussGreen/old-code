#ifndef _MARKETDATA_H_
#define _MARKETDATA_H_

#include <list>

#include "ito33/sharedptr.h"

using namespace ito33;
using namespace ito33::finance;
using namespace ito33::hg;


class ito33::finance::CDS;
class ito33::finance::EDS;
class ito33::finance::Option;
class ito33::finance::SessionData;
class ito33::hg::UnderlyingProcess;


/// The companies for which we currently have market data for testing
enum Company
{
  Accor,
  Carrefour,
  FranceTelecom,
  Lafarge,
  Valeo,
  BNBParibas,
  Endesea,
  //Iberdrola,
  Renault
};


/// Read market data from the given file
std::list< shared_ptr<Option> > ReadOptionData(std::string sFilename,
                                              shared_ptr<SessionData> pSessionData);

/// Read market data from the given file using format from 2nd (Elie) data set
std::list< shared_ptr<Option> > ReadOptionData2(std::string sFilename,
                                               shared_ptr<SessionData> pSessionData);

/// Get cds data for the given company
std::list< shared_ptr<CDS> > GetCDSData(Company company, 
                                       shared_ptr<SessionData> pSessionData);

/// Get eds data for the given company
std::list< shared_ptr<EDS> > GetEDSData(Company company,
                                       shared_ptr<SessionData> pSessionData);

/// Setup the sessio for the specified company
shared_ptr<SessionData> InitSessionDataMarket(Company company);

/// Return which file to read
std::string GetInputFile(Company company, bool bIsAdi);

/// Return a starting underlying process for calibration
shared_ptr<hg::UnderlyingProcess> GetUnderlyingProcess(Company company);

#endif // #ifndef _MARKETDATA_H_
