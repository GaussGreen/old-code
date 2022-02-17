#ifndef ARM_LOCAL_XXXProject_H
#define ARM_LOCAL_XXXProject_H

#include "ARM_result.h"


#include <string>
#include <vector>

using std::string;
using std::vector;

#include <GP_Infra\gpinfra\gramfunctorarg.h>
#include "ARM_local_gp_genericaddin.h"
using ARM::ARM_GramFctorArg;


extern long ARMLOCAL_MktDatas_Create(	const double&			asOfDate,
										const vector<string>&	keys,
										const vector<long>&		mktDatas,
										ARM_result&				result,
										long					objId = -1);

extern long ARMLOCAL_ARM_XXX_Price(		long					secId,
										long					mktDt,
										ARM_result&				result);

extern long ARMLOCAL_Hedge_Create(		const long&				secId, 
										const long&				scenId,
										const long&				mktDt,
 										ARM_result&				result, 
										long					objId=-1);

extern long ARMLOCAL_Hedge_GetData(		long					hedgeId,
										const string&			key,
										ARM_GramFctorArg&		argResult,
										ARM_result&				result );


extern long ARMLOCAL_ARM_XXX_Scenario(	const double&			shift,
										const string&			currency,
										const string&			type_Scenario,
										const string&			subType_Scenario,
										const string&			stress_Order,
										const long&				relatif,
										const long&				cumulInv,
										const long&				perturbative,
										ARM_result&				result,
										long					objId = -1);

extern long ARMLOCAL_Scenari_Compose(	const long&					scenariId1,
										const long&					scenariId2,
										ARM_result&				result,
										long					objId=-1);

class ARM_GP_XXXMktDataFromMktDataMgerFunctor : public ARM_GenericAddinFunctor
{
public:
	ARM_GP_XXXMktDataFromMktDataMgerFunctor() {}

	virtual	long operator()( ARM_result& result, long objId );
};

class XXX_MktData_ApplyScenarioFunctor : public ARM_GenericAddinFunctor
{
public:
	XXX_MktData_ApplyScenarioFunctor() {}

	virtual	long operator()( ARM_result& result, long objId );
};

#endif