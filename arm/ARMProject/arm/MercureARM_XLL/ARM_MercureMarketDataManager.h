#ifndef ARM_MERCURE_MARKETDATAMANAGER_H
#define ARM_MERCURE_MARKETDATAMANAGER_H

#include "armglob.h"
#include "MarketDataManager.h"
#include <string>

using namespace std;
using namespace mercure;


class ARM_MercureMarketDataManager : public ARM_Object
{
	private:
		string	itsAsOfDate;
		MarketDataManager*		itsMarketDataManager;
		map<string, ARM_Object*>	itsMarketDataList;

	public:
		ARM_MercureMarketDataManager(void)
		{
			itsMarketDataManager = NULL;
			SetName(ARM_MERCURE_MARKETDATAMANAGER);
		}

		ARM_MercureMarketDataManager(MarketDataManager* aMarketDataManager, string aDate)
		{
			itsMarketDataManager = aMarketDataManager;
			itsAsOfDate = aDate;
			SetName(ARM_MERCURE_MARKETDATAMANAGER);
		}

		~ARM_MercureMarketDataManager(void)
		{
			delete	itsMarketDataManager;
			itsMarketDataManager = NULL;	
		}

		MarketDataManager*	GetMarketDataManager(void)
		{
			return	itsMarketDataManager;
		}

		string	GetAsOfDate(void)
		{
			return	itsAsOfDate;
		}

		map<string, ARM_Object*>&	GetMarketDataList()
		{
			return	itsMarketDataList;
		}

		virtual void View(char* id = NULL, FILE* ficOut = NULL);
};


#endif