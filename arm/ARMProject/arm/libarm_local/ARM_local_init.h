#ifndef _ARM_LOCAL_INIT_H
#define _ARM_LOCAL_INIT_H


#include <CInitReader.h>
#include <vector>



class ARMLOCAL_Init : public CInitReader 
{
   public:
      
    ARMLOCAL_Init(const char* pFileName);

	   string default_currency;
	   string data_folder;

	   string wsetkprod;
	   string wsetkrec;
	   string wspricingprod;
	   string wspricingrec;

	   vector<string>	currencyList;
	   vector<int>		indirectList;

	virtual void Init ();

	string ToString ();
};


#endif
/*---- End Of File ----*/