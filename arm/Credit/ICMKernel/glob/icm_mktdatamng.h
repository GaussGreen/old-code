/********************************************************************************/
/*! \file icm_mktdatamng.h
 * 
 *  \brief Market data manager withstd::string key
 *  \author 
 *	\version 1.0
 *	\date   febuary 2005 */
/*
 *********************************************************************************/

#ifndef _ICM_MARKET_DATA_MANAGER_H_
#define _ICM_MARKET_DATA_MANAGER_H_


#include <string>
#include "ARMKernel/glob/armglob.h"
#include "ARMKernel/glob/dates.h"
#include "ICMKernel/glob/icm_sensi_parameters.h"
#include <map>
// #include "gpinflation/infcurv.h"
// #include "ICMKernel/crv/icm_implied_loss_tree.h"

//using namespace ARM;
// using namespace ARM;
namespace ARM
{
	class ARM_InfCurv; 
} ; 

class ARM_ZeroCurve;
class ICM_DefaultCurve;
class ARM_VolCurve;
class ICM_Correlation;
class ICM_Fixing_Curve; 
class ICM_ImpLossTree; 

// temporaire pour model multicurve
#define IRCURVE "ircurve"
#define INFCURV "infcurve"
#define CORREL  "correlation"

/*********************************************************************************/
/*! \class  ICM_MktDataMng icm_mktdatamng.h "icm_mktdatamng.h"
 *  \author D Pouponneau
 *	\version 1.0
 *	\date   febuary 2005
 *	\brief  Creates a Market data manager withstd::string key
/***********************************************************************************/

std::string GetZeroCurveName(const const std::string& ccy,const ARM_Date& AsOf);
std::string GetInflaCurveName(const const std::string& label,const ARM_Date& AsOf);
std::string GetDefaultCurveName(const std::string& label,const ARM_Date& AsOf);
std::string GetVolCurveName(const std::string& label,const std::string& ccy,const ARM_Date& AsOf);
std::string GetCorrelationName(const std::string& structName,const std::string& index1, const std::string& index2,const ARM_Date& AsOf);
std::string GetLossTreeName(const std::string& label,const ARM_Date& AsOf);
std::string GetCorrelationIndex(const string& structname, const string& index1, const string& index2);
std::string GetFixingCurveName(const std::string&indexName,const std::string&ccy); 

void deducelabelforobject(const ARM_Object& ,std::string& label);

// typedef map<string,ARM_Object*> MapMgn;
//	the bool is "isOwner"
typedef std::map<std::string, std::pair<bool,ARM_Object*> > MapMgn;



class ICM_MktDataMng : public ARM_Object
{
private:

	MapMgn its_Manager;

public: 

	inline void Init(void) 
	{ 
		SetName(ICM_MKTDATAMNG);
		clear(); 
	}

	ICM_MktDataMng(void) {Init();};

	~ICM_MktDataMng() 
	{
		clear() ; 
	};
	
public:
	//	clears the manager : delete those items that the manager owns.
	//
	void clear()
	{
		MapMgn::iterator it = its_Manager.begin() ;
		while (it != its_Manager.end()) 
		{
			if (it->second.first && it->second.second) {
				delete it->second.second;
				it->second.second = NULL;
			}
			++it; 
		}
		its_Manager.clear() ;
	}

public:
	// this will insert a copy of the argument. null argument does nothing
	// the manager will be owner of the object and will delete it
	//		return true  : insertion OK (or null arg)
	//		return false : already an element
	//
	inline bool insert(ARM_Object*  object,std::string input = "AUTO")
	{
		//JLA 
		if (!object) return true; 
		string label;

		if (input != "AUTO") label = input;
		else deducelabelforobject(*object,label);

		if (its_Manager.find(label) == its_Manager.end())
		{
			its_Manager[label]= std::make_pair(true,object->Clone());
			return true;
		}
		return false;
	}
	void insert(const ARM_Object&object,std::string input="AUTO")
	{
		std::string label;
		if (input != "AUTO") label = input;
		else deducelabelforobject(object,label);
		erase(label); 
		its_Manager[label]= std::make_pair(true,unconst(object).Clone());
	}
	// this will insert the argument. null argument does nothing
	// the manager will take ownership of the object and will delete it
	//
	void adopt(ARM_Object*  object,std::string input = "AUTO")
	{
		//JLA 
		if (!object) return ; 
		std::string label;
		if (input != "AUTO") label = input;
		else deducelabelforobject(*object,label);
		erase(label); 
		its_Manager[label]= std::make_pair(true,object);
	}
	// this will insert the argument without copy, null argument does nothing
	// the manager is not owner of the item and will not delete it. 
	//	
	void associate(ARM_Object* object,std::string input="AUTO") 
	{
		//JLA 
		if (!object) return ;
		std::string label;
		if (input != "AUTO") label = input;
		else deducelabelforobject(*object,label);
		erase(label); 
		its_Manager[label]= std::make_pair(false,object);
	}
	//	find the specified item; returns null if not found.
	//
	inline ARM_Object* find(const std::string& input) const 
	{
		MapMgn::const_iterator	it=its_Manager.find(input) ;
		if (it == its_Manager.end()) return NULL;
		return it->second.second ;
	}
	//	detachs the specified item dives the address of the detached pair object
	std::pair<bool,ARM_Object*>* detach(const std::string& input) 
	{
		MapMgn::iterator	it=its_Manager.find(input) ;
		if (it == its_Manager.end()) return NULL; 
		if (it->second.first) its_Manager.erase(it) ;
		return &(it->second);
	}
	//	removes the specified item. 
	void erase(const std::string& input) 
	{
		MapMgn::iterator	it=its_Manager.find(input) ;
		if (it == its_Manager.end()) return ; 
		if (it->second.first && it->second.second) { delete it->second.second; it->second.second= NULL;} 
		its_Manager.erase(it) ;
	}

	//	returns true if all mkdata keys are present
	//
	bool CheckMktData(const std::vector<std::string>& mktdata) const 
	{
		for (int i=0; i<mktdata.size(); i++)
		{ 
			if
				(find(mktdata[i])==NULL)
				return false;
		}
		return true;
	}

	void GetMktDataVectotKey(std::vector<std::string>& mktdata) const ;
	// ----------------------------
	//	Copy of members data
	// ----------------------------
	void BitwiseCopy(const ARM_Object* src)
	{
		ICM_MktDataMng* Mbg = (ICM_MktDataMng*) src;
		// its_Manager.clear();
		clear() ;
		MapMgn::iterator ite;
		for (ite=Mbg->its_Manager.begin();ite != Mbg->its_Manager.end();ite++)
		{
			its_Manager[(*ite).first]= std::make_pair(true,ite->second.second->Clone());
		}
	}

	// -------------
	//	Copy Method 
	// -------------
	void Copy(const ARM_Object* src)
	{
		ARM_Object::Copy(src);
 		BitwiseCopy(src);
	}

	// --------------
	//	Clone Method
	// --------------
	ARM_Object* Clone(void)
	{
     ICM_MktDataMng* theClone = new ICM_MktDataMng();
     theClone->Copy(this);
     return(theClone);
	}

	void View(char* id, FILE* ficOut)
	{
	FILE* fOut;
	char  fOutName[200];
	int i = 0;

	 if ( ficOut == NULL )
		{
			ARM_GetViewFile(ARM_VIEW_FILE, id, fOutName);

			(void) unlink(fOutName);

			fOut = fopen(fOutName, "w"); 
		}
		else
		{
			fOut = ficOut;
		} 

		fprintf(fOut, "\n-----------------------------------------------------\n");
		fprintf(fOut, "---------------- Market Data Manager ----------------\n");
		fprintf(fOut, "-----------------------------------------------------\n\n");
		
		if (its_Manager.size()>0)
		{

		MapMgn::iterator ite;

		for (ite=its_Manager.begin(),i=0;ite != its_Manager.end();ite++,i++)
		{ fprintf(fOut, "DATAMARKET N°%i:%s:%s\n",i,
			ite->second.first?"owner":"not owner",
			((*ite).first.c_str())); }

		fprintf(fOut,"\n");

		for (ite=its_Manager.begin();ite != its_Manager.end();ite++)
		{
			fprintf(fOut, "\n\nDATAMARKET : %s\n",((*ite).first.c_str()));
			ite->second.second->View(id,fOut);
		}

		}

		if ( ficOut == NULL )
		{
			fclose(fOut);
		}
	}

	// All these functions will return 0 if not found. 
	//
	ARM_ZeroCurve* GetZeroCurve(const std::string& ccy,const ARM_Date& AsOf) const ;
	ARM::ARM_InfCurv*   GetInflaCurve(const const std::string& label,const ARM_Date& AsOf) const ;
	ICM_DefaultCurve* GetDefaultCurve(const std::string& label,const ARM_Date& AsOf) const ;
	ARM_VolCurve* GetVolCurve(const std::string& label,const std::string& ccy,const ARM_Date& AsOf) const ;
	ICM_Correlation* GetCorrelation(const std::string& label,const std::string& index1, const std::string& index2,const ARM_Date& AsOf) const ;
	ICM_Correlation* GetCorrelation(const std::string& label,const ARM_Date& AsOf) const ;
	ICM_ImpLossTree* GetLossTree(const std::string& label,const ARM_Date& AsOf) const ;
	ICM_Fixing_Curve* GetFixingCurve(const std::string&indexName,const std::string&ccy) const; 	

	// temporaire pour le model multicurve
	ARM_ZeroCurve* GetUnicZeroCurve() const ;
	ARM::ARM_InfCurv*   GetUnicInflaCurve() const ;
	ICM_Correlation* GetUnicCorrelation() const ;
	//	Apply a Shift the the MarketDataManager. 
	//	
	ICM_MktDataMng* GenerateShiftMng(ICM_BUMP<ICM_BASIS_BUMP>& parameters);
	ICM_MktDataMng* GenerateShiftMng(qSENSITIVITY_TYPE typesensi,
										const std::string& plot, 
										const std::string& label,
										double epsilon );

private:
	ICM_MktDataMng& operator=(const ICM_MktDataMng&); // NA
	ICM_MktDataMng(const ICM_MktDataMng&); 

	// ARM_ZeroCurve* GetZeroCurveFromKey(const std::string& Key) const ;
	// ARM::ARM_InfCurv*   GetInflaCurveFromKey(const std::string& Key) const ;
	// ICM_DefaultCurve* GetDefaultCurveFromKey(const std::string& Key) const ;
	// ARM_VolCurve* GetVolCurveFromKey(const std::string& Key) const ;
	// ICM_Correlation* GetCorrelationFromKey(const std::string& Key) const ;
	// ICM_ImpLossTree* GetLossTreeFromKey(const std::string& Key) const ;
	// ICM_Fixing_Curve* GetFixingCurveFromKey(const std::string& key) const; 
};

#endif
