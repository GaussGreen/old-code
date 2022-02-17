#include "ARMKernel\glob\firsttoinc.h"
#include "ICMKernel\util\icm_macro.h"
#include "ICMKernel\util\icm_utils.h"
#include "ICMKernel/mod/icm_lossunits.h"
#include <iomanip>
// local helper
//-----------------------------------------------------------------------------------------
static inline ICM_LossUnits::cont_t::const_iterator 
find(const ICM_LossUnits::cont_t& cont,const ARM_Date&date) 
{
	if (cont.empty()) return cont.end(); 
	ICM_LossUnits::cont_t::const_iterator it = cont.begin(); 
	if (date<it->first) return it ;
	while  ( (++it)  != cont.end()) 
	{
		if (it->first>date) return --it; 
	}
	return --it; 
}

//-----------------------------------------------------------------------------------------
ICM_LossUnits::ICM_LossUnits()
{
}
//-----------------------------------------------------------------------------------------
ICM_LossUnits::ICM_LossUnits(const ICM_LossUnits&ref)
{
	itsLossUnits=ref.itsLossUnits ;
}
//-----------------------------------------------------------------------------------------
ICM_LossUnits& ICM_LossUnits::operator=(const ICM_LossUnits&ref)
{
	if (this!=&ref)
	{
		this->~ICM_LossUnits() ;
		new(this)ICM_LossUnits(ref); 
	}
	return *this; 
}
//-----------------------------------------------------------------------------------------
ICM_LossUnits::~ICM_LossUnits()
{}

//-----------------------------------------------------------------------------------------
void 
ICM_LossUnits::insert(const ARM_Date&date,const lossunit_data_t&item)
{
	itsLossUnits[date]=item ;
}
//-----------------------------------------------------------------------------------------
double 
ICM_LossUnits::getLossUnit(const ARM_Date&date) const
{
	cont_t::const_iterator it = find(itsLossUnits,date); 
	if (it==itsLossUnits.end()) 
		ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_LossUnits::getLossUnit: empty"); 
	return it->second.LU; 
}

bool ICM_LossUnits::IsLossUnitNull(const ARM_Date&date) const
{
	cont_t::const_iterator it = find(itsLossUnits,date); 
	if (it==itsLossUnits.end()) 
		ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_LossUnits::getLossUnit: empty"); 
	return it->second.IsNull; 
}

//-----------------------------------------------------------------------------------------
const ARM_Vector& 
ICM_LossUnits::getLossRates(const ARM_Date&date) const
{
	cont_t::const_iterator it = find(itsLossUnits,date); 
	if (it==itsLossUnits.end()) 
		ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_LossUnits::getLossUnit: empty"); 
	return it->second.LossRates; 
}
//-----------------------------------------------------------------------------------------
const ARM_Vector& 
ICM_LossUnits::getLSLossRates(const ARM_Date&date) const
{
	cont_t::const_iterator it = find(itsLossUnits,date); 
	if (it==itsLossUnits.end()) 
		ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_LossUnits::getLossUnit: empty"); 
	return it->second.LSLossRates; 
}
//-----------------------------------------------------------------------------------------
const std::vector<int>& 
ICM_LossUnits::getCollatRank(const ARM_Date&date) const
{
	cont_t::const_iterator it = find(itsLossUnits,date); 
	if (it==itsLossUnits.end()) 
		ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_LossUnits::getLossUnit: empty"); 
	return it->second.CollatRank; 
}
//-----------------------------------------------------------------------------------------
const ARM_Date& 
ICM_LossUnits::getDate(unsigned int i) const
{
	if (i>=itsLossUnits.size()) 
		ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_LossUnits::getDate: out of bounds "<<i); 
	cont_t::const_iterator it = itsLossUnits.begin(); 
	for(unsigned int k=0;k<i;k++) { ++it; }
	return it->first ;
}
//-----------------------------------------------------------------------------------------
void 
ICM_LossUnits::getDates(ICM_Vector<ARM_Date>&ret) const
{
	ret.resize(itsLossUnits.size()); 
	cont_t::const_iterator it = itsLossUnits.begin(); 
	unsigned int i=0; 
	while (it!=itsLossUnits.end())
	{
		ret[i]=it->first; 
		++it; ++i;
	}
}
//-----------------------------------------------------------------------------------------
// virtual 
void 
ICM_LossUnits::View(char* id , FILE* ficOut )
{
	FILE* fOut;
	char  fOutName[200];
	if ( ficOut == NULL )
	{
		ARM_GetViewFile(ARM_VIEW_FILE, id, fOutName);
		(void) unlink(fOutName);
		fOut = fopen(fOutName, "w"); 
	}
	else  {fOut = ficOut;} 

	int size = itsLossUnits.size() ;

	fprintf(fOut,"--------------------- LossUnits View ----------------------------- \n"); 
	fprintf(fOut,"LossUnit size: %d\n",size); 
	if (size!=0)
	{
		int nbIssuers = itsLossUnits.begin()->second.LossRates.size(); 
		std::stringstream sstr1,sstr2; 
		cont_t::const_iterator it = itsLossUnits.begin(); 	
		while (it != itsLossUnits.end()) 
		{
			sstr1<<it->first<<"\t"; 
			++it; 
		}
		fprintf(fOut,"%s\n",sstr1.str().c_str()); 
		// sstr2<<std::
		it = itsLossUnits.begin();
		while (it != itsLossUnits.end()) 
		{
			sstr2<<"LU="<<std::setw(8)<<it->second.LU<<"\t"; 
			++it; 
		}
		fprintf(fOut,"%s\n",sstr2.str().c_str()); 
		for(int i = 0;i<nbIssuers;i++) 
		{
			cont_t::const_iterator it = itsLossUnits.begin(); 	
			std::stringstream sstr; 
			while (it != itsLossUnits.end()) 
			{
				sstr<<"LR="<<std::setw(8)<<it->second.LossRates.Elt(i) ;
				sstr<<",RNK="<<it->second.CollatRank[i] ;
				sstr<<"\t"; 
				++it ;
			}
			fprintf(fOut,"%s\n",sstr.str().c_str()); 
		}
	}
	fprintf(fOut, "\n");
	if ( ficOut == NULL )
	{fclose(fOut);}
}

