

#include "ICMKernel\util\icm_matrix.h"


//	-----------------------------------------------------------------------------------
bool 
ICM_Parameters::getParam(const std::string&paramName,double&ret,bool throwOnError) const
{
	ARM_Vector local; 
	if (!getParam(paramName,local,throwOnError)) return false ; 
	if (local.GetNumLines()==0) 
	{
		if (throwOnError) ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_Parameters::getParam "<<paramName<<": empty value! "); 
		return false ; 
	}
	ret=local.Elt(0); 
	return true; 
}
//	-----------------------------------------------------------------------------------
bool 
ICM_Parameters::getParam(const std::string&paramName,long&ret,bool throwOnError) const
{
	ARM_Vector local; 
	if (!getParam(paramName,local,throwOnError)) return false ; 
	if (local.GetNumLines()==0) 
	{
		if (throwOnError) ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_Parameters::getParam "<<paramName<<": empty value! "); 
		return false ; 
	}
	ret=local.Elt(0); 
	return true; 
}
//	-----------------------------------------------------------------------------------
bool 
ICM_Parameters::getParam(const std::string&paramName,ARM_Vector&ret,bool throwOnError) const
{
	ARM_Vector* item = GetColVect(paramName) ;
	if (item==0) 
	{
		ret.Resize(0); 
		if (throwOnError) ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_Parameters::getParam "<<paramName<<": not found !"); 
		return false ;
	}
	ret=*item; 
	return true; 
}
//	--------------------------------------------------------------------
// virtual 
ICM_PropertyList::~ICM_PropertyList()
{
	cont_t::iterator it = itsValues.begin(); 
	while (it!=itsValues.end()) { delete it->second; ++it; }
}
//	--------------------------------------------------------------------
ICM_PropertyList::ICM_PropertyList(const ICM_PropertyList&ref) : ARM_Object(ref)
{
	cont_t::const_iterator it = ref.itsValues.begin(); 
	while (it!=ref.itsValues.end()) 
	{ 
		itsValues[it->first]=it->second->clone() ;
		++it; 
	}
}
//	--------------------------------------------------------------------
// virtual 
ARM_Object*
ICM_PropertyList::Clone() 
{
	return new ICM_PropertyList(*this); 
}
//	--------------------------------------------------------------------
// virtual 
void ICM_PropertyList::View(char* id, FILE* ficOut )
{
	std::stringstream sstr; 
	sstr<<" ----------- ICM_PropertyList -----------------"<<std::endl; 
	sstr<<"Size: "<<itsValues.size()<<std::endl; 
	cont_t::const_iterator it = itsValues.begin(); 
	while (it!=itsValues.end()) 
	{
		sstr<<it->first<<":"<<*it->second<<std::endl; 
		++it; 
	}

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

	fprintf(fOut, "%s\n",sstr.str().c_str() );

	if ( ficOut == NULL )
	{
		fclose(fOut);
	}
}
