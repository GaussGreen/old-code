/****************************************************************
 * Module:	PenGuin
 * Submodule:	
 * File:	vpatom.cxx
 * Function:	
 * test
 * Author:	David Liu
 ****************************************************************/
#include "kvpatom.h"
#include "kutilios.h"


//---------------------------------------------------------------
// 
void
KVPAtom::SetName(const char *name)
{
	if (name == NULL) {
		mName = String(format("%p", this));
	} else {
		mName = String(name);
	}
}


//---------------------------------------------------------------
// 
ostream& 
KVPAtom::PutRecursive(ostream& os) const
{
	int     idx;

	Put(os);
	os << GetName() << " DEPENDENCIES: " << NumDep() << endl;
	for (idx=0; idx<NumDep(); idx++) {
	    os << GetName() << format(" DEPENDENCY NO %d/%d ",
                               idx, NumDep()) << Dep(idx)->GetName() << endl;
	    Dep(idx)->PutRecursive(os);
	}

	return(os);
}
	


//---------------------------------------------------------------
 
static  inline  void    osIdent(ostream& os, int ident)
{
	int     k;
	for(k=0; k<ident; k++) os << " ";
}


//---------------------------------------------------------------
// 
ostream&
KVPAtom::PutTree(ostream& os, int ident)
{
	int	idx;
 
	osIdent(os, ident);
	os << GetName() << endl;
 
	for (idx=0; idx<NumDep(); idx++) {
		osIdent(os, ident);
		os << format("%d-> ", idx+1) << endl;
		Dep(idx)->PutTree(os, ident+4);
	}
 
	return(os);
}



//---------------------------------------------------------------
// 
ostream& 
KVPAtom::YacctionWriteRecursive(ostream& os)
{

	int     idx;

	for (idx=0; idx< NumDep(); idx++)
		Dep(idx)->YacctionWriteRecursive(os);

	YacctionWrite(os);

	return(os);
}

	
//---------------------------------------------------------------
// 
void	
KVPAtom::SetWriteFlag(bool flag)
{
	int	idx;
	mWriteFlag = flag;
	for (idx=0; idx< NumDep(); idx++)
		Dep(idx)->SetWriteFlag(flag);
	
}


//---------------------------------------------------------------
// 
int
KVPAtom::IsType(const char *name) const
{
	return(!strcmp(name, TypeName()));
}
 

