/************************************************************************
 * Module:	PenGuin
 * File:
 * Function:    
 * Author:      C. Daher
 * Revision:	$Header$
 ************************************************************************/
#include "vpbundle.h"
#include "kutilios.h"

//---------------------------------------------------------------

istream&
KVPWBundle::Get(istream& is, int drw)
{

	return(is);
}


//---------------------------------------------------------------

ostream& 
KVPWBundle::Put(ostream& os, int indent) const
{
	int	idx;

	os << "NAME: " << GetName() << endl;
	os << "DEPENDENCIES: " << NumDep() << endl;
	if (NumDep() != 0) {
		os << "     IDX      WEIGHT   NAME" << endl;
		for (idx=0; idx<NumDep(); idx++) {
			os << format("  %3d/%3d   %9.5f  ",
				idx, NumDep(), mWeights[idx]);
			os << Dep(idx)->GetName();
			os << endl;
		}
	}



	return(os);
}





ostream& 
KVPWBundle::YacctionWrite(ostream& os, int indent)
{
	int	idx;

	if (GetWriteFlag())
	{	
	    os << GetName() << "=SUM(" << endl;
	    if (NumDep() != 0) {
		for (idx=0; idx<NumDep()-1; idx++) {
			os << format("%15s,  %9.5f,\n",
				Dep(idx)->GetName(), mWeights[idx]);
		}

		//
		// The last dependent. No "," at the end
		os << format("%15s,  %9.5f);\n",
				Dep(idx)->GetName(), mWeights[idx]);
	    }
	    os << endl << endl;

	    WriteDone();
	}

	return(os);

}



KVPWBundle&
KVPWBundle::AddDepWeight(
	const SharedPointer<KVPInstr> instr,
	double weight)
{

	SharedPointer<KVPAtom> instrAtom;

	SharedPointerConvertTo(instr, instrAtom);
	KVPAtom::AddDep(instrAtom);
	mWeights.insert(mWeights.end(), weight);


	return(*this);

}


//
// Return the first dependent discount curve name
//
const String&
KVPWBundle::GetDiscName() const
{
	SharedPointer<KVPInstr> instr;

	SharedPointerDownCastTo(Dep(0), instr);

    return instr->GetDiscName();
}
