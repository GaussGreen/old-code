/************************************************************************
 * Module:	PenGuin
 * File:
 * Function:    
 * Author:      C. Daher, David Liu
 * Revision:	$Header$
 ************************************************************************/
#include "vtlbase.h"

#include "vtlbundle.h"
#include "vtlcashfl.h"
#include "vtloption.h"
#include "vtlkio.h"
#include "vtlkio2idx.h"
#include "vtlfleg.h"
#include "vtlfleg2Idx.h"
#include "vtlfleg3Idx.h"
#include "vtlrate.h"
#include "vtlprotleg.h"
#include "vtlfleg.h"
#include "vtlfleg2Idx.h"
#include "vtlfleg3Idx.h"
#include "vtldefprotect.h"

#include "kutilios.h"


//===============================================================
//
//===============================================================



//---------------------------------------------------------------


void
KVPToolAtom::Initialize()
{
	//
	// We reset tehe counter at #tp + 1, so that
	// we can decrease it by one at first step.
	//
	SetFlag(mVTree->TPNum() + 1);

}

//---------------------------------------------------------------



int
KVPToolAtom::NeedsUpdate()
{
	int	tpIdx = mVTree->TPIdxCurrent();
	int	idx;

	if (GetFlag() == tpIdx+1) {
		// We last updated the tree at the previous node.
		// Need to update at current node.
		for (idx=0; idx<NumDep(); idx++) {
			Dep(idx)->Update();
		}

		return (TRUE);
	} else if (GetFlag() == tpIdx) {
		// already done, just skip
		//
		return (FALSE);
	} else {
		throw KFailure("%s: inconsistent update: "
			"(current tpIdx %d != last update tpIdx %d - 1).\n",
			GetName(), tpIdx, GetFlag());
	}
	return (FALSE);

}


//---------------------------------------------------------------

void
KVPToolAtom::UpdateDone()
{
	// We set the flag to the last update tpIdx
	//
	SetFlag(mVTree->TPIdxCurrent());
}

//---------------------------------------------------------------

 
KMap(String,double)
KVPToolAtom::GetResults()
{
static	char	routine[] = "KVPToolAtom::GetResults";
	KMap(String,double)	results;

	if (mVTree->TPIdxCurrent() != 0) {
		throw KFailure("%s: current tpIdx (%d) != 0 on %s.\n",
			routine, mVTree->TPIdxCurrent(), GetName());
	}

	// Compute pv and add to result map
	//
	double	pv = GetValue().GetCenter();
	results.insert(KMap(String, double)::value_type("PV", pv));

	return (results);
}



//---------------------------------------------------------------


static	SharedPointer<KVPToolAtom>
KVPToolAtomNewRoot(
  const	SharedPointer<KVPAtom> ins,  // (I) instrument from which tool created
	KVTree &vt,		     // (I) virtual tree for evaluation
	SharedPointer<KVPToolAtom> preTool,     // (B) NULL if newly created
	SharedPointer<KVPToolAtom> rootTool)	// (B) NULL if newly created
{
static	char	routine[] = "KVPToolAtomNewRoot";
	SharedPointer<KVPToolAtom>	tool;
	int				idx;


    try {
	if(debugLevel)
		dppLog << "Processing " << ins->GetName() << endl;

	//
	//
	//
	if (rootTool) {
		tool = rootTool->FindAtom(ins);
	}


	if (tool == NULL) {
		//
		// Not found, create new tool
		//
		if (ins->IsType("KVPCashFlows")) {
		    //
		    // It is a KVPCashFlows
		    //
		    SharedPointer<KVPCashFlows> insA;
		    SharedPointerDownCastTo(ins, insA);

		    ASSERT_OR_THROW(insA != NULL);	// just in case !

		    SharedPointerConvertTo(
			Raw2SharedPointer(new KVPToolCashFlows(insA, vt)),
			tool);

		    ASSERT_OR_THROW(tool != NULL);

		} else if (ins->IsType("KVPOption")) {
		    //
		    // It is a KVPOption
		    //
		    SharedPointer<KVPOption> insA;
		    SharedPointerDownCastTo(ins, insA);

		    ASSERT_OR_THROW(insA != NULL);	// just in case !

		    SharedPointerConvertTo(
			Raw2SharedPointer(new KVPToolOption(insA, vt)),
			tool);

		    ASSERT_OR_THROW(tool != NULL);

		} else if (ins->IsType("KVPKnockIO")) {
		    //
		    // It is a KVPKnockIO
		    //
		    SharedPointer<KVPKnockIO> insA;
		    SharedPointerDownCastTo(ins, insA);

		    ASSERT_OR_THROW(insA != NULL);	// just in case !

		    SharedPointerConvertTo(
			Raw2SharedPointer(new KVPToolKnockIO(insA, vt)),
			tool);

		    ASSERT_OR_THROW(tool != NULL);

		} else if (ins->IsType("KVPKnockIO2Idx")) {
		    //
		    // It is a KVPKnockIO2Idx
		    //
		    SharedPointer<KVPKnockIO2Idx> insA;
		    SharedPointerDownCastTo(ins, insA);

		    ASSERT_OR_THROW(insA != NULL);	// just in case !

		    SharedPointerConvertTo(
			Raw2SharedPointer(new KVPToolKnockIO2Idx(insA, vt)),
			tool);

		    ASSERT_OR_THROW(tool != NULL);

		} else if (ins->IsType("KVPWBundle")) {
		    //
		    // It is a KVPWBundle
		    //
		    SharedPointer<KVPWBundle> insA;
		    SharedPointerDownCastTo(ins, insA);

		    ASSERT_OR_THROW(insA != NULL);	// just in case !

		    SharedPointerConvertTo(
			Raw2SharedPointer(new KVPToolWBundle(insA, vt)),
			tool);

		    ASSERT_OR_THROW(tool != NULL);

		} else if (ins->IsType("KVPFloatLeg")) {
		    //
		    // It is a KVPFloatLeg
		    //
		    SharedPointer<KVPFloatLeg> insA;
		    SharedPointerDownCastTo(ins, insA);

		    ASSERT_OR_THROW(insA != NULL);	// just in case !

		    SharedPointerConvertTo(
			Raw2SharedPointer(new KVPToolFloatLeg(insA, vt)),
			tool);

		    ASSERT_OR_THROW(tool != NULL);

		} else if (ins->IsType("KVPFloatLeg2Idx")) {
		    //
		    // It is a KVPFloatLeg2Idx
		    //
		    SharedPointer<KVPFloatLeg2Idx> insA;
		    SharedPointerDownCastTo(ins, insA);

		    ASSERT_OR_THROW(insA != NULL);	// just in case !

		    SharedPointerConvertTo(
			Raw2SharedPointer(new KVPToolFloatLeg2Idx(insA, vt)),
			tool);

		    ASSERT_OR_THROW(tool != NULL);

		} else if (ins->IsType("KVPFloatLeg3Idx")) {
		    //
		    // It is a KVPFloatLeg3Idx
		    //
		    SharedPointer<KVPFloatLeg3Idx> insA;
		    SharedPointerDownCastTo(ins, insA);

		    ASSERT_OR_THROW(insA != NULL);	// just in case !

		    SharedPointerConvertTo(
			Raw2SharedPointer(new KVPToolFloatLeg3Idx(insA, vt)),
			tool);

		    ASSERT_OR_THROW(tool != NULL);

		} else if (ins->IsType("KRate")) {
		    //
		    // It is a KRate
		    //
		    SharedPointer<KRate> insA;
		    SharedPointerDownCastTo(ins, insA);

		    ASSERT_OR_THROW(insA != NULL);	// just in case !

		    SharedPointerConvertTo(
			Raw2SharedPointer(new KVPToolRate(insA, vt)),
			tool);

		    ASSERT_OR_THROW(tool != NULL);

		} else if (ins->IsType("KVPProtLeg")) {
		    //
		    // It is a KVPProtLeg
		    //
		    SharedPointer<KVPProtLeg> insA;
		    SharedPointerDownCastTo(ins, insA);

		    ASSERT_OR_THROW(insA != NULL);	// just in case !

		    SharedPointerConvertTo(
			Raw2SharedPointer(new KVPToolProtLeg(insA, vt)),
			tool);

		    ASSERT_OR_THROW(tool != NULL);

		} else if (ins->IsType("KVPDefProtect")) {
		    //
		    // It is a KVPDefProtect
		    //
		    SharedPointer<KVPDefProtect> insA;
		    SharedPointerDownCastTo(ins, insA);

		    ASSERT_OR_THROW(insA != NULL);	// just in case !

		    SharedPointerConvertTo(
			Raw2SharedPointer(new KVPToolDefProtect(insA, vt)),
			tool);

		    ASSERT_OR_THROW(tool != NULL);

		} else if (ins->IsType("KVPFloatLeg")) {
		    //
		    // It is a KVPFloatLeg
		    //
		    SharedPointer<KVPFloatLeg> insA;
		    SharedPointerDownCastTo(ins, insA);

		    ASSERT_OR_THROW(insA != NULL);	// just in case !

		    SharedPointerConvertTo(
			Raw2SharedPointer(new KVPToolFloatLeg(insA, vt)),
			tool);

		    ASSERT_OR_THROW(tool != NULL);

		} else if (ins->IsType("KVPFloatLeg2Idx")) {
		    //
		    // It is a KVPFloatLeg2Idx
		    //
		    SharedPointer<KVPFloatLeg2Idx> insA;
		    SharedPointerDownCastTo(ins, insA);

		    ASSERT_OR_THROW(insA != NULL);	// just in case !

		    SharedPointerConvertTo(
			Raw2SharedPointer(new KVPToolFloatLeg2Idx(insA, vt)),
			tool);

		    ASSERT_OR_THROW(tool != NULL);

		} else if (ins->IsType("KVPFloatLeg3Idx")) {
		    //
		    // It is a KVPFloatLeg3Idx
		    //
		    SharedPointer<KVPFloatLeg3Idx> insA;
		    SharedPointerDownCastTo(ins, insA);

		    ASSERT_OR_THROW(insA != NULL);	// just in case !

		    SharedPointerConvertTo(
			Raw2SharedPointer(new KVPToolFloatLeg3Idx(insA, vt)),
			tool);

		    ASSERT_OR_THROW(tool != NULL);


		} else {
		    throw KFailure("%s: unable to process instrument "
				"%s of type %s.\n",
				routine,
				ins->GetName(),
				ins->TypeName());
		}

		//
		// If no rootTool passed, newly created becomes rootTool.
		//
		if (rootTool == NULL)
			rootTool = tool;

		if(debugLevel)
			dppLog << "Created " << tool->GetName() << endl;


		// Link
		if (preTool) {
			preTool->AddDep(tool);
			if (debugLevel)
				dppLog << "Linking " << preTool->GetName() 
					<< " to " << tool->GetName() << endl;
		}


		//
		// Process all depedencies
		//
		for (idx=0; idx<ins->NumDep(); idx++) {
			SharedPointer<KVPAtom>      depIns = ins->Dep(idx);
			SharedPointer<KVPToolAtom>  depTool;


			depTool = KVPToolAtomNewRoot(
				depIns,
				vt,
				tool,
				rootTool);
			ASSERT_OR_THROW(depTool != NULL);
		}


	} else {
		//
		// tool already exists
		//
		if(debugLevel)
			dppLog << "Found " << tool->GetName() << endl;


		// link
		if (preTool) {
			if (debugLevel)
				dppLog << "Linking " << preTool->GetName() 
					<< " to " << tool->GetName() << endl;
			preTool->AddDep(tool);
		}



		return(tool);
	}

	return(tool);
    }
    catch (...) {
	throw KFailure("%s: failed.\n", routine);
    }
}



/*---------------------------------------------------------------
 * Creates recursively a new KVPToolAtom and all corresponding
 * dependencies given a instrument and
 * a virtual tree.
 */


SharedPointer<KVPToolAtom>
NewToolRecursive(
  const	SharedPointer<KVPAtom> ins,	// (I) instrument to value
	KVTree &vt)			// (I) virtual tree
{
	SharedPointer<KVPToolAtom> rootTool;
	SharedPointer<KVPToolAtom> nullTool;

	if(debugLevel)
		dppLog << "Start creating root tool " << ins->GetName() << endl;

	rootTool = KVPToolAtomNewRoot(ins, vt, nullTool, nullTool);

	if(debugLevel)
		dppLog << "Done creating root tool " << ins->GetName() << endl;

	return rootTool;
}



//---------------------------------------------------------------



const char*
KVPToolAtom::GetName() 
{
	SharedPointer<KVPAtom>	p = Atom();

	mName = String(format("%s[%s]",
			TypeName(),
			(p ? p->GetName() : "????")
			));

	return(mName.c_str());
}


//---------------------------------------------------------------
// Find if an VPAtom is associated to a tool or its successors.
// The routine uses the Atom() method to retrieve the KVPAtom
// associated to a KVPToolAtom.
// Returns the corresponding tool or NULL.


SharedPointer<KVPToolAtom>
KVPToolAtom::FindAtom(SharedPointer<KVPAtom> ins)
{
	int				idx;
	SharedPointer<KVPToolAtom>	ptr;



	if (Atom() == ins) {
		/*dppLog << format("    KVPToolAtom::Find: "
			"Checking %s atom %p on %s-%p TRUE.\n",
			GetName(), Atom(), ins->GetName(), ins);*/
		return(Raw2SharedPointer(this));
	}

	for (idx=0; idx<NumDep(); idx++) {
		ptr = Dep(idx)->FindAtom(ins);
		if (ptr != NULL)
			return(ptr);
	}

	/*dppLog << format("    KVPToolAtom::Find: "
		"Checking %s atom %p on %s-%p TRUE.\n",
		GetName(), Atom(), ins->GetName(), ins);*/
	return(ptr);
}



//---------------------------------------------------------------


static	inline	void	osIdent(ostream& os, int ident)
{
	int	k;
	for(k=0; k<ident; k++) os << " ";
}

ostream& 
KVPToolAtom::Put(ostream& os, int ident) const
{
	int	idx;

	osIdent(os, ident);
	//os << GetName() << endl;

	for (idx=0; idx<NumDep(); idx++) {
		osIdent(os, ident);
		os << format("%d-> ", idx+1) << endl;
		Dep(idx)->Put(os, ident+4);
	}

	return(os);
}




ostream& 
PutPV(SharedPointer<KVPToolAtom> atom, ostream& os, int ident)
{
	int	idx;

	osIdent(os, ident);
	os << atom->GetName();
#ifdef	_SKIP
	os << format("  (PV: %.6f)", atom->PresentValue());
	os << endl;
#endif

	KMap(String,double) results = atom->GetResults();

	os << "\t\t";
	for(KMap(String,double)::iterator p=results.begin();
     	    p != results.end();
            ++p) {
		os << "  " << (*p).first << ":" << results[(*p).first];
	}
	os << "" << endl;




	for (idx=0; idx<atom->NumDep(); idx++) {
	    if (!(atom->Dep(idx)->IsType("KVPToolRate")))
	    {	
		osIdent(os, ident);
		os << format("%d-> ", idx+1) << endl;
		PutPV(atom->Dep(idx), os, ident+4);
	    }
	}

	return(os);
}



//---------------------------------------------------------------


void
KVPAtomKVTreeEvalNpv(
	SharedPointer<KVPAtom> vpRoot,	// (I) instrument to value
	KVTree &vt,			// (I) virtual tree for evaluation
	double *npv)			// (O) value

{
static	char	routine[] = "KVPAtomKVTreeEvalNpv";
	SharedPointer<KVPToolAtom>	vpToolRoot;

    try {


	dppLog << "====================================="
			"======================================" << endl;

	//
        // Perform Pricing
        //
 
        vpToolRoot = NewToolRecursive(vpRoot, vt);
        ASSERT_OR_THROW(vpToolRoot != NULL);
 
 
	dppLog << "====================================="
			"======================================" << endl;

	//
	// Calibrate tree
	//
	vt.Calibrate();

	//
	// Initialize tool
	//
	vpToolRoot->Initialize();

	//
	// Rollback
	//
	for (int tpIdx=vt.TPNum(); tpIdx >= 0; tpIdx--) {
		vt.Update(tpIdx);

		vpToolRoot->Update();
	}


	//
        //
        //
        PutPV(vpToolRoot, dppLog, 0);


	// Output price.
	// Forward value at value date.
	//
	*npv = vpToolRoot->GetResults()["PV"];

    }
    catch (KFailure) {
	throw KFailure("%s: failed.\n", routine);
    }
}




