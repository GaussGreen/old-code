/*
 * $log$
 */


#include "mpricer.h"
#include "frmtree.h"
#include "miter.h"

#include "frmutils.h"


/**********************************************/
/*         Version 1D                         */
/**********************************************/





ARM_MTreePricer1::ARM_MTreePricer1(ARM_Security *sec, ARM_Model *mod)
    :ARM_TreePricer(sec,mod)
{
    Init();
}    




// 0 : exercise value
// 1 ... NbHold : payoffs
// N}bHold+1 .. end : other


double ARM_MTreePricer1::Price(void)
{

    //typedef MArray<ARM_Vector*, 1> TSlice;
    typedef M1Array<ARM_Vector*> TSlice;

    
    ARM_FrmTree1 *frmtree  = dynamic_cast<ARM_FrmTree1 *> (itsModel);   
    ARM_Epochs *epochs = frmtree->GetEpochs();

    ARM_Date calcDate = frmtree->GetStartDate();


    // If we are pricing a european option at exercise date,
    // return the exercise value

    if ( (itsSecurity->GetExercisingDates()->GetSize() == 1) &&
         (itsSecurity->GetExercisingDates()->Elt(0) == calcDate.GetJulian()))
    {
        return itsSecurity->ComputeNPVForGrid(calcDate, calcDate, itsModel);
    }


    itsSecurity->PrepareToPrice(calcDate);					

    BeFittedTo(itsSecurity);
    frmtree->BeFittedTo(itsSecurity, frmtree->GetZeroCurve()->GetCurrencyUnit()->GetCcyName(), frmtree->GetMfine());						

    // YAN 02/2002
	itsSecurity->ComputeUnderlyingPrice(calcDate, frmtree);

    bool CurveMemory = TREE_MANAGEMENT_MODE::IsPrecomputeMode(frmtree->GetMemoryMode());
    int currentIndex = frmtree->GetEpochs()->get_i_exp();
    int valueIndex   = 0;       
    int payoffMode;


    TSlice *slice;
    slice = frmtree->GetLattice()->CreateSlice(currentIndex, 1 + itsNbHoldValue + itsNbExtraInfo);
    frmtree->SetTimeStep(currentIndex);


    
    frmtree->GetIter()->InitLayer(currentIndex);
    frmtree->GetIter()->InitLayer(currentIndex-1);





    if (frmtree->IsPathDepCor())
    {
        payoffMode = 0;
    }
    else
    {
        payoffMode = PAYOFF_MODE::Payoff;
    }

    //DEBUG
	payoffMode  |= PAYOFF_MODE::Probas;						


    // Could be out in the big loop if the memory can change in the tree
    if (!CurveMemory)
    {
        payoffMode  |= PAYOFF_MODE::Probas;						
    }


    //DEBUG
    #ifdef DEBUGMODE
    payoffMode  |= PAYOFF_MODE::Payoff;						
    frmtree->GetLattice()->Print("lattice.txt");
    #endif

    // Last slice : ambiant payoffs must be computed to initialise exercise boundary

    itsSecurity->PayoffSlice( currentIndex, slice, frmtree,payoffMode | PAYOFF_MODE::Payoff);						

    frmtree->InitPropagation(slice, currentIndex);

    frmtree->PropagateOneStepBackward(currentIndex, slice, 0);

    currentIndex--;

    

	while ( currentIndex >  valueIndex ) 
	{

        itsSecurity->PayoffSlice( currentIndex, slice, frmtree,payoffMode);						


    	frmtree->PropagateOneStepBackward(currentIndex, slice, 0);



        currentIndex--;
				
	}




    double price = 	slice->rbegin()->Elt(1);

    slice->Cleanit();
	frmtree->EndPricing();
	
	cleanit(slice);
	
				
	return price;
}






void ARM_MTreePricer1::BeFittedTo(ARM_Security *sec)
{
    if (! sec->IsSimplePayoff())
    {
    	ARM_SecurityExt* ext = (ARM_SecurityExt*) sec->GetExt();
    
        itsSimplePayoff    = false;
        itsNbHoldValue     = ext->GetNbHoldValue();
        itsNbExtraInfo     = ext->GetNbExtraInfo();
    }
}




/**********************************************/
/*         Version Template                   */
/**********************************************/


/*
template <int N>
ARM_MTreePricer<N>::ARM_MTreePricer(ARM_Security *sec, ARM_Model *mod)
    :ARM_TreePricer(sec,mod)
{
    Init();
}    




// 0 : exercise value
// 1 ... NbHold : payoffs
// N}bHold+1 .. end : other

template <int N>
double ARM_MTreePricer<N>::Price(void)
{

    //typedef MArray<ARM_Vector*, 1> TSlice;
    typedef M1Array<ARM_Vector*> TSlice;

    //ARM_TreeSolver *treeModel= dynamic_cast<ARM_TreeSolver *> (itsModel);
    ARM_FrmTree<N> *frmtree  = dynamic_cast<ARM_FrmTree<N> *> (itsModel);   
    ARM_Epochs *epochs = frmtree->GetEpochs();

    ARM_Date calcDate = frmtree->GetStartDate();

    itsSecurity->PrepareToPrice(calcDate);					

    BeFittedTo(itsSecurity);
    frmtree->BeFittedTo(itsSecurity, frmtree->GetZeroCurve()->GetCurrencyUnit()->GetCcyName(), frmtree->GetMfine());						



    bool CurveMemory = TREE_MANAGEMENT_MODE::IsPrecomputeMode(frmtree->GetMemoryMode());
    int currentIndex = frmtree->GetEpochs()->get_i_exp();
    int valueIndex   = 0;       
    int payoffMode;


    TSlice *slice;
    slice = frmtree->GetLattice()->CreateSlice(currentIndex, 1 + itsNbHoldValue + itsNbExtraInfo);
    frmtree->SetTimeStep(currentIndex);


    
    frmtree->GetIter()->InitLayer(currentIndex);
    frmtree->GetIter()->InitLayer(currentIndex-1);

    //FOUT
    ofstream fout("payoff_slice.txt", ios::out);
    fout.close();




    if (frmtree->IsPathDepCor())
    {
        payoffMode = 0;
    }
    else
    {
        payoffMode = PAYOFF_MODE::Payoff;
    }

    


    // Could be out in the big loop if the memory can change in the tree
    if (!CurveMemory)
    {
        payoffMode  |= PAYOFF_MODE::Probas;						
    }


    //DEBUG
    payoffMode  |= PAYOFF_MODE::Payoff;						

    // Last slice : ambiant payoffs must be computed to initialise exercise boundary

    itsSecurity->PayoffSlice( currentIndex, slice, frmtree,payoffMode | PAYOFF_MODE::Payoff);						

    frmtree->InitPropagation(slice, currentIndex);

    frmtree->PropagateOneStepBackward(currentIndex, slice, 0);

    currentIndex--;

    

	while ( currentIndex >  valueIndex ) 
	{

        itsSecurity->PayoffSlice( currentIndex, slice, frmtree,payoffMode);						


    	frmtree->PropagateOneStepBackward(currentIndex, slice, 0);



        currentIndex--;
				
	}




    double price = 	slice->rbegin()->Elt(1);

    slice->Cleanit();
	frmtree->EndPricing();
	
	cleanit(slice);
	
				
	return price;
}





template <int N>
void ARM_MTreePricer<N>::BeFittedTo(ARM_Security *sec)
{
    if (! sec->IsSimplePayoff())
    {
        itsSimplePayoff    = false;
        itsNbHoldValue     = sec->GetExt()->GetNbHoldValue();
        itsNbExtraInfo     = sec->GetExt()->GetNbExtraInfo();
    }
}




*/
