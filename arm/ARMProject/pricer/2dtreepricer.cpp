/*
 * $Log: 2dtreepricer.cpp,v $
 * Revision 1.10  2001/11/06 17:42:35  vberger
 * retour en arrière sur la gestion des exceptions : on gère
 * désormais le pb par la nouvelle méthode SummitPrice
 * écrite dans la classe base ARM_Pricer
 *
 * Revision 1.9  2001/10/24 08:49:46  vberger
 * ajout/modif de commentaire
 * gestion des exeptions générées lors d'un pricing
 *
 * Revision 1.8  2001/01/29 17:18:50  vberger
 * modif de price : simplification
 * une partie de la gestion des dates est
 * desormais feite dans le reverse
 *
 * Revision 1.7  2000/12/26 13:57:08  vberger
 * modifs en liaison avec la modification de la fct UndrelyingSwapValue
 * de ARM_Reverse
 *
 * Revision 1.6  2000/11/20 15:28:45  vberger
 * appel de fonctions du model pour la remontee dans
 * itsParameters d'output : vol fx, correl ...
 *
 * Revision 1.5  2000/09/25 15:39:59  vberger
 * modif de l'algorithme de pricing : on prend en compte lexistence
 * d un pricer de sj analytique
 *
 * Revision 1.4  2000/02/29 13:17:34  vberger
 * modif due au chgt de nom de la fonction LastResetDate en LastFlowDate
 *
 * Revision 1.3  2000/02/29 08:59:47  vberger
 * correction de la fonction Price
 *
 */ 
#include <iostream>
#include "2dtreepricer.h"
#include "2dlattice.h"
#include "treesolv.h"
#include "constree.h"
#include "reverse.h"


ARM_2DTreePricer::ARM_2DTreePricer(ARM_Security *sec, ARM_Model *mod)
:ARM_TreePricer( sec, mod)

{
	;
}
double ARM_2DTreePricer::Price(void)
{
	/**********************************************************************
	on prepare le modele en fonction des specificites de la security;
	entre autre on r&up&-e la derniere date de reset du produit qui est la date
	correspondant a la proba forward neutre sous laquelle est construif l'arbre 
	************************************************************************/
	try
	{
		itsModel->BeFittedTo(itsSecurity);
	}
	catch(Exception& m)
    {
        m.DebugPrint();
 
        throw Exception(__LINE__, __FILE__, ERR_SWAP_CALC_PB,
                      "Problem in  function BeFittedTo() launched in ARM_2DTreePricer::Price");
    }
    
    
    /**********************************************************************
	on prepare la security en fonction des specificites du modele;
	************************************************************************/
	
	try
	{
		itsSecurity->PrepareToPrice(itsModel->GetStartDate());
	}
	catch(Exception& m)
    {
        m.DebugPrint();
 
        throw Exception(__LINE__, __FILE__, ERR_SWAP_CALC_PB,
                      "Problem in  function ARM_2DTreePricer::Price");
    }	
	/********les dates locales sont toutes en format Julian ***********/
	double valuedate = itsModel->GetStartDate().GetJulian();
	double previousreset = itsSecurity->LastFlowDate();
	double previousexe = itsSecurity->LastExerciseDate() ;//sera egal a -1 s'il n'y a oas d'exercice
	
	/********construction de la slice de r‰tropropagation*****************************/
	int sliceIndex ;
	ARM_Matrix ** slice = new ARM_Matrix*[1] ; 
	double currentdate, oldcurrentdate ;

	//on cherche la première date d'exercice strictement inférieure à la dernière date de reset
	while (previousexe > previousreset)
	{
		previousexe = ((ARM_Reverse*)itsSecurity)->PreviousExerciseDate(previousexe);
	}
	
	/*******	si les toutes les dates d'exercices sont inférieures ou égales à la date d'‰valuation le produit pric‰ ne comporte
	que le sous-jacent qui est ‰valu‰ analytiquement *****/
	if (previousexe <= valuedate)
	{
		double result =  ((ARM_Reverse*)itsSecurity)->UnderlyingSwapValue
			(  valuedate, valuedate,  itsSecurity->LastFlowDate() +1,
			(ARM_TreeSolver *)itsModel,  0 )->Elt(0,0) ;
		((ARM_DFHWSigVarTree *)itsModel)->UpdateOutputParameters();
		itsModel->SetParameters(((ARM_DFHWSigVarTree *)itsModel)->GetOutputParameters());			
		return result;
	}

	
	/*******	‰valuation par arbre : cas g‰n‰ral******/

	/***********localisation de la premiere slice*******************/
	/*la premiˆre slice se localise au niveau de la derniˆre date d'execice*/
	currentdate = previousexe ;



	// introduction du pricing analytique du sous-jacent
	
	/***********construction de la premiere slice*******************/
			
	sliceIndex = ((ARM_TreeSolver *)itsModel)->DateToIndex((ARM_Date)currentdate);
	try
	{
	
	slice[0] =(ARM_Matrix*)(((ARM_Reverse*)itsSecurity)->UnderlyingSwapValue
							(  currentdate,currentdate + 1.0,  itsSecurity->LastFlowDate() +1,
								(ARM_TreeSolver *)itsModel,1 )); 
	}
	catch(Exception& m)
    {
        m.DebugPrint();
	}
	/********r‰tropropagation generale*****************************/
	while ( currentdate >=  valuedate ) 
	{
		
		itsSecurity->Payoff2DTree((ARM_Date)currentdate, (ARM_GenMatrix* *)slice, itsModel);						
		oldcurrentdate = currentdate ; 	
		try
		{
		currentdate = ((ARM_Reverse*)itsSecurity)->PreviousFixingDate(currentdate);
		}
		catch(Exception& m)
		{
        m.DebugPrint();
		}
		if (( currentdate > -1.0 )&& ( currentdate >= valuedate ))
			((ARM_TreeSolver *)itsModel)->PropagateBackward((ARM_Date)currentdate,
				(ARM_Date)oldcurrentdate, (ARM_GenMatrix**)slice , 0 );
				
	}
	
	((ARM_TreeSolver *)itsModel)->PropagateBackward((ARM_Date)valuedate,
				(ARM_Date)oldcurrentdate, (ARM_GenMatrix**)slice, 1);

	//pour r‰cup‰rer les ouputs : vol fx , correls...d'ou le downcast qui n'est pas tres beau...
	((ARM_DFHWSigVarTree *)itsModel)->UpdateOutputParameters();
	itsModel->SetParameters(((ARM_DFHWSigVarTree *)itsModel)->GetOutputParameters());			
	return (*slice)->Elt(0,0);
	

	


 
}
