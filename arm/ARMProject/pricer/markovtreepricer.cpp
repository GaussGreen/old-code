/*
 * $Log: markovtreepricer.cpp,v $
 * Revision 1.23  2003/08/06 10:19:52  jmprie
 * ComputeSliceSize() est remplace par GetTotalStateSize()
 *
 * Revision 1.22  2003/07/31 16:06:47  jmprie
 * desactivation de la retro analytique car pas assez testee
 *
 * Revision 1.21  2003/07/23 08:48:36  jmprie
 * remplacement de l'appel a ComputeRoundIndexLeft du modele par SearchIdx()
 *
 * Revision 1.20  2003/07/21 16:00:24  jmprie
 * modifs pour calcul des probas d'exercice
 *
 * Revision 1.19  2003/07/04 15:02:22  jmprie
 * appel aux probas de la slice + tests de perf en temps (via clock() remis en //)
 * on reste tjrs en esperances purement numeriques sans lissage de la frontiere
 *
 * Revision 1.18  2003/06/30 08:39:33  jmprie
 * nlle correction fuite ds ComputeNodeAnalyticalExpectation
 *
 * Revision 1.17  2003/06/30 07:40:38  jmprie
 * correction de fuites memoire avec lissage de la frontiere d'exercice
 *
 * Revision 1.14  2003/06/24 14:46:48  jmprie
 *  ajout des fcts pour le calcul analytique des esperances via regression
 *  polynomiale
 *
 * Revision 1.13  2003/06/06 15:53:58  jmprie
 * le Price appelle le ComputePrices + ds cette fct boucle sur le nombre
 * cde dates de notification et pas de correction avec le gap sur le prix
 *
 * Revision 1.12  2003/05/30 14:46:26  jmprie
 * ds le Price() merge fait avec les Exet et non plus les Reset
 *
 * Revision 1.11  2003/05/15 09:16:42  ykhlif
 *  Add computePrices for pricing option and underlyings
 *
 * Revision 1.10  2003/04/16 14:36:12  ykhlif
 * mis a jour avec markovtree.h et .cpp
 *
 * Revision 1.9  2003/04/04 12:58:03  jmprie
 * plus de reference au schedule de MarkovTree via les yearterm
 *
 * Revision 1.8  2003/03/27 09:45:46  ykhlif
 * adapter aux dates des resets.
 *
 * Revision 1.7  2003/03/20 16:55:45  ykhlif
 * mis a jour
 *
 * Revision 1.6  2003/03/19 11:35:42  ykhlif
 * adapter aux changements markovtree
 *
 * Revision 1.5  2003/03/14 11:15:13  jmprie
 * oubli du dos2unix
 *
 * Revision 1.4  2003/03/14 11:08:44  jmprie
 * dans price delete du slice
 *
 * Revision 1.3  2003/03/13 10:40:55  jmprie
 * dans Price(), on set la date de la slice + appel a StateEval() avec la slice
 * recastee en ARM_Object *
 *
 * Revision 1.2  2003/03/12 08:22:35  jmprie
 * dans price(), inversion des appels SetModelVariable() et PrepareToPrice()
 *
 * Revision 1.1  2003/03/07 17:40:40  mab
 * Initial revision
 *
 * Revision 1.1  2003/02/10  08:52:12  mab
 * Initial revision
 *
 */ 




/*----------------------------------------------------------------------------*
 
  markovtreepricer.cpp 
 
  Markov Tree Pricer Methods
 
*----------------------------------------------------------------------------*/



#include "markovtreepricer.h"
#include "security.h"
#include "markovtree.h"

#include "interpol.h"
#include "gaussian.h"

#include "time.h"

ARM_MarkovTreePricer::ARM_MarkovTreePricer(ARM_Security* sec, ARM_Model* mod)
                                    :ARM_TreePricer( sec, mod)
{
}

// Recuperation du prix seul de l'option
double ARM_MarkovTreePricer::Price(void)
{

    double optionPrice=UPPER_INFINITE_BOUND;

//    clock_t ts=clock();

    bool isWithExerProba=false;
    ARM_Vector *prices=ComputePrices(isWithExerProba);
    optionPrice=prices->Elt(0);

//    clock_t te=clock();
//    FILE* fOut=fopen("c:\\temp\\dumpTreeKernel.txt","a+");
//    fprintf(fOut,"%lf sec\n",((double)(te-ts))/CLOCKS_PER_SEC);
//    fclose(fOut);

    delete prices;

    return optionPrice;
}

// Recuperation des prix de l'option et du ss-jacent + les probas d'exercice
ARM_Vector* ARM_MarkovTreePricer::ComputePrices(void)
{
    bool isWithExerProba=true;
    return ComputePrices(isWithExerProba);
}

ARM_Vector* ARM_MarkovTreePricer::ComputePrices(bool isWithExerProba)
{
    int nbPrice = itsSecurity->GetPayoffNumber();
    int i,j,k;

    // Initialisation de l'arbre
    itsSecurity->SetModelVariable(itsModel);

    itsModel->BeFittedTo(itsSecurity);

    itsSecurity->PrepareToPrice(itsModel->GetStartDate());


    ARM_FRMMarkovTree* frmModel=(ARM_FRMMarkovTree*)itsModel;
    double asOfDate = frmModel->GetStartDate().GetJulian();

    ARM_Vector *exerSched=frmModel->GetJulianExerciseSched();
	int nbExer = exerSched->GetSize(); // asOfDate + nb notif dates

	ARM_Vector* treeSched =frmModel->GetJulianTreeSched();
    int nbTreeStep=treeSched->GetSize();
    if(treeSched->Elt(nbTreeStep-1) != exerSched->Elt(nbExer-1) ||
        exerSched->Elt(0) != asOfDate || treeSched->Elt(0) != asOfDate)
    {
       throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
        "Inconsistency between notification and tree schedules");
    }
    

    // Augmentation de la taille pour avoir les valeurs d'exercice et
    // de continuation de l'option en vue du calcul analytique
    int analyticExpectationStep=frmModel->GetAnalyticExpectationStep();
    bool isAnalyticExpectation=(analyticExpectationStep>0);
    isAnalyticExpectation=false; // pas active pour le moment

    int nbTreePrice;
    int nbNode = frmModel->GetTotalStateSize(nbTreeStep-1);
    int nbExerProba=0;
    if(isWithExerProba)
        nbExerProba=1;
    ARM_Matrix *payoff,*tmpPayOff;
    if(isAnalyticExpectation)
    {
        if(nbPrice < 2)
        {
           throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
            "Underlying price must be used to go on with analytical expectation");
        }
        nbTreePrice = nbPrice + 2; // ajout pour exercice & continuation
        payoff = new ARM_Matrix(nbNode,nbTreePrice+nbExerProba,0.0);

        // Init pour identification de la zone exercice & continuation
        for(i=0;i<nbNode;++i)
            for(j=nbPrice;j<nbTreePrice;++j)
                payoff->Elt(i,j)=UPPER_INFINITE_BOUND;
    }
    else
    {
        nbTreePrice=nbPrice;
        payoff = new ARM_Matrix(nbNode,nbTreePrice+nbExerProba,0.0);
    }

	int nbFwd = frmModel->GetJulianStartFWDSched()->GetSize()-1; //[0]=asOfDate
	double terminalDate = frmModel->GetJulianEndFWDSched()->Elt(nbFwd);

    ARM_FRMSlice* slice;

    double nextExerDate=terminalDate;
    double exerDate=exerSched->Elt(nbExer-1);
    double prevExerDate;

    int exerTreeIdx=nbTreeStep-1;
    int prevExerTreeIdx,prevAnaTreeIdx;

    for(int exerIdx=nbExer-1; exerIdx>0; --exerIdx)
    {
        // Evaluation de la security entre la date d'exercice
        // courante et la suivante (ou la date terminale)
        slice = frmModel->GenerateSliceConstanteSpacing(exerTreeIdx);
        slice->SetDate(exerDate);
			
		itsSecurity->StateEval(slice,exerDate,nextExerDate,
            (ARM_GenMatrix*)payoff);


        // Retropropagation jusqu'a la date d'exercice precedente
        // (l'indice [0] correspond a l'asOfDate)
        prevExerDate = exerSched->Elt(exerIdx-1);
		prevExerTreeIdx = SearchIdx(treeSched,prevExerDate,K_NEW_DOUBLE_TOL);

        if(isAnalyticExpectation)
        {
            // Lissage via regression sur les valeurs d'exercice et de
            // continuation de l'option puis calcul analytique des esperances
            prevAnaTreeIdx = MAX(prevExerTreeIdx,exerTreeIdx -
                analyticExpectationStep);
            ARM_Matrix* exerProba;
            if(isWithExerProba)
            {
                // Separation des probas et des prix
                nbNode = payoff->GetNumLines();
                tmpPayOff = new ARM_Matrix(nbNode,nbTreePrice);
                exerProba = new ARM_Matrix(nbNode,nbExerProba);
                for(i=0;i<nbNode;++i)
                {
                    for(j=0;j<nbTreePrice;++j)
                        tmpPayOff->Elt(i,j)=payoff->Elt(i,j);
                    for(j=0,k=nbTreePrice;j<nbExerProba;++j,++k)
                        exerProba->Elt(i,j)=payoff->Elt(i,k);
                }

                // Traitement des prix avec lissage analytique
                ComputeSliceAnalyticalExpectation(exerTreeIdx,slice,
                    prevAnaTreeIdx,(ARM_GenMatrix**)&tmpPayOff);

                // Traitement des probas par retropropagation numerique
		        frmModel->PropagateBackward(prevAnaTreeIdx,exerTreeIdx,
                       (ARM_GenMatrix**)&exerProba);

                // Regroupement des probas et des prix
                delete payoff;
                nbNode = tmpPayOff->GetNumLines();
                payoff = new ARM_Matrix(nbNode,nbPrice+nbExerProba);
                for(i=0;i<nbNode;++i)
                {
                    for(j=0;j<nbPrice;++j)
                        payoff->Elt(i,j)=tmpPayOff->Elt(i,j);
                    for(j=0,k=nbPrice;j<nbExerProba;++j,++k)
                        payoff->Elt(i,k)=exerProba->Elt(i,j);
                }
                
                delete tmpPayOff;
                delete exerProba;
            }
            else
            {
                ComputeSliceAnalyticalExpectation(exerTreeIdx,slice,
                    prevAnaTreeIdx,(ARM_GenMatrix**)&payoff);
            }


            // Attention apres le calcul analytique, le nombre de colonnes
            // de "payoff" a diminue de 2 car celles de la valeur d'exercice
            // et de continuation de l'option ont disparu (pas utile de faire
            // la retropropagation numerique qui suit avec)
            if(prevAnaTreeIdx > prevExerTreeIdx)
            {
                // Calcul numerique jusqu'a la precedente notif
		        frmModel->PropagateBackward(prevExerTreeIdx,prevAnaTreeIdx,
                       (ARM_GenMatrix**)&payoff);
            }
        }
        else
        {
             // Calcul uniquement numerique des esperances
		    frmModel->PropagateBackward(prevExerTreeIdx,exerTreeIdx,
                   (ARM_GenMatrix**)&payoff);
        }

        nextExerDate=exerDate;
        exerDate=prevExerDate;
        exerTreeIdx=prevExerTreeIdx;

        if(exerIdx>1 && (isAnalyticExpectation || isWithExerProba))
        {
            if(isWithExerProba)
                nbExerProba += 1;

            // Retour a la taille initiale pour le prochain calcul
            // analytique : les prix + 2 (exerice & continuation) +
            // les probas (si demandees)
            nbNode = payoff->GetNumLines();
            tmpPayOff = new ARM_Matrix(nbNode,nbTreePrice+nbExerProba);
            for(i=0;i<nbNode;++i)
            {
                // Recuparation des prix
                for(j=0;j<nbPrice;++j)
                    tmpPayOff->Elt(i,j)=payoff->Elt(i,j);

                // Indentification de la zone exercice & continuation
                for(j=nbPrice;j<nbTreePrice;++j)
                    tmpPayOff->Elt(i,j)=UPPER_INFINITE_BOUND;

                // Recuparation des probas d'exercice en decalant d'un cran
                // vers la drte, zone pour la proba a l'exercice precedent
                int l=nbPrice;
                for(j=0,k=nbTreePrice+1;j<nbExerProba-1;++j,++k,++l)
                    tmpPayOff->Elt(i,k)=payoff->Elt(i,l);
            }
            delete payoff;
            payoff=tmpPayOff;
        }

        delete slice;
    }


    // Calcul du ss jacent eventuel entre asOfDate et 1ere date de notif
    slice = frmModel->GenerateSliceConstanteSpacing(0);
    slice->SetDate(asOfDate);

    itsSecurity->StateEval(slice,asOfDate,exerDate,(ARM_GenMatrix*)payoff);

	double terminalZc = frmModel->GetZeroCurve()->
        DiscountPrice((terminalDate-asOfDate)/K_YEAR_LEN);

    // Actualisation en 0 des prix capitalises en date terminale
    ARM_Vector* result = new ARM_Vector(nbPrice+nbExerProba+1);
    for(i=0;i<nbPrice;i++)
        result->Elt(i) = terminalZc * payoff->Elt(0,i);

    if(isWithExerProba)
    {
        // Recuperation et traitement final des probas d'exercice
        double lastProba=0.0;
        result->Elt(nbPrice)=payoff->Elt(0,nbPrice+nbExerProba-1);
        for(i=nbPrice;i<nbPrice+nbExerProba;++i)
        {
            result->Elt(i+1) = payoff->Elt(0,i)-lastProba;
            lastProba = payoff->Elt(0,i);
        }
    }

    // Liberation memoire
    delete slice;
    delete payoff;

	ARM_Vector tmpResult(*result);
	tmpResult *= itsSecurity->GetPorS();
	*result = tmpResult;

    return(result);
}

/*--------------------------------------------------------------------------*/
/* A partir d'un noeud de la slice courante, calcul de la regression        */
/* polynomiale d'ordre n passant par une selection de noeuds de la slice    */
/* suivante                                                                 */
/*--------------------------------------------------------------------------*/

ARM_Vector* ARM_MarkovTreePricer::NodePolynomialInterpolation(int nodeIdx,
                                        bool isTrunc,int nbPt,
                                        ARM_Vector* state,
                                        int payOffIdx,ARM_Matrix* payOff)
{
    // Les points de la slice courante servant au calcul du polynome
    // sont deduits du rang du noeud sur la slice precedente

    int nbNode=state->GetSize();
    if(nbNode<2)
    {
       throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
            "Not enought node to compute polynomial interpolation");
    }

    int mink;
    if(isTrunc)
    {
        // Troncature de l'arbre
        if(nodeIdx==0)
        {
            mink=0;
        }
        else if(nodeIdx==nbNode-1)
        {
            mink=nbNode-nbPt;
        }
        else
        {
            mink=nodeIdx-nbPt/2;
        }
    }
    else
    {
        mink=nodeIdx -nbPt/2 + 1;
    }
    mink=MAX(mink,0);
    int maxk=MIN(mink+nbPt-1,nbNode-1);
    nbPt=maxk-mink+1;

    ARM_Vector x(nbPt),fx(nbPt);
    int k,kk;
    for(k=0;k<nbPt;k++)
    {
        kk=mink+k;
        x.Elt(k)=state->Elt(kk);
        fx.Elt(k)=payOff->Elt(kk,payOffIdx);
    }

    return PolynomialInterpolation(&x,&fx);
}

ARM_Vector* ARM_MarkovTreePricer::NodePolynomialInterpolation(int nodeIdx,
                                        double stateMin,double stateMax,
                                        int nbPt,ARM_Vector* state,
                                        int payOffIdx,ARM_Matrix* payOff)
{
    // Les points de la slice courante servant au calcul du polynome
    // sont ceux dont l'etat est ds l'intervalle donne en input

    int nbNode=state->GetSize();
    if(nbNode<2)
    {
       throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
            "Not enought node to compute polynomial interpolation");
    }

    int k;
    for(k=nodeIdx;k>=0;k--)
    {
        if(state->Elt(k) >= stateMax - K_NEW_DOUBLE_TOL)
            break;
    }
    int mink=MAX(k,0);
    for(k=nodeIdx;k<nbNode;k++)
    {
        if(state->Elt(k) <= stateMin + K_NEW_DOUBLE_TOL)
            break;
    }
    int maxk=MIN(k,nbNode-1);
    int nbPtRange=maxk-mink+1;

    int step=1;
    if(nbPtRange>nbPt)
// FIXMEFRED: mig.vc8 (22/05/2007 15:49:47): (int)(floor(int/int))+1 = int/int+1
		step=nbPtRange/nbPt+1;
    else
        nbPt=MAX(nbPtRange,2);

    ARM_Vector idx(nbPt);
    idx.Elt(0)=nodeIdx;
    mink=nodeIdx;
    maxk=nodeIdx;
    bool isMin=(nodeIdx==0);
    bool isMax=(nodeIdx==nbNode-1);
    for(k=1;k<nbPt;k++)
    {
        if(isMin)
        {
            maxk += step;
            if(maxk<nbNode-1)
                idx.Elt(k)=maxk;
            else
            {
                idx.Elt(k)=nbNode-1;
                isMax=true;
            }
        }
        else if(isMax)
        {
            mink -= step;
            if(mink>0)
                idx.Elt(k)=mink;
            else
            {
                idx.Elt(k)=0;
                isMin=true;
            }
        }
        else
        {
            if(k%2 == 0)
            {
                mink -= step;
                if(mink>0)
                    idx.Elt(k)=mink;
                else
                {
                    idx.Elt(k)=0;
                    isMin=true;
                }
            }
            else
            {
                maxk += step;
                if(maxk<nbNode-1)
                    idx.Elt(k)=maxk;
                else
                {
                    idx.Elt(k)=nbNode-1;
                    isMax=true;
                }
            }
        }
    }
    ARM_Vector sIdx=idx.Sort();
    ARM_Vector x(nbPt),fx(nbPt);
    int kk;
    for(k=0;k<nbPt;k++)
    {
        kk=sIdx.Elt(k);
        x.Elt(k)=state->Elt(kk);
        fx.Elt(k)=payOff->Elt(kk,payOffIdx);
    }

    return PolynomialInterpolation(&x,&fx);
}


/*--------------------------------------------------------------------------*/
/*  Calcul analytique de l'esperance conditionnelle en un noeud             */
/*--------------------------------------------------------------------------*/
double ARM_MarkovTreePricer::ComputeNodeAnalyticalExpectation(double mean,
                                    double stdDev,ARM_Vector* polyInf,
                                    double xstar,ARM_Vector* polySup)
{
    // Le polynome "inf" est integre de -infini a xstart et le polynome "sup"
    // de xstar a +infini

    double value=0;

    int nbInf=0,nbSup=0;
    double *momentInf=NULL,*momentSup=NULL;
    double lastMoment;
    if(polyInf)
    {
        nbInf=polyInf->GetSize();
        momentInf = new double [nbInf];
        lastMoment=cumulNormalMoment(xstar,nbInf-1,mean,stdDev,momentInf);
    }

    if(polySup)
    {
        nbSup=polySup->GetSize();
        momentSup = new double[nbSup];
        lastMoment=cumulNormalMoment(-xstar,nbSup-1,-mean,stdDev,momentSup);
    }

    int i;
    double sign=1.0;
    for(i=0;i<MIN(nbInf,nbSup);i++)
    {
        value += polyInf->Elt(i)*momentInf[i] +
            (sign>0 ? polySup->Elt(i) : -polySup->Elt(i))*momentSup[i];
        sign = -sign;
    }

    if(nbInf<nbSup)
    {
        for(i=nbInf;i<nbSup;i++)
        {
            value += (sign>0 ? polySup->Elt(i) : -polySup->Elt(i))*momentSup[i];
            sign = -sign;
        }
    }
    else if(nbInf>nbSup)
    {
        for(i=nbSup;i<nbInf;i++)
        {
            value += polyInf->Elt(i)*momentInf[i];
        }
    }

    if(momentInf)
        delete momentInf;
    if(momentSup)
        delete momentSup;

    return value;
}

double ARM_MarkovTreePricer::ComputeNodeAnalyticalExpectation(double mean,
                                    double stdDev,ARM_Vector* poly)
{
    // Le polynome est integre de -infini a +infini

    double value=0;
    int nb=poly->GetSize();

    double* moment = new double[nb];
    double lastMoment=normalMoment(nb-1,mean,stdDev,moment);
    for(int i=0;i<nb;i++)
    {
        value += poly->Elt(i)*moment[i];
    }

    delete moment;

    return value;
}

double ARM_MarkovTreePricer::ComputeNodeNumericalExpectation(int nodeIdx,
                                bool isTrunc,ARM_Vector* proba,
                                int payOffIdx,ARM_Matrix* payOff)
{
    double value;
    if(isTrunc)
    {
        if(nodeIdx==0 || nodeIdx==payOff->GetNumLines()-1)
        {
            value=payOff->Elt(nodeIdx,payOffIdx);
        }
        else
        {
            if(proba->GetSize()==T_TRINOMIAL_CONSTANT_TREE)
            {
                value = proba->Elt(T_UPDOWN)*
                    ( payOff->Elt(nodeIdx+1,payOffIdx) +
                    payOff->Elt(nodeIdx-1,payOffIdx) ) +
                    proba->Elt(T_MID)*payOff->Elt(nodeIdx,payOffIdx);
            }
            else
            {
                value=proba->Elt(T_DOWN)*payOff->Elt(nodeIdx+1,payOffIdx) +
                        proba->Elt(T_MID)*payOff->Elt(nodeIdx,payOffIdx) +
                        proba->Elt(T_UP)*payOff->Elt(nodeIdx-1,payOffIdx);
            }
        }
    }
    else
    {
        if(proba->GetSize()==T_TRINOMIAL_CONSTANT_TREE)
        {
            value = proba->Elt(T_UPDOWN)*
                    ( payOff->Elt(nodeIdx+2,payOffIdx) +
                    payOff->Elt(nodeIdx,payOffIdx) ) +
                    proba->Elt(T_MID)*payOff->Elt(nodeIdx+1,payOffIdx);
        }
        else
        {
            value=proba->Elt(T_DOWN)*payOff->Elt(nodeIdx+2,payOffIdx) +
                    proba->Elt(T_MID)*payOff->Elt(nodeIdx+1,payOffIdx) +
                    proba->Elt(T_UP)*payOff->Elt(nodeIdx,payOffIdx);
        }
    }

    return value;
}


/*--------------------------------------------------------------------------*/
/*  Calcul des esperances conditionnelles a la slice indexee par            */
/*  "prevTreeIdx". Les fcts d'exercice et de conservation de l'option sont  */
/*  discretisees a la slice de notification indexee par "treeIdx" et        */
/*  pointee par "treeState"                                                 */
/*--------------------------------------------------------------------------*/
void ARM_MarkovTreePricer::ComputeSliceAnalyticalExpectation(int treeIdx,
                                              ARM_Object* treeState,
                                              int prevTreeIdx,
                                              ARM_GenMatrix** treePrice)
{
    int OPT_IDX=0;
    int UNDL_IDX=1;
    int UNEXER_IDX=2;
    int EXER_IDX=3;

    int MAX_ITER_EXER=20;
    double TOL_EXER_COEF=10000.0;

    // Nbre de pts utilises pour les regressions polynomiales
    // si isSigmaArea=false : on utilise NBPT_OUT_INTERP pts consecutifs
    // si isSigmaArea=true : on utilise ts les pts ds + ou - coef*stdDev
    // sans depasser NBPT_OUT_INTERP
    bool isSigmaArea=false;
    int NBPT_OUT_INTERP=3;
    int NBPT_IN_INTERP=4;

    // + ou - EXER_AREA_COEF*stdDev pour delimiter la zone d'exercice
    double EXER_AREA_COEF=5.0;


    ARM_FRMSlice *slice = (ARM_FRMSlice *)treeState;
    ARM_Vector* sliceStates = slice->BuiltSliceStateVariables();
    int nbNode = sliceStates->GetSize();

    ARM_Matrix* payOff =(ARM_Matrix *)(*treePrice);
    int nbValue=payOff->GetNumLines();
    int nbPrice=payOff->GetNumCols();
    if(nbPrice < 4 || nbValue != nbNode)
    {
       throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
            "Inconsistency in analytical slice expectation computation");
    }


    // Localisation de la frontiere
    bool isUnexerFirst=(payOff->Elt(0,UNEXER_IDX) > payOff->Elt(0,EXER_IDX));
    int i,j,i0;
    for(i0=1;i0<nbValue;i0++)
    {
        if((payOff->Elt(i0,UNEXER_IDX) > payOff->Elt(i0,EXER_IDX)) !=
            isUnexerFirst)
            break;
    }

    double y;
    ARM_Vector* areaPoly;
    int areaIdx;

    ARM_FRMMarkovTree* frmModel = (ARM_FRMMarkovTree*)itsModel;
    int nbPrevPrice=nbPrice-2;
    ARM_Matrix* prevPayOff;

    bool isAnaEqNum=(prevTreeIdx==treeIdx-1 && NBPT_OUT_INTERP==3 && !isSigmaArea);

    if(i0==nbValue && isAnaEqNum)
    {
        // Le calcul est alors equivalent a la methode numerique
        // donc on fait une retropropagation classique + efficace
        // Resize pour ne pas retropropager les fcts d'exercice ni
        // de continuation de l'option
        
        prevPayOff = new ARM_Matrix(nbNode,nbPrevPrice);

        for(i=0;i<nbNode;i++)
            for(j=0;j<nbPrevPrice;j++)
                prevPayOff->Elt(i,j)=payOff->Elt(i,j);

        delete payOff;

        *treePrice=(ARM_GenMatrix *)prevPayOff;
        frmModel->PropagateOneStepBackward(prevTreeIdx,treePrice);

        return;
    }


    // Creation de la slice anterieure et de la matrice de prix
    ARM_FRMSlice* prevSlice=frmModel->GenerateSliceConstanteSpacing(prevTreeIdx);
    ARM_Vector* prevSliceStates = prevSlice->BuiltSliceStateVariables();
    int nbPrevNode = prevSliceStates->GetSize();
    bool isTrunc=(nbPrevNode==nbNode);

    prevPayOff = new ARM_Matrix(nbPrevNode,nbPrevPrice);

    double stdDev=sqrt(frmModel->GetTotalDELTASlices(treeIdx)-
        frmModel->GetTotalDELTASlices(prevTreeIdx));
    double dy=EXER_AREA_COEF*stdDev;

    if(i0==nbValue)
    {
        // Exercice ou conservation de l'option sur tte la slice
        areaIdx=(isUnexerFirst ? UNEXER_IDX : EXER_IDX);
        for(i=0;i<nbPrevNode;i++)
        {
            y=prevSliceStates->Elt(i);

            if(isSigmaArea)
                areaPoly=NodePolynomialInterpolation(i,y-dy,y+dy,
                    NBPT_OUT_INTERP,sliceStates,areaIdx,payOff);
            else
                areaPoly=NodePolynomialInterpolation(i,isTrunc,
                    NBPT_OUT_INTERP,sliceStates,areaIdx,payOff);

            // Une seule fct prise en compte de -infini -> +infini
            // pour le calcul de l'option
            prevPayOff->Elt(i,OPT_IDX)=ComputeNodeAnalyticalExpectation(y,
                stdDev,areaPoly);
            delete areaPoly;

            // Esperance du sous-jacent
            if(isSigmaArea)
                areaPoly=NodePolynomialInterpolation(i,y-dy,y+dy,
                    NBPT_OUT_INTERP,sliceStates,UNDL_IDX,payOff);
            else
                areaPoly=NodePolynomialInterpolation(i,isTrunc,
                    NBPT_OUT_INTERP,sliceStates,UNDL_IDX,payOff);

            prevPayOff->Elt(i,UNDL_IDX)=ComputeNodeAnalyticalExpectation(y,
                stdDev,areaPoly);
            delete areaPoly;

        }// for node prevSlice
    }
    else
    {
        // x0Prev > x0 car on a trace les etats par valeurs decroissantes
        double x0Prev=sliceStates->Elt(i0-1);
        double x0=sliceStates->Elt(i0);

        // Interpolation des fcts d'exercice et de conservation de l'option
        // autour de la frontiere reperee
        int nbPt=MIN(nbNode,NBPT_IN_INTERP);
        if(nbPt<2)
        {
           throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
                "Not enought points to regress exercise frontier");
         }
        int imin=i0-nbPt/2;
        imin=MAX(imin,0);
        int imax=imin+nbPt-1;
        imax=MIN(imax,nbNode-1);
        nbPt=imax-imin+1;
        ARM_Vector x(nbPt);
        ARM_Vector unexer(nbPt),exer(nbPt);
        int ii=MIN(imin+1,i0);
        double y0up=sliceStates->Elt(ii);
        double y0down=sliceStates->Elt(MAX(imin+nbPt-2,ii));
        for(i=0;i<nbPt;i++)
        {
            ii=imin+i;
            x.Elt(i)=sliceStates->Elt(ii);
            unexer.Elt(i)=payOff->Elt(ii,UNEXER_IDX);
            exer.Elt(i)=payOff->Elt(ii,EXER_IDX);
        }
        ARM_Vector* unexerPoly=PolynomialInterpolation(&x,&unexer);
        ARM_Vector* exerPoly=PolynomialInterpolation(&x,&exer);
        ARM_Vector* downPoly;
        ARM_Vector* upPoly;
        if(isUnexerFirst)
        {
            upPoly=unexerPoly;
            downPoly=exerPoly;
        }
        else
        {
            upPoly=exerPoly;
            downPoly=unexerPoly;
        }

        // Extimation du seuil d'exercice : exer(x*)=unexer(x*)
        double dx=x0Prev-x0;
        double tol=fabs(dx)/TOL_EXER_COEF;
        int nbIter=0;
        double u,fxstar,dfxstar;
        double xstar=0.5*(x0+x0Prev);
        while(fabs(dx) > tol && nbIter<MAX_ITER_EXER)
        {
            u=1.0;
            fxstar=unexerPoly->Elt(0)-exerPoly->Elt(0);
            dfxstar=unexerPoly->Elt(1)-exerPoly->Elt(1);
            for(i=2;i<nbPt;i++)
            {
                u *= xstar;
                fxstar += (unexerPoly->Elt(i-1)-exerPoly->Elt(i-1))*u;
                dfxstar += i*(unexerPoly->Elt(i)-exerPoly->Elt(i))*u;
            }
            u *= xstar;
            fxstar += (unexerPoly->Elt(nbPt-1)-exerPoly->Elt(nbPt-1))*u;
            dx=-fxstar/dfxstar;
            if(dx>0.0)
                xstar=MIN(xstar+dx,x0Prev);
            else
                xstar=MAX(xstar+dx,x0);
            nbIter++;
        }

/***
        FILE* fOut=fopen("c:\\temp\\dump.txt","a+");
        fprintf(fOut,"x0 xstar x0Prev nbIter\t%lf\t%lf\t%lf\t%d\n",x0,xstar,x0Prev,nbIter);
        fprintf(fOut,"y0down-dy y0down y0up y0up+dy\t%lf\t%lf\t%lf\t%lf\n",y0down-dy,y0down,y0up,y0up+dy);
        fprintf(fOut,"Unexer\tExer\n");
        for(int ij=0;ij<nbPt;ij++)
        {
            fprintf(fOut,"%lf\t%lf\n",unexerPoly->Elt(ij),exerPoly->Elt(ij));
        }
        fprintf(fOut,"\n\n");
        fclose(fOut);
***/

        if(nbIter==MAX_ITER_EXER)
        {
           throw Exception(__LINE__, __FILE__, ERR_MAX_ITER_NUM_EXD,
                "Max loop reached in searching exercise frontier value");
        }

        // Calcul analytiques des esperances conditionnelles
        double value;
        ARM_Vector* proba=prevSlice->GetConstantProbabilities();
        for(i=0;i<nbPrevNode;i++)
        {
            // On trace les etats par valeurs decroissantes
            y=prevSliceStates->Elt(i);

            // Esperance du sous-jacent
            if(isAnaEqNum)
            {
                // Calcul numerique identique au calcul analytique
                prevPayOff->Elt(i,UNDL_IDX)=ComputeNodeNumericalExpectation(i,
                    isTrunc,proba,UNDL_IDX,payOff);
            }
            else
            {
                // Calcul purement analytique
                if(isSigmaArea)
                    areaPoly=NodePolynomialInterpolation(i,y-dy,y+dy,
                        NBPT_OUT_INTERP,sliceStates,UNDL_IDX,payOff);
                else
                    areaPoly=NodePolynomialInterpolation(i,isTrunc,
                        NBPT_OUT_INTERP,sliceStates,UNDL_IDX,payOff);

                prevPayOff->Elt(i,UNDL_IDX)=ComputeNodeAnalyticalExpectation(y,
                                                stdDev,areaPoly);
                delete areaPoly;
            }

            // Esperance de l'option
            if(isAnaEqNum && (y - dy > y0up || y + dy < y0down))
            {
                // Calcul numerique identique au calcul analytique
                value=ComputeNodeNumericalExpectation(i,isTrunc,proba,
                            OPT_IDX,payOff);
            }
            else if(y > y0up + K_NEW_DOUBLE_TOL)
            {
                // Noeud menant a une zone de conservation de l'option
                areaIdx=(isUnexerFirst ? UNEXER_IDX : EXER_IDX);

                if(isSigmaArea)
                    areaPoly=NodePolynomialInterpolation(i,y-dy,y+dy,
                        NBPT_OUT_INTERP,sliceStates,areaIdx,payOff);
                else
                    areaPoly=NodePolynomialInterpolation(i,isTrunc,
                        NBPT_OUT_INTERP,sliceStates,areaIdx,payOff);

                if(y - dy > y0up)
                {
                    // Une seule fct prise en compte de -infini -> +infini
                    value=ComputeNodeAnalyticalExpectation(y,stdDev,areaPoly);
                }
                else
                {
                    // Exercice et conservation pris en compte :
                    // -infini -> xstar : downPoly
                    // xstar -> y0up : upPoly
                    // y0up -> +infini : areaPoly (regression locale)
                    value=ComputeNodeAnalyticalExpectation(y,stdDev,
                        downPoly,xstar,upPoly);
                    value -= ComputeNodeAnalyticalExpectation(y,stdDev,
                        NULL,y0up,upPoly);
                    value += ComputeNodeAnalyticalExpectation(y,stdDev,
                        NULL,y0up,areaPoly);
                }
                delete areaPoly;

            }// if node au dessus de la frontiere a la slice

            else if(y < y0down - K_NEW_DOUBLE_TOL)
            {
                // Noeud menant a une zone d'exercice de l'option
                areaIdx=(isUnexerFirst ? EXER_IDX : UNEXER_IDX);

                if(isSigmaArea)
                    areaPoly=NodePolynomialInterpolation(i,y-dy,y+dy,
                        NBPT_OUT_INTERP,sliceStates,areaIdx,payOff);
                else
                    areaPoly=NodePolynomialInterpolation(i,isTrunc,
                        NBPT_OUT_INTERP,sliceStates,areaIdx,payOff);

                if(y + dy < y0down)
                {
                    // Une seule fct prise en compte de -infini -> +infini
                    value=ComputeNodeAnalyticalExpectation(y,stdDev,areaPoly);
                }
                else
                {
                    // Exercice et conservation pris en compte :
                    // -infini -> y0down : areaPoly (regression locale)
                    // y0down -> xstar : downPoly
                    // xstar -> +infini : upPoly
                    value=ComputeNodeAnalyticalExpectation(y,stdDev,
                        downPoly,xstar,upPoly);
                    value -= ComputeNodeAnalyticalExpectation(y,stdDev,
                        downPoly,y0down,NULL);
                    value += ComputeNodeAnalyticalExpectation(y,stdDev,
                        areaPoly,y0down,NULL);
                }
                delete areaPoly;

            }// if node en dessous de la frontiere a la slice

            else
            {
                // Le noeud mene dans la zone d'exercice
                // -infini -> xstar : downPoly
                // xstar -> +infini : upPoly
                value=ComputeNodeAnalyticalExpectation(y,stdDev,
                    downPoly,xstar,upPoly);
            }

            prevPayOff->Elt(i,OPT_IDX)=value;

        } // for node prevSlice

        delete downPoly;
        delete upPoly;

    }// if frontiere d'exerice

    delete sliceStates;
    delete prevSliceStates;
    delete prevSlice;
    delete payOff;

    *treePrice=(ARM_GenMatrix *)prevPayOff;
}

/*---------------------------------------------------------------------------*/
/*---- End Of File ----*/
