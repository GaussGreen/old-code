
#include "DKMaille.h"
#include "DKMaille2D.h"
#include "Tree.h"
#include "utils.h"

#include <math.h>

Node::Node()
{
    P=NULL;
    ZC = NULL;
    SWOi = 0.;
}

Node::~Node()
{
    if( P )
        delete []P;
    if( ZC )
        delete []ZC;
}


void Node::Init()
{
    P = new double[4];
    P[0] = 0.;
    j = 0;
    ArrowDebreu = 0.;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////


Slice::Slice()
{
    nodes = NULL;
    isNotice = false;
    Multiplicateur = 1.0;
}

Slice::~Slice()
{
    if (nodes)
        delete []nodes;
}


void Slice::Init(int i)
{
    nodes=new Node[i];
    for( int j = 0; j < i; j++ )
        nodes[j].Init();

    nbNodes=i;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////


Tree::Tree()
{
    slices=NULL;
}



Tree::~Tree()
{

    T.clear();
    dT.clear();
    Sigma.clear();
    Beta.clear();

    if(slices)
        delete []slices;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Tree::Clean()
{
    T.clear();
    dT.clear();
    Sigma.clear();
    Beta.clear();

    if(slices)
        delete []slices;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void Tree::Init( double StartPoint,
            int dNbSteps,
            DKMaille<double> MaturityDates,
            DKMaille<double> dSigma,
            DKMaille<double> dMeanReversion)
{
    int i,j, l, nodeMin, nbNodes, nodeMin1, nbNodes1;
    Node* nodes;
    Node* nodes1;
    double k, X, dX, dX1, Mu, Eta;

    X0=StartPoint;
    TreeLimite = 1.e-10;

    // On suppose que les T correspondent au dates des slices donc T[0]=0;
    T = MaturityDates;
    Sigma = dSigma;
    Beta = dMeanReversion;

    nbSteps=dNbSteps;
    nbSlices=nbSteps+1;


    if( Sigma.entries()!=T.entries() || Beta.entries()!=T.entries() )
        throw("Tree constructor : T, Sigma et Beta doivent avoir le meme nombre d'elements");

    slices=new Slice[nbSlices];

    dT.resize(T.entries()-1);
    for(i=0;i<T.entries()-1;i++)
        dT.at(i)=T.at(i+1)-T.at(i);

    // on set la taille des ecarts entre les differents noeuds sur une slice
    slices[0].dX=0.;
    for(i=1;i<nbSlices;i++)
    {
        slices[i].dX=Sigma.at(i-1)*sqrt(3.*dT.at(i-1));
    }


    // Creation du noeud de depart
    slices[0].nbNodes=1;
    slices[0].nodeMax=0;
    slices[0].nodeMin=0;
    slices[0].nodes=new Node[1];
    slices[0].nodes[0].Init();

    for( i = 0; i < nbSlices - 1; i++ )
    {
        nodeMin=slices[i].nodeMin;
        nbNodes=slices[i].nbNodes;
        nodes=slices[i].nodes;
        dX=slices[i].dX;
        dX1=slices[i+1].dX;

        for( j = 0; j < nbNodes; j++)
        {
            X=X0+nodes[j].j*dX;
            Mu=-Beta.at(i)*X;

            k=(X+Mu*dT.at(i)-X0)/dX1;

            nodes[j].k=(int) k;
            if( k-nodes[j].k>0.5 ) nodes[j].k++;
            if( nodes[j].k-k>0.5 ) nodes[j].k--;

            Eta=Mu*dT.at(i)+nodes[j].j*dX-nodes[j].k*dX1;

            nodes[j].P[1]=0.5*pow(Sigma.at(i)/dX1,2)*dT.at(i)+0.5*pow(Eta/dX1,2)+0.5*Eta/dX1;
            nodes[j].P[2]=1.-pow(Sigma.at(i)/dX1,2)*dT.at(i)-pow(Eta/dX1,2);
            nodes[j].P[3]=0.5*pow(Sigma.at(i)/dX1,2)*dT.at(i)+0.5*pow(Eta/dX1,2)-0.5*Eta/dX1;

        }

        slices[i+1].nodeMin=nodes[0].k-1;
        slices[i+1].nodeMax=nodes[nbNodes-1].k+1;

        if( fabs(Beta.at(i)*dT.at(i)) > 1 )
        {
            for( l = 0 ; l < nbNodes;l++)
            {
                if( nodes[l].k-1 < slices[i+1].nodeMin )
                    slices[i+1].nodeMin = nodes[l].k-1;
                if( nodes[l].k+1 > slices[i+1].nodeMax )
                    slices[i+1].nodeMax = nodes[l].k+1;
            }
        }

        slices[i+1].nbNodes = slices[i+1].nodeMax - slices[i+1].nodeMin + 1;
        slices[i+1].nodes=new Node[slices[i+1].nbNodes];
        for( j = 0; j < slices[i+1].nbNodes; j++ )
            slices[i+1].nodes[j].Init();

        nodeMin1=slices[i+1].nodeMin;
        nbNodes1=slices[i+1].nbNodes;

        nodes1=slices[i+1].nodes;

        for( j = 0; j < nbNodes1; j++ )
            nodes1[j].j=nodeMin1+j;
    }
}
// Tree::Init


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Tree::SetMultiplicateurInConstVolCase(DKMaille<double> TimeDependantSigma)
{
    if( Sigma.entries() != TimeDependantSigma.entries() )
    {
        throw("Tree.SetMultiplicateurInConstVolCase : Les vecteurs de vol dependantes du temps et de vol constantes doivent avoir la meme taille");
    }

    for( int i = 0; i < nbSlices; i++)
    {
        slices[i].Multiplicateur = TimeDependantSigma.at(i) / Sigma.at(i);
    }
}
// Tree::SetMultiplicateurInConstVolCase

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Tree::SetArrowDebreu(   double (*dFromXtoRate)(double, double, double),
                        double (*dFromXtoRate_deriv)(double, double, double),
                        CcyMarketData ccyMarketData,
                        DKMaille<double> q,
                        DKMaille<double> K)
{
    int i, j, nbNodes, nbNodes1, nodeMin1, iCount;
    double dX, dX1, discountFactor, Multiplicateur1, T1, dT1;;
    double dGuess, dGuessMin, dGuessMax, ZCtoFit, ZC_Calculated, ZC_Calculated_prime, MaxError;
    Node* nodes;
    Node* nodes1;

    // set the Arrow Debreu Price of the first node at 1
    slices[0].nodes[0].ArrowDebreu=1;

    slices[0].X_Adjustment = FromRateToX( - log(ccyMarketData.BasisZC(T.at(1))) / dT.at(0), q.at(0), K.at(0)) / slices[0].Multiplicateur;
    slices[0].nodes[0].X = slices[0].X_Adjustment * slices[0].Multiplicateur;
    slices[0].nodes[0].r = dFromXtoRate( slices[0].nodes[0].X, q.at(0), K.at(0));



    //for( i = 0; i < nbSlices - 1 && i + 2 < T.entries(); i++ )
    for( i = 0; i < nbSlices - 1; i++ )
    {
        nodes = slices[i].nodes;
        nbNodes = slices[i].nbNodes;
        dX = slices[i].dX;
        nodes1 = slices[i+1].nodes;
        nodeMin1 = slices[i+1].nodeMin;
        dX1 = slices[i+1].dX;
        nbNodes1 = slices[i+1].nbNodes;
        Multiplicateur1 = slices[i+1].Multiplicateur;

        // ayant l'ajustement et les prix d'Arrow Debreu sur la slice i
        // on calcule les prix d'Arrow Debreu sur la slice i+1
        for( j = 0; j < nbNodes; j++)
        {
            discountFactor=exp( -nodes[j].r * dT.at(i) ) * nodes[j].ArrowDebreu;

            nodes1[nodes[j].k+1-nodeMin1].ArrowDebreu+=nodes[j].P[1]*discountFactor;
            nodes1[nodes[j].k-nodeMin1].ArrowDebreu+=nodes[j].P[2]*discountFactor;
            nodes1[nodes[j].k-1-nodeMin1].ArrowDebreu+=nodes[j].P[3]*discountFactor;
        }

        // On calcule l'ajustement sur la slice i+1
        if( i + 2 < T.entries() )
        {
            ZCtoFit = ccyMarketData.BasisZC(T.at(i+2));
            dT1 = dT.at(i+1);
        }
        else
        {
            T1 = T.at(T.entries()-1) + (T.at(T.entries()-1) -T.at(T.entries()-2)) * (i + 3 - T.entries());
            dT1 = T.at(T.entries()-1) -T.at(T.entries()-2) * (i + 3 - T.entries());
            ZCtoFit = ccyMarketData.BasisZC( T1 ) ;
        }

        // on prend comme dGuess l'ajutement de la slice precedente
        dGuess = slices[i].X_Adjustment;
        dGuessMin = -1.e6;
        dGuessMax = 1.e6;

        ZC_Calculated=0.;
        ZC_Calculated_prime=0.;
        MaxError=1.e-10;

        for( j = 0; j < nbNodes1; j++)
        {
            if(nodes1[j].ArrowDebreu > TreeLimite)
            {
                ZC_Calculated += nodes1[j].ArrowDebreu
                                 * exp(-dFromXtoRate((X0+nodes1[j].j*dX1+dGuess)*Multiplicateur1, q.at(i+1), K.at(i+1))*dT1);
            }
        }

        iCount=0;
        while( fabs(ZC_Calculated-ZCtoFit) > MaxError )
        {
            if( iCount > 100 )
            {
                throw("Tree.SetArrowDebreu : Calcul de l'ajustement impossible");
            }

            if( (ZC_Calculated -ZCtoFit) > 0 )
            {
                dGuessMin = dGuess;
            }
            else
            {
                dGuessMax = dGuess;
            }

            ZC_Calculated_prime = 0.;
            for( j = 0; j < nbNodes1; j++)
            {
                if(nodes1[j].ArrowDebreu > TreeLimite)
                {
                    ZC_Calculated_prime -= nodes1[j].ArrowDebreu
                                           * exp(-dFromXtoRate((X0+nodes1[j].j*dX1+dGuess)*Multiplicateur1, q.at(i+1), K.at(i+1))*dT1)
                                           * dFromXtoRate_deriv((X0+nodes1[j].j*dX1+dGuess)*Multiplicateur1, q.at(i+1), K.at(i+1))*Multiplicateur1*dT1;
                }
            }

            if( fabs( ZC_Calculated_prime ) > 1.e-20 )
            {
                if(fabs(( ZC_Calculated - ZCtoFit ) / ZC_Calculated_prime ) < 30 )
                {
                    dGuess = - ( ZC_Calculated - ZCtoFit ) / ZC_Calculated_prime + dGuess;
                }
                else
                {
                    dGuess -=  0.5 * fabs(( ZC_Calculated - ZCtoFit ) / ZC_Calculated_prime) / (( ZC_Calculated - ZCtoFit ) / ZC_Calculated_prime) * fabs(dGuess);
                }
            }
            else
            {
                throw("Tree.SetArrowDebreu : Calcul de l'ajustement impossible, Derivee nulle");
            }


            if( dGuess <= dGuessMin )
            {
                dGuess = FMIN( dGuessMin + 10., (dGuessMin + dGuessMax) / 2.);
            }
            else if( dGuess >= dGuessMax )
            {
                dGuess = FMAX( dGuessMax - 10., (dGuessMin + dGuessMax) / 2.);
            }



            ZC_Calculated = 0.;
            for( j = 0; j < nbNodes1; j++)
            {
                if(nodes1[j].ArrowDebreu > TreeLimite)
                {
                    ZC_Calculated += nodes1[j].ArrowDebreu
                                     * exp(-dFromXtoRate((X0+nodes1[j].j*dX1+dGuess)*Multiplicateur1, q.at(i+1), K.at(i+1))*dT1);
                }
            }

            iCount++;
        }

        slices[i+1].X_Adjustment = dGuess;

        // Ayant l'ajustement on calcule le taux court a chaque noeud
        for( j = 0; j < nbNodes1; j++)
        {
            nodes1[j].X = (X0+nodes1[j].j*dX1+dGuess)*Multiplicateur1;
            nodes1[j].r = dFromXtoRate(nodes1[j].X, q.at(i+1), K.at(i+1));
            nodes1[j].X -= dGuess * Multiplicateur1;
        }
    }
}
// Tree::SetArrowDebreu


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void Tree::SetZCinTree()
{
    int i, j, k, nbSlicesAfter, nbNodes, nodeMin1;
    Node* nodes;
    Node* nodes1;
    Node nodeJ;

    for( i = nbSlices - 2; i >= 0; i-- )
    {
        nbSlicesAfter = nbSlices - i -1;
        nbNodes = slices[i].nbNodes;
        nodes = slices[i].nodes;
        nodes1 = slices[i+1].nodes;
        nodeMin1=slices[i+1].nodeMin;

        for( j = 0; j < nbNodes; j++)
        {
            nodes[j].ZC = new double[nbSlicesAfter];
            nodes[j].ZC[0] = exp( - nodes[j].r * dT[i]);
            for( k = 1; k < nbSlicesAfter; k++ )
            {
                nodes[j].ZC[k] = exp( - nodes[j].r * dT[i]) * (
                                     nodes[j].P[1]*nodes1[nodes[j].k+1-nodeMin1].ZC[k-1]
                                     +nodes[j].P[2]*nodes1[nodes[j].k-nodeMin1].ZC[k-1]
                                     +nodes[j].P[3]*nodes1[nodes[j].k-1-nodeMin1].ZC[k-1]);
            }
        }
    }
}
// Tree::SetZCinTree

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double Tree::CalculateZCinTree(int slice,
                               int node,
                               double t)
{
    if( t < T.at(slice))
    {
        throw("Tree.CalculateZCinTree : l'expiration du ZC doit etre au dela de la date de la slice");
    }

    if( t == T.at(slice))
    {
        return 1.;
    }

    int i = slice + 1;

    while( i < nbSlices && T.at(i) < t )
    {
        i++;
    }

    if( i == nbSlices )
    {
        throw("Tree.CalculateZCinTree : l'expiration du ZC est au dela de la derniere date de l'arbre!");
    }

    double ZC1, ZC2;

    if( i - 1 == slice )
    {
        ZC1 = 1.;
    }
    else
    {
        ZC1 = slices[slice].nodes[node].ZC[i-2-slice];
    }

    ZC2 = slices[slice].nodes[node].ZC[i-1-slice];

    return ZC_interpole(T.at(i-1), T.at(i), ZC1, ZC2, t);

}
// Tree::CalculateZCinTree

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
