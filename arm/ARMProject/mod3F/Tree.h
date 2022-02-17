#include "DKMaille.h"
#include "DKMaille2D.h"
#include "MarketData.h"

// FIXMEFRED: mig.vc8 (25/05/2007 11:10:37): missing return types
class Node  
{
    public:
	    Node();
	    virtual ~Node();
        void Init();

	    int j;						// situation du noeuds
	    int k;						// noeud fils central
	    double* P;					// probabilités d'atteindre les noeuds k+1, k, k-1 en P[1], 
	    double r;					// taux court sur le noeud
        double X;                   // X egal X_tild multiplie par le mutiplicateur de la slice
        double GreenTaux;           // fonction de green multiplie par exp(-r.dt)
        double ArrowDebreu;         // prix de l'actif qui paye 1 si ce noeud est atteint     
		double* ZC;
        double SWOi;

};




class Slice
{
    public:
	    Slice();
	    virtual ~Slice();

        void Init(int);

	    int nbNodes;		    // Nb de noeuds
        int nodeMin;		    // Noeud minimum
	    int nodeMax;		    // Noeud maximum
	    Node* nodes;		    // Tableau des noeuds de l'étape
        double dX;              // taille des ecarts entre les noeuds de l'arbre
                                // proportionnel a Sigma*SQRT(3.dT)
        double X_Adjustment;    // ajustement a appliquer a tous les X de la slice 
                                // pour ajuster E(f(Xi+dX)) au prix du zero-coupon
		double Multiplicateur;  // valeur par laquelle on doit multiplier X 
								// quand on a cree un arbre en Sigma constant
        bool isNotice;
};




class Tree
{
    public:
        Tree();
        virtual ~Tree();

        void Clean();
        
        void Init(	double StartPoint, 
				int dNbSteps,
				DKMaille<double> MaturityDates, 
				DKMaille<double> TimeDependantSigma,
				DKMaille<double> TimeDependantMeanReversion);
        
		void SetArrowDebreu(	double (*func)(double, double, double), 
                        double (*dFromXtoRate_deriv)(double, double, double),
						CcyMarketData ccyMarketData,
                        DKMaille<double> q, 
                        DKMaille<double> K);
        
        void SetZCinTree();

		void SetMultiplicateurInConstVolCase(DKMaille<double> TimeDependantSigma);

		double CalculateZCinTree(	int slice, 
									int node, 
									double t);
		
		double X0;
        DKMaille<double> T;
        DKMaille<double> dT;
        DKMaille<double> Sigma;
		DKMaille<double> Beta;
        int nbSlices;
        int nbSteps;                // nbSteps=nbSlices-1 puisque les slices sont
                                    // les piquets et les steps sont les intervalles entre les piquets
        Slice* slices;
        double TreeLimite;

};