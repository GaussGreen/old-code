#include "xxxproject/scenario.h"
#include "xxxproject/argconvdefault.h"

#include <mod/bssmiled.h>
#include <crv/zerocurv.h>
#include <crv/volcurv.h>
#include <crv/zeroint.h>
#include <crv/zerointspreaded.h>



CC_BEGIN_NAMESPACE( ARM )



ARM_Scenario::ARM_Scenario	(	const double & shift, const string & ccy, const int & isRelative, const int & isCumulInv, const int & isPerturbative ):	ARM_RootObject	(		), 
	itsShift		( shift	),
	itsCurrency		( ccy	),
	isRelative		( isRelative ),
	isPerturbative	( isPerturbative ),	
	isCumulInv		( isCumulInv ){
	ARM_Scenario::ARM_Scenario(); 
}


ARM_Scenario::ARM_Scenario	(	const ARM_Scenario & rhs):	
	ARM_RootObject	(	rhs					),
	itsDim			(	rhs.itsDim			),
	itsShift		(	rhs.itsShift		),
	itsCurrency		(	rhs.itsCurrency		),
	isRelative		(	rhs.isRelative		),
	isCumulInv		(	rhs.isCumulInv		),
	isPerturbative	(	rhs.isPerturbative	),	
	itsNbBump		(	rhs.itsNbBump		),
	itsPosition		(	rhs.itsPosition		),
	itsPermOrder	(	rhs.itsPermOrder	),
	itsPlotsPos		(	rhs.itsPlotsPos		),
	itsPlotsLab		(	rhs.itsPlotsLab		){ }


void	ARM_Scenario::InitDimScenario	(	const string & scen, const int & n ){

	itsPlotsPos.resize(n);
	itsPlotsLab.resize(n);
	itsPosition.resize(n);
	itsPermOrder.resize(n);

	ConvertString( scen );

	InitPosition( );
}

int		ARM_Scenario::GetNbEffShift ( )		const{

	map<string,ScenarioType>::const_iterator	it;

	int Index=0;
	for(it = itsPriceCritera.begin(); it !=itsPriceCritera.begin(); it++){
		if ( it->second != NONE)
			Index++;
	}
	return Index; 
}


void	ARM_Scenario::ShiftNextPlot	(){
	
	int nbDims(itsPosition.size());
 
	for(int j=0;j<nbDims;j++){
		itsPosition[itsPermOrder[j]]++;
		if (	itsPosition[itsPermOrder[j]] == itsPlotsLab[itsPermOrder[j]].size()	)
				itsPosition[itsPermOrder[j]] = 0;
		else
				break;
	}
}

void	ARM_Scenario::InitPosition			( ) {

	for ( int i=0; i< itsPosition.size(); i++)
		itsPosition[i] = 0;
	itsPosition[itsPermOrder[0]] = -1;
}


vector<string>	ARM_Scenario::ReduceString	( const string & str ){

	int					pos;
	string				tmp;
	string				tmpStr;
	vector< string >	tmpV_Str;

	pos = str.find("=");
	if ( pos>=0 && pos<= str.size() )
			tmpStr.insert(0, str.begin()+pos+1,str.end() );
	else
			tmpStr.resize(0);

	pos = tmpStr.find("=");
	if( pos>=0 && pos<= tmpStr.size() ) 
		throw Exception(__LINE__, __FILE__, ERR_DATE_NOT_VALID, "the string should contain at most one sign =" );


	while ( tmpStr.size() != 0 ){

		pos = tmpStr.find(",");

		if ( pos>=0 && pos<= tmpStr.size() ) {
			tmp.erase();
			tmp.insert(0,tmpStr.begin(),tmpStr.begin()+pos );
			tmpV_Str.push_back( tmp );
			tmp.erase();
			tmp.insert(0,tmpStr.begin()+pos+1,tmpStr.end() );
			tmpStr = tmp;
		}
		else{
			tmpV_Str.push_back( tmpStr );
			break;
		}
	}
	return tmpV_Str;
}
 
void		ARM_Scenario::OrderData( ARM_GP_VectorPtr data){

	if ( isCumulInv) {
		int					dim	=	data ->size();	
		ARM_GP_VectorPtr	tmp (  (ARM_GP_T_Vector<double>*) (&*data)->Clone() );

		for( int i=0; i<  dim; i++)	data->Elt(i)= tmp->Elt(dim-i-1);
	}
}


void		ARM_Scenario::OrderData( ARM_GP_MatrixPtr data){

	if ( isCumulInv) {
		int					row	=	data -> rows();	
		int					col =	data -> cols();	
		ARM_GP_MatrixPtr	tmp (  (ARM_GP_T_Matrix<double>*) (&*data)->Clone() );

		for( int i=0; i<  row; i++)	{
			for( int j=0; j< col; j++)
					data->Elt(i,j)= tmp->Elt(row-i-1, col-j-1);
		}
	}
}

string  ARM_Scenario::toString	(	const string & indent, 	const string & nextIndent)const {

	CC_Ostringstream	os;

	os	<< "\n";
	os	<< CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(30)<<" ";
	os	<< CC_NS(std,setfill('=')) << CC_NS(std,fixed) << CC_NS(std,setw)(30)<<"\n";
	os	<< CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(50)<<"SCENARIO"<<"\n";
	os	<< CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(30)<<" ";
	os	<< CC_NS(std,setfill('=')) << CC_NS(std,fixed) << CC_NS(std,setw)(30)<<"\n";

	os	<< "\n\n";
			
	os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(30)<< CC_NS(std,setprecision)(3);
	os <<CC_NS(std,right)<<"Shift";
	os <<"  :  "<< itsShift<<"\n";

	os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(30);
	os <<CC_NS(std,right)<<"Currency";
	os <<"  :  "<< itsCurrency<<"\n";

	os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(30);
	os <<CC_NS(std,right)<<"is Cumulative";
	if ( isCumulInv==K_YES)
		os <<"  :  "<< string("YES")<<"\n";
	else
		os <<"  :  "<< string("NO")<<"\n";

	os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(30);
	os <<CC_NS(std,right)<<"is Relative";
	if ( isRelative==K_YES)
		os <<"  :  "<< string("YES")<<"\n";
	else
		os <<"  :  "<< string("NO")<<"\n";

	os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(30);
	os <<CC_NS(std,right)<<"is Perturbative";
	if ( isPerturbative==K_YES)
		os <<"  :  "<< string("YES")<<"\n";
	else
		os <<"  :  "<< string("NO")<<"\n";

	os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(30);
	os <<CC_NS(std,right)<<"Scenario type";
	os <<"  :  "<< GetScenarioKey()<<"\n";

	if ( isInitialized ){
		os	<< "\n" <<"=========> Stress Table :"<<"\n\n";
		os<< ViewPosition();
	}

	return os.str();
}

			/*	*	*	*	*	*	*	*	*	*	*	*	*	*	*/
			/*														*/
			/*						ND SCENARIO						*/
			/*														*/
			/*	*	*	*	*	*	*	*	*	*	*	*	*	*	*/



string	ARM_ND_Scenario::toString(	const string& indent, const string& nextIndent)  const {


	CC_Ostringstream		os;
	string					key;
	string					tmp("NO");
	map<string,ScenarioType>::const_iterator	it;
	vector<ARM_StringVector	>		lab;
	ARM_IntVector					pos;

	for ( int i =0; i<itsDim; i++)
		os<<itsScenario[i]->toString(indent,nextIndent);

	os	<< "\n" <<"=========> Composite Stress Table :"<<"\n\n";


	const_cast<ARM_ND_Scenario*> (this)->InitPosition();

	for( i=0; i<itsNbBump; i++ ){
		const_cast<ARM_ND_Scenario*> (this)->ShiftNextPlot( );
		const_cast<ARM_ND_Scenario*> (this)->NotifySubScenario( );
		key.resize(0);
		for(int j=0; j<itsDim; j++){
			if ( key.size()!=0 )	key = key + string(":");
			key = key + itsPlotsLab[j][itsPosition[j] ];
		}
			
		it = itsPriceCritera.find(key);
		switch ( it->second){
			case VALUE:{
				os << CC_NS(std,setfill(' ')) << CC_NS(std,fixed) << CC_NS(std,setw)(50);
				os <<CC_NS(std,right)<<key;
				os <<"  -  ";
				os <<CC_NS(std,right)<<tmp;
				os << "\n";
				tmp = key;
				break;
		   }
			case REF:
				tmp = key;
				break;
		}
	}
	const_cast<ARM_ND_Scenario*> (this)->InitPosition();

	return os.str( );
}


ARM_ND_Scenario::~ARM_ND_Scenario() { 
	for ( int i =0; i< itsDim; i++) { delete itsScenario[i]; itsScenario[i]=NULL;} 
}

ARM_ND_Scenario::ARM_ND_Scenario( vector<ARM_Scenario*> & sce ) : ARM_Scenario(){ 
	itsDim =sce.size();
	itsScenario.resize( itsDim) ;

	for ( int i=0 ; i<itsDim; i++){
		itsScenario[i] = dynamic_cast<ARM_Scenario*> ( sce[i] ->Clone() );
	}
	InitScenario();

	isRelative		= K_NO;
	isCumulInv		= K_NO;
	isPerturbative	= K_NO;
}

void ARM_ND_Scenario::InitScenario(const string& theScenario ) {	InitDimScenario	(theScenario,itsDim);	}

void ARM_ND_Scenario::ConvertString	(   const string & ){ 

	itsPosition.clear();
	for ( int i =0 ; i<itsDim; i++)
		itsPermOrder[i] = i;
}

ARM_ND_Scenario::ARM_ND_Scenario(	const ARM_ND_Scenario & rhs):ARM_Scenario(rhs){
	
	for ( int i=0 ; i<itsDim; i++)
		itsScenario[i] = dynamic_cast<ARM_Scenario *>  (rhs.itsScenario[i]->Clone() );

	isRelative		= rhs.isRelative;
	isCumulInv		= rhs.isCumulInv;
	isPerturbative	= rhs.isPerturbative;
}

void ARM_ND_Scenario::InitMktData	( ARM_MktData *	mkt		){

	int									isAcceptable = 1;
	string								key;
	string								tmp;
	int									dim;

	map<string,ScenarioType>			tmpMap;
	vector<ARM_StringVector	>			lab;
	ARM_IntVector						pos;

	vector<string>						tmpStr;
	vector<int>							tmpSce;
	
	ARM_StringVector::iterator			it;

	itsNbBump = 1;
	for( int i=0; i<itsDim;i++){
		itsScenario[i]->InitMktData(mkt);
		itsNbBump *= (1+ itsScenario[i]->GetNbShift());	

		tmpMap = itsScenario[i]->GetPriceCritera();
		itsPlotsLab[i].clear();
		itsPlotsPos[i].clear();
		itsPlotsLab[i].push_back(string("NO"));
		itsPlotsPos[i].push_back(string("NO"));

		lab		= itsScenario[i]->GetPlotsLabel();
		tmpMap	= itsScenario[i]->GetPriceCritera();

		itsScenario[i]->InitPosition();
		dim= itsScenario[i]->GetDim()>0?itsScenario[i]->GetDim():1;

		for( int j = 0; j< itsScenario[i]->GetNbShift(); j++){
	
			itsScenario[i]->ShiftNextPlot( );
			pos=itsScenario[i]->GetPosition();
			key.resize(0);
			
			for( int k = 0; k<dim; k++)	key = key + lab[k][pos[k]];
			itsPlotsLab[i].push_back(key);
			switch( tmpMap[key] ){
			
			case VALUE:
				itsPlotsPos[i].push_back(key);
				break;	
			case REF_VALUE:
				itsPlotsPos[i].push_back(key);
				break;
			default:
				break;
			} 
		}
		itsScenario[i]->InitPosition();
	}
	
	InitPosition();

	for( i=0; i<itsNbBump; i++ ){
		ShiftNextPlot( );
		key.resize(0);
		isAcceptable = 1;
		for(int j=0; j<itsDim; j++){
			tmp = itsPlotsLab[j][itsPosition[j] ];
			if ( key.size()!=0 ) key = key + string(":");
			key = key + tmp;
			it=find( itsPlotsPos[j].begin(), itsPlotsPos[j].end(), tmp );

			if( j==0)
				tmpSce.push_back(tmpMap[tmp]);

		}
		tmpStr.push_back(key);
	}
	
	for( i=0 ; i < tmpSce.size();i++){
		itsPriceCritera.insert(pair<string,ScenarioType>(tmpStr[i],(ScenarioType)tmpSce[i]) );
	}

	InitPosition();

}

void ARM_ND_Scenario::ShiftMktData( ARM_MktData *	mkt){
	NotifySubScenario();
	for( int i=0; i<itsDim;i++){
		if( itsPlotsLab[i][itsPosition[i] ] != string("NO") )
			itsScenario[i]->ShiftMktData(mkt);
	}
}

void ARM_ND_Scenario::UnShiftMktData( ARM_MktData *	mkt){
	NotifySubScenario();
	for( int i=0; i<itsDim;i++){
		if( itsPlotsLab[i][itsPosition[i] ] != string("NO") )
			itsScenario[i]->UnShiftMktData(mkt);
	}
}

void ARM_ND_Scenario::Finalize	( ARM_MktData *	mkt	)	const{
	const_cast<ARM_ND_Scenario*> (this)->NotifySubScenario();
	for( int i=0; i<itsDim;i++){
		itsScenario[i]->Finalize(mkt);
	}
}

string	ARM_ND_Scenario::GetScenarioKey	( )	const{
	string tmp("");

	for( int i=0; i<itsDim-1;i++)
		tmp += itsScenario[i] -> GetScenarioKey	( ) +"/";
	tmp += itsScenario[itsDim-1] -> GetScenarioKey	( );
	return tmp;
}


void	ARM_ND_Scenario::InitPosition( ) {

	ARM_Scenario::InitPosition( );

	for ( int i=0; i< itsDim; i++)
		itsScenario[i]->InitPosition();

	itsPosition.resize(itsDim);

	for (    i=0; i< itsDim; i++)
		itsPosition[i] = 0;
	itsPosition[itsPermOrder[0]] = -1;
}

void	ARM_ND_Scenario::NotifySubScenario( ){

	ARM_IntVector	pos;
	string			tmp;
	int				dim;
	
	for(int i=0;i<itsDim;i++){
		dim = itsScenario[i]->GetDim();
		itsScenario[i]->InitPosition();

		if( itsPlotsLab[i][itsPosition[i] ] != string("NO") ){
			
			vector<ARM_StringVector>	lab = itsScenario[i]->GetPlotsLabel();
			for( int k= 0; k<itsScenario[i]->GetNbShift(); k++ ){
				
				itsScenario[i]->ShiftNextPlot();
				tmp.resize(0);
				pos = itsScenario[i]->GetPosition();
				for(int  j=0; j<dim ;j++) tmp = tmp + lab[j][pos[j]]; 
				if ( tmp == itsPlotsLab[i][itsPosition[i] ] )	break;
			}
		}
	}
}


vector<string>	ARM_ND_Scenario::GetSubScenario	( )	const { 
	
	vector<string> tmp;

	for( int i = 0; i<itsDim; i++)
		tmp.push_back(itsScenario[i]->GetScenarioKey() );
	
	return tmp; 
}


CC_END_NAMESPACE()

