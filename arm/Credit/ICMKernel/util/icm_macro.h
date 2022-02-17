//	IXIS CIB
//
//

#ifndef ICM_MACRO_H_
#define ICM_MACRO_H_

#include <expt.h>
#include <crtdbg.h>
#include <sstream>

//	--	Helper Class for Logging 
//
//	LogLevel 
//
//		SILENT	-	message is sent in the debugger & localhost (default for ICMLOG) 
//		INFO	-	message is sent in debugger, localhost, remotehost 
//		WARN	-	message is sent in debugger, localhost, remotehost 
//		FATAL	-	exception 
// 

typedef enum { SILENT,INFO,WARN,FATAL} qIcmLoggerLevel; 
class	IcmLogger
{ 
public:	static void logDebugger(qIcmLoggerLevel ,const std::string& file,long line,const std::string& msg); 
public:	static void logStream(std::ostream& o,const std::string& file,long line,const std::string&msg); 
} ;

// --	
extern std::ostream nout ;


//	--	Helper function to output memory report
//
static 
inline std::ostream& operator<<(std::ostream& o,const _CrtMemState& state)
{
	o<<" : NORMAL_BLOCK : " <<state.lSizes[_NORMAL_BLOCK]<<" : CRT_BLOCK : "<< state.lSizes[_CRT_BLOCK]<<" : CLIENT_BLOCK : " <<state.lSizes[_CLIENT_BLOCK] ; 
	o<<" : TotalCount : "<<state.lTotalCount<<" : HighWaterCount : "<< state.lHighWaterCount;
	return o ;
}

//
//		ICMLOG
//
//		ICMLOG("Dumping i="<<i); 			
//
#ifndef	ICMLOG
	#define ICMLOG(msgICMLOG) \
	{ \
		std::stringstream _sstr; _sstr<<msgICMLOG<<std::endl; \
		IcmLogger::logDebugger(SILENT ,__FILE__,__LINE__,_sstr.str()); \
	} 
#endif	// MYLOG

//
//		ICMMSG
//
//		ICMMSG(FATAL,"Dumping i="<<i); 			
//
#ifndef	ICMMSG
	#define ICMMSG(level,msgICMLOG) \
	{ \
		std::stringstream _sstr; _sstr<<msgICMLOG<<std::endl; \
		IcmLogger::logDebugger(level,__FILE__,__LINE__,_sstr.str()); \
	} 
#endif	// MYLOG

//
//		ICMSTREAM
//
//		ICMSTREAM(cout,"Dumping i="<<i); 
//
#ifndef	ICMSTREAM
	#define ICMSTREAM(streamICMSTREAM, msgICMSTREAM) \
	{ \
		std::stringstream _sstr; _sstr<<msgICMSTREAM<<std::endl; \
		IcmLogger::logStream(streamICMSTREAM,__FILE__,__LINE__,_sstr.str()); \
	} 
#endif	// ICMSTREAM

//
//		ICMTHROW
//
//		ICMTHROW(ERR_INVALID_MODEL,"Invalid Model Choice "<<modelName); 
//
#ifndef	ICMTHROW
	#define ICMTHROW(errICMTHROW,msgICMTHROW) \
	{ \
	std::stringstream _sstr; _sstr<<msgICMTHROW<<std::endl; \
		IcmLogger::logDebugger(FATAL,__FILE__,__LINE__,_sstr.str()); \
		throw Exception(__LINE__, __FILE__, errICMTHROW, _sstr.str()); \
	} 
#endif	// ICMTHROW


//
//		ICMTEST
//
//		bool ret(true)
//		ICMTEST(cout,i==1) 
//
#ifndef	ICMTEST
	#define ICMTEST(osICMTEST,theTestICMTEST) \
	{ \
	bool b = ( theTestICMTEST) ; \
	std::string str ; \
	if (b) str="NONREG:   SUCCESS:"; else str="NONREG:###FAILED:"; \
	ICMSTREAM(osICMTEST, str+ #theTestICMTEST ) \
	ret&=b ; \
	}
#endif // ICMTEST

//
//		ICM_BEGIN_MEMCHECK
//
#ifndef ICM_BEGIN_MEMCHECK
	#define ICM_BEGIN_MEMCHECK \
	{ _CrtMemState s1, s2, s3; \
	 _CrtMemCheckpoint( &s1 ); {
#endif	//BEGIN_MEMCHECK


//
//		ICM_END_MEMCHECK
//
#ifndef	ICM_END_MEMCHECK
	#define ICM_END_MEMCHECK(osICM_END_MEMCHECK) \
	} _CrtMemCheckpoint( &s2 ); \
	bool _memok = !_CrtMemDifference( &s3, &s1, &s2) ; \
	ICMTEST(osICM_END_MEMCHECK, _memok ) ; \
	if (!_memok) ICMSTREAM(osICM_END_MEMCHECK,s3); \
	} 
#endif	// ICM_END_MEMCHECK

//	--	Helper Class for Precise Chrono (ms level) 
//
class PreciseChronoImpl; // forward declaration of the (win32) implementation
class PreciseChrono
{
private:
	PreciseChronoImpl* impl; 
	PreciseChrono(const PreciseChrono&); 
	PreciseChrono& operator=(const PreciseChrono&); 
public:
	PreciseChrono(void) ; 
	~PreciseChrono(void) ; 
	void Start(void) ;
	void Stop(void) ;
	double GetDurationInSeconds(void) const ;
	double GetDurationInTicks(void) const ;
	double GetDurationInMilliseconds(void) const ;
	double GetDurationInMicroseconds(void) const ;
};

//
//		ICMTIMER
//
//		ICMTIMER(cout,"new+delete of T",1000, T* t=new T ;delete t;)
//
#ifndef ICMTIMER
	#ifdef _DEBUG
static std::string __timerdbg("  DEBUG"); 
	#else 
static std::string __timerdbg("RELEASE"); 
	#endif	
	#define ICMTIMER(osTIMER,__title,__i,__expr) \
	{	PreciseChrono timer; timer.Start() ; \
		for(unsigned int __j=0;__j<__i;__j++) { __expr ; } ; \
		timer.Stop(); \
		if ((__i/timer.GetDurationInMilliseconds())<1.) \
{ ICMSTREAM(osTIMER,"TIMER:"<<__title<<":"<<timer.GetDurationInMilliseconds()/__i<<" :ms : "<<__i/timer.GetDurationInMilliseconds()<<" :op/ms :"<<__i<<":"<<__timerdbg) ;} \
		else \
{ ICMSTREAM(osTIMER,"TIMER:"<<__title<<":"<<timer.GetDurationInMilliseconds()/__i<<" :ms : "<<(int)(__i/timer.GetDurationInMilliseconds())<<" :op/ms :"<<__i<<":"<<__timerdbg) ; }\
	}
#endif // ICMTIMER


class	TimeElapser
{
	unsigned long line ; 
	std::string msg,file; 
	PreciseChrono chrono; 
public:
	// TimeElapser(){}; 
	void set(const std::string& arg,const std::string& file_,long line_) 
	{
		line=line_; 
		file=file_;
		msg=arg;  
		chrono.Start(); 
	}
	~TimeElapser() 
	{
		chrono.Stop(); 
		std::stringstream sstr; 
		sstr<<"TimeElapser:"<<msg<<":elapsedTime:"<<chrono.GetDurationInMilliseconds()<<" ms"<<std::endl ; 
		IcmLogger::logDebugger(SILENT,file,line,sstr.str()); 
	}
}; 

//
//		ICMELAPSER
//
#ifndef ICMELAPSER
	#define ICMELAPSER(arg) \
	TimeElapser __elapser; { std::stringstream sstr; sstr<<arg; __elapser.set(sstr.str(),__FILE__,__LINE__); }
#endif // ICMQUANTIFYER



//	--	Helper Class for Quantifyer 
class QuantifyerInfo
{
	friend class Quantifyer; 
	unsigned long	callCount,line ;
	double			avgTime,maxTime,minTime,elapsedTime;
	std::string		msg,file; 
public:
	QuantifyerInfo(const std::string&arg,const std::string& _file,long _line) ; 
	~QuantifyerInfo() ; 
} ; 

//	--	Helper Class for Quantifyer
//
class Quantifyer
{
	QuantifyerInfo	&quantifyerInfo; 
	PreciseChrono chrono; 
public:
	Quantifyer(QuantifyerInfo&	ref) ; 
	~Quantifyer() ; 
}; 
 

//
//	QUANTIFYER
//
//		QUANTIFYER("MyPosition") 
//
#ifndef ICMQUANTIFYER
#define ICMQUANTIFYER(arg) \
	static QuantifyerInfo __quantifyerInfo(arg,__FILE__,__LINE__); \
	Quantifyer __quantifyer(__quantifyerInfo);
#endif // ICMQUANTIFYER

//
//	Temporary ARM Helper 
//	
//	remove the const thing for those ARM functions that should accept const objects
template <class T> 
static inline 
T& 
unconst(const T&t) 
{ 
	return const_cast<T&>(t); 
}

//	
//	NAG HELPERS
// 
#ifndef ICM_NAGSILENT 
#define ICM_NAGSILENT(argICM_NAGSILENT) \
	{	e04xxc(&argICM_NAGSILENT) ; \
		argICM_NAGSILENT.print_level  = Nag_NoPrint; \
		argICM_NAGSILENT.output_level = Nag_NoOutput; \
		argICM_NAGSILENT.minor_print_level = Nag_NoPrint ; \
		argICM_NAGSILENT.list = FALSE ;  \
		argICM_NAGSILENT.print_iter=-1 ; \
		argICM_NAGSILENT.print_gcheck=FALSE ; \
		argICM_NAGSILENT.print_deriv = Nag_D_NoPrint; \
	}
#endif // ICM_NAGSILENT
#ifndef ICM_NAGFREE 
#define ICM_NAGFREE(argICM_NAGFREE) \
	{	NagError fail; \
		INIT_FAIL(fail) ; \
		e04xzc(&argICM_NAGFREE,"all",&fail) ; \
	}
#endif // ICM_NAGSILENT

class VersionInfoSingleton 
{
public:
	static std::string version() ; 
	static std::string longversion() ;
	VersionInfoSingleton() ; 
	~VersionInfoSingleton() ;
}; 


#endif	// ICM_MACRO_H_

