/*

Obtianed from : 

 
http://support.microsoft.com/default.aspx?scid=kb;EN-US;Q201321



Debugging by David Cuin



HOWTO: Alternative Implementation of ATL Singleton



Applies To
This article was previously published under Q201321 
SUMMARY
This article shows an implementation of singletons 
that differs from the default ATL implementation in 
the following ways: 


The lifetime of an ATL singleton is tied to the class 
factory; the singleton gets destroyed only when the class 
factory is destroyed, which is when the EXE/DLL unloads 
[CComModule::RevokeClassObject() and CComModule::Term()]. 
The singleton implementation in this article 
maintains a refcount and deletes itself when it goes 
to zero. 

Implementing singletons in a DLL can be problematic. 
One scenario is if your singleton object marked "Apartment",
is in a DLL. The client creates two STA threads, 
each one creates your singleton object. Both STA threads have 
the same raw pointer to your singleton, 
allowing them to simultaneously call into the singleton. 
You are probably not synchronizing access to instance data 
in your singleton object, because according to COM rules, 
you don't need to synchronize access to instance data 
for objects marked "Apartment". 

It's only a matter of time before data is corrupted.

The singleton implementation in this article works 
around this problem by returning a marshaled pointer 
in IClassFactory::CreateInstance(). 
But marshaling the pointer introduces another 
potential problem. All calls are now marshaled to the 
first Apartment that calls IClassFactory::CreateInstance. 
If the first Apartment goes away, calls to the 
singleton from other Apartments will fail. 
If you use the singleton implementation in this article in 
a DLL marked "Apartment", you must ensure that the first 
Apartment that uses the component stays alive. Because of all 
these potential problems, you should only use singletons 
in a DLL when absolutely necessary

If you just want to share data in the DLL within the same process, 
an alternative is to create global/static variables to hold the data. 
You just need to synchronize access to this data.


*/


/********** Start of singleton code **********/ 
template <class Base>
class  CMyComObject: public CComObject<Base>
{
public:
   // Variable that lets class factory know COM object is deleted.
   CMyComObject** m_ppObj;

   STDMETHOD_(ULONG, Release)()
   {
      ULONG l = InternalRelease();
      if (l == 0)
      {
         // Set to NULL so class factory knows to recreate us.
         *m_ppObj = NULL;
         delete this;
      }
      return l;
   }

   static HRESULT WINAPI CreateInstance(CMyComObject<Base>** pp)
   {
      ATLASSERT(pp != NULL);
      HRESULT hRes = E_OUTOFMEMORY;
      CMyComObject<Base>* p = NULL;
      ATLTRY(p = new CMyComObject<Base>())
      if (p != NULL)
      {
         p->SetVoid(NULL);
         p->InternalFinalConstructAddRef();
         hRes = p->FinalConstruct();
         p->InternalFinalConstructRelease();
         if (hRes != S_OK)
         {
            delete p;
            p = NULL;
         }
      }
      *pp = p;
      // Initialize m_ppObj to CComObject<> pointer in class factory.
      if (p) p->m_ppObj = pp;
      return hRes;
   }
};

template <class T>
class CMyComClassFactorySingleton : public CComClassFactory
{
public:
   CMyComObject<T>* m_pObj;   // actual COM object
   CComAutoCriticalSection m_sec;  // A critical section is needed
   // because class factory is stored in CComObjectCached. All clients
   // get the same raw pointer to the class factory. Any instance
   // data in this class factory must be synchronized.
   // IClassFactory
   STDMETHOD(CreateInstance)(LPUNKNOWN pUnkOuter, REFIID riid,
                             void** ppvObj)
   {
      if (NULL == ppvObj)
         return E_POINTER;
      *ppvObj = NULL;
      // Aggregation is not supported in Singletons.
      ATLASSERT(pUnkOuter == NULL);
      if (pUnkOuter != NULL)
         return CLASS_E_NOAGGREGATION;
      HRESULT hRes = E_FAIL;
      // Lock critical section.
      m_sec.Lock();
      // Create COM object.
      if (NULL == m_pObj)
      {
         hRes = CMyComObject<T>::CreateInstance(&m_pObj);
         if (SUCCEEDED(hRes))
         {
            // Free marshaled interface and release stream.
            if (m_pStream)
            {
               AtlFreeMarshalStream(m_pStream);
               m_pStream = NULL;
            }
            // Marshal interface into a stream.
            // Can't use AtlMarshalPtrInProc since we need to use
            // MSHLFLAGS_TABLEWEAK.
            hRes = CreateStreamOnHGlobal(NULL, TRUE, &m_pStream);
            if (SUCCEEDED(hRes))
            {
				// modified to use QI to get Iunknown* instead of casting pointer to CMyComObject<T> 
				CComPtr<IUnknown> spUnknown ; 
				IUnknown* pUnk = NULL; 
					
				 hRes = m_pObj->QueryInterface(&pUnk);
				 if (FAILED(hRes))
				 {
					 return E_FAIL;
				 }
	
               hRes = CoMarshalInterface(m_pStream, IID_IUnknown, pUnk, MSHCTX_INPROC, NULL, MSHLFLAGS_TABLEWEAK);
			   //hRes = CoMarshalInterface(m_pStream, IID_IUnknown, pUnk, MSHCTX_NOSHAREDMEM , NULL, MSHLFLAGS_TABLEWEAK);
               if (FAILED(hRes))
               {
                  m_pStream->Release();
                  m_pStream = NULL;
               }
            }
         }
         if (FAILED(hRes))
         {
            // Unlock critical section.
            m_sec.Unlock();
            return hRes;
         }
      }
      CComPtr<IUnknown> spUnk;
      // Unmarshal interface from stream.
      hRes = AtlUnmarshalPtr(m_pStream, IID_IUnknown, &spUnk);
      // unlock critical section
      m_sec.Unlock();
      if (FAILED (hRes))
         return hRes;
      // Call QI on unmarshaled interface.
      hRes = spUnk->QueryInterface(riid, ppvObj);
         return hRes;
   }
   LPSTREAM m_pStream; // IStream* used for storing marshaled interface.

   void FinalRelease()
   {
      // Free marshaled interface and release stream.
      if (m_pStream)
         AtlFreeMarshalStream (m_pStream);
   }

   CMyComClassFactorySingleton() : m_pObj(NULL), m_pStream(NULL)
   {
   }
};

#define DECLARE_CLASSFACTORY_SAFEPROC_SINGLETON(obj)    DECLARE_CLASSFACTORY_EX(CMyComClassFactorySingleton<obj>)

/********** End of singleton code **********/ 