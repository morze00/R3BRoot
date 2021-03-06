#ifndef TCACONNECTOR
#define TCACONNECTOR

#include "FairRootManager.h"
#include "TClonesArray.h"
#include "TString.h"
#include <exception>
#include <utility>
#include <vector>

template <typename T>
class TCAInputConnector
{
  private:
    TString fBranchName;
    TString fClassName;
    TClonesArray* fTCA; // non-owning

  public:
    TCAInputConnector(TString b)
        : fBranchName(std::move(b))
        , fClassName(T().ClassName())
        , fTCA(nullptr)
    {
    }

    void Init()
    {
        auto ioman = FairRootManager::Instance();
        if (ioman == nullptr)
        {
            throw std::runtime_error("TCAInputConnector: No FairRootManager");
        }
        fTCA = (TClonesArray*)ioman->GetObject(fBranchName);
        if (fTCA == nullptr)
        {
            throw std::runtime_error(TString("TCAInputConnector: No TClonesArray called " + fBranchName +
                                             " could be obtained from the FairRootManager")
                                         .Data());
        }
        if (!TString(fTCA->GetClass()->GetName()).EqualTo(fClassName))
        {
            throw std::runtime_error(TString("TCAInputConnector: TClonesArray " + fBranchName +
                                             " does not contain elements of type " + fClassName)
                                         .Data());
        }
    }

    std::vector<T*> Retrieve() const
    {
        std::vector<T*> fV;
        if (fTCA == nullptr)
        {
            throw std::runtime_error("TCAInputConnector: TClonesArray not available");
        }

        const Int_t n = fTCA->GetEntries();
        fV.reserve(n);
        for (Int_t i = 0; i < n; i++)
        {
            fV.emplace_back((T*)fTCA->At(i));
        }
        return fV;
    }

    std::vector<T> RetrieveObjects() const
    {
        std::vector<T> fV;
        if (fTCA == nullptr)
        {
            throw std::runtime_error("TCAInputConnector: TClonesArray not available");
        }

        const Int_t n = fTCA->GetEntries();
        fV.reserve(n);
        for (Int_t i = 0; i < n; i++)
        {
            fV.emplace_back(*(T*)fTCA->At(i));
        }
        return fV;
    }
};

template <typename T>
class TCAOptionalInputConnector
{
  private:
    TString fBranchName;
    TString fClassName;
    TClonesArray* fTCA; // non-owning

  public:
    TCAOptionalInputConnector(TString b)
        : fBranchName(std::move(b))
        , fClassName(T().ClassName())
        , fTCA(nullptr)
    {
    }

    void Init()
    {
        auto ioman = FairRootManager::Instance();
        if (ioman == nullptr)
        {
            throw std::runtime_error("TCAInputConnector: No FairRootManager");
        }
        fTCA = (TClonesArray*)ioman->GetObject(fBranchName);
        if (fTCA != nullptr && !TString(fTCA->GetClass()->GetName()).EqualTo(fClassName))
        {
            throw std::runtime_error(TString("TCAInputConnector: TClonesArray " + fBranchName +
                                             " does not contain elements of type " + fClassName)
                                         .Data());
        }
    }

    std::vector<T*> Retrieve() const
    {
        std::vector<T*> fV;
        if (fTCA == nullptr)
        {
            return fV;
        }

        const Int_t n = fTCA->GetEntries();
        fV.reserve(n);
        for (Int_t i = 0; i < n; i++)
        {
            fV.emplace_back((T*)fTCA->At(i));
        }
        return fV;
    }

    std::vector<T> RetrieveObjects() const
    {
        std::vector<T> fV;
        if (fTCA == nullptr)
        {
            return fV;
        }

        const Int_t n = fTCA->GetEntries();
        fV.reserve(n);
        for (Int_t i = 0; i < n; i++)
        {
            fV.emplace_back(*(T*)fTCA->At(i));
        }
        return fV;
    }
};

template <typename T>
class TCAOutputConnector
{
  private:
    TString fBranchName;
    TString fClassName;
    TClonesArray* fTCA; // non-owning

  public:
    TCAOutputConnector(TString b)
        : fBranchName(std::move(b))
        , fClassName(T().ClassName())
        , fTCA(nullptr)
    {
    }

    void Init()
    {
        auto ioman = FairRootManager::Instance();
        if (ioman == nullptr)
        {
            throw std::runtime_error("TCAOutputConnector: No FairRootManager");
        }
        fTCA = ioman->Register(fBranchName, fClassName, "", kTRUE);
    }

    void Store(std::vector<T>& v)
    {
        if (fTCA == nullptr)
        {
            throw std::runtime_error("TCAOutputConnector: TClonesArray not available");
        }

        fTCA->Clear();

        for (auto& o : v)
        {
            new ((*fTCA)[fTCA->GetEntries()]) T(std::move(o));
        }

        v.clear();
    }
};

#endif // TCACONNECTOR
