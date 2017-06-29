#include <algorithm>
#include <iostream>

#include "BinningData.h"

BinningData::BinningData()
{
    SetPtId(0);
    SetEtaId(0);
}

BinningData::~BinningData()
{

}

bool BinningData::SetPtId(int id)
{
    if ((id < 0) || (id > _MAX_PT_BINS))
    {
        std::cout << "Error! BinningData::SetPtId - Wrong bin number" << endl;
        return false;
    }
    
    NptBins=id;
    return true;
}

bool BinningData::SetEtaId(int id)
{
    if ((id < 0) || (id > _MAX_ETA_BINS))
    {
        std::cout << "Error! BinningData::SetEtaId - Wrong bin number" << endl;
        return false;
    }
    
    NetaBins=id;
    return true;
}

bool BinningData::SetRapidityId(int id)
{
    if ((id < 0) || (id > _MAX_RAPIDITY_BINS))
    {
        std::cout << "Error! BinningData::SetRapidityId - Wrong bin number" << endl;
        return false;
    }
    
    NrapidityBins=id;
    return true;
}

int BinningData::GetPtBin(Float_t fValue) const
{
    const Float_t *endBins = ptBins + NptBins;
    const Float_t *f = std::lower_bound(ptBins, endBins, fValue);
	return (f == ptBins || f == endBins) ? -1 : (f - ptBins) - 1;
}

int BinningData::GetEtaBin(Float_t fValue) const
{
    const Float_t *endBins = etaBins + NetaBins;
    const Float_t *f = std::lower_bound(etaBins, endBins, fValue);
	return (f == etaBins || f == endBins) ? -1 : (f - etaBins) - 1;
}

int BinningData::GetRapidityBin(Float_t fValue) const
{
    const Float_t *endBins = etaBins + NrapidityBins;
    const Float_t *f = std::lower_bound(rapidityBins, endBins, fValue);
	return (f == rapidityBins || f == endBins) ? -1 : (f - rapidityBins) - 1;
}

Float_t BinningData::GetPtBinContent(int iValue) const
{
    return ptBins[iValue];
}

Float_t BinningData::GetEtaBinContent(int iValue) const
{
    return etaBins[iValue];
}

Float_t BinningData::GetRapidityBinContent(int iValue) const
{
    return rapidityBins[iValue];
}

void BinningData::SetPtBinContent(Float_t fValue)
{
    ptBins[NptBins] = fValue;
    NextPt();
}

void BinningData::SetEtaBinContent(Float_t fValue)
{
    etaBins[NetaBins] = fValue;
    NextEta();
}

void BinningData::SetRapidityBinContent(Float_t fValue)
{
    rapidityBins[NetaBins] = fValue;
    NextRapidity();
}

void BinningData::NextPt()
{
    NptBins++;
}

void BinningData::NextEta()
{
    NetaBins++;
}

void BinningData::NextRapidity()
{
    NrapidityBins++;
}

void BinningData::SetPtBins(Float_t fValue[], int dim)
{
    if (SetPtId(dim) && checkBinsArray(fValue, dim))
    {
        std::copy(fValue, fValue + dim, ptBins);
    }
}

void BinningData::SetPtBins(const Float_t fValue[], const int dim)
{
    if (SetPtId(dim) && checkBinsArray(fValue, dim))
    {
        std::copy(fValue, fValue + dim, ptBins);
    }
}
 

void BinningData::SetEtaBins(Float_t fValue[], int dim)
{
    if (SetEtaId(dim) && checkBinsArray(fValue, dim))
    {
        std::copy(fValue, fValue + dim, etaBins); 
    }
}

void BinningData::SetEtaBins(const Float_t fValue[], const int dim)
{
    if (SetEtaId(dim) && checkBinsArray(fValue, dim))
    {
        std::copy(fValue, fValue + dim, etaBins); 
    }
}

void BinningData::SetRapidityBins(Float_t fValue[], int dim)
{
    if (SetRapidityId(dim) && checkBinsArray(fValue, dim))
    {
        std::copy(fValue, fValue + dim, rapidityBins); 
    }
}

void BinningData::SetRapidityBins(const Float_t fValue[], const int dim)
{
    if (SetRapidityId(dim) && checkBinsArray(fValue, dim))
    {
        std::copy(fValue, fValue + dim, rapidityBins); 
    }
}

int BinningData::GetPtBinSize() const
{
    return NptBins;
}

int BinningData::GetEtaBinSize() const
{
    return NetaBins;
}

int BinningData::GetRapidityBinSize() const
{
    return NrapidityBins;
}

bool BinningData::checkBinsArray(Float_t *fValue, int dim)
{
    for (int i = 0; i < (dim - 1); i++)
    {
        if (fValue[i] >= fValue[i + 1])
        {
            std::cout << "Error! BinningData::checkBinsArray - wrong borders of " << i << "'s bin!" << endl;
            return false;
        }
        
    }
    
    return true;
}

bool BinningData::checkBinsArray(const Float_t *fValue, const int dim)
{
    for (int i = 0; i < (dim - 1); i++)
    {
        if (fValue[i] >= fValue[i + 1])
        {
            std::cout << "Error! BinningData::checkBinsArray - wrong borders of " << i << "'s bin!" << endl;
            return false;
        }
        
    }
    
    return true;
}