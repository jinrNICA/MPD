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

Float_t BinningData::GetPtBinContent(int iValue) const
{
    return ptBins[iValue];
}

Float_t BinningData::GetEtaBinContent(int iValue) const
{
    return etaBins[iValue];
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

void BinningData::NextPt()
{
    NptBins++;
}

void BinningData::NextEta()
{
    NetaBins++;
}

void BinningData::SetPtBins(Float_t fValue[], int dim)
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

int BinningData::GetPtBinSize() const
{
    return NptBins;
}

int BinningData::GetEtaBinSize() const
{
    return NetaBins;
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