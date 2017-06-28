#include <algorithm>

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
        cout << "Error! BinningData::SetPtId - Wrong bin number";
        return false;
    }
    
    NptBins=id;
    return true;
}

bool BinningData::SetEtaId(int id)
{
    if ((id < 0) || (id > _MAX_ETA_BINS))
    {
        cout << "Error! BinningData::SetEtaId - Wrong bin number";
        return false;
    }
    
    NetaBins=id;
    return true;
}

int BinningData::GetPtBin(Float_t fValue)
{
	int pt_bin = -1;
	for (int i = 0; i < NptBins; ++i)
	if ((fValue > ptBins[i]) && (fValue <= ptBins[i + 1])) 
		pt_bin = i;
	return pt_bin;
}

int BinningData::GetEtaBin(Float_t fValue)
{
	int eta_bin = -1;
	for (int i = 0; i < NetaBins; ++i)
	if ((fValue > etaBins[i]) && (fValue <= etaBins[i + 1])) 
		eta_bin = i;
	return eta_bin;
}

Float_t BinningData::GetPtBinContent(int iValue)
{
    NptBins++;
    return ptBins[iValue];
}

Float_t BinningData::GetEtaBinContent(int iValue)
{
    NetaBins++;
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
    if (SetPtId(dim))
    {
        std::copy(fValue, fValue+dim+1, ptBins);
    }
}
 

void BinningData::SetEtaBins(Float_t fValue[], int dim)
{
    if (SetEtaId(dim))
    {
        std::copy(fValue, fValue+dim+1, etaBins);
    }
}

int BinningData::GetPtBinSize()
{
    return NptBins;
}

int BinningData::GetEtaBinSize()
{
    return NetaBins;
}