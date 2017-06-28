#include "BinningData.h"

BinningData::BinningData()
{
    SetPtId(0);
    SetEtaId(0);
}

BinningData::~BinningData()
{

}

void BinningData::SetPtId(int id)
{
    NptBins=id;
}

void BinningData::SetEtaId(int id)
{
    NetaBins=id;
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
    if (dim > 0 && dim < _MAX_PT_BINS){
        SetPtId(dim);
        for (int i=0;i<dim;i++){
            ptBins[i] = fValue[i];
        }
    }
}

void BinningData::SetEtaBins(Float_t fValue[], int dim)
{
    if (dim > 0 && dim <_MAX_ETA_BINS){
        SetEtaId(dim);
        for (int i=0;i<dim;i++){
            etaBins[i] = fValue[i];
        }
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