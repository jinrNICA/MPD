#ifndef BINNING_DATA_H

#define BINNING_DATA_H

#define _MAX_PT_BINS 100
#define _MAX_ETA_BINS 100

class BinningData {
    public:
        BinningData();
        ~BinningData();
        int GetPtBin(Float_t fValue) const;
        int GetEtaBin(Float_t fValue) const;
        Float_t GetPtBinContent(int iValue) const;
        Float_t GetEtaBinContent(int iValue) const;
        int GetPtBinSize() const;
        int GetEtaBinSize() const;

        void SetPtBinContent(Float_t fValue);
        void SetPtBins(Float_t *fValue, int dim);
        void SetEtaBinContent(Float_t fValue);
        void SetEtaBins(Float_t *fValues, int dim);

    private:
        bool SetPtId(int id);
        bool SetEtaId(int id);
        void NextPt();
        void NextEta();
        
        bool checkBinsArray(Float_t *fValue, int dim);

        int NptBins;
        int NetaBins;
        Float_t ptBins[_MAX_PT_BINS];
        Float_t etaBins[_MAX_ETA_BINS];
        
        

};

#endif