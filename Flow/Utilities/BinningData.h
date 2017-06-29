#ifndef BINNING_DATA_H

#define BINNING_DATA_H

#define _MAX_PT_BINS 100
#define _MAX_ETA_BINS 100
#define _MAX_RAPIDITY_BINS 100

class BinningData {
    public:
        BinningData();
        ~BinningData();
        int GetPtBin(Float_t fValue) const;
        int GetEtaBin(Float_t fValue) const;
        int GetRapidityBin(Float_t fValue) const;
        Float_t GetPtBinContent(int iValue) const;
        Float_t GetEtaBinContent(int iValue) const;
        Float_t GetRapidityBinContent(int iValue) const;
        int GetPtBinSize() const;
        int GetEtaBinSize() const;
        int GetRapidityBinSize() const;

        void SetPtBinContent(Float_t fValue);
        void SetPtBins(Float_t *fValue, int dim);
        void SetPtBins(const Float_t *fValue, const int dim);
        void SetEtaBinContent(Float_t fValue);
        void SetEtaBins(Float_t *fValues, int dim);
        void SetEtaBins(const Float_t *fValues, const int dim);
        void SetRapidityBinContent(Float_t fValue);
        void SetRapidityBins(Float_t *fValues, int dim);
        void SetRapidityBins(const Float_t *fValues, const int dim);

    private:
        bool SetPtId(int id);
        bool SetEtaId(int id);
        bool SetRapidityId(int id);
        void NextPt();
        void NextEta();
        void NextRapidity();
        
        static bool checkBinsArray(Float_t *fValue, int dim);
        static bool checkBinsArray(const Float_t *fValue, const int dim);

        int NptBins;
        int NetaBins;
        int NrapidityBins;
        Float_t ptBins[_MAX_PT_BINS];
        Float_t etaBins[_MAX_ETA_BINS];
        Float_t rapidityBins[_MAX_RAPIDITY_BINS];
        

};

#endif