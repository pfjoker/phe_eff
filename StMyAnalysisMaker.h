#ifndef StMyAnalysisMaker_h
#define StMyAnalysisMaker_h

#include "StMaker.h"
#include "StarClassLibrary/SystemOfUnits.h"
#include "StThreeVectorF.hh"
#include "TString.h"
#include "TLorentzVector.h"
#include "string"
#include "THnSparse.h"

#define Mpion 0.1396
#define Mkaon 0.4937
#define Melectron 0.000511
const int nVzBin = 20;
const int nCenBin = 10;
const int nMaxEventsInBuffer = 20;
const int nMaxPions = 1000;
const int nMaxElectronP = 2000;
const int nMaxElectronG = 200;
class StPicoDst;
class StPicoDstMaker;
class TString;
class StRefMultCorr;
class TFile;
class TTree;
class TH1F;
class TH2F;
class TH3F;
class TH2D;
class TH3D;
class THnSparse;
class TLorentzVector;
class TProfile2D;
class StMyAnalysisMaker : public StMaker
{
public:
    StMyAnalysisMaker(const char *name, StPicoDstMaker *picoMaker, const char *outName);
    virtual ~StMyAnalysisMaker();

    virtual Int_t Init();
    virtual Int_t Make();
    virtual void  Clear(Option_t *opt = "");
    virtual Int_t Finish();
    TLorentzVector boost(const TLorentzVector& x, const TLorentzVector& y);

    void setHistoFileName(TString name) {mHistogramOutputFileName = name;} // Make name available to member functions

    void    WriteHistograms();
    Int_t   GetCentrality();

private:
    StPicoDstMaker *mPicoDstMaker;
    StPicoDst      *mPicoDst;
    StRefMultCorr* refmultCorrUtil;
    TString    mOutName;
    Int_t      mCent1;    //bin9
    Int_t      mCent2;    //bin16
    Double_t   mCentWeight;
    int Current_nCanE;
    int Current_nCanGE;
    int Current_CanEId[nMaxElectronP];
    int Current_CanGEId[nMaxElectronG];
    int Current_CanECharge[nMaxElectronP];
    int Current_CanGECharge[nMaxElectronG];
    float Current_CanEPt[nMaxElectronP];
    float Current_CanEP[nMaxElectronG];        float Current_CanEEta[nMaxElectronP];
    float Current_CanEPhi[nMaxElectronP];
    float Current_CanGEPt[nMaxElectronG];
    float Current_CanGEEta[nMaxElectronG];
    float Current_CanGEPhi[nMaxElectronG];
    StThreeVectorF Current_CanEMom[nMaxElectronP];
    StThreeVectorF Current_CanEOrigin[nMaxElectronP];
    StThreeVectorF Current_CanGEMom[nMaxElectronG];
    StThreeVectorF Current_CanGEOrigin[nMaxElectronG];
    //float Current_CanEHelixPx[nMaxElectronP];
    //float Current_CanEHelixPy[nMaxElectronP];
    //float Current_CanEHelixPz[nMaxElectronP];
    // float Current_CanEHelixX[nMaxElectronP];
    //float Current_CanEHelixY[nMaxElectronP];
    //float Current_CanEHelixZ[nMaxElectronP];
    //float Current_CanGEHelixPx[nMaxElectronG];
    //float Current_CanGEHelixPy[nMaxElectronG];
    //float Current_CanGEHelixPz[nMaxElectronG];
    //float Current_CanGEHelixX[nMaxElectronG];
    //float Current_CanGEHelixY[nMaxElectronG];
    //float Current_CanGEHelixZ[nMaxElectronG];
    float Current_CanEDca[nMaxElectronP];
    float Current_CanEnHitsFit[nMaxElectronP];
    float Current_CanEnHitsDedx[nMaxElectronP];
    float Current_CanEnHitsMax[nMaxElectronP];
    float Current_CanEnSigE[nMaxElectronP];
    float Current_CanEBeta[nMaxElectronP];
    float Current_CanEE0[nMaxElectronP];
    float Current_CanEDistZ[nMaxElectronP];
    float Current_CanEDistPhi[nMaxElectronP];
    float Current_CanEnEta[nMaxElectronP];
    float Current_CanEnPhi[nMaxElectronP];
    double Current_Vx;
    double Current_Vy;
    double Current_Vz;
    double Current_bField;

    int current_q[2000];
    TLorentzVector current_pion[2000];
    int current_kq[2000];
    TLorentzVector current_kaon[2000];
    float PI;
    float twoPI;

    float  cutAbsVertexZ;
    int    cutTriggerId;

    float mDcaCut[2];
    int mHitsFitCut[2];
    float mFitRatioCut[2];
    float mEtaCut[2];

    TH1F*         hVertexZAll;
    TH1F*         hVzDiff;
    TH3F*         hVertexXYZ;
    TH2F*         hVzvsVpd;
    TH2F*         hCent_ref;
    TH2F*         hCent_ref_un;
    TH3F*         hDenP;
    TH3F*         hDenM;
    TH3F*         hNum;

    TH2F*         hOrigin_xy;
    TH2F*         hOrigin_xz;
    TH2F*         hOrigin_yz;
    TH3F*         hDecayL;
    TH3F*         hDcaV0;
    TH3F*         hpairmass;
    TH3F*         hpairdca;
    THnSparse   *hnHitsFit;
    THnSparse   *hnHitsDedx;
    THnSparse   *hnSigE;
    THnSparse   *htofmatch;
    THnSparse   *hbemcmatch;
    THnSparse   *hneta;
    THnSparse   *hnphi;
    THnSparse   *hphiDist;
    THnSparse   *hzDist;

    TString       mHistogramOutputFileName ;
    TFile*        histogram_output;

    UInt_t        mEventsProcessed ;

    Int_t   mRefMult;
    Int_t   cent;
    void  initCutParams ( ) ;
    void  initConst ( ) ;
    void  initHisto ( ) ;
    int makephe();

    ClassDef(StMyAnalysisMaker, 1)
};

#endif
