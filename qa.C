#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TFile.h"
#include "TObject.h"
#include "TTree.h"
#include "TBranch.h"
#include "TSystem.h"
#include "TChain.h"
#include "TVector2.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TProfile.h"
#include "TString.h"
#include "THnSparse.h"
#include <iostream>
#define Pi 3.14159265
#define twoPi 6.28318530


class StMyElectronEvent;
class StRefMultCorr;

StMyElectronEvent *mElectronEvent;
class StMyElectron;
StMyElectron *mElectron;
//event info
TH1F *hMcVertexZ = new TH1F("hmcVertexZ", "mcVertexZ;Vz^{mc} (cm)", 400, -200, 200);
TH2F *hMcVertexXY = new TH2F("hmcVertexXY", "mcVertexXY;Vx^{mc} (cm);Vy^{mc} (cm)", 40, -2, 2, 40, -2, 2);
TH1F *hRcVertexZ = new TH1F("hrcVertexZ", "rcVertexZ;Vz^{rc} (cm)", 400, -200, 200);
TH1F *hVertexZdiff = new TH1F("hvertexZdiff", "vertexZdiff;(Vz^{rc}-Vz^{mc} (cm)", 100, -5, 5);
TH1F *hRefMult = new TH1F ("hrefMult", "refMult;Reference Multiplicity", 1000, 0, 1000);
TH3F *hnMce = new TH3F("hnMce", "", 10, -0.5, 9.5, 5000, 0, 5000., 100, 0., 1.);
//mc electron info
Int_t hMcptetaphi_bin[4] = {10, 200, 40, 100};
Double_t hMcptetaphi_xmin[4] = { -0.5, 0., -2, -Pi};
Double_t hMcptetaphi_xmax[4] = {9.5, 20., 2, Pi};
THnSparse *hMcptetaphi = new THnSparseF("hMcptetaphi", "sign:cent:pt:eta:phi", 4, hMcptetaphi_bin, hMcptetaphi_xmin, hMcptetaphi_xmax);
hMcptetaphi->Sumw2();

//reconstruct electron info without cut
Int_t hRcptetaphi_bin[4] = {10, 200, 40, 100};
Double_t hRcptetaphi_xmin[4] = { -0.5, 0., -2, -Pi};
Double_t hRcptetaphi_xmax[4] = {9.5, 20., 2, Pi};
THnSparse *hRcptetaphi = new THnSparseF("hRcptetaphi", "sign:cent:pt:eta:phi", 4, hRcptetaphi_bin, hRcptetaphi_xmin, hRcptetaphi_xmax);
hRcptetaphi->Sumw2();

//reconstruct electron info with TPC cut
Int_t htpcptetaphi_bin[4] = {10, 200, 40, 100};
Double_t htpcptetaphi_xmin[4] = { -0.5, 0., -2, -Pi};
Double_t htpcptetaphi_xmax[4] = {9.5, 20., 2, Pi};
THnSparse *htpcptetaphi = new THnSparseF("htpcptetaphi", "sign:cent:pt:eta:phi", 4, htpcptetaphi_bin, htpcptetaphi_xmin, htpcptetaphi_xmax);
htpcptetaphi->Sumw2();

//reconstruct electron info with TPC cut && BEMC match
Int_t hbemcptetaphi_bin[4] = {10, 200, 40, 100};
Double_t hbemcptetaphi_xmin[4] = { -0.5, 0., -2, -Pi};
Double_t hbemcptetaphi_xmax[4] = {9.5, 20., 2, Pi};
THnSparse *hbemcptetaphi = new THnSparseF("hbemcptetaphi", "sign:cent:pt:eta:phi", 4, hbemcptetaphi_bin, hbemcptetaphi_xmin, hbemcptetaphi_xmax);
hbemcptetaphi->Sumw2();

//reconstruct electron info with TPC cut && BEMC match && poE cut
Int_t hpoEptetaphi_bin[4] = {10, 200, 40, 100};
Double_t hpoEptetaphi_xmin[4] = { -0.5, 0., -2, -Pi};
Double_t hpoEptetaphi_xmax[4] = {9.5, 20., 2, Pi};
THnSparse *hpoEptetaphi = new THnSparseF("hpoEptetaphi", "sign:cent:pt:eta:phi", 4, hpoEptetaphi_bin, hpoEptetaphi_xmin, hpoEptetaphi_xmax);
hpoEptetaphi->Sumw2();

//reconstruct electron info with TPC cut && BEMC match && poE cut && SMD cut
Int_t hsmdptetaphi_bin[4] = {10, 200, 40, 100};
Double_t hsmdptetaphi_xmin[4] = { -0.5, 0., -2, -Pi};
Double_t hsmdptetaphi_xmax[4] = {9.5, 20., 2, Pi};
THnSparse *hsmdptetaphi = new THnSparseF("hsmdptetaphi", "sign:cent:pt:eta:phi", 4, hsmdptetaphi_bin, hsmdptetaphi_xmin, hsmdptetaphi_xmax);
hsmdptetaphi->Sumw2();


//

Int_t hcutptetaphi_bin[4] = {10, 200, 40, 100};
Double_t hcutptetaphi_xmin[4] = { -0.5, 0., -2, -Pi};
Double_t hcutptetaphi_xmax[4] = {9.5, 20., 2, Pi};
THnSparse *hcutptetaphi = new THnSparseF("hcutptetaphi", "sign:cent:pt:eta:phi", 4, hcutptetaphi_bin, hcutptetaphi_xmin, hcutptetaphi_xmax);
hcutptetaphi->Sumw2();


TH3F *hRcdca          = new TH3F("hRcdca"        , "", 10, -0.5, 9.5, 200, 0, 20, 50,  0.0,  5.0);
TH3F *hRcnHitsRatio     = new TH3F("hRcnHitsRatio"   , "", 10, -0.5, 9.5, 200, 0, 20, 100, 0, 1);
TH3F *hRcnHitsFit     = new TH3F("hRcnHitsFit"   , "", 10, -0.5, 9.5, 200, 0, 20, 60, -0.5, 59.5);
TH3F *hRcnHitsDedx    = new TH3F("hRcnHitsDedx"  , "", 10, -0.5, 9.5, 200, 0, 20, 60, -0.5, 59.5);
TH3F *hRcnSigE        = new TH3F("hRcnSigE"      , "", 10, -0.5, 9.5, 200, 0, 20, 60, -3.0,  3.0);
//TH3F *hRctpce         = new TH3F("hRctpce"       , "", 10, -0.5, 9.5, 200, 0, 20, );
//TH3F *hRctofmatch_e   = new TH3F("hRctofmatch_e" , "", 10, -0.5, 9.5, 200, 0, 20,);
TH3F *hRcbemcmatch    = new TH3F("hRcbemcmatch"  , "", 10, -0.5, 9.5, 200, 0, 20, 100, 0.0, 10.0);
TH3F *hRcneta         = new TH3F("hRcneta"       , "", 10, -0.5, 9.5, 200, 0, 20, 20, -0.5, 19.5);
TH3F *hRcnphi         = new TH3F("hRcnphi"       , "", 10, -0.5, 9.5, 200, 0, 20, 20, -0.5, 19.5);
TH3F *hRcphiDist      = new TH3F("hRcphiDist"    , "", 10, -0.5, 9.5, 200, 0, 20, 200, -0.1, 0.1);
TH3F *hRczDist        = new TH3F("hRczDist"      , "", 10, -0.5, 9.5, 200, 0, 20, 200, -10.0, 10.0);

//reconstruct electron info with cut
TH3F *hdca          = new TH3F("hdca"        , "", 10, -0.5, 9.5, 200, 0, 20, 50,  0.0,  5.0);
TH3F *hnHitsRatio   = new TH3F("hnHitsRatio"   , "", 10, -0.5, 9.5, 200, 0, 20, 100, 0, 1);
TH3F *hnHitsFit     = new TH3F("hnHitsFit"   , "", 10, -0.5, 9.5, 200, 0, 20, 60, -0.5, 59.5);
TH3F *hnHitsDedx    = new TH3F("hnHitsDedx"  , "", 10, -0.5, 9.5, 200, 0, 20, 60, -0.5, 59.5);
TH3F *hnSigE        = new TH3F("hnSigE"      , "", 10, -0.5, 9.5, 200, 0, 20, 60, -3.0,  3.0);
//TH3F *htpce         = new TH3F("htpce"       , "", 10, -0.5, 9.5, 200, 0, 20, );
//TH3F *htofmatch_e   = new TH3F("htofmatch_e" , "", 10, -0.5, 9.5, 200, 0, 20,);
TH3F *hbemcmatch    = new TH3F("hbemcmatch"  , "", 10, -0.5, 9.5, 200, 0, 20, 100, 0.0, 10.0);
TH3F *hneta         = new TH3F("hneta"       , "", 10, -0.5, 9.5, 200, 0, 20, 20, -0.5, 19.5);
TH3F *hnphi         = new TH3F("hnphi"       , "", 10, -0.5, 9.5, 200, 0, 20, 20, -0.5, 19.5);
TH3F *hphiDist      = new TH3F("hphiDist"    , "", 10, -0.5, 9.5, 200, 0, 20, 200, -0.1, 0.1);
TH3F *hzDist        = new TH3F("hzDist"      , "", 10, -0.5, 9.5, 200, 0, 20, 200, -10.0, 10.0);


double mCentWeight = -999.0;


void qa(const char *fileInDir = "./outMB_Electron_100", const char *OutFile = "jpsi_qa_eletron_100")
{
    gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
    loadSharedLibraries();
    gSystem->Load("StMyElectronMaker");
    gSystem->Load("StRefMultCorr");
    const Int_t nFilesMax = 10000;
    const Double_t emcScale = 1.00;
    gRandom = new TRandom3();
    mElectronEvent = new StMyElectronEvent();
    TChain chMc("mcT");
    chMc.SetBranchAddress("mcE", &mElectronEvent);
    void *dir = gSystem->OpenDirectory(gSystem->ExpandPathName(fileInDir));
    int nruns = 0;
    char *file_name;
    TString Tname;
    char file_list[nFilesMax][256];

    do {
        file_name = (char *)gSystem->GetDirEntry(dir);
        Tname = file_name;

        if (file_name && Tname.Contains("myminimc.root") && Tname.Contains("emb")) {
            sprintf(file_list[nruns], "%s/%s", fileInDir, file_name);
            chMc.Add(file_list[nruns]);
            cout << " read in " << file_list[nruns] << endl;
            nruns++;
        }
    }
    while (file_name && nruns <= nFilesMax);

    int iret = 0;
    int nb = 0;
    cout << chMc.GetEntries() << " events in chain" << endl;
    int nevents = chMc.GetEntries();

    for (int i = 0; i < nevents; i++) {
        nb += chMc.GetEvent(i);
        chMc.GetEvent(i);

        if (i % 10000 == 0) cout << "event# " << i << endl;

        int n;

        if (!mElectronEvent) continue;

        if (mElectronEvent->eventID() <= 0) continue;

        int Mult = mElectronEvent->refMultPos() + mElectronEvent->refMultNeg();

        if (Mult <= 0) continue;

        hRefMult->Fill(Mult);
        //cout<<"vertexZ = "<<mElectronEvent->vertexZ()<<endl;
        int mCent = getCentrality(Mult, mElectronEvent->vertexZ(), mElectronEvent->RunId());

        if ((mCent < 9) && (mCent > -1)) {

            hMcVertexZ->Fill(mElectronEvent->mcVertexZ());
            hMcVertexXY->Fill(mElectronEvent->mcVertexX(), mElectronEvent->mcVertexY());
            hRcVertexZ->Fill(mElectronEvent->vertexZ());
            hVertexZdiff->Fill(mElectronEvent->vertexZ() - mElectronEvent->mcVertexZ());

        }

        int nMCelectron = mElectronEvent->nReal();
        double nMCe = (double)nMCelectron / (Double_t)Mult;
        hnMce->Fill(mCent, Mult, nMCe);
        double MCelectron[4] = {0., 0., 0., 0.};
        double RCelectron[4] = {0., 0., 0., 0.};
        Double_t cutelectron[4] = {0., 0., 0.0, 0.};

        for (int j = 0; j < mElectronEvent->nReal(); j++) {
            mElectron = (StMyElectron) mElectronEvent->real()->UncheckedAt(j);

            //cout<<j<<"  "<<mElectron->geantId<<"  "<<mElectron->pGeantId<<endl;
            if (mElectron->pGeantId > 0) continue;  //not original electron

            if (mElectron->mcId < 0) continue;

            Double_t dsmAdc0 = mElectron->dsmAdc0;
            Double_t adc0 = mElectron->adc0;
            Double_t energy = mElectron->energy;
            dsmAdc0 = int((dsmAdc0 + gRandom->Uniform(0, 1. - 1.e-9)) * emcScale);
            adc0 *= emcScale;
            energy *= emcScale;


            if (mElectron->mcId >= 0) {
                MCelectron[0] = mCent;
                MCelectron[1] = mElectron->mcPt;
                MCelectron[2] = mElectron->mcEta;
                MCelectron[3] = mElectron->mcPhi;

                if (MCelectron[3] > Pi) MCelectron[3] = MCelectron[3] - twoPi;

                hMcptetaphi->Fill(MCelectron);
            }

            if (mElectron->mcId >= 0 && mElectron->id >= 0) {
                //cout<<mElectron->id<<endl;
                RCelectron[0] = mCent;
                RCelectron[1] = mElectron->mcPt;
                RCelectron[2] = mElectron->mcEta;
                RCelectron[3] = mElectron->mcPhi;

                if (RCelectron[3] > Pi) RCelectron[3] = RCelectron[3] - twoPi;

                hRcptetaphi->Fill(RCelectron);
                double pt = mElectron->pt;

                if (fabs(mElectron->eta) > 1)  continue;

                //                if (pt < 0.2) continue;

                //hRcdca->Fill(mCent, pt, mElectron->dca);
                hRcnHitsFit->Fill(mCent, pt, mElectron->nFitPts);
                hRcnHitsDedx->Fill(mCent, pt, mElectron->nDedxPts);
                Float_t ratio = (Double_t)(mElectron->nFitPts) / (Double_t)(mElectron->nMaxPts);
                hRcnHitsRatio->Fill(mCent, pt, ratio);
                Float_t p = mElectron->p;
                Float_t poE = (energy > 0) ? p / energy : -999;
                hRcnSigE->Fill(mCent, pt, mElectron->nSigE);
                hRcbemcmatch->Fill(mCent, pt, poE);
                hRcneta->Fill(mCent, pt, mElectron->nEta);
                hRcnphi->Fill(mCent, pt, mElectron->nPhi);
                hRcphiDist->Fill(mCent, pt, mElectron->phiDist);
                hRczDist ->Fill(mCent, pt, mElectron->zDist);

                hdca->Fill(mCent, pt, mElectron->dca);

                if (mElectron->dca > 2) continue;

                if ((mElectron->nFitPts) > 15) {
                    hnHitsRatio->Fill(mCent, pt, ratio);
                }

                if (ratio < 0.52) continue;

                hnHitsFit->Fill(mCent, pt, mElectron->nFitPts);

                if (mCent < 6) {
                    if ((mElectron->nFitPts ) < 20) continue;

                    hnHitsDedx->Fill(mCent, pt, mElectron->nDedxPts );
                    htpcptetaphi->Fill(RCelectron);

                    //                   if ((mElectron->nDedxPts ) < 12) continue;
                }

                else {
                    if ((mElectron->nFitPts )  < 25) continue;

                    hnHitsDedx->Fill(mCent, pt, mElectron->nDedxPts );
                    htpcptetaphi->Fill(RCelectron);

                    //                 if ((mElectron->nDedxPts ) < 15) continue;
                }

                hRcnSigE->Fill(mCent, pt, mElectron->nSigE);

                // double nse_p = p < 1 ? (1.5 * p - 2.5) : -1;
                //
                // if ((mElectron->nSigE < nse_p) || (mElectron->nSigE > 2)) continue;
                if (poE < 0) continue;

                hbemcptetaphi->Fill(RCelectron);
                hbemcmatch->Fill(mCent, pt, poE);

                if ((poE > 0.3 ) && (poE < 1.5)) {
                    hpoEptetaphi->Fill(RCelectron);
                    int nEta = mElectron->nEta;
                    int nPhi = mElectron->nPhi;
                    nEta = nEta > 0 ? nEta : 0;
                    nPhi = nPhi > 0 ? nPhi : 0;
                    hneta->Fill(mCent, pt, nEta);
                    hnphi->Fill(mCent, pt, nPhi);

                    if (((fabs(mElectron->zDist) < 2) && (nEta > 0)) || ((fabs(mElectron->phiDist) < 0.01) && (nPhi > 0))) {
                        hcutptetaphi->Fill(RCelectron);
                    }

                    if (nEta >= 1) {
                        hzDist ->Fill(mCent, pt, mElectron->zDist);

//                        if (fabs(mElectron->zDist) < 2) {
//
//                        }
                    }

                    if (nPhi >= 1) {
                        hphiDist->Fill(mCent, pt, mElectron->phiDist);

//                        if (fabs(mElectron->phiDist) < 0.01) {
//                        }
                    }//#SMD strips
                }//PE
            }
        }//end of RC electron loop
    }

    char buf[1024];
    sprintf(buf, "%s.root", OutFile);
    TFile *f = new TFile(buf, "recreate");
    f->cd();
    hMcVertexZ->Write();
    hMcVertexXY->Write();
    hRcVertexZ->Write();
    hVertexZdiff->Write();
    hRefMult->Write();
    hnMce->Write();
    hRcdca->Write();
    hRcnHitsRatio->Write();
    hRcnHitsFit->Write();
    hRcnHitsDedx->Write();
    hRcnSigE->Write();
    hRcbemcmatch->Write();
    hRcneta->Write();
    hRcnphi->Write();
    hRcphiDist->Write();
    hRczDist->Write();
    hdca->Write();
    hnHitsRatio->Write();
    hnHitsFit->Write();
    hnHitsDedx->Write();
    hnSigE->Write();
    hbemcmatch->Write();
    hneta->Write();
    hnphi->Write();
    hphiDist->Write();
    hzDist->Write();
    hMcptetaphi->Write();
    hRcptetaphi->Write();
    htpcptetaphi->Write();
    hbemcptetaphi->Write();
    hpoEptetaphi->Write();

    hcutptetaphi->Write();
    f->Close();
}

//=======================================
int getCentrality(int refmult, float vz, int runnumber)
{
    int mCentrality = -1;
    StRefMultCorr *refmultCorrUtil = new StRefMultCorr();
    refmultCorrUtil->init(runnumber);
    refmultCorrUtil->initEvent(refmult, vz) ;
    int cent = refmultCorrUtil->getCentralityBin9() ;
    mCentWeight = refmultCorrUtil->getWeight();

    if (mCentrality == -1) mCentrality = 9;

    delete refmultCorrUtil;
    return cent;
}

