#include "StMyAnalysisMaker.h"
#include "StRoot/StPicoDstMaker/StPicoDst.h"
#include "StRoot/StPicoDstMaker/StPicoEvent.h"
#include "StRoot/StPicoDstMaker/StPicoTrack.h"
#include "StRoot/StPicoDstMaker/StPicoV0.h"
#include "StRoot/StPicoDstMaker/StPicoDstMaker.h"
#include "StThreeVectorF.hh"
#include "StThreeVectorD.hh"
#include "StLorentzVectorD.hh"
#include "TH1F.h"
#include "TH2F.h"
#include "TRandom.h"
#include "TRandom3.h"

#include "TFile.h"
#include "StRoot/StRefMultCorr/StRefMultCorr.h"
#include <iostream>

#include "StMessMgr.h"
#include "PhysicalConstants.h"
#include "TProfile2D.h"
#include "TVector3.h"
#include "TMath.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TFile.h"
#include "TString.h"
#include "string"
#include "TTree.h"
#include "StPhysicalHelixD.hh"
#include "SystemOfUnits.h"

ClassImp(StMyAnalysisMaker)

//-----------------------------------------------------------------------------
StMyAnalysisMaker::StMyAnalysisMaker(const char *name, StPicoDstMaker *picoMaker, const char *outName)
	: StMaker(name)
{
	mPicoDstMaker = picoMaker;
	mPicoDst = 0;
	refmultCorrUtil = new StRefMultCorr();
	histogram_output = NULL  ;
	mEventsProcessed = 0     ;
	mHistogramOutputFileName = outName;
}
TLorentzVector StMyAnalysisMaker::boost(const TLorentzVector &x, const TLorentzVector &y)
{
	Double_t mass = y.M();
	TVector3 eta = (-1. / mass) * y.Vect();
	Double_t gamma = y.E() / mass;
	TVector3 pl = (x.Vect() * eta) / (eta * eta) * eta;
	TLorentzVector c(x.Vect() + (gamma - 1.)*pl - x.E()*eta, gamma * x.E() - x.Vect()*eta);
	return c;
};

//-----------------------------------------------------------------------------
StMyAnalysisMaker::~StMyAnalysisMaker()
{
	delete refmultCorrUtil;
}

//-----------------------------------------------------------------------------
void StMyAnalysisMaker::initCutParams()
{

	cutAbsVertexZ  = 30.;
	cutTriggerId = 0;

	mDcaCut[0]      = 0.0;    mDcaCut[1]      = 2.0;
	mHitsFitCut[0]  = 15;     mHitsFitCut[1]  = 50;
	mFitRatioCut[0] = 0.52;   mFitRatioCut[1] = 1.05;
	mEtaCut[0]      = -1.0;   mEtaCut[1]      = 1.0;

	return;
}

//-----------------------------------------------------------------------------
void StMyAnalysisMaker::initHisto()
{
	TH1::SetDefaultSumw2();
	hVertexZAll  = new TH1F( "hVertexZAll", "All Event Vertex Z Position", 500, -250.0, 250.0);
	hVzDiff  = new TH1F( "hVzDiff", "VpdVz-TPCVz", 2000, -10.0, 10.0);
	hVzvsVpd  = new TH2F( "hVzvsVpd", "VpdVz:TPCVz", 2000, -50.0, 50.0, 2000, -50.0, 50.0);
	hVertexXYZ = new TH3F( "hVertexXYZ", "Vertex Position", 100, -5.0, 5.0 , 100, -5., 5, 1000, -50, 50);
	//hVertexXY2 = new TH2F( "VertexXY2", "Vertex XY Position(Rxy<2)", 200, -10.0, 10.0, 200, -10., 10);
	//refmult for correction
	hCent_ref_un = new TH2F("hCent_ref_un", "Cent:RefMult", 10, -0.5, 9.5, 3000, 0, 3000);
	//refult after correction
	//hCent_ref = new TH2F("hCent_ref", "Cent:RefMult", 9, 0, 9., 3000, 0, 3000);

	hNum =  new TH3F("hNum", "Centrality:Pt:M", 10, -0.5, 9.5, 500, 0, 5, 200, 0, 0.2);
	hDenP =  new TH3F("hDenP", "Centrality:Pt:M", 10, -0.5, 9.5, 500, 0, 5, 200, 0.0, 0.2);
	hDenM =  new TH3F("hDenM", "Centrality:Pt:M", 10, -0.5, 9.5, 500, 0, 5, 200, 0.0, 0.2);



	hOrigin_xy = new TH2F("hOrigin_xy", "", 1000, -100, 100, 1000, -100, 100);
	hOrigin_yz = new TH2F("hOrigin_yz", "", 1000, -100, 100, 1000, -100, 100);
	hOrigin_xz = new TH2F("hOrigin_xz", "", 1000, -100, 100, 1000, -100, 100);
	//hDecayL = new TH3F("hDecayL", "Centrality:P:decayL", 10, -0.5, 9.5, 100, 0, 10, 1000, -100, 100);
	hpairmass = new TH3F("hpairmass", "cent:pt:mass", 10, -0.5, 9.5, 100, 0, 10, 500, 0, 5);
	hpairdca = new TH3F("hpairdca", "cent:pt:dca", 10, -0.5, 9.5, 100, 0, 10, 100, 0, 10);

	Int_t hnHitsFit_bin[4] = {3, 10, 100, 60};
	Double_t hnHitsFit_xmin[4] = { -1.5, -0.5, 0., -0.5};
	Double_t hnHitsFit_xmax[4] = {1.5, 9.5, 10., 59.5};
	hnHitsFit = new THnSparseF("hnHitsFit", "sign:cent:pt:nHitsFit", 4, hnHitsFit_bin, hnHitsFit_xmin, hnHitsFit_xmax);
	hnHitsFit->Sumw2();


	Int_t hnHitsDedx_bin[4] = {3, 10, 100, 60};
	Double_t hnHitsDedx_xmin[4] = { -1.5, -0.5, 0., -0.5};
	Double_t hnHitsDedx_xmax[4] = {1.5, 9.5, 10., 59.5};
	hnHitsDedx = new THnSparseF("hnHitsDedx", "sign:cent:pt:nHitsDedx", 4, hnHitsDedx_bin, hnHitsDedx_xmin, hnHitsDedx_xmax);
	hnHitsDedx->Sumw2();

	Int_t hnSigE_bin[4] = {3, 10, 100, 80};
	Double_t hnSigE_xmin[4] = { -1.5, -0.5, 0., -4};
	Double_t hnSigE_xmax[4] = {1.5, 9.5, 10., 4};
	hnSigE = new THnSparseF("hnSigE", "sign:cent:pt:nSigE", 4, hnSigE_bin, hnSigE_xmin, hnSigE_xmax);
	hnSigE->Sumw2();


	Int_t htofmatch_bin[4] = {3, 10, 100, 1000};
	Double_t htofmatch_xmin[4] = { -1.5, -0.5, 0., 0.};
	Double_t htofmatch_xmax[4] = {1.5, 9.5, 10., 10.};
	htofmatch = new THnSparseF("htofmatch", "sign:cent:pt:ovb", 4, htofmatch_bin, htofmatch_xmin, htofmatch_xmax);
	htofmatch->Sumw2();
	Int_t hbemcmatch_bin[4] = {3, 10, 100, 1000};
	Double_t hbemcmatch_xmin[4] = { -1.5, -0.5, 0., 0.};
	Double_t hbemcmatch_xmax[4] = {1.5, 9.5, 10., 10.};
	hbemcmatch = new THnSparseF("hbemcmatch", "sign:cent:pt:pve", 4, hbemcmatch_bin, hbemcmatch_xmin, hbemcmatch_xmax);
	hbemcmatch->Sumw2();

	Int_t hnphi_bin[4] = {3, 10, 100, 20};
	Double_t hnphi_xmin[4] = { -1.5, -0.5, 0., -0.5};
	Double_t hnphi_xmax[4] = {1.5, 9.5, 10., 19.5};
	hnphi = new THnSparseF("hnphi", "sign:cent:pt:nphi", 4, hnphi_bin, hnphi_xmin, hnphi_xmax);
	hnphi->Sumw2();
	Int_t hneta_bin[4] = {3, 10, 100, 20};
	Double_t hneta_xmin[4] = { -1.5, -0.5, 0., -0.5};
	Double_t hneta_xmax[4] = {1.5, 9.5, 10., 19.5};
	hneta = new THnSparseF("hneta", "sign:cent:pt:neta", 4, hneta_bin, hneta_xmin, hneta_xmax);
	hneta->Sumw2();

	Int_t hphiDist_bin[4] = {3, 10, 100, 200};
	Double_t hphiDist_xmin[4] = { -1.5, -0.5, 0., -0.1};
	Double_t hphiDist_xmax[4] = {1.5, 9.5, 10., 0.1};
	hphiDist = new THnSparseF("hphiDist", "sign:cent:pt:phiDist", 4, hphiDist_bin, hphiDist_xmin, hphiDist_xmax);
	hphiDist->Sumw2();
	Int_t hzDist_bin[4] = {3, 10, 100, 200};
	Double_t hzDist_xmin[4] = { -1.5, -0.5, 0., 10.};
	Double_t hzDist_xmax[4] = {1.5, 9.5, 10., 10.};
	hzDist = new THnSparseF("hzDist", "sign:cent:pt:zDist", 4, hzDist_bin, hzDist_xmin, hzDist_xmax);
	hzDist->Sumw2();

	return;
}
//-----------------------------------------------------------------------------
void StMyAnalysisMaker::initConst()
{

	PI = 3.1415926;
	twoPI = 6.2831852;

	return;
}
//-----------------------------------------------------------------------------


Int_t StMyAnalysisMaker::Init()
{
	initConst();

	initCutParams();  // initialize cuts  parameters

	if (mHistogramOutputFileName == "")        //FIXME: ALWAYS USE { } HERE!!! LOG_XXX is a if()xxx macro!!! it sucks!
	{
		LOG_ERROR << "Please input output histrogram file" << endm;
	}

	else
	{
		histogram_output = new TFile( mHistogramOutputFileName, "recreate" ) ;
	}

	initHisto();

	return kStOK;
}
//-----------------------------------------------------------------------------
Int_t StMyAnalysisMaker::Finish()
{
	if (histogram_output != NULL)
	{
		histogram_output ->cd();
		histogram_output -> Write() ; // Write all histograms to disk
		hnHitsFit->Write();
		hnHitsDedx->Write();
		hnSigE->Write();
		htofmatch->Write();
		hbemcmatch->Write();
		hnphi->Write();
		hneta->Write();
		hphiDist->Write();
		hzDist->Write();

	}

	cout << "Total Events Processed in DstMaker " << mEventsProcessed << endl ;

	return kStOK;
}


//-----------------------------------------------------------------------------
void StMyAnalysisMaker::WriteHistograms()
{
}
//----------------------------------------------------------------------------
Int_t StMyAnalysisMaker::GetCentrality()
{
	StPicoEvent *event = (StPicoEvent *)mPicoDst->event();

	if (!event) return 0;

	Int_t RunNum = event->runId();
	int badrun = refmultCorrUtil->isBadRun(RunNum);

	if (badrun)  return -1;

	refmultCorrUtil->init(RunNum);
	Int_t refmul_for_corr  = event->refMult();
	Float_t vertexZ_for_corr = event->primaryVertex().z();
	Float_t zdcx_for_corr    = event->ZDCx();
	refmultCorrUtil->initEvent( refmul_for_corr, vertexZ_for_corr, zdcx_for_corr);
	mCent1             = refmultCorrUtil->getCentralityBin9();
	mCent2             = refmultCorrUtil->getCentralityBin16();
	mCentWeight             = refmultCorrUtil->getWeight();
	return 1;
}
//-----------------------------------------------------------------------------
void StMyAnalysisMaker::Clear(Option_t *opt)
{
}

//-----------------------------------------------------------------------------
Int_t StMyAnalysisMaker::Make()
{
	if (!mPicoDstMaker)
	{
		LOG_WARN << " No PicoDstMaker! Skip! " << endm;
		return kStWarn;
	}

	mPicoDst = mPicoDstMaker->picoDst();

	if (!mPicoDst)
	{
		LOG_WARN << " No PicoDst! Skip! " << endm;
		return kStWarn;
	}

	if (GetCentrality() != 1)
	{
		LOG_WARN << " Bad Run! Skip! " << endm;
		return kStWarn;
	}

	StPicoEvent *event = (StPicoEvent *)mPicoDst->event();

	if (!(event->isMinBias()))
	{
		return kStOK;
	}

	int    cent             = mCent1;

	if ( cent < 0 || cent > 8) return kStOk;

	double refmul_for_corr  = event->refMult();
	double vpdVz = event->vzVpd();
	Current_Vx = event->primaryVertex().x();
	Current_Vy = event->primaryVertex().y();
	Current_Vz = event->primaryVertex().z();
	Current_bField = event->bField() * kilogauss;
	hVzDiff->Fill(Current_Vz - vpdVz);
	hVzvsVpd->Fill(Current_Vz, vpdVz);
	hVertexXYZ->Fill(Current_Vx, Current_Vy, Current_Vz);
	hVertexZAll->Fill(Current_Vz);


	if ( fabs(Current_Vx) < 1e-4 && fabs(Current_Vy) < 1e-4 && fabs(Current_Vz) < 1e-4 )  return kStOK ;

	if (fabs(Current_Vz) > cutAbsVertexZ) return kStOK ;

	if (refmul_for_corr < 0) return kStOK ;

	if (fabs(Current_Vz) > 30) return kStOK;

	hCent_ref_un->Fill(mCent1, refmul_for_corr, mCentWeight);
	int Ntracks = (int)mPicoDst->numberOfTracks();
	float px = -999., py = -999., pz = -999.;
	Current_nCanE = Current_nCanGE = 0;

	for (Int_t i = 0; i < Ntracks; i++)
	{
		StPicoTrack *ptrack = (StPicoTrack *)mPicoDst->track(i);

		if (!ptrack)continue;

		if (fabs(ptrack->charge()) != 1)continue;

		Int_t nFitPts = ptrack->nHitsFit();
		Int_t nHitsMax = ptrack->nHitsMax();
		Int_t nHitsDedx = ptrack->nHitsDedx();
		double dca = ptrack->dca();
		double hitsRatio = ((double)nFitPts) / ((double)nHitsMax);

		if (nFitPts < 15) continue;

		if (hitsRatio < 0.52) continue;

		//   if (dca > 3) continue;

		double nSigE = ptrack->nSigmaElectron();
		double beta = ptrack->btofBeta();
		double trackE0 = ptrack->e0();
		double trackPhiDist = ptrack->phiDist();
		double trackZDist = ptrack->zDist();
		Int_t nEta = ptrack->nEta();
		Int_t nPhi = ptrack->nPhi();


		px = ptrack->pMom().x();
		py = ptrack->pMom().y();
		pz = ptrack->pMom().z();
		TVector3 mom(px, py, pz);
		double p = mom.Mag();
		double pt = mom.Perp();

		if (pt < 0.2)continue;

		double eta = mom.PseudoRapidity();
		double phi = mom.Phi();

		if (fabs(eta) > 1) continue;

		if (fabs(nSigE) > 4) continue;

		Current_CanEId[Current_nCanE] = i;
		Current_CanECharge[Current_nCanE] = ptrack->charge();
		Current_CanEPt[Current_nCanE] = pt;
		Current_CanEP[Current_nCanE] = p;
		Current_CanEEta[Current_nCanE] = eta;
		Current_CanEPhi[Current_nCanE] = phi;
		Current_CanEMom[Current_nCanE] = ptrack->gMom();
		Current_CanEOrigin[Current_nCanE] = ptrack->origin();
		Current_CanEDca[Current_nCanE] = dca;
		Current_CanEnHitsFit[Current_nCanE] = nFitPts;
		Current_CanEnHitsDedx[Current_nCanE] = nHitsDedx;
		Current_CanEnHitsMax[Current_nCanE] = nHitsMax;
		Current_CanEnSigE[Current_nCanE] = nSigE;
		Current_CanEBeta[Current_nCanE] = beta;
		Current_CanEE0[Current_nCanE] = trackE0;
		Current_CanEDistZ[Current_nCanE] = trackZDist;
		Current_CanEDistPhi[Current_nCanE] = trackPhiDist;
		Current_CanEnEta[Current_nCanE] = nEta;
		Current_CanEnPhi[Current_nCanE] = nPhi;
		Current_nCanE++;

		bool tof_cut = fabs(1 - 1 / beta) < 0.03;

		bool btowe_cut = (p / trackE0 > 0) && (p / trackE0 < 2);

		if (p < 1.5)
		{
			if (!tof_cut) continue;
		}

		else
		{
			if (!btowe_cut) continue;
		}

		if ( fabs(nSigE) > 2) continue;

		Current_CanGEId[Current_nCanGE] = i;
		Current_CanGECharge[Current_nCanGE] = ptrack->charge();
		Current_CanGEPt[Current_nCanGE] = pt;
		Current_CanGEEta[Current_nCanGE] = eta;
		Current_CanGEEta[Current_nCanGE] = phi;
		Current_CanGEMom[Current_nCanGE] = ptrack->gMom();
		Current_CanGEOrigin[Current_nCanGE] = ptrack->origin();
		Current_nCanGE++;


	}


	makephe();
	mEventsProcessed++ ;

	return kStOK;
}
Int_t StMyAnalysisMaker::makephe()
{
	TVector3 pheMom1(0, 0, 0);
	TVector3 pheMom2(0, 0, 0);
	TLorentzVector FourpheMom(0, 0, 0, 0);
	TLorentzVector FourpheMom1(0, 0, 0, 0);
	TLorentzVector FourpheMom2(0, 0, 0, 0);
	TVector3 PrimaryVertex(Current_Vx, Current_Vy, Current_Vz);
	Double_t pairdca = 0.0;
	Double_t pairmass = 0.0;
	Double_t pairpt = 0;
	int sign = 0;
	Double_t fill_content[4] = {0., 0., 0., 0.};
	Double_t p = 0.0;
	Double_t nSigE = 0.;
	Double_t nse_p = 0.;
	Bool_t nse_cut = 0;

	for (int i = 0; i < Current_nCanE; i++)
	{
		int id1 = Current_CanEId[i];
		int charge1 = Current_CanECharge[i];


		StThreeVectorD Helix1P(Current_CanEMom[i].x(), Current_CanEMom[i].y(), Current_CanEMom[i].z());
		StThreeVectorD Helix1O(Current_CanEOrigin[i].x(), Current_CanEOrigin[i].y(), Current_CanEOrigin[i].z());

		StPhysicalHelixD Helix1(Helix1P, Helix1O, Current_bField, charge1);

		for (int j = 0; j < Current_nCanGE; j++)
		{
			int id2 = Current_CanGEId[j];

			if (id1 == id2) continue;


			int charge2 = Current_CanGECharge[j];
			StThreeVectorD Helix2P(Current_CanGEMom[j].x(), Current_CanGEMom[j].y(), Current_CanGEMom[j].z());
			StThreeVectorD Helix2O(Current_CanGEOrigin[j].x(), Current_CanGEOrigin[j].y(), Current_CanGEOrigin[j].z());
			StPhysicalHelixD Helix2(Helix2P, Helix2O, Current_bField, charge2);

			pairD pair = Helix2.pathLengths(Helix1);
			double aPathLength = pair.first;
			double ePathLength = pair.second;
			StThreeVectorD aOrigin = Helix2.at(aPathLength);
			StThreeVectorD eOrigin = Helix1.at(ePathLength);
			StThreeVectorD V0Origin = (aOrigin + eOrigin) * 0.5;
			StThreeVectorD diff = aOrigin - eOrigin;
			StThreeVectorD V0pheMom1 = Helix1.momentumAt(ePathLength, Current_bField);
			FourpheMom1.SetXYZM(V0pheMom1.x(), V0pheMom1.y(), V0pheMom1.z(), Melectron);
			StThreeVectorD V0pheMom2 = Helix2.momentumAt(aPathLength, Current_bField);
			FourpheMom2.SetXYZM(V0pheMom2.x(), V0pheMom2.y(), V0pheMom2.z(), Melectron);
			FourpheMom = FourpheMom1 + FourpheMom2;
			pairpt = FourpheMom.Pt();
			pairmass = FourpheMom.M();

			if (charge1 * charge2 == -1)
			{
				hNum->Fill(mCent1, pairpt, pairmass, mCentWeight);
			}

			else if (charge1 == -1 && charge2 == -1)
			{
				hDenM->Fill(mCent1, pairpt, pairmass, mCentWeight);
			}

			else
			{
				hDenP->Fill(mCent1, pairpt, pairmass, mCentWeight);
			}

			if (pairmass < 0.05 && (charge1 * charge2 == -1))
			{
				pairdca = diff.mag();
				hpairdca->Fill(mCent1, pairpt, pairdca, mCentWeight);

				if (pairdca < 3)
				{
					hOrigin_xy->Fill(V0Origin.x(), V0Origin.y(), mCentWeight);
					hOrigin_yz->Fill(V0Origin.y(), V0Origin.z(), mCentWeight);
					hOrigin_xz->Fill(V0Origin.x(), V0Origin.z(), mCentWeight);

				}
			}

			if (pairmass < 0.05 && pairdca < 3)
			{
				if (charge1 * charge2 == -1) sign = 0;
				else if (charge1 == -1 && charge2 == -1) sign = -1;
				else sign = 1;

				fill_content[0] = sign;
				fill_content[1] = mCent1;
				fill_content[2] = Current_CanEPt[i];
				fill_content[3] = Current_CanEnHitsFit[i];
				hnHitsFit->Fill(fill_content, mCentWeight);

				if (cent < 6)
				{
					if (fill_content[3] < 20) continue;

					fill_content[3] = Current_CanEnHitsDedx[i];
					hnHitsDedx->Fill(fill_content, mCentWeight);

					if (fill_content[3] < 12) continue;
				}

				else
				{
					if (fill_content[3] < 25) continue;

					fill_content[3] = Current_CanEnHitsDedx[i];
					hnHitsDedx->Fill(fill_content, mCentWeight);

					if (fill_content[3]  < 15) continue;
				}

				nSigE = Current_CanEnSigE[i];
				p = Current_CanEP[i];
				fill_content[3] = nSigE;

				hnSigE->Fill(fill_content, mCentWeight);
				nse_p = p < 1 ? (1.5 * p - 2.5) : -1;
				nse_cut = nSigE > nse_p && nSigE < 2;

				if (!nse_cut) continue;

				if (Current_CanEBeta[i] > (0 + 1e-4))
				{
					fill_content[3] = 1 / Current_CanEBeta[i];
					htofmatch->Fill(fill_content, mCentWeight);
				}

				if (Current_CanEE0[i] > 1E-4)
				{
					fill_content[3] = Current_CanEP[i] / Current_CanEE0[i];
					hbemcmatch->Fill(fill_content, mCentWeight);

					if (Current_CanEP[i] / Current_CanEE0[i] > 0 && Current_CanEP[i] / Current_CanEE0[i] < 2)
					{
						fill_content[3] = Current_CanEnEta[i];
						hneta->Fill(fill_content, mCentWeight);

						if (fill_content[3] > 0)
						{
							fill_content[3] = Current_CanEDistZ[i];
							hzDist->Fill(fill_content, mCentWeight);

						}

						fill_content[3] = Current_CanEnPhi[i];
						hnphi->Fill(fill_content, mCentWeight);

						if (fill_content[3] > 0)
						{
							fill_content[3] = Current_CanEDistPhi[i];
							hphiDist->Fill(fill_content, mCentWeight);

						}
					}
				}
			}
		}
	}

	return 1;
}

