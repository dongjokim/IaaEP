
#include "include/Filipad.h"
#include "include/rootcommon.h"
#include "TVirtualFitter.h"

TGraphErrors* get_ratio( TGraphErrors * l, TGraphErrors *r );
double FitGeneralizedGaus(double *x, double *par);
double FitGeneralizedGausPlusBG(double *x, double *par);
double FitKaplan(double *x, double *par);
void DrawAfterFlip(int iPTT, int iPTA);
void LoadData(TString inFile="sysErrors/Signal_GG_LHC15o_pass1_CentralBarrelTracking_hadronPID_FieldConfigs_5146_JCIAA_GlobalSDD_LHC17p_pass1_CENT_woSDD_Iaa_R0.2_1.0_1.60_Near_Wing0.root");
void run2();

const int kMAXD       = 20; //maximal number of pT trigger bins
const int kCENT       = 10; //maximal number of pT trigger bins
const int kZvtx       = 15; //maximal number of pT trigger bins
enum dataType { AA, pp };
int NC;

TH1D *hTriggPtBin[2][kCENT][kMAXD]; 
TH1D *hAssocPtBin[2][kCENT][kMAXD][kMAXD]; 

TFile *fin;

TH1D *hDeltaEtaFlip[2][kCENT][kMAXD][kMAXD]; // Flip Deta around 0 to positive side
TH1D *hDeltaEtaSig[2][kCENT][kMAXD][kMAXD]; // background substracted signal based on fit
// IAA before substraction and after
// IAA before substraction and after

// Fits only after Fip at this moment
TF1 *fKaplan[2][kCENT][kMAXD][kMAXD]; // 
TF1 *fGG[2][kCENT][kMAXD][kMAXD]; //

Bool_t saveRoot = kTRUE;
double lowx=-0.8;
double highx=0.8;
double ly = -0.1;
double hy = 0.4;
double lowIAA = 0.0;
double highIAA = 4.;

TLatex latexRun;
TString strRun = "Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV";
Bool_t useGG = kTRUE;//kFALSE; // for background sub

void LoadData(TString inFile="sysErrors/_AA_moon1_pp_moon1_Iaa_R0.2_1.0_1.60_Near_Wing0.root"){
	cout <<"+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<endl;
	cout <<"PbPb : "<< inFile << endl;

	fin = TFile::Open(inFile);

	TVector *TriggPtBorders;
	TVector *AssocPtBorders;
	TVector *CentBinBorders;

	TriggPtBorders             = (TVector*) fin->Get("TriggPtBorders");
	AssocPtBorders             = (TVector*) fin->Get("AssocPtBorders");
	CentBinBorders             = (TVector*) fin->Get("CentBinBorders");

	//int NumCent[2] = {1,1};// for jewel
	int NumCent[2]    =  { CentBinBorders->GetNoElements()-2, 1}; //{1,1};// for jewel
	int NumPtt     = TriggPtBorders->GetNoElements()-1;
	int NumPta     = AssocPtBorders->GetNoElements()-1;
	cout <<"PbPb"<<endl;
	cout <<"bins:  c="<<  NumCent[AA] <<" ptt="<< NumPtt <<" pta="<< NumPta  << endl; 
	cout <<"+++++++++++++++++++++++++++++++++++++++++++++++"<<endl;

	NC = NumCent[AA];
	int NPTT = NumPtt;
	int NPTA = NumPta;
	double MeanPta[NC][NPTT][NPTA];


	//------------ R e a d    D a t a ------------    
	cout <<"Reading data...."<<endl;
	for(int idtyp=0; idtyp<2; idtyp++){ // 0 = AA, 1 = pp  //(*fCentralityBinBorders)[i+1]
		for(int ic=0; ic<NumCent[idtyp]; ic++){
			for(int iptt=0; iptt<NPTT; iptt++){
				for(int ipta=0;ipta<NPTA;ipta++) {
					hDeltaEtaFlip[idtyp][ic][iptt][ipta] = (TH1D *)fin->Get(Form("hDeltaEtaFlip%02dC%02dT%02dA%02d",idtyp,ic,iptt,ipta));
					fKaplan[idtyp][ic][iptt][ipta] = (TF1*)fin->Get(Form("fKaplanDeltaEtaSig%02dC%02dT%02dA%02d",idtyp,ic,iptt,ipta));
					fGG[idtyp][ic][iptt][ipta]  = (TF1*)fin->Get(Form("fGGDeltaEtaSig%02dC%02dT%02dA%02d",idtyp,ic,iptt,ipta));
					hDeltaEtaSig[idtyp][ic][iptt][ipta] = (TH1D *)fin->Get(Form("hDeltaEtaSig%02dC%02dT%02dA%02d",idtyp,ic,iptt,ipta));
				} // ipta
			} // iptt 
		} // ic
	} // pp or AA

}


//------------------------------------------------------------------------------------------------
void DrawAfterFlip(int iPTT, int iPTA) {
	lowx = -0.01;
	for(int ic=0;ic<NC;ic++) {
		Filipad *fpad = new Filipad(ic+1, 1.1, 0.5, 100, 100, 0.7,NC);
		fpad->Draw();
		//==== Upper pad
		TPad *p = fpad->GetPad(1); //upper pad
		p->SetTickx(); p->SetLogx(0); p->SetLogy(1); p->cd();
		//hy = hDeltaEtaFlip[AA][ic][iPTT][iPTA]->GetMaximum()*2.0;
		TH2F *hfr = new TH2F("hfr"," ", 100,lowx, highx, 10, 1e-2, 2); // numbers: tics x, low limit x, upper limit x, tics y, low limit y, upper limit y
		hset( *hfr, "|#Delta#eta|", "1/N_{trigg} dN/d#Delta#eta",1.1,1.0, 0.09,0.09, 0.01,0.01, 0.04,0.05, 510,505);//settings of the upper pad: x-axis, y-axis
		hfr->Draw();
		//Legend definition
		TLegend *leg = new TLegend(0.45,0.6,0.85,0.82,"","brNDC");
		leg->SetTextSize(0.037);leg->SetBorderSize(0);leg->SetFillStyle(0);//legend settings;

		latexRun.DrawLatexNDC( 0.47, 0.85 ,strRun);

		double zeropar_gamma = 1.2;       // width, exponential part, alpha
        double firstpar_w = 0.3;     // division part, beta
        TString opt = "QRN";
		TString fitname = Form("fGG%02d%02d%02d", ic, iPTT, iPTA);
		TF1 *fFitGG     = new TF1(fitname, FitGeneralizedGausPlusBG, 0,0.8, 4); // 4 Parameters
		fFitGG->SetParameter(0,zeropar_gamma);
		fFitGG->SetParameter(1,firstpar_w);
		fFitGG->SetParLimits(1, 0.05, 3.3);

		hDeltaEtaFlip[AA][ic][iPTT][iPTA]->SetMarkerStyle(20);
		hDeltaEtaFlip[AA][ic][iPTT][iPTA]->Draw("p,same");
		hDeltaEtaFlip[pp][0][iPTT][iPTA]->SetMarkerStyle(21);
		hDeltaEtaFlip[pp][0][iPTT][iPTA]->Draw("p,same");
		fKaplan[AA][ic][iPTT][iPTA]->Draw("same");
		fKaplan[pp][0][iPTT][iPTA]->SetLineColor(kBlue);
		fKaplan[pp][0][iPTT][iPTA]->Draw("same");

		hDeltaEtaFlip[AA][ic][iPTT][iPTA]->Fit((TF1*) fFitGG,opt);
    	/*Create a histogram to hold the confidence intervals*/
   		TH1D *hint = new TH1D("hint","Fitted gaussian with .95 conf.band", 100, 0, 0.8);
   		(TVirtualFitter::GetFitter())->GetConfidenceIntervals(hint);
   		//Now the "hint" histogram has the fitted function values as the
   		//bin contents and the confidence intervals as bin errors
   		hint->SetStats(kFALSE);
   		hint->SetFillColor(kGray);
   		hint->Draw("e3 same");
   
		fGG[AA][ic][iPTT][iPTA]->SetLineStyle(2);
		fGG[pp][0][iPTT][iPTA]->SetLineStyle(2);
		fGG[AA][ic][iPTT][iPTA]->Draw("same");
		fGG[pp][0][iPTT][iPTA]->SetLineColor(kBlue);
		fGG[pp][0][iPTT][iPTA]->Draw("same");

		//substracted
		hDeltaEtaSig[AA][ic][iPTT][iPTA]->SetMarkerStyle(24);
		hDeltaEtaSig[AA][ic][iPTT][iPTA]->Draw("p,same");
		hDeltaEtaSig[pp][0][iPTT][iPTA]->SetMarkerStyle(25);
		hDeltaEtaSig[pp][0][iPTT][iPTA]->Draw("p,same");

		leg->AddEntry(hDeltaEtaFlip[AA][ic][iPTT][iPTA],hDeltaEtaFlip[AA][ic][iPTT][iPTA]->GetTitle(),"p");
		leg->AddEntry(hDeltaEtaFlip[pp][0][iPTT][iPTA],hDeltaEtaFlip[pp][0][iPTT][iPTA]->GetTitle(),"p");
		

		leg->Draw();

		//==== Lower pad
		p = fpad->GetPad(2);
		p->SetTickx(); p->SetGridy(1); p->SetLogx(0), p->SetLogy(0); p->cd();
		TH2F *hfr1 = new TH2F("hfr1"," ", 100, lowx, highx, 10, 0.50, 1.5);
		hset( *hfr1, "|#Delta#eta|", "GG/Kaplan",1.1,1.0, 0.09,0.09, 0.01,0.01, 0.08,0.08, 510,505);

		TH1D *htemp[4];
		for(int i=0;i<4;i++ ){ htemp[i] = (TH1D*)hDeltaEtaFlip[AA][ic][iPTT][iPTA]->Clone(); htemp[i]->Reset();}
		htemp[0]->Add(fKaplan[AA][ic][iPTT][iPTA]);
		htemp[1]->Add(fGG[AA][ic][iPTT][iPTA]);
		htemp[2]->Add(fKaplan[pp][0][iPTT][iPTA]);
		htemp[3]->Add(fGG[pp][0][iPTT][iPTA]);
		
		htemp[1]->Divide(htemp[0]); // GG/Kap AA
		htemp[3]->Divide(htemp[2]); // GG/Kap pp

		htemp[1]->SetLineColor(2);
		htemp[1]->SetMarkerStyle(20);
		htemp[1]->SetMarkerColor(2);
		
		htemp[3]->SetLineStyle(1);
		htemp[3]->SetLineColor(1);

		hfr1->Draw();
		htemp[1]->Draw("p,same");
		htemp[3]->Draw("same");
	}
}


//------------------------------------------------------------------------------------------------
double FitGeneralizedGaus(double *x, double *par) {
	return par[2]*par[0]/(2*par[1]*TMath::Gamma(1/par[0]))*TMath::Exp(-TMath::Power(TMath::Abs(x[0])/par[1], par[0]));
}

//------------------------------------------------------------------------------------------------
double FitGeneralizedGausPlusBG(double *x, double *par) {
	return FitGeneralizedGaus(x, par) + par[3];
}

//------------------------------------------------------------------------------------------------
double FitKaplan(double *x, double *par){  //Fit Deta near with const + Kaplan
	double deta = x[0];
	double bg   = par[0];
	double ampl = par[1];
	double k    = par[2];
	double n    = par[3];

	return bg + ampl*pow( 1 + k*deta*deta,-n);
}


//------------------------------------------------------------------------------------------------
TH1D* SubtractBg(TH1D *h, double bg, double ebg){
	int nb =h->GetNbinsX();
	TString hname = h->GetName();

	TH1D *hsig  = (TH1D*) h->Clone(Form("%s_sig",hname.Data()));
	hsig->Reset();

	for(int ib=1; ib<=nb; ib++){
		double val = h->GetBinContent(ib);
		double err = h->GetBinError(ib);
		hsig->SetBinContent(ib,val-bg);
		hsig->SetBinError(ib,sqrt(err*err +  ebg * ebg));
	}

	return hsig;
}