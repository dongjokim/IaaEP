
#include "include/Filipad.h"
#include "include/rootcommon.h"

TH1D *Flip(TH1D* hin, int idtyp);
TGraphErrors* get_ratio( TGraphErrors * l, TGraphErrors *r );
double FitGeneralizedGaus(double *x, double *par);
double FitGeneralizedGausPlusBG(double *x, double *par);
double FitKaplan(double *x, double *par);
TH1D* SubtractBg(TH1D *h, double bg, double ebg);
void DrawAfterFlip(int iPTT, int iPTA);
void DrawBeforeFlip(int iPTT, int iPTA);
void DrawSignal(int iPTT, int iPTA);
void DoAnalysis(TString inFile="sysErrors/_AA_moon1_pp_moon1_Iaa_R0.2_1.0_1.60_Near_Wing0.root",  TString oname="");
void run2();

const int kMAXD       = 20; //maximal number of pT trigger bins
const int kCENT       = 10; //maximal number of pT trigger bins
const int kZvtx       = 15; //maximal number of pT trigger bins
enum dataType { AA, pp };
int NC;

TH1D *hTriggPtBin[2][kCENT][kMAXD]; 
TH1D *hAssocPtBin[2][kCENT][kMAXD][kMAXD]; 

TH1D *hIAAEta[kCENT][kMAXD][kMAXD]; // in eta,phi 

TFile *fin;

// with $\Delta\phi$ < 0.2 ??? check
TH1D *hDeltaEta[2][kCENT][kMAXD][kMAXD]; // summed DeltaEta AA-0 pp-1
TH1D *hDeltaEtaFlip[2][kCENT][kMAXD][kMAXD]; // Flip Deta around 0 to positive side
TH1D *hDeltaEtaSig[2][kCENT][kMAXD][kMAXD]; // background substracted signal based on fit
// IAA before substraction and after
TH1D *hIAADeltaEta[kCENT][kMAXD][kMAXD];
TH1D *hIAADeltaEtaFlip[kCENT][kMAXD][kMAXD];
TH1D *hIAADeltaEtaSig[kCENT][kMAXD][kMAXD]; // background substracted signal IAA

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
//Bool_t useGG = kTRUE; // for background sub
Bool_t useGG = kFALSE; // for background sub

void run1() {

	const int Nsets = 2;
	TString infiles[Nsets] = {
		"sysErrors/_AMPT_LHC13f3c_JCIaa_KineOnly_pythia8230_pp5.02TeV_GF0_SoftQCD_Iaa_R0.2_1.0_1.60_Near_Wing0.root",
		"sysErrors/_LHC10h_AOD86_MgFpMgFm_5217_JCIAA_TPCOnly_H0_T0_LHC11a_p4_AOD113_noSDD_Iaa_R0.2_1.0_1.60_Near_Wing0.root"
	};

	TObjArray *outString[Nsets];
	TString outrootname[Nsets];
	for(int i=0;i<Nsets;i++) { 
		outString[i] = infiles[i].Tokenize("/");
		TString sDir = ((TObjString *)outString[i]->At(0))->GetString();
		TString sName = ((TObjString *)outString[i]->At(1))->GetString();
		if(useGG) {
			outrootname[i] = Form("%s/Signal_GG%s",sDir.Data(),sName.Data());
		} else {
		outrootname[i] = Form("%s/Signal%s",sDir.Data(),sName.Data());
		}
		//cout << outrootname[i] << endl;
	}
	for(int i=0;i<Nsets;i++) { 
		DoAnalysis(infiles[i],outrootname[i]);
	}
//	DoAnalysis ("sysErrors/_AA_moon1_pp_moon1_Iaa_R0.2_1.0_1.60_Near_Wing0.root","sysErrors/_Signal_AA_moon1_pp_moon1_Iaa_R0.2_1.0_1.60_Near_Wing0.root");
}

void run2systematics() {

	const int Nsets = 4;
	TString infiles[Nsets] = {
		"sysErrors/_LHC15o_CentralBarrelTracking_hadronPID-8356_Hybrid_vtx08_JCIAA_Hybrid_LHC17p_pass1_CENT_woSDD_Iaa_R0.2_1.0_1.60_Near_Wing0.root"
		,"sysErrors/_LHC15o_CentralBarrelTracking_hadronPID_FieldConfigs-8322_JCIAA_GlobalSDD_VTX08_LHC17p_pass1_CENT_woSDD_Iaa_R0.2_1.0_1.60_Near_Wing0.root"
	    ,"sysErrors/_LHC15o_pass1_CentralBarrelTracking_hadronPID_FieldConfigs_5146_JCIAA_GlobalSDD_LHC17p_pass1_CENT_woSDD_Iaa_R0.2_1.0_1.60_Near_Wing0.root" 
        ,"sysErrors/_LHC15o_pass1_CentralBarrelTracking_hadronPID_FieldConfigs_829_Hybrid_JCIAA_TPCOnly_LHC17p_pass1_CENT_woSDD_Iaa_R0.2_1.0_1.60_Near_Wing0.root"
	};

	TObjArray *outString[Nsets];
	TString outrootname[Nsets];
	for(int i=0;i<Nsets;i++) { 
		outString[i] = infiles[i].Tokenize("/");
		TString sDir = ((TObjString *)outString[i]->At(0))->GetString();
		TString sName = ((TObjString *)outString[i]->At(1))->GetString();
		if(useGG) {
			outrootname[i] = Form("%s/Signal_GG%s",sDir.Data(),sName.Data());
		} else {
		outrootname[i] = Form("%s/Signal%s",sDir.Data(),sName.Data());
		}
		//cout << outrootname[i] << endl;
	}
	for(int i=0;i<1;i++) { 
		DoAnalysis(infiles[i],outrootname[i]);
	}
}

void runOnflyModel() {

	const int Nsets = 2;
	TString infiles[Nsets] = {
		"sysErrors/_JEWEL_v2.0.2_JCIaa_KineOnly_JEWEL_vacuum_Iaa_R0.2_1.0_1.60_Near_Wing0.root",
		"sysErrors/_JEWEL_v2.0.2_KeepRecoil_JCIaa_KineOnly_JEWEL_vacuum_Iaa_R0.2_1.0_1.60_Near_Wing0.root"
		//"sysErrors/_JEWEL_JCIaa_KineOnly_JEWEL_vacuum_Iaa_R0.2_1.0_1.60_Near_Wing0.root"
	};

	TObjArray *outString[Nsets];
	TString outrootname[Nsets];
	for(int i=0;i<Nsets;i++) { 
		outString[i] = infiles[i].Tokenize("/");
		TString sDir = ((TObjString *)outString[i]->At(0))->GetString();
		TString sName = ((TObjString *)outString[i]->At(1))->GetString();
		if(useGG) {
			outrootname[i] = Form("%s/Signal_GG%s",sDir.Data(),sName.Data());
		} else {
			outrootname[i] = Form("%s/Signal%s",sDir.Data(),sName.Data());
		}
		//cout << outrootname[i] << endl;
	}
	for(int i=0;i<Nsets;i++) { 
		DoAnalysis(infiles[i],outrootname[i]);
	}
//	DoAnalysis ("sysErrors/_AA_moon1_pp_moon1_Iaa_R0.2_1.0_1.60_Near_Wing0.root","sysErrors/_Signal_AA_moon1_pp_moon1_Iaa_R0.2_1.0_1.60_Near_Wing0.root");
}

void DoAnalysis(TString inFile="sysErrors/_AA_moon1_pp_moon1_Iaa_R0.2_1.0_1.60_Near_Wing0.root",  TString oname=""){
	cout <<"+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<endl;
	cout <<"PbPb : "<< inFile << endl;
	cout <<"Discription :"<< oname << endl;

	fin = TFile::Open(inFile);

	TVector *TriggPtBorders;
	TVector *AssocPtBorders;
	TVector *CentBinBorders;

	TriggPtBorders             = (TVector*) fin->Get("TriggPtBorders");
	AssocPtBorders             = (TVector*) fin->Get("AssocPtBorders");
	CentBinBorders             = (TVector*) fin->Get("CentBinBorders");

	int NumCent[2] = {1,1};// for jewel
	//int NumCent[2]    =  { CentBinBorders->GetNoElements()-2, 1}; //{1,1};// for jewel
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
					hDeltaEta[idtyp][ic][iptt][ipta] = (TH1D *)fin->Get(Form("hDeltaEtaType%02dC%02dT%02dA%02d",idtyp,ic,iptt,ipta));
					//hDeltaEta[idtyp][ic][iptt][ipta]->Print();
					hDeltaEta[idtyp][ic][iptt][ipta]->SetXTitle("#Delta#eta");
					hDeltaEta[idtyp][ic][iptt][ipta]->GetXaxis()->CenterTitle(kTRUE);
					hDeltaEta[idtyp][ic][iptt][ipta]->GetXaxis()->SetTitleOffset(2);
					hDeltaEta[idtyp][ic][iptt][ipta]->SetYTitle("1/N_{trigg} dN/d#Delta#eta"); // done in z01
					hDeltaEtaFlip[idtyp][ic][iptt][ipta] = Flip((TH1D*) hDeltaEta[idtyp][ic][iptt][ipta], idtyp);
				} // ipta
			} // iptt 
		} // ic
	} // pp or AA

	// Need to substract background from Generalized Gaussians..
	// Fit functions
	double fitRange = -1;
	double MaxEta = 0.8;
	// GG fit
	double zeropar_gamma = 1.2;       // width, exponential part, alpha
    double firstpar_w = 0.3;     // division part, beta
	for(int idtyp=0; idtyp<2; idtyp++){ // 0 = AA, 1 = pp  //(*fCentralityBinBorders)[i+1]
		for(int ic=0; ic<NumCent[idtyp]; ic++){
			for(int iptt=0; iptt<NPTT; iptt++){
				for(int ipta=0;ipta<NPTA;ipta++) {
					// Kaplan
					double maxRange = (fitRange > 0) ? fitRange : 2*MaxEta;
					TString fitname = Form("Kaplan%02d%02d%02d%02d", idtyp, ic, iptt, ipta);
					fKaplan[idtyp][ic][iptt][ipta] = new TF1(fitname, FitKaplan, 0, maxRange,4); 
					// GG
					fitname = Form("GG%02d%02d%02d%02d", idtyp, ic, iptt, ipta);
					fGG[idtyp][ic][iptt][ipta]     = new TF1(fitname, FitGeneralizedGausPlusBG, 0,maxRange, 4); // 4 Parameters
					fGG[idtyp][ic][iptt][ipta]->SetParameter(0,zeropar_gamma);
					fGG[idtyp][ic][iptt][ipta]->SetParameter(1,firstpar_w);
				    fGG[idtyp][ic][iptt][ipta]->SetParLimits(1, 0.05, 3.3);

					//estimate fit parametes
					TH1D *hFlipDeta = (TH1D*) hDeltaEtaFlip[idtyp][ic][iptt][ipta];
					double bg       = hFlipDeta->Integral(hFlipDeta->FindBin(0.4), hFlipDeta->GetNbinsX())/(hFlipDeta->GetNbinsX()-hFlipDeta->FindBin(0.4));
					double peakAmpl = hFlipDeta->GetBinContent(hFlipDeta->FindBin(0))-bg;
					fKaplan[idtyp][ic][iptt][ipta]->SetParameters(bg, peakAmpl, 20.0, 1.0);
					//fKaplan[idtyp][ic][iptt][ipta]->SetParLimits(0,bg/10.,bg*4.);
					//fKaplan[idtyp][ic][iptt][ipta]->SetParLimits(1,2e-3,100);
					//fKaplan[idtyp][ic][iptt][ipta]->SetParLimits(3,2e-3,10);
					TString opt = "QRN";
					hFlipDeta->Fit((TF1*) fKaplan[idtyp][ic][iptt][ipta],opt);
					hFlipDeta->Fit((TF1*) fGG[idtyp][ic][iptt][ipta],opt);
					bg         = fKaplan[idtyp][ic][iptt][ipta]->GetParameter(0);
					double ebg = fKaplan[idtyp][ic][iptt][ipta]->GetParError(0);
					double bg_GG = fGG[idtyp][ic][iptt][ipta]->GetParameter(3);
					double ebg_GG = fGG[idtyp][ic][iptt][ipta]->GetParError(3);
					//cout << Form("bg = %.3f:%.3f, err= %.3f:%.3f",bg,bg_GG,ebg,ebg_GG) << endl;
					if(useGG) {
						hDeltaEtaSig[idtyp][ic][iptt][ipta] = SubtractBg(hFlipDeta,bg_GG,ebg_GG); //subtract bg
					} else {
          			    hDeltaEtaSig[idtyp][ic][iptt][ipta] = SubtractBg(hFlipDeta,bg,ebg); //subtract bg
          			}
				} // ipta
			} // iptt 
		} // ic
	} // pp or AA

	cout <<"Calculationg IAA..."<<endl;

	for(int ic=0; ic<NumCent[AA]; ic++){
		for(int iptt=0; iptt<NPTT; iptt++){
			for(int ipta=0;ipta<NPTA;ipta++) {
				hIAADeltaEta[ic][iptt][ipta] = (TH1D*)hDeltaEta[AA][ic][iptt][ipta]->Clone();
				hIAADeltaEta[ic][iptt][ipta]->Divide(hDeltaEta[pp][0][iptt][ipta]);

				hIAADeltaEtaFlip[ic][iptt][ipta] = (TH1D*)hDeltaEtaFlip[AA][ic][iptt][ipta]->Clone();
				hIAADeltaEtaFlip[ic][iptt][ipta]->Divide(hDeltaEtaFlip[pp][0][iptt][ipta]);

				hIAADeltaEtaSig[ic][iptt][ipta]= (TH1D*)hDeltaEtaSig[AA][ic][iptt][ipta]->Clone();
				hIAADeltaEtaSig[ic][iptt][ipta]->Divide(hDeltaEtaSig[pp][0][iptt][ipta]);

			} // iptt 
		} // ic
	}


	int iPTT=3;
	int iPTA=4;

	//DrawBeforeFlip(iPTT, iPTA);
	DrawAfterFlip(iPTT, iPTA);
	//DrawSignal(iPTT, iPTA);

	if(saveRoot) {
		cout <<"Writing the results into a file..."<< endl;
		TFile *fout = new TFile(Form("%s",oname.Data()),"recreate");
		fout->cd();
		// Deltaeta
		for(int idtyp=0; idtyp<2; idtyp++){ // 0 = AA, 1 = pp
			for(int ic=0; ic<NumCent[idtyp]; ic++){
				for(int iptt=0; iptt<NPTT; iptt++){
					for(int ipta=0;ipta<NPTA;ipta++) {
						hDeltaEtaSig[idtyp][ic][iptt][ipta]->Write(Form("hDeltaEtaSig%02dC%02dT%02dA%02d",idtyp,ic,iptt,ipta));
						if(idtyp==AA) hIAADeltaEtaSig[ic][iptt][ipta]->Write(Form("hIAADeltaEtaSigC%02dT%02dA%02d",ic,iptt,ipta));
						hDeltaEtaFlip[idtyp][ic][iptt][ipta]->Write(Form("hDeltaEtaFlip%02dC%02dT%02dA%02d",idtyp,ic,iptt,ipta));
						fKaplan[idtyp][ic][iptt][ipta]->Write(Form("fKaplanDeltaEtaSig%02dC%02dT%02dA%02d",idtyp,ic,iptt,ipta));
						fGG[idtyp][ic][iptt][ipta]->Write(Form("fGGDeltaEtaSig%02dC%02dT%02dA%02d",idtyp,ic,iptt,ipta));
					} // pta
				} // ptt 
			} // cent
		} // type 
		fout->cd();
		CentBinBorders->Write("CentBinBorders");
		TriggPtBorders->Write("TriggPtBorders");
		AssocPtBorders->Write("AssocPtBorders");
		fout->Close();
	}

}

//------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------
void DrawSignal(int iPTT, int iPTA) {
	lowx = -0.01;
	for(int ic=0;ic<NC;ic++) {
		Filipad *fpad = new Filipad(ic+1, 1.1, 0.5, 100, 100, 0.7,NC);
		fpad->Draw();
		//==== Upper pad
		TPad *p = fpad->GetPad(1); //upper pad
		p->SetTickx(); p->SetLogx(0); p->SetLogy(0); p->cd();
		hy = hDeltaEtaSig[AA][ic][iPTT][iPTA]->GetMaximum()*1.2;
		TH2F *hfr = new TH2F("hfr"," ", 100,lowx, highx, 10, ly, hy); // numbers: tics x, low limit x, upper limit x, tics y, low limit y, upper limit y
		hset( *hfr, "|#Delta#eta|", "1/N_{trigg} dN/d#Delta#eta",1.1,1.0, 0.09,0.09, 0.01,0.01, 0.04,0.05, 510,505);//settings of the upper pad: x-axis, y-axis
		hfr->Draw();
		//Legend definition
		TLegend *leg = new TLegend(0.45,0.6,0.85,0.82,"","brNDC");
		leg->SetTextSize(0.037);leg->SetBorderSize(0);leg->SetFillStyle(0);//legend settings;

		latexRun.DrawLatexNDC( 0.47, 0.85 ,strRun);

		hDeltaEtaSig[AA][ic][iPTT][iPTA]->SetMarkerStyle(20);
		hDeltaEtaSig[AA][ic][iPTT][iPTA]->Draw("p,same");
		hDeltaEtaSig[pp][0][iPTT][iPTA]->SetMarkerStyle(24);
		hDeltaEtaSig[pp][0][iPTT][iPTA]->Draw("p,same");


		leg->AddEntry(hDeltaEtaSig[AA][ic][iPTT][iPTA],hDeltaEtaSig[AA][ic][iPTT][iPTA]->GetTitle(),"p");
		leg->AddEntry(hDeltaEtaSig[pp][0][iPTT][iPTA],hDeltaEtaSig[pp][0][iPTT][iPTA]->GetTitle(),"p");
		
		leg->Draw();

		//==== Lower pad
		p = fpad->GetPad(2);
		p->SetTickx(); p->SetGridy(1); p->SetLogx(0), p->SetLogy(0); p->cd();
		TH2F *hfr1 = new TH2F("hfr1"," ", 100, lowx, highx, 10, lowIAA, highIAA);
		hset( *hfr1, "|#Delta#eta|", "AA/pp",1.1,1.0, 0.09,0.09, 0.01,0.01, 0.08,0.08, 510,505);
		hfr1->Draw();
		hIAADeltaEtaSig[ic][iPTT][iPTA]->SetMarkerStyle(21);
		hIAADeltaEtaSig[ic][iPTT][iPTA]->Draw("p,same");
		//gPad->GetCanvas()->SaveAs(Form("figs_svn/FigA4_v%d_modelcomparisonBest.eps",i+2));
	}
}

//------------------------------------------------------------------------------------------------
void DrawAfterFlip(int iPTT, int iPTA) {
	lowx = -0.01;
	for(int ic=0;ic<NC;ic++) {
		Filipad *fpad = new Filipad(ic+1, 1.1, 0.5, 100, 100, 0.7,NC);
		fpad->Draw();
		//==== Upper pad
		TPad *p = fpad->GetPad(1); //upper pad
		p->SetTickx(); p->SetLogx(0); p->SetLogy(0); p->cd();
		hy = hDeltaEtaFlip[AA][ic][iPTT][iPTA]->GetMaximum()*1.2;
		TH2F *hfr = new TH2F("hfr"," ", 100,lowx, highx, 10, ly, hy); // numbers: tics x, low limit x, upper limit x, tics y, low limit y, upper limit y
		hset( *hfr, "|#Delta#eta|", "1/N_{trigg} dN/d#Delta#eta",1.1,1.0, 0.09,0.09, 0.01,0.01, 0.04,0.05, 510,505);//settings of the upper pad: x-axis, y-axis
		hfr->Draw();
		//Legend definition
		TLegend *leg = new TLegend(0.45,0.6,0.85,0.82,"","brNDC");
		leg->SetTextSize(0.037);leg->SetBorderSize(0);leg->SetFillStyle(0);//legend settings;

		latexRun.DrawLatexNDC( 0.47, 0.85 ,strRun);

		hDeltaEtaFlip[AA][ic][iPTT][iPTA]->SetMarkerStyle(20);
		hDeltaEtaFlip[AA][ic][iPTT][iPTA]->Draw("p,same");
		hDeltaEtaFlip[pp][0][iPTT][iPTA]->SetMarkerStyle(24);
		hDeltaEtaFlip[pp][0][iPTT][iPTA]->Draw("p,same");
		fKaplan[AA][ic][iPTT][iPTA]->Draw("same");
		fKaplan[pp][0][iPTT][iPTA]->SetLineColor(kBlue);
		fKaplan[pp][0][iPTT][iPTA]->Draw("same");

		fGG[AA][ic][iPTT][iPTA]->SetLineStyle(2);
		fGG[pp][0][iPTT][iPTA]->SetLineStyle(2);
		fGG[AA][ic][iPTT][iPTA]->Draw("same");
		fGG[pp][0][iPTT][iPTA]->SetLineColor(kBlue);
		fGG[pp][0][iPTT][iPTA]->Draw("same");

		leg->AddEntry(hDeltaEtaFlip[AA][ic][iPTT][iPTA],hDeltaEtaFlip[AA][ic][iPTT][iPTA]->GetTitle(),"p");
		leg->AddEntry(hDeltaEtaFlip[pp][0][iPTT][iPTA],hDeltaEtaFlip[pp][0][iPTT][iPTA]->GetTitle(),"p");
		
		leg->Draw();

		//==== Lower pad
		p = fpad->GetPad(2);
		p->SetTickx(); p->SetGridy(1); p->SetLogx(0), p->SetLogy(0); p->cd();
		TH2F *hfr1 = new TH2F("hfr1"," ", 100, lowx, highx, 10, lowIAA, highIAA);
		hset( *hfr1, "|#Delta#eta|", "AA/pp",1.1,1.0, 0.09,0.09, 0.01,0.01, 0.08,0.08, 510,505);
		hfr1->Draw();
		hIAADeltaEtaFlip[ic][iPTT][iPTA]->SetMarkerStyle(21);
		hIAADeltaEtaFlip[ic][iPTT][iPTA]->Draw("p,same");
		//gPad->GetCanvas()->SaveAs(Form("figs_svn/FigA4_v%d_modelcomparisonBest.eps",i+2));
	}
}

//------------------------------------------------------------------------------------------------
void DrawBeforeFlip(int iPTT, int iPTA) {
	for(int ic=0;ic<NC;ic++) {
		Filipad *fpad = new Filipad(ic+1, 1.1, 0.5, 100, 100, 0.7,NC);
		fpad->Draw();
		//==== Upper pad
		TPad *p = fpad->GetPad(1); //upper pad
		p->SetTickx(); p->SetLogx(0); p->SetLogy(0); p->cd();
		hy = hDeltaEta[AA][ic][iPTT][iPTA]->GetMaximum()*1.2;
		TH2F *hfr = new TH2F("hfr"," ", 100,lowx, highx, 10, ly, hy); // numbers: tics x, low limit x, upper limit x, tics y, low limit y, upper limit y
		hset( *hfr, "#Delta#eta", "1/N_{trigg} dN/d#Delta#eta",1.1,1.0, 0.09,0.09, 0.01,0.01, 0.04,0.05, 510,505);//settings of the upper pad: x-axis, y-axis
		hfr->Draw();
		//Legend definition
		TLegend *leg = new TLegend(0.45,0.6,0.85,0.82,"","brNDC");
		leg->SetTextSize(0.037);leg->SetBorderSize(0);leg->SetFillStyle(0);//legend settings;

		latexRun.DrawLatexNDC( 0.47, 0.85 ,strRun);
		hDeltaEta[AA][ic][iPTT][iPTA]->SetMarkerStyle(20);
		hDeltaEta[AA][ic][iPTT][iPTA]->Draw("p,same");
		hDeltaEta[pp][0][iPTT][iPTA]->SetMarkerStyle(24);
		hDeltaEta[pp][0][iPTT][iPTA]->Draw("p,same");
		hDeltaEtaFlip[AA][ic][iPTT][iPTA]->SetMarkerStyle(25);
		hDeltaEtaFlip[AA][ic][iPTT][iPTA]->Draw("p,same");

		leg->AddEntry(hDeltaEta[AA][ic][iPTT][iPTA],hDeltaEta[AA][ic][iPTT][iPTA]->GetTitle(),"p");
		leg->AddEntry(hDeltaEta[pp][0][iPTT][iPTA],hDeltaEta[pp][0][iPTT][iPTA]->GetTitle(),"p");
		
		leg->Draw();

		//==== Lower pad
		p = fpad->GetPad(2);
		p->SetTickx(); p->SetGridy(1); p->SetLogx(0), p->SetLogy(0); p->cd();
		TH2F *hfr1 = new TH2F("hfr1"," ", 100, lowx, highx, 10, lowIAA, highIAA);
		hset( *hfr1, "#Delta#eta", "AA/pp",1.1,1.0, 0.09,0.09, 0.01,0.01, 0.08,0.08, 510,505);
		hfr1->Draw();
		hIAADeltaEta[ic][iPTT][iPTA]->SetMarkerStyle(20);
		hIAADeltaEta[ic][iPTT][iPTA]->Draw("p,same");
		hIAADeltaEtaFlip[ic][iPTT][iPTA]->SetMarkerStyle(24);
		hIAADeltaEtaFlip[ic][iPTT][iPTA]->Draw("p,same");
		//gPad->GetCanvas()->SaveAs(Form("figs_svn/FigA4_v%d_modelcomparisonBest.eps",i+2));
	}
}

//------------------------------------------------------------------------------------------------
TH1D *Flip(TH1D* hin, int idtyp){
	int nb  = hin->GetNbinsX();
	double scale = 2.0; // MUST Have dN/deta -> dN/d|eta|
	double max = hin->GetBinLowEdge(nb+1);
	TString hname = hin->GetName();
	TString newName = Form("%s_flip%d",hname.Data(),idtyp);

	TH1D *hFlip = new TH1D(newName.Data(), newName.Data(), (int) nb/2., 0, max);
	hFlip->SetTitle(hin->GetTitle());
	int zero = hin->FindBin(0.00001);
	for(int ib=zero; ib<=nb; ib++){
		double valPos = hin->GetBinContent(ib);
		double errPos = hin->GetBinError(ib);
		double valNeg = hin->GetBinContent(nb - ib+1);
		double errNeg = hin->GetBinError(nb - ib+1);

		hFlip->SetBinContent(ib-zero+1, (valPos+valNeg)/scale);
		hFlip->SetBinError(ib-zero+1, sqrt( errPos*errPos + errNeg * errNeg )/scale);
	}

	return hFlip;
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