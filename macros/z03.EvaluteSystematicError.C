#include "include/Filipad.h"
#include "include/rootcommon.h"

int qColors[] = {kAzure, kRed+1, kBlack, kGreen+3};
int fillStylesTheory[] = {3001,3004,3345,3354,1001,3008}; 
int fillColors[] = {kAzure-4, kRed-9, 12, kGreen+2};
int fillsty = 3001;
double mSize = 1.8;
double lwidth = 1.5;


void Barlow( TH1D* hd, TH1D* hv, TH1D* hbar, int corrtype ){
 for(int j=0;j<hd->GetNbinsX();j++){
	hbar->Fill( ( hd->GetBinContent(j+1) - hv->GetBinContent(j+1) ) /
	sqrt( fabs( pow( hd->GetBinError(j+1),2 ) + corrtype*pow( hv->GetBinError(j+1),2 ) ) ) );
 }
}


TH1D* Smooth( TH1D* hT ){
 TH1D* hSmooth = (TH1D*)hT->Clone(0); hSmooth->Reset();
 for(int i=1;i<hSmooth->GetNbinsX()-1;i++){
        hSmooth->SetBinContent( i+1, ( hT->GetBinContent(i) + hT->GetBinContent(i+1) + hT->GetBinContent(i+2) )/3.0 );
 }
 hSmooth->SetBinContent( 1, hT->GetBinContent(1) );
 hSmooth->SetBinContent( hSmooth->GetNbinsX(), hT->GetBinContent(hSmooth->GetNbinsX()) );
 return hSmooth;
}

TH1D* GetPlotWithSyst( TH1D* hCntl, TH1D* hFracSyst ){
 TH1D* hReturn = (TH1D*)hCntl->Clone(0);
 for(int i=0;i<hReturn->GetNbinsX();i++){
	hReturn->SetBinError( i+1, hReturn->GetBinContent(i+1)*hFracSyst->GetBinContent(i+1) );
 }
 return hReturn;
}

TH1D* GetRelativeSyst(TH1D* hFracSyst ){
 TH1D* hReturn = (TH1D*)hFracSyst->Clone(0);
 for(int i=0;i<hReturn->GetNbinsX();i++){
	hReturn->SetBinContent( i+1, 1.+ hFracSyst->GetBinContent(i+1) );
 }
 return hReturn;
}

TGraphErrors* h2g(TH1D *hid);
void LoadData();
void Compare();
void DrawSignal(int padid, int iPTT, int iPTA);
void DrawIAA(int padID, int iPTT, int iPTA);
void ObtainSyst();
void BarlowTest();
void ObtainTotalSyst();
void SaveFinalResults();

double lowx=-0.1;
double highx=0.3;
double ly = -0.05;
double hy = 3.0;
double rlow = 0.5;  // ratio
double rhigh = 1.5;

TLatex latexRun;
TString strRun = "Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV";

// mainly for track cut 
// vtx change..
const int Nsets = 3;
TString infiles[Nsets] = {
	"sysErrors/Signal_LHC15o_pass1_CentralBarrelTracking_hadronPID_FieldConfigs_5146_JCIAA_GlobalSDD_LHC17p_pass1_CENT_woSDD_Iaa_R0.2_1.0_1.60_Near_Wing0.root",
	"sysErrors/Signal_LHC15o_pass1_CentralBarrelTracking_hadronPID_FieldConfigs_829_Hybrid_JCIAA_TPCOnly_LHC17p_pass1_CENT_woSDD_Iaa_R0.2_1.0_1.60_Near_Wing0.root",
	"sysErrors/Signal_LHC15o_CentralBarrelTracking_hadronPID_FieldConfigs-8322_JCIAA_GlobalSDD_VTX08_LHC17p_pass1_CENT_woSDD_Iaa_R0.2_1.0_1.60_Near_Wing0.root"
};
TFile *fin[Nsets];

TString sLeg[Nsets] = {
	"GlobalSDD default",
	"Hybrid",
	"zvtx08"
};

int gMarkers[] = {20,24,21,25,23,27,29,30};
int gColors[]={kBlack, kRed, kBlue, kDeepSea, kPink, kGray, kRed, kBlack};

const int kMAXD       = 20; //maximal number of pT trigger bins
const int kCENT       = 10; //maximal number of pT trigger bins
enum dataType { AA, pp };


TH1D *hDeltaEtaSig[Nsets][2][kCENT][kMAXD][kMAXD]; // background substracted signal based on fit
TH1D *hIAADeltaEtaSig[Nsets][kCENT][kMAXD][kMAXD]; // background substracted signal IAA
TGraphErrors* grIAADeltaEtaSig[kCENT][kMAXD][kMAXD]; // final result to graph (default)

TH1D *hRatio_DeltaEtaSig[Nsets][kCENT][kMAXD][kMAXD]; //Data TPC V0A V0P to AMPT inclusive
TH1D *hRatio_IAA[Nsets][kCENT][kMAXD][kMAXD];
TH1D *hRatio_IAA_fordrawing[Nsets][kCENT][kMAXD][kMAXD];

TGraphAsymmErrors *grFinalIAAStat[Nsets][kCENT][kMAXD][kMAXD];      //For Filip's QM, IAA with statstistical 
TGraphAsymmErrors *grFinalIAASyst[Nsets][kCENT][kMAXD][kMAXD];      //For Filip's QM, IAA with Systematics
TGraphAsymmErrors *grFinalIAAScale[Nsets][kCENT][kMAXD][kMAXD];      //For Filip's QM, Converting box to TGraph
TBox *sysbox[Nsets][kCENT][kMAXD][kMAXD];                           //For Filip's QM, Scaling Box


TVector *TriggPtBorders;
TVector *AssocPtBorders;
TVector *CentBinBorders;
int NumCent[2];
int NPTT;
int NPTA;
int iRef=0; // default


TH1D* hBarlowDist[Nsets];
// individual 0-1 relative 
TH1D* hSystSmooth[Nsets][kCENT][kMAXD][kMAXD]; 
TH1D* hSystSmoothTotal[kCENT][kMAXD][kMAXD];  // total

TH1D* hSystSmoothRelativeSyst[Nsets][kCENT][kMAXD][kMAXD];

// individual systematics mean is taken from the defaults
TH1D* hIAADeltaEtaSig_syst[Nsets][kCENT][kMAXD][kMAXD]; 
TH1D* hIAADeltaEtaSig_Totsyst[kCENT][kMAXD][kMAXD]; // total systenatics
TGraphErrors* grIAADeltaEtaSig_syst[Nsets][kCENT][kMAXD][kMAXD];
TGraphErrors* grIAADeltaEtaSig_Totsyst[kCENT][kMAXD][kMAXD];
void RemovePoints(TGraphErrors *ge, double xmin, double xmax);

//------------------------------------------------------------------------------------------------
void LoadData() {
	
	for(int i=0;i<Nsets;i++){
		fin[i] = TFile::Open(infiles[i]);
	}

	int irefD = 0;
	TriggPtBorders             = (TVector*) fin[irefD]->Get("TriggPtBorders");
	AssocPtBorders             = (TVector*) fin[irefD]->Get("AssocPtBorders");
	CentBinBorders             = (TVector*) fin[irefD]->Get("CentBinBorders");
	NumCent[AA]    = CentBinBorders->GetNoElements()-2;  // 5TeV 1 less cent
	NumCent[pp]    = 1; 
	NPTT     = TriggPtBorders->GetNoElements()-1;
	NPTA     = AssocPtBorders->GetNoElements()-1;
	cout <<"PbPb"<<endl;
	cout <<"bins:  c="<<  NumCent[AA] <<" ptt="<< NPTT <<" pta="<< NPTA  << endl; 
	cout <<"+++++++++++++++++++++++++++++++++++++++++++++++"<<endl;
	//------------ R e a d    D a t a ------------    
	cout <<"Reading data...."<<endl;
	for(int i=0;i<Nsets;i++){
		for(int idtyp=0; idtyp<2; idtyp++){ // 0 = AA, 1 = pp  //(*fCentralityBinBorders)[i+1]
			for(int ic=0; ic<NumCent[idtyp]; ic++){
				for(int iptt=0; iptt<NPTT; iptt++){
					for(int ipta=0;ipta<NPTA;ipta++) {
						hDeltaEtaSig[i][idtyp][ic][iptt][ipta] = (TH1D *)fin[i]->Get(Form("hDeltaEtaSig%02dC%02dT%02dA%02d",idtyp,ic,iptt,ipta));
						if(idtyp==AA) hIAADeltaEtaSig[i][ic][iptt][ipta] = (TH1D *)fin[i]->Get(Form("hIAADeltaEtaSigC%02dT%02dA%02d",ic,iptt,ipta));
					} // ipta
				} // iptt 
			} // ic
		} // pp or AA
	} // iset
	// Calculation Various Ratios
	// 1. Out/In ratio for each data setfor(int i=0;i<Nsets;i++){
	for(int ic=0; ic<NumCent[AA]; ic++){
		for(int iptt=0; iptt<NPTT; iptt++){
			for(int ipta=0;ipta<NPTA;ipta++) {
				for(int i=0;i<Nsets;i++){
					hRatio_DeltaEtaSig[i][ic][iptt][ipta] = (TH1D*)hDeltaEtaSig[i][AA][ic][iptt][ipta]->Clone();
					hRatio_DeltaEtaSig[i][ic][iptt][ipta]->Divide(hDeltaEtaSig[iRef][AA][ic][iptt][ipta]);

					// for systematics
					hRatio_IAA[i][ic][iptt][ipta] = (TH1D*)hIAADeltaEtaSig[i][ic][iptt][ipta]->Clone();
					hRatio_IAA[i][ic][iptt][ipta]->Divide(hIAADeltaEtaSig[iRef][ic][iptt][ipta]);
					// for drawing purpose in this macro
					hRatio_IAA_fordrawing[i][ic][iptt][ipta] = (TH1D*)hIAADeltaEtaSig[i][ic][iptt][ipta]->Clone();
					hRatio_IAA_fordrawing[i][ic][iptt][ipta]->Divide(hIAADeltaEtaSig[iRef][ic][iptt][ipta]);
				}
			} // ipta
		} // iptt 
	} // ic
}

//------------------------------------------------------------------------------------------------
void Compare(){
    ObtainSyst();
	int ic=0;
	for(int iptt=3; iptt<NPTT; iptt++){
		for(int ipta=3;ipta<NPTA-1;ipta++) {
				DrawIAA(ic++,iptt, ipta);
	 	}
	}
	SaveFinalResults();
}


//------------------------------------------------------------------------------------------------
void DrawIAA(int padID, int iPTT, int iPTA) {
	Filipad *fpad[NumCent[AA]];
	lowx = -0.01;
	for(int ic=0;ic<NumCent[AA];ic++) {
		fpad[ic] = new Filipad(padID+ic+1, 1.1, 0.5, 100, 100, 0.7,NumCent[AA]);
		fpad[ic]->Draw();
		//==== Upper pad
		TPad *p = fpad[ic]->GetPad(1); //upper pad
		p->SetTickx(); p->SetLogx(0); p->SetLogy(0); p->cd();
		TH2F *hfr = new TH2F("hfr"," ", 100,lowx, highx, 10, ly, hy); // numbers: tics x, low limit x, upper limit x, tics y, low limit y, upper limit y
		hset( *hfr, "|#Delta#eta|", "IAA(|#Delta#eta|)",1.1,1.0, 0.09,0.09, 0.01,0.01, 0.04,0.05, 510,505);//settings of the upper pad: x-axis, y-axis
		hfr->Draw();
		//Legend definition
		TLegend *leg = new TLegend(0.45,0.4,0.85,0.78,"","brNDC");
		leg->SetTextSize(0.037);leg->SetBorderSize(0);leg->SetFillStyle(0);//legend settings;

		latexRun.DrawLatexNDC( 0.25, 0.85 ,strRun);

		leg->AddEntry((TObject*)NULL,hIAADeltaEtaSig[0][ic][iPTT][iPTA]->GetTitle(),"");

		// add systematics
		grIAADeltaEtaSig_Totsyst[ic][iPTT][iPTA]->SetFillStyle(fillStylesTheory[0]);
		grIAADeltaEtaSig_Totsyst[ic][iPTT][iPTA]->SetFillColor(fillColors[0]);
		grIAADeltaEtaSig_Totsyst[ic][iPTT][iPTA]->Draw("same2");
		
		for(int iS=0;iS<Nsets;iS++) {
	

			hIAADeltaEtaSig[iS][ic][iPTT][iPTA]->SetMarkerStyle(gMarkers[iS]);
			hIAADeltaEtaSig[iS][ic][iPTT][iPTA]->SetMarkerColor(gColors[iS]);
			hIAADeltaEtaSig[iS][ic][iPTT][iPTA]->SetLineColor(gColors[iS]);
			hIAADeltaEtaSig[iS][ic][iPTT][iPTA]->Draw("p,same");
			leg->AddEntry(hIAADeltaEtaSig[iS][ic][iPTT][iPTA],sLeg[iS],"pl");

		}
		
		leg->Draw();

		//==== Lower pad
		p = fpad[ic]->GetPad(2);
		p->SetTickx(); p->SetGridy(1); p->SetLogx(0), p->SetLogy(0); p->cd();
		TH2F *hfr1 = new TH2F("hfr1"," ", 100, lowx, highx, 10, rlow, rhigh);
		hset( *hfr1, "|#Delta#eta|", "Ratio to Default",1.1,1.0, 0.09,0.09, 0.01,0.01, 0.08,0.08, 510,505);
		hfr1->Draw();
		for(int i=0;i<Nsets;i++) {
			hSystSmoothRelativeSyst[i][ic][iPTT][iPTA]->SetLineColor(gColors[i]);
			hSystSmoothRelativeSyst[i][ic][iPTT][iPTA]->Draw("histosame");
			hRatio_IAA_fordrawing[i][ic][iPTT][iPTA]->SetMarkerStyle(gMarkers[i]);
			hRatio_IAA_fordrawing[i][ic][iPTT][iPTA]->SetMarkerColor(gColors[i]);
			hRatio_IAA_fordrawing[i][ic][iPTT][iPTA]->SetLineColor(gColors[i]);
			hRatio_IAA_fordrawing[i][ic][iPTT][iPTA]->Draw("p,same");		
		}
		//gPad->GetCanvas()->SaveAs(Form("figs_syst/IAA_C%02dT%02dA%02d.pdf",ic,iPTT,iPTA));
	}
	//for(int ic=0;ic<NumCent[AA];ic++) delete fpad[ic];
}

//------------------------------------------------------------------------------------------------
void BarlowTest(){
 LoadData();
 for(int i=0;i<Nsets;i++){
	hBarlowDist[i] = new TH1D(Form("hBarlowDist_%d",i),Form("hBarlowDist_%d",i),200,-10,10);
 }
 for(int ic=0; ic<NumCent[AA]; ic++){
	for(int iptt=0; iptt<NPTT; iptt++){
        	for(int ipta=0;ipta<NPTA;ipta++) {
                	for(int i=0;i<Nsets;i++){
				Barlow(  hIAADeltaEtaSig[iRef][ic][iptt][ipta],  hIAADeltaEtaSig[i][ic][iptt][ipta], hBarlowDist[i], 1 ); 			
			}
		}
	}
 }
}

void ObtainSyst(){
 LoadData();
 BarlowTest();
 for(int ic=0; ic<NumCent[AA]; ic++){
    for(int iptt=0; iptt<NPTT; iptt++){
        for(int ipta=0;ipta<NPTA;ipta++) {
            for(int i=0;i<Nsets;i++){
				for(int l=0;l<hRatio_IAA[i][ic][iptt][ipta]->GetNbinsX();l++){
					hRatio_IAA[i][ic][iptt][ipta]->AddBinContent( l+1, -1.0 );
					hRatio_IAA[i][ic][iptt][ipta]->SetBinContent( l+1, fabs( hRatio_IAA[i][ic][iptt][ipta]->GetBinContent(l+1) ) );
				}
				hSystSmooth[i][ic][iptt][ipta] = (TH1D*)Smooth( hRatio_IAA[i][ic][iptt][ipta] );
				hIAADeltaEtaSig_syst[i][ic][iptt][ipta] = (TH1D*)GetPlotWithSyst( hIAADeltaEtaSig[iRef][ic][iptt][ipta], hSystSmooth[i][ic][iptt][ipta] );
				grIAADeltaEtaSig_syst[i][ic][iptt][ipta] = (TGraphErrors*)h2g(hIAADeltaEtaSig_syst[i][ic][iptt][ipta]);
				RemovePoints(grIAADeltaEtaSig_syst[i][ic][iptt][ipta],0.01,0.27);
				hSystSmoothRelativeSyst[i][ic][iptt][ipta]= (TH1D*)GetRelativeSyst(hSystSmooth[i][ic][iptt][ipta]);
			}
		}
	}
 }
 ObtainTotalSyst(); // calculate total systematics
}
// calculate total systematics
void ObtainTotalSyst(){
	for(int ic=0; ic<NumCent[AA]; ic++){
    	for(int iptt=0; iptt<NPTT; iptt++){
        	for(int ipta=0;ipta<NPTA;ipta++) {
        		hSystSmoothTotal[ic][iptt][ipta] = (TH1D*)hSystSmooth[0][ic][iptt][ipta]->Clone();
            	hSystSmoothTotal[ic][iptt][ipta]->Reset();
 				for(int i=1;i<=hSystSmoothTotal[ic][iptt][ipta]->GetNbinsX();i++){
 					double totsys = 0.;
 					for(int j=1;j<Nsets;j++){
						if( hBarlowDist[j]->GetRMS() < 1.0 ) continue;
 						totsys =+ hSystSmooth[j][ic][iptt][ipta]->GetBinContent(i)*hSystSmooth[j][ic][iptt][ipta]->GetBinContent(i);
 					}
        			hSystSmoothTotal[ic][iptt][ipta]->SetBinContent(i,TMath::Sqrt(totsys));
        		} // end of histo binx
        		hIAADeltaEtaSig_Totsyst[ic][iptt][ipta] = (TH1D*)GetPlotWithSyst( hIAADeltaEtaSig[iRef][ic][iptt][ipta], hSystSmoothTotal[ic][iptt][ipta] );
        		grIAADeltaEtaSig_Totsyst[ic][iptt][ipta] = (TGraphErrors*)h2g(hIAADeltaEtaSig_Totsyst[ic][iptt][ipta]);
        		RemovePoints(grIAADeltaEtaSig_Totsyst[ic][iptt][ipta],0.01,0.27);
        		grIAADeltaEtaSig[ic][iptt][ipta] = (TGraphErrors*)h2g( hIAADeltaEtaSig[iRef][ic][iptt][ipta]);
        		RemovePoints(grIAADeltaEtaSig[ic][iptt][ipta],0.01,0.27);
			}
		}
 	}
}

void SaveFinalResults(){

	TFile *fout = new TFile("results/Iaa_PbPb5.02TeV_results.root","recreate");
	for(int ic=0; ic<NumCent[AA]; ic++){
    	for(int iptt=0; iptt<NPTT; iptt++){
        	for(int ipta=0;ipta<NPTA;ipta++) {
        		grIAADeltaEtaSig[ic][iptt][ipta]->Write(Form("grIAADeltaEtaSigC%02dT%02dA%02d",ic,iptt,ipta));
        		grIAADeltaEtaSig_Totsyst[ic][iptt][ipta]->Write(Form("grIAADeltaEtaSigC%02dT%02dA%02d_syst",ic,iptt,ipta));
			}
		}
 	}
}

TGraphErrors* h2g(TH1D *hid){
	int offstart=1;
    int    NC =  hid->GetNbinsX()+offstart;
    double binx[NC];
    double ebinx[NC];
    double co[NC];
    double eco[NC];

    for(int i=offstart; i<NC; i++) {
    	binx[i] = hid->GetBinCenter(i);
    	ebinx[i] = 0.01;
    	co[i] = hid->GetBinContent(i);
    	//cout << i<<"\t"<< hid->GetXaxis()->GetBinLabel(i) <<endl;
    	eco[i]=hid->GetBinError(i);
    }
    TGraphErrors *grr = new TGraphErrors(NC,binx,co,ebinx,eco);
    grr->SetTitle(hid->GetTitle());

return grr;
}

void RemovePoints(TGraphErrors *ge, double xmin, double xmax)
{
  // Remove zero points from TGraphErrors.

  if(!ge){return;}

  Int_t nPoints = ge->GetN();
  Double_t x = 0.;
  Double_t y = 0.;
  int p =0;
  while(p<nPoints) {
    ge->GetPoint(p,x,y);
    if( x < xmin || x > xmax ) // Npart cent < 60%
    {
      ge->RemovePoint(p);
      //            cout<<Form(" WARNING (%s): point %d is < 1.e-10 and it was removed from the plot !!!!",ge->GetName(),p+1)<<endl;
      nPoints = ge->GetN();
    } else {
      p++;
    }
  } // end of for(Int_t p=0;p<nPoints;p++)

  //    cout<<endl;
  return;

} // end of void RemoveZeroPoints(TGraphErrors *ge)

