const int NC = 6; 
const int NPTT = 2;
const int NPTA = 2;
int ptStart = 2; // from 6-8
int ptaBins[NPTT][NPTA] = {{3,4},{4,5}};
TH1D *hIAADeltaEta[NC][NPTT][NPTA]; 
TGraphErrors *grIAADeltaEta[NC][NPTT][NPTA];

TGraphAsymmErrors *grAsymmIAADeltaEtaSystPointByPoint[NC][NPTT][NPTA];
TGraphAsymmErrors *grAsymmIAADeltaEtaSystScaling[NC][NPTT][NPTA];
TBox *boxIAADeltaEtaSystScaling[NC][NPTT][NPTA];
void LoadData();

TGraphErrors* h2g(TH1D *hid);
void RemovePoints(TGraphErrors *ge, double xmin, double xmax);

void LoadData(){

    cout <<"+++++++++++++++++++++++++++++++++++++++++++++++"<<endl;
	//------------ R e a d    D a t a ------------    
	cout <<"Reading data...."<<endl;
	TFile *fin = TFile::Open("Final_Marton.root");
  for(int ic=0;ic<NC;ic++) {
    for(int it=0;it<NPTT;it++) {
      for(int ia=0;ia<NPTA;ia++) {
        hIAADeltaEta[ic][it][ia] = (TH1D*)fin->Get(Form("hIAADeltaEtaC%02dT%02dA%02d",ic,it+ptStart,ptaBins[it][ia]));
        grAsymmIAADeltaEtaSystPointByPoint[ic][it][ia] = (TGraphAsymmErrors*)fin->Get(Form("hIAADeltaEtaSystPointByPointC%02dT%02dA%02d",ic,it+ptStart,ptaBins[it][ia]));
        boxIAADeltaEtaSystScaling[ic][it][ia] = (TBox*)fin->Get(Form("hIAADeltaEtaSystScalingC%02dT%02dA%02d",ic,it+ptStart,ptaBins[it][ia]));
      }
    }
  }
		
}

void WriteData(){

	cout <<"+++++++++++++++++++++++++++++++++++++++++++++++"<<endl;
	TFile *fout = new TFile("Fianl_Marton_graphs.root","recreate");
	//------------ R e a d    D a t a ------------    
	cout <<"Writing data...."<<endl;
	fout->cd();

  for(int ic=0;ic<NC;ic++) {
    for(int it=0;it<NPTT;it++) {
      for(int ia=0;ia<NPTA;ia++) {
        cout << Form("C%02dT%02dA%02d",ic,it+ptStart,ptaBins[it][ia]) << endl;
        grIAADeltaEta[ic][it][ia] = (TGraphErrors*)h2g(hIAADeltaEta[ic][it][ia]);
        RemovePoints(grIAADeltaEta[ic][it][ia],0.01,0.27);
        //grIAADeltaEta[ic][it][ia]->Print();
        grIAADeltaEta[ic][it][ia]->Write(Form("grIAADeltaEtaC%02dT%02dA%02d",ic,it+ptStart,ptaBins[it][ia]));
        grAsymmIAADeltaEtaSystPointByPoint[ic][it][ia]->Write(Form("grAsymmIAADeltaEtaSystPointByPointC%02dT%02dA%02d",ic,it+ptStart,ptaBins[it][ia]));
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
    	ebinx[i] = 0.50;
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