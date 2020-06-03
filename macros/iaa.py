
import numpy as np
import ROOT

import scipy
from scipy import interpolate

import sys
sys.path.append("JPyPlotRatio");
#sys.path.append("/home/jasper/Asiakirjat/projects/JPyPlotRatio");


import JPyPlotRatio



#fData    = ROOT.TFile("sysErrors/Signal_LHC15o_GlobalSDD_JCIAA_GlobalSDD_LHC17p_pass1_CENT_woSDD_Iaa_R0.2_1.0_1.60_Near_Wing0.root","read");
fData    = ROOT.TFile("sysErrors/Signal_LHC15o_pass1_CentralBarrelTracking_hadronPID_FieldConfigs_5146_JCIAA_GlobalSDD_LHC17p_pass1_CENT_woSDD_Iaa_R0.2_1.0_1.60_Near_Wing0.root","read");
fMarton    = ROOT.TFile("results/Fianl_Marton_graphs.root","read");

Modelfiles = [
			  "sysErrors/Signal_JEWEL_JCIaa_KineOnly_JEWEL_vacuum_Iaa_R0.2_1.0_1.60_Near_Wing0.root",
			  "sysErrors/Signal_AMPT_LHC13f3c_JCIAA_EPInclusive_LHC12f1a_Pythia_2760GeV_KineOnly_Iaa_R0.2_1.0_1.60_Near_Wing0.root",
#    		  "sysErrors/Signal_AMPT_LHC13f3a_JCIAA_EPInclusive_LHC12f1a_Pythia_2760GeV_KineOnly_Iaa_R0.2_1.0_1.60_Near_Wing0.root",
#			  "sysErrors/Signal_MCGen_PbPb_AMPT_5TeV_modPars2_JCIaa_KineOnly_MCGen_pp_amptpp_Iaa_R0.2_1.0_1.60_Near_Wing0.root"
			];


fModel = [ROOT.TFile(elm) for elm in Modelfiles];

ModelLabel = [
"~~~JEWEL",
"~~~AMPT String Melting",
#"~~~AMPT default",
#"AMPT ST Fly"
];


#fJewel   = ROOT.TFile("sysErrors/Signal_JEWEL_JCIaa_KineOnly_JEWEL_vacuum_Iaa_R0.2_1.0_1.60_Near_Wing0.root","read");
#fAmpt    = ROOT.TFile("sysErrors/Signal_AMPT_LHC13f3c_JCIAA_EPInclusive_pythia8230_pp2.76TeV_GF0_SoftQCD_Iaa_R0.2_1.0_1.60_Near_Wing0.root","read");
#fAmpt    = ROOT.TFile("sysErrors/sysErrors/Signal_AMPT_LHC13f3c_JCIAA_EPInclusive_pythia8230_pp2.76TeV_GF0_SoftQCD_Iaa_R0.2_1.0_1.60_Near_Wing0.root","read");
dataTypePlotParams = [
	{'plotType':'data','color':'r','fmt':'o','markersize':5.0},
	{'plotType':'data','color':'k','fmt':'s','markersize':5.0},
	{'plotType':'theory','facecolor':'C0','edgecolor':'C0','alpha':0.5,'linestyle':'solid','linecolor':'C0'},
	{'plotType':'theory','facecolor':'C1','edgecolor':'C1','alpha':0.5,'linestyle':'dotted','linecolor':'C0'},
	{'plotType':'theory','facecolor':'C2','edgecolor':'C2','alpha':0.5,'linestyle':'dashed','linecolor':'C0'},
	{'plotType':'theory','facecolor':'C3','edgecolor':'C3','alpha':0.5,'linestyle':'dashdot','linecolor':'C0'},
	{'plotType':'data','color':'k','fmt':'o','fillstyle':'none','markersize':5.0} #PP
];


# define panel/xaxis limits/titles
nrow = 1;
ncol = 2;
xlimits = [(0.005,0.23),(0.01,0.23)];
ylimits = [(0.4,2.8)];
rlimits = [(0.0,3.0)];

centrality =["0-5\%","5-10\%"];
TriggPtBorders = [8,15];
AssocPtBorders = [4,6,8];
pttN = len(TriggPtBorders);
ptaN = len(AssocPtBorders);
startPttBin = 3; #{3,4,6,8,15} 
ptaBins = [4,5];#{0.6,1,2,3,4,6,8,10}
xtitle = ["$|\\Delta\\eta|$"];
ytitle = ["$I_{AA}$"];
plabelptt = {i: "{}$<p_{{Tt}}<$ {} GeV/$c$".format(TriggPtBorders[i],TriggPtBorders[i+1]) for i in range(0,pttN-1)};
plabelpta = {i: "{}$<p_{{Tt}}<$ {} $\\otimes$ {}$<p_{{Ta}}<$ {}".format(TriggPtBorders[0],TriggPtBorders[1],AssocPtBorders[i],AssocPtBorders[i+1]) for i in range(0,ptaN-1)};
# Following two must be added
toptitle = "PbPb $\\sqrt{s_{NN}}$ = 5.02 TeV"; # need to add on the top
xlong = 0;

dataDetail = "$|\\eta| < 0.8$";


plot = JPyPlotRatio.JPyPlotRatio(panels=(nrow,ncol),
	rowBounds=ylimits,  # for nrow
	colBounds=xlimits,  # for ncol
	panelLabel=plabelpta,  # nrowxncol
	ratioBounds=rlimits,# for nrow
	panelLabelLoc=(0.12,0.92),panelLabelSize=10,panelLabelAlign="left",
	legendPanel=1,
	legendLoc=(0.40,0.65),legendSize=9,xlabel=xtitle[0],ylabel=ytitle[0]);



plot.EnableLatex(True);

plotMatrix = np.empty((pttN,ptaN),dtype=int);

for it in range(0,pttN-1):
	for ia in range(0,ptaN-1):
		print("{d}",ptaBins[ia]);
		grData = fData.Get("hIAADeltaEtaSigC{:02d}T{:02d}A{:02d}".format(0,startPttBin+it,ptaBins[ia]));
		plotMatrix[it,ia] = plot.AddTH1(ia,grData,**dataTypePlotParams[0],label="ALICE, 0-5\%");
		grMarton = fMarton.Get("grIAADeltaEtaC{:02d}T{:02d}A{:02d}".format(0,startPttBin+it,ptaBins[ia]));
		grMarton.Print();
		plotMatrixMarton = plot.Add(ia,grMarton,**dataTypePlotParams[1],label="2.76TeV");
		gr_sys = fMarton.Get("grAsymmIAADeltaEtaSystPointByPointC{:02d}T{:02d}A{:02d}".format(0,startPttBin+it,ptaBins[ia]));
		gr_sys.Print();
		plot.AddSyst(plotMatrixMarton,gr_sys);
		for im in range(len(Modelfiles)):
			fModel[im].Print();	
			grm = fModel[im].Get("hIAADeltaEtaSigC{:02d}T{:02d}A{:02d}".format(0,startPttBin+it,ptaBins[ia]));
			pm = plot.AddTH1(ia,grm,**dataTypePlotParams[im+2],label=ModelLabel[im]);	
			plot.Ratio(pm,plotMatrix[it,ia],style="default"); #Calculate and plot ratio between data and theory

fData.Close();

plot.GetPlot().text(0.2,0.75,toptitle,fontsize=9);
#plot.GetPlot().text(0.25,0.70,dataDetail,fontsize=9);
#plot.GetPlot().text(0.23,0.77,strXlong[xlong],fontsize=9);
#plot.GetAxes(3).text(0.1,0.1,dataDetail,fontsize=9);

plot.Plot();

#plot.GetRatioAxes(3).remove();

#plot.Save("figs/Fig1_jtdist.pdf");
plot.Show();

