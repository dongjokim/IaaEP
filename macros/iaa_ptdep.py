
import numpy as np
import ROOT

import scipy
from scipy import interpolate

import sys
sys.path.append("JPyPlotRatio");
#sys.path.append("/home/jasper/Asiakirjat/projects/JPyPlotRatio");


import JPyPlotRatio



#fData    = ROOT.TFile("sysErrors/Signal_LHC15o_GlobalSDD_JCIAA_GlobalSDD_LHC17p_pass1_CENT_woSDD_Iaa_R0.2_1.0_1.60_Near_Wing0.root","read");
fData    = ROOT.TFile("results/Iaa_PbPb5.02TeV_results.root","read");
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
nrow = 2;
ncol = 2;
xlimits = [(0.005,0.23),(0.01,0.23)];
ylimits = [(0.4,2.8),(0.4,2.8)];
rlimits = [(0.0,3.0)];

histnames = [
			 ["C00T02","C01T02"],  # check it with ROOT file Title
			 ["C00T03","C01T03"]#,"0_6_11"]
		    ];

iaStart = [3,4];
ptaLabel = [
			["$3.0 <p_{{Ta}}< 4.0$","$4.0 <p_{{Ta}}< 6.0$"],
			["$4.0 <p_{{Ta}}< 6.0$","$6.0 <p_{{Ta}}< 8.0$"]
			];
plables = [
			"$6.0 <p_{{Tt}}< 8.0$ 0-5\%","$6.0 <p_{{Tt}}< 8.0$ 5-10\%",
			"$8.0 <p_{{Tt}}< 15.0$ 0-5\%","$8.0 <p_{{Tt}}< 15.0$ 5-10\%",
		 ];

centrality =["0-5\%","5-10\%","10-20\%","20-40\%","40-60\%"];

xtitle = ["$|\\Delta\\eta|$"];
ytitle = ["$I_{AA}$"];
# Following two must be added
toptitle = "PbPb $\\sqrt{s_{NN}}$ = 5.02 TeV"; # need to add on the top

dataDetail = "$|\\eta| < 0.8$";


plot = JPyPlotRatio.JPyPlotRatio(panels=(nrow,ncol),
	rowBounds=ylimits,  # for nrow
	colBounds=xlimits,  # for ncol
	panelLabel=plables,  # nrowxncol
	ratioBounds=rlimits,# for nrow
	disableRatio=[0,1], # disable ratio..
	panelLabelLoc=(0.12,0.92),panelLabelSize=10,panelLabelAlign="left",
	legendPanel=1,
	legendLoc=(0.40,0.65),legendSize=9,xlabel=xtitle[0],ylabel=ytitle[0]);



plot.EnableLatex(True);

plotMatrix = np.empty((nrow,ncol),dtype=int);

for i in range(0,nrow):
	for j in range(0,ncol):
		index = i*ncol+j; # for each panel 
		for ia in range(2):
			grData = fData.Get("grIAADeltaEtaSig{}A{:02d}".format(histnames[i][j],iaStart[i]+ia));
			plotMatrix[i,j] = plot.Add(index,grData,**dataTypePlotParams[ia],label=ptaLabel[i][ia]);
			grData_sys = fData.Get("grIAADeltaEtaSig{}A{:02d}_syst".format(histnames[i][j],iaStart[i]+ia));
			_,_,_,syst = JPyPlotRatio.TGraphErrorsToNumpy(ROOT.TGraphErrors(grData_sys));
			plot.AddSyst(plotMatrix[i,j],syst);
		

fData.Close();

plot.GetPlot().text(0.2,0.75,toptitle,fontsize=9);
#plot.GetPlot().text(0.25,0.70,dataDetail,fontsize=9);
#plot.GetPlot().text(0.23,0.77,strXlong[xlong],fontsize=9);
#plot.GetAxes(3).text(0.1,0.1,dataDetail,fontsize=9);

plot.Plot();

#plot.GetRatioAxes(3).remove();

#plot.Save("figs/Fig1_jtdist.pdf");
plot.Show();

