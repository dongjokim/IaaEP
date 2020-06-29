
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
fMarton    = ROOT.TFile("results/Final_Marton_graphs.root","read");

Modelfiles = [
			  "sysErrors/Signal_JEWEL_JCIaa_KineOnly_JEWEL_vacuum_Iaa_R0.2_1.0_1.60_Near_Wing0.root",
			  "sysErrors/Signal_AMPT_LHC13f3c_JCIaa_KineOnly_pythia8230_pp5.02TeV_GF0_SoftQCD_Iaa_R0.2_1.0_1.60_Near_Wing0.root"
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
	{'plotType':'data','color':'k','fmt':'s','markersize':5.0},
	{'plotType':'data','color':'r','fmt':'o','markersize':5.0},
	{'plotType':'data','color':'b','fmt':'P','markersize':5.0},
	{'plotType':'data','color':'m','fmt':'X','markersize':5.0},
	{'plotType':'data','color':'k','fmt':'*','markersize':5.0},
	{'plotType':'data','color':'k','fmt':'s','fillstyle':'none','markersize':5.0},
	{'plotType':'data','color':'r','fmt':'o','fillstyle':'none','markersize':5.0},
	{'plotType':'data','color':'b','fmt':'P','fillstyle':'none','markersize':5.0},
	{'plotType':'data','color':'m','fmt':'X','fillstyle':'none','markersize':5.0},
	{'plotType':'data','color':'k','fmt':'o','fillstyle':'none','markersize':5.0} 
];


# define panel/xaxis limits/titles
nrow = 2;
ncol = 2;
xlimits = [(0.005,0.28),(0.005,0.28)];
ylimits = [(0.3,1.8),(0.3,1.8)];
rlimits = [(0.0,3.0)];

histnames = [
			 ["T02A03","T02A04"],  # check it with ROOT file Title
			 ["T03A04","T03A05"]#,"0_6_11"]
		    ];

plables = [
			"$6.0 <p_{{Tt}}< 8.0 \\otimes 3.0<p_{{Ta}}< 4.0$","$6.0 <p_{{Tt}}< 8.0 \\otimes 4.0<p_{{Ta}}< 6.0$",
			"$8.0 <p_{{Tt}}< 15.0 \\otimes 4.0<p_{{Ta}}< 6.0$","$8.0 <p_{{Tt}}< 15.0 \\otimes 6.0<p_{{Ta}}< 8.0$",
		 ];

centrality =["0-5\%","5-10\%","10-20\%","20-40\%","40-60\%"];#,"60-80\%"];
xtitle = ["$|\\Delta\\eta|$"];
ytitle = ["$I_{AA}$"];
# Following two must be added
toptitle = "PbPb $\\sqrt{s_{NN}}$ = 2.76 TeV"; # need to add on the top

dataDetail = "$|\\eta| < 0.8$";


plot = JPyPlotRatio.JPyPlotRatio(panels=(nrow,ncol),
	rowBounds=ylimits,  # for nrow
	colBounds=xlimits,  # for ncol
	panelLabel=plables,  # nrowxncol
	ratioBounds=rlimits,# for nrow
	disableRatio=[0,1], # disable ratio..
	panelLabelLoc=(0.06,0.90),panelLabelSize=10,panelLabelAlign="left",
	legendPanel=0,
	legendLoc=(0.20,0.30),legendSize=9,xlabel=xtitle[0],ylabel=ytitle[0]);



plot.EnableLatex(True);

for i in range(0,nrow):
	for j in range(0,ncol):
		index = i*ncol+j; # for each panel
		plot.GetAxes(index).set_xticks([0,0.1,0.2]); 
		for ic in range(len(centrality)):
			#grData = fData.Get("grIAADeltaEtaSigC{:02d}{}".format(ic,histnames[i][j]));
			#plotMatrix = plot.Add(index,grData,**dataTypePlotParams[ic],label=centrality[ic]);
			#grData_sys = fData.Get("grIAADeltaEtaSigC{:02d}{}_syst".format(ic,histnames[i][j]));
			#_,_,_,syst = JPyPlotRatio.TGraphErrorsToNumpy(ROOT.TGraphErrors(grData_sys));
			#plot.AddSyst(plotMatrix,syst);
			grMarton = fMarton.Get("grIAADeltaEtaC{:02d}{}".format(ic,histnames[i][j]));
			#grMarton.Print();
			plotMatrixMarton = plot.Add(index,grMarton,**dataTypePlotParams[ic],label=centrality[ic]);
			gr_sys = fMarton.Get("grAsymmIAADeltaEtaSystPointByPointC{:02d}{}".format(ic,histnames[i][j]));
			#gr_sys.Print();
			plot.AddSyst(plotMatrixMarton,gr_sys);
			#plot.Ratio(plotMatrixMarton,plotMatrix,style="default"); #Calculate and plot ratio between data and theory
		
fData.Close();

plot.GetPlot().text(0.29,0.53,toptitle,fontsize=9);
#plot.GetPlot().text(0.25,0.70,dataDetail,fontsize=9);
#plot.GetPlot().text(0.23,0.77,strXlong[xlong],fontsize=9);
#plot.GetAxes(3).text(0.1,0.1,dataDetail,fontsize=9);

plot.Plot();

#plot.GetRatioAxes(3).remove();

#plot.Save("figs/Fig1_jtdist.pdf");
plot.Show();

