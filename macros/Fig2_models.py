
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
			  "sysErrors/Signal_JEWEL_v2.0.2_sqrts2760_0_5_JCIaa_KineOnly_JEWEL_vacuum_Iaa_R0.2_1.0_1.60_Near_Wing0.root",
			  "sysErrors/Signal_JEWEL_v2.0.2_KeepRecoil_sqrts2760_0_5_JCIaa_KineOnly_JEWEL_vacuum_Iaa_R0.2_1.0_1.60_Near_Wing0.root",
			  "sysErrors/Signal_AMPT_LHC13f3c_JCIaa_KineOnly_pythia8230_pp5.02TeV_GF0_SoftQCD_Iaa_R0.2_1.0_1.60_Near_Wing0.root"
			  #"sysErrors/Signal_JEWEL_v2.0.2_sqrts2760_JCIaa_KineOnly_JEWEL_vacuum_Iaa_R0.2_1.0_1.60_Near_Wing0.root",
			  #"sysErrors/Signal_JEWEL_v2.0.2_KeepRecoil_sqrts2760_JCIaa_KineOnly_JEWEL_vacuum_Iaa_R0.2_1.0_1.60_Near_Wing0.root",
#			  "sysErrors/Signal_GG_JEWEL_v2.0.2_KeepRecoil_JCIaa_KineOnly_JEWEL_vacuum_Iaa_R0.2_1.0_1.60_Near_Wing0.root", # no difference GG
			  #"sysErrors/Signal_GG_JEWEL_JCIaa_KineOnly_JEWEL_vacuum_Iaa_R0.2_1.0_1.60_Near_Wing0.root",

			];


fModel = [ROOT.TFile(elm) for elm in Modelfiles];

ModelLabel = [
"~~~JEWEL",
"~~~JEWEL with Recoils",
"~~~AMPT String Melting",
#"~~~JEWEL 2.76 0-10\%",
#"~~~JEWEL with Recoils 2,76 0-10\%",
#"~~~AMPT String Melting GG",
#"~~~AMPT default"
#,"AMPT ST Fly"
];


#fJewel   = ROOT.TFile("sysErrors/Signal_JEWEL_JCIaa_KineOnly_JEWEL_vacuum_Iaa_R0.2_1.0_1.60_Near_Wing0.root","read");
#fAmpt    = ROOT.TFile("sysErrors/Signal_AMPT_LHC13f3c_JCIAA_EPInclusive_pythia8230_pp2.76TeV_GF0_SoftQCD_Iaa_R0.2_1.0_1.60_Near_Wing0.root","read");
#fAmpt    = ROOT.TFile("sysErrors/sysErrors/Signal_AMPT_LHC13f3c_JCIAA_EPInclusive_pythia8230_pp2.76TeV_GF0_SoftQCD_Iaa_R0.2_1.0_1.60_Near_Wing0.root","read");
dataTypePlotParams = [
	{'plotType':'data','color':'r','fmt':'s','markersize':5.0},
	{'plotType':'data','color':'k','fmt':'o','markersize':5.0},
	{'plotType':'theory','facecolor':'C0','edgecolor':'C0','alpha':0.5,'linestyle':'solid','linecolor':'C0'},
	{'plotType':'theory','facecolor':'C1','edgecolor':'C1','alpha':0.5,'linestyle':'dotted','linecolor':'C1'},
	{'plotType':'theory','facecolor':'C2','edgecolor':'C2','alpha':0.5,'linestyle':'dashed','linecolor':'C2'},
	{'plotType':'theory','facecolor':'C3','edgecolor':'C3','alpha':0.5,'linestyle':'dashdot','linecolor':'C3'},
	{'plotType':'theory','facecolor':'xkcd:silver','edgecolor':'xkcd:silver','alpha':0.5,'linestyle':'solid','linecolor':'xkcd:silver'},
	{'plotType':'theory','facecolor':'C5','edgecolor':'C5','alpha':0.5,'linestyle':'dashdot','linecolor':'C5'},
	{'plotType':'data','color':'k','fmt':'o','fillstyle':'none','markersize':5.0} #PP
];


# define panel/xaxis limits/titles
nrow = 2;
ncol = 2;
xlimits = [(0.005,0.28),(0.005,0.28)];
ylimits = [(0.1,2.1),(0.1,2.1)];
rlimits = [(0.0,1.5),(0.0,3.0)];

histnames = [
			 ["C00T02A03","C00T02A04"],  # check it with ROOT file Title
			 ["C00T03A04","C00T03A05"]#,"0_6_11"]
			];

iaStart = [3,4];

plables = [
			"$6.0 <p_{{Tt}}< 8.0 \\otimes 3.0<p_{{Ta}}< 4.0$","$6.0 <p_{{Tt}}< 8.0 \\otimes 4.0<p_{{Ta}}< 6.0$",
			"$8.0 <p_{{Tt}}< 15.0 \\otimes 4.0<p_{{Ta}}< 6.0$","$8.0 <p_{{Tt}}< 15.0 \\otimes 6.0<p_{{Ta}}< 8.0$",
		 ];

xtitle = ["$|\\Delta\\eta|$"];
ytitle = ["$I_{AA}$"];

# Following two must be added
#toptitle = "PbPb $\\sqrt{s_{NN}}$ = 5.02 TeV"; # need to add on the top
toptitle = "PbPb 0--5\% $|\\eta| < 0.8$"
dataDetail = "$|\\eta| < 0.8$";


plot = JPyPlotRatio.JPyPlotRatio(panels=(nrow,ncol),
	rowBounds=ylimits,  # for nrow
	colBounds=xlimits,  # for ncol
	panelLabel=plables,  # nrowxncol
	ratioBounds=rlimits,# for nrow
	disableRatio=[0,1], # disable ratio..
	#ratioSystPlot=True,
	panelLabelLoc=(0.06,0.90),panelLabelSize=10,panelLabelAlign="left",
	#legendPanel=0,
	#legendLoc=(0.45,0.25),
	legendPanel={0:0,1:1},legendLoc={0:(0.27,0.18),1:(0.36,0.18)},
	legendSize=9,xlabel=xtitle[0],ylabel=ytitle[0]);



plot.EnableLatex(True);

plotMatrix = np.empty((nrow,ncol),dtype=int);

for i in range(0,nrow):
	for j in range(0,ncol):
		index = i*ncol+j; # for each panel 
		plot.GetAxes(index).set_xticks([0,0.1,0.2]);
		grData = fData.Get("grIAADeltaEtaSig{}".format(histnames[i][j]));
		#plotMatrix[i,j] = plot.Add(index,grData,**dataTypePlotParams[0],labelLegendId=0,label="5.02 TeV");
		#grData_sys = fData.Get("grIAADeltaEtaSig{}_syst".format(histnames[i][j]));
		#_,_,_,syst = JPyPlotRatio.TGraphErrorsToNumpy(ROOT.TGraphErrors(grData_sys));
		#plot.AddSyst(plotMatrix[i,j],syst);		
		grMarton = fMarton.Get("grIAADeltaEta{}".format(histnames[i][j]));
		plotMatrixMarton = plot.Add(index,grMarton,**dataTypePlotParams[1],labelLegendId=0,label="ALICE");
		gr_sys = fMarton.Get("grAsymmIAADeltaEtaSystPointByPoint{}".format(histnames[i][j]));
		#gr_sys.Print();
		plot.AddSyst(plotMatrixMarton,gr_sys);
		gr_sysSC = fMarton.Get("grIAADeltaEtaSystScaling{}".format(histnames[i][j]));
		plotMatrixMartonSCsys = plot.Add(index,gr_sysSC,**dataTypePlotParams[6],labelLegendId=0,label="Scale error");
		for im in range(len(Modelfiles)):
			fModel[im].Print();	
			grm = fModel[im].Get("hIAADeltaEtaSig{}".format(histnames[i][j]));
			pm = plot.AddTH1(index,grm,**dataTypePlotParams[im+2],labelLegendId=1,label=ModelLabel[im]);	
			plot.Ratio(pm,plotMatrixMarton,style="default"); #Calculate and plot ratio between data and theory

fData.Close();

plot.GetPlot().text(0.18,0.80,toptitle,fontsize=9);
#plot.GetPlot().text(0.35,0.80,dataDetail,fontsize=9);
#plot.GetPlot().text(0.23,0.77,strXlong[xlong],fontsize=9);
#plot.GetAxes(3).text(0.1,0.1,dataDetail,fontsize=9);
for i in range(4):
	plot.GetAxes(i).plot([0,2.6], [1.0,1.0], linestyle='--',color='k');
	
plot.Plot();

#plot.GetRatioAxes(3).remove();

plot.Save("figs/Fig2_models.pdf");
plot.Show();

