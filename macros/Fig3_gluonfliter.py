
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
	"sysErrors/Signal_AMPT_LHC13f3c_JCIAA_EPInclusive_pythia8230_pp2.76TeV_GF0_SoftQCD_Iaa_R0.2_1.0_1.60_Near_Wing0.root",
	"sysErrors/Signal_AMPT_LHC13f3c_JCIAA_EPInclusive_pythia8230_pp2.76TeV_GF1_SoftQCD_Iaa_R0.2_1.0_1.60_Near_Wing0.root",
	"sysErrors/Signal_AMPT_LHC13f3c_JCIAA_EPInclusive_pythia8230_pp2.76TeV_QF1_SoftQCD_Iaa_R0.2_1.0_1.60_Near_Wing0.root"
			];


fModel = [ROOT.TFile(elm) for elm in Modelfiles];

ModelLabel = [
	"SoftQCD",
	"SoftQCD, Gluon Filtering",
	"SoftQCD, Quark Filtering"
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
xlimits = [(0.005,0.6),(0.005,0.6)];
ylimits = [(-0.1,0.85),(-0.1,0.85)];
rlimits = [(0.0,2.0),(0.0,2.0)];

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
ytitle = ["$\\frac{1}{N_{trigg}} \\frac{dN}{d|\\Delta\\eta|}$"];

# Following two must be added
toptitle = "PYTHIA8 pp $\\sqrt{s}$ = 2.76 TeV"; # need to add on the top
#toptitle = "PbPb ALICE"
dataDetail = "$|\\eta| < 0.8$";


plot = JPyPlotRatio.JPyPlotRatio(panels=(nrow,ncol),
	rowBounds=ylimits,  # for nrow
	colBounds=xlimits,  # for ncol
	panelLabel=plables,  # nrowxncol
	ratioBounds=rlimits,# for nrow
	#disableRatio=[0,1], # disable ratio..
	#ratioSystPlot=True,
	panelLabelLoc=(0.06,0.90),panelLabelSize=10,panelLabelAlign="left",
	legendPanel=0,
	legendLoc=(0.5,0.65),
	legendSize=9,xlabel=xtitle[0],ylabel=ytitle[0]);



plot.EnableLatex(True);

plotMatrix = np.empty((nrow,ncol),dtype=int);
plotMatrixModel = np.empty((len(Modelfiles)),dtype=int);

for i in range(0,nrow):
	for j in range(0,ncol):
		index = i*ncol+j; # for each panel 
		plot.GetAxes(index).set_xticks([0,0.1,0.2,0.3,0.4,0.5]);
		for im in range(len(Modelfiles)):
			grm = fModel[im].Get("hDeltaEtaSig01{}".format(histnames[i][j]));
			plotMatrixModel[im] = plot.AddTH1(index,grm,**dataTypePlotParams[im+2],label=ModelLabel[im]);
			if(im>0):	
				plot.Ratio(plotMatrixModel[im],plotMatrixModel[0],style="default"); #Calculate and plot ratio between data and theory

fData.Close();

plot.GetPlot().text(0.6,0.80,toptitle,fontsize=9);
plot.GetPlot().text(0.6,0.78,dataDetail,fontsize=9);
#plot.GetPlot().text(0.23,0.77,strXlong[xlong],fontsize=9);
#plot.GetAxes(3).text(0.1,0.1,dataDetail,fontsize=9);
#for i in range(4):
#	plot.GetAxes(i).plot([0,2.6], [1.0,1.0], linestyle='--',color='k');
	
plot.Plot();

#plot.GetRatioAxes(3).remove();

plot.Save("figs/Fig3_gluonfilter.pdf");
plot.Save("figs/Fig3_gluonfilter.png");
plot.Show();

