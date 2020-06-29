
import numpy as np
import ROOT

import scipy
from scipy import interpolate

import sys
sys.path.append("JPyPlotRatio");
#sys.path.append("/home/jasper/Asiakirjat/projects/JPyPlotRatio");


import JPyPlotRatio



fData    = ROOT.TFile("sysErrors/Signal_GG_LHC15o_pass1_CentralBarrelTracking_hadronPID_FieldConfigs_829_Hybrid_JCIAA_TPCOnly_LHC17p_pass1_CENT_woSDD_Iaa_R0.2_1.0_1.60_Near_Wing0.root","read");
fMarton    = ROOT.TFile("results/Final_Marton_graphs.root","read");

Modelfiles = {"sysErrors/Signal_JEWEL_v2.0.2_JCIaa_KineOnly_JEWEL_vacuum_Iaa_R0.2_1.0_1.60_Near_Wing0.root",
			  "sysErrors/Signal_JEWEL_v2.0.2_sqrts2760_JCIaa_KineOnly_JEWEL_vacuum_Iaa_R0.2_1.0_1.60_Near_Wing0.root",
			  "sysErrors/Signal_AMPT_LHC13f3c_JCIaa_KineOnly_pythia8230_pp5.02TeV_GF0_SoftQCD_Iaa_R0.2_1.0_1.60_Near_Wing0.root"
			};



fModel = [ROOT.TFile(elm) for elm in Modelfiles];

ModelLabel = ["~~~JEWEL","~~~JEWEL 2.76","~~~AMPT String Melting"];

dataTypePlotParams = [
	{'plotType':'data','color':'r','fmt':'o','markersize':5.0},
	{'plotType':'theory','facecolor':'C0','edgecolor':'C0','alpha':0.5,'linestyle':'solid','linecolor':'C0'},
	{'plotType':'theory','facecolor':'C1','edgecolor':'C1','alpha':0.5,'linestyle':'dotted','linecolor':'C1'},
	{'plotType':'theory','facecolor':'C2','edgecolor':'C2','alpha':0.5,'linestyle':'dashed','linecolor':'C2'},
	{'plotType':'theory','facecolor':'C3','edgecolor':'C3','alpha':0.5,'linestyle':'dashdot','linecolor':'C3'},
	{'plotType':'theory','facecolor':'xkcd:silver','edgecolor':'xkcd:silver','alpha':0.5,'linestyle':'solid','linecolor':'xkcd:silver'},
	{'plotType':'theory','facecolor':'C5','edgecolor':'C5','alpha':0.5,'linestyle':'dashdot','linecolor':'C5'},
];


# define panel/xaxis limits/titles
nrow = 2;
ncol = 2;
xlimits = [(0,0.3),(0,0.3),(0,0.3),(0,0.3)];
ylimits = [(-0.5,7.1),(-0.5,7.1)];
rlimits = [(0.1,20.0),(0.1,20.0)];

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
ytitle = ["$\\frac{1}{N_{trigg}}\\frac{dN}{|\\Delta\\eta|}$"];

# Following two must be added
toptitle = "PbPb $\\sqrt{s_{NN}}$ = 5.02 TeV"; # need to add on the top

dataType =['AA', 'pp' ];
print(list(enumerate(dataType, start = 0)));
dataDetail = "$|\\eta| < 0.8$";

plot = JPyPlotRatio.JPyPlotRatio(panels=(nrow,ncol),
	rowBounds=ylimits,  # for nrow
	colBounds=xlimits,  # for ncol
	panelLabel=plables,  # nrowxncol
	ratioBounds=rlimits,# for nrow
	panelLabelLoc=(0.12,0.92),panelLabelSize=10,panelLabelAlign="left",
	legendPanel=1,
	legendLoc=(0.40,0.65),legendSize=9,xlabel=xtitle[0],ylabel=ytitle[0]);



plot.EnableLatex(True);

plotMatrix = np.empty((nrow,ncol),dtype=int);

isy = 0; #AA
for i in range(0,nrow):
	for j in range(0,ncol):
		index = i*ncol+j; # for each panel 
		grData = fData.Get("hDeltaEtaSig{:02d}{}".format(isy,histnames[i][j]));
		grData.Print();
		plotMatrix[i,j] = plot.AddTH1(index,grData,**dataTypePlotParams[0],label="ALICE, 0-5\%");
		for im in range(len(Modelfiles)):	
			grm = fModel[im].Get("hDeltaEtaSig{:02d}{}".format(isy,histnames[i][j]));
			grm.Print();
			pm = plot.AddTH1(index,grm,**dataTypePlotParams[im+1],label=ModelLabel[im]);	
			plot.Ratio(pm,plotMatrix[i,j],style="default"); #Calculate and plot ratio between data and theory

fData.Close();

plot.GetPlot().text(0.2,0.75,toptitle,fontsize=9);
#plot.GetPlot().text(0.25,0.70,dataDetail,fontsize=9);
#plot.GetPlot().text(0.23,0.77,strXlong[xlong],fontsize=9);
#plot.GetAxes(3).text(0.1,0.1,dataDetail,fontsize=9);

plot.Plot();

#plot.GetRatioAxes(3).remove();

#plot.Save("figs/Fig1_jtdist.pdf");
plot.Show();

