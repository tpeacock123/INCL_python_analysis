from incl_analysis_functions import *
import numpy as np
import sys 

ROOT.gSystem.Load("libNEUTROOTClass.so")
ROOT.gSystem.Load("libNEUTOutput.so")
ROOT.gSystem.Load("libNEUTReWeight.so")
ROOT.TH1.AddDirectory(False)

filename = sys.argv[1]
filename_chunk = filename.split(".")[0].split("out_")[1]
f = ROOT.TFile(filename)
t = f.Get("neuttree")

h_src = ROOT.TH2D("h_src", "SRC Events;Missing Energy (MeV);Neutron Momentum (MeV/c)", 
                  100, 0, 500, 100, 0, 300)
h_nonsrc = ROOT.TH2D("h_nonsrc", "Non-SRC Events;Missing Energy (MeV);Neutron Momentum (MeV/c)", 
                     100, 0, 500, 100, 0, 300)

for event in t:
    nvect = event.vectorbranch
    nvect_class = nvect_reader(nvect)
    
    miss_E = nvect_class.missing_E_calc()
    miss_p = nvect_class.neutron_mom()


    if nvect_class.eventType == EventType.SRC:
        h_src.Fill(miss_p, miss_E)
    else:
        h_nonsrc.Fill(miss_p, miss_E)


c1 = ROOT.TCanvas("c1", "SRC Comparison", 1200, 600)
c1.Divide(2, 1)

c1.cd(1)
h_src.Draw("COLZ") 

c1.cd(2)
h_nonsrc.Draw("COLZ")

c1.SaveAs("src_comparison_2d_{}.pdf".format(filename_chunk))