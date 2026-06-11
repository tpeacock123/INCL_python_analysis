from incl_analysis_functions import *
from default_processing import *
from kdar_processing import *
from paper_plots import *
import numpy as np
from array import array
import sys 
import argparse

ROOT.gSystem.Load("libNEUTROOTClass.so")
ROOT.gSystem.Load("libNEUTOutput.so")
ROOT.gSystem.Load("libNEUTReWeight.so")
ROOT.TH1.AddDirectory(False)

import ROOT
import numpy as np


def reader(filename):
    
    filename_chunk = filename.split(".")[0].split("out_")[1]

    f = ROOT.TFile(filename)
    t = f.Get("neuttree")   
    channel = []

    for event in t:
        
        nvect = event.vectorbranch
        nvect_class = nvect_reader(nvect)
        
        string_val = nvect_class.Print() 
        channel.append(string_val)

    file = ROOT.TFile("channel_{}.root".format(filename_chunk),"RECREATE")
    hist2 = ROOT.TH1D("h_channel", "Channel; Category; Counts", 3, 0, 3)
    hist2.SetCanExtend(ROOT.TH1.kAllAxes) 

    for value in channel:
        hist2.Fill(value, 1.0)

    # Optional: Shrink the histogram to remove empty bins if you enabled CanExtend

    hist2.LabelsDeflate()

    hist2.Write()
    file.Close()
 
def excitation_energy(filename):
    filename_chunk = filename.split(".")[0].split("out_")[1]

    f = ROOT.TFile(filename)
    t = f.Get("neuttree")

    noCascadeFSI = []
    QE_deex = []
    multipleNucleon_noCluster = []
    multipleNucleon_noCluster_src = []
    nuclearCluster = []
    nuclearCluster_src = []
    oneProton = []
    protonPion = []
    no_protons = []
    other = []
    src = []


    for event in t:
        nvect = event.vectorbranch
        nvect_class = nvect_reader(nvect)
        
        prefsi_proton_mom = nvect_class.PreFSIProtMom
        eventType = nvect_class.eventType
        intChannel = nvect_class.intChannel
        miss_E = nvect_class.E_miss
        miss_p = nvect_class.P_miss

        

        if eventType == EventType.MF:
            excitation_E = nvect_class.excitation_E_CCQE()
            if (intChannel == intChannel_CCQE.muOnly) or (intChannel == intChannel_CCQE.neutronPion):
                no_protons.append(prefsi_proton_mom)
                continue
                
            if intChannel == intChannel_CCQE.noCascadeFSI:
                noCascadeFSI.append(excitation_E)


            elif intChannel == intChannel_CCQE.qeDeEX:
                QE_deex.append(excitation_E)

            
            elif intChannel == intChannel_CCQE.oneProton:
                oneProton.append(excitation_E)

            
            elif intChannel == intChannel_CCQE.multipleNucleon:
                multipleNucleon_noCluster.append(excitation_E)

            elif intChannel == intChannel_CCQE.nuclearCluster:
                nuclearCluster.append(excitation_E)

            elif intChannel == intChannel_CCQE.protonPion: 
                protonPion.append(excitation_E)


            elif intChannel == intChannel_CCQE.other:
                other_counter +=1 
                nvect_class.Print()
                other.append(excitation_E)

        elif eventType == EventType.SRC:
            src.append(nvect_class.excitation_E_SRC())


    file = ROOT.TFile("excitation_energy_{}.root".format(filename_chunk),"RECREATE")

    c1 = ROOT.TCanvas("c1", "Excitation_energy_distribution", 800, 600)
    c1.cd() 
    pad1 = ROOT.TPad("pad1", "pad1", 0, 0.3, 1, 1.0)
    pad1.SetBottomMargin(0.02) 
    pad1.SetLeftMargin(0.12)
    pad1.SetTicks(1, 1)        
    pad1.Draw()

    n_bins = 300
    x_min = -60
    x_max = 60

    h_nucCluster = create_histo("h_nuc", "Nuclear Clusters",       ROOT.kViolet-4, 3001, nuclearCluster,n_bins,x_min,x_max)
    h_noCascade  = create_histo("h_noC", "No Cascade FSI",         ROOT.kOrange-4, 3001, noCascadeFSI,n_bins,x_min,x_max)
    h_QE         = create_histo("h_qe",  "QE Proton + De-ex",      ROOT.kOrange+7, 3001, QE_deex,n_bins,x_min,x_max)
    h_multiNuc   = create_histo("h_mul", "Multiple Nucleons",      ROOT.kRed+2,    3001, multipleNucleon_noCluster,n_bins,x_min,x_max)
    h_oneProton  = create_histo("h_one", "One Proton",             ROOT.kGreen,  3001, oneProton,n_bins,x_min,x_max)
    h_protonPion = create_histo("h_pion","Proton + Pion",          ROOT.kYellow-9, 3001, protonPion,n_bins,x_min,x_max)
    h_other      = create_histo("h_oth", "No protons",             ROOT.kGray,     3001, no_protons,n_bins,x_min,x_max)
    h_src        = create_histo("h_src", "SRC",                    ROOT.kBlue,     3001, src,n_bins,x_min,x_max)


    hs1 = ROOT.THStack("hs_ExcitationE", "excitationE_distribution")
    hs1.Add(h_protonPion)
    hs1.Add(h_other)
    hs1.Add(h_oneProton)
    hs1.Add(h_multiNuc)
    hs1.Add(h_nucCluster)
    hs1.Add(h_noCascade)
    hs1.Add(h_QE)
    hs1.Add(h_src)

    hs1.Draw("hist") 
    hs1.SetTitle("; excitation E [MeV];Number of Events") 

    hs1.GetYaxis().SetLabelSize(0.03)
    hs1.GetYaxis().SetTitleSize(0.04)
    hs1.GetYaxis().SetTitleOffset(1.1)

    hs1.GetXaxis().SetTitle("Excitation Energy [MeV]")

    max_h = hs1.GetMaximum()
    hs1.SetMaximum(max_h * 1.2) 

    legend1 = ROOT.TLegend(0.55, 0.55, 0.88, 0.88)
    legend1.SetBorderSize(0)
    legend1.SetFillStyle(0) # Transparent
    legend1.AddEntry(h_noCascade, "no cascade FSI", "f")
    legend1.AddEntry(h_QE, "QE proton + de-excitation", "f")
    legend1.AddEntry(h_multiNuc, "multiple nucleons", "f")
    legend1.AddEntry(h_nucCluster, "nuclear clusters", "f")
    legend1.AddEntry(h_oneProton, "one proton", "f")
    legend1.AddEntry(h_protonPion, "proton + pion", "f")
    legend1.AddEntry(h_other, "no protons", "f")
    legend1.AddEntry(h_src, "src", "f")

    legend1.Draw()

    file.cd()
    c1.Write()
    file.Close()

def missing_E_2p2h(filename):
    filename_chunk = filename.split(".")[0].split("out_")[1]

    f = ROOT.TFile(filename)
    t = f.Get("neuttree")
    missing_E2p2h = []
    excitation_E = []
    for event in t:
        nvect = event.vectorbranch
        nvect_class = nvect_reader(nvect)
        missing_E2p2h.append(nvect_class.missing_energy_2p2h())
        excitation_E.append(nvect_class.excitation_E_2p2h())
        print(nvect_class.excitation_E_2p2h())

    file = ROOT.TFile("missing_E_2p2h_{}.root".format(filename_chunk),"RECREATE")
    hist = ROOT.TH1D("h_missing_E2p2h", "Missing Energy 2p2h;Missing E [MeV];Counts", 100, 0, 100)

    for value in missing_E2p2h:
        hist.Fill(value)

    hist2 = ROOT.TH1D("h_excitation_E2p2h", "Excitation E 2p2h; Excitation E [MeV];Counts", 100, -10, 100)
    for value in excitation_E:
        hist2.Fill(value)

    hist.Write()
    hist2.Write()
    file.Close()


def main():
    additional_file = None
    input_path = []
    parser = argparse.ArgumentParser()
    parser.add_argument("-t", "--type", help="type of processing")
    parser.add_argument('-i', '--input', nargs='+', help="Input file(s)")
    parser.add_argument("-f", "--file",  nargs='+', help="file/s")
    args = parser.parse_args()

    input_files = args.input if args.input else []
    additional_files = args.file if args.input else []

    for file in input_files:
        print(f"Logic: Input file queued: '{file}'")
    if additional_files:
        for additional_file in additional_files:
            print(f"Logic: Additional file detected: '{additional_file}'")


    if args.type:
        function_to_run = args.type
        if function_to_run == "reader":
            for file in input_files:
                reader(file)

        elif function_to_run == "kdar":
            for file in input_files:
                missing_energy_kdar(file, additional_files)

        elif function_to_run == "excitation":
            for file in input_files:
                excitation_energy(file)
        
        elif function_to_run == "spectral_function":
            for file in input_files:
                for add_file in additional_files:
                    spectral_function2d(file, add_file)
        
        elif function_to_run == "src":
            for file in input_files:
                for add_file in additional_files:
                    SRC_plot(file, add_file)

        elif function_to_run == "processing":
            for file in input_files:
                INCL_processing(file)

        elif function_to_run == "2p2h":
            for file in input_files:
                missing_E_2p2h(file)
        elif function_to_run == "transparency":
            for file in input_files:
                transparency(file, additional_files)


        elif function_to_run == "CCQE_comps":
            for file in input_files:
                ccqe_combined_plots(file,additional_files)

        else:
            print(f"Error: '{function_to_run}' is not a recognized option.")

    else:
        # Default behavior if no type is specified
        for file in input_files:
            INCL_processing(file)


main()