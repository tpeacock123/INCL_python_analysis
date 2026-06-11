from incl_analysis_functions import *
import ROOT
import numpy as np
import sys 
import argparse



ROOT.gSystem.Load("libNEUTROOTClass.so")
ROOT.gSystem.Load("libNEUTOutput.so")
ROOT.gSystem.Load("libNEUTReWeight.so")
ROOT.TH1.AddDirectory(False)




def process_events(t):

    output_dic = {
        "noCascadeFSI_pt": [], "QE_deex_pt": [], "multipleNucleon_noCluster_pt": [],
        "nuclearCluster_pt": [], "oneProton_pt": [], "protonPion_pt": [], "other_pt": [],

        "noCascadeFSI_at": [], "QE_deex_at": [], "multipleNucleon_noCluster_at": [],
        "nuclearCluster_at": [], "oneProton_at": [], "protonPion_at": [], "other_at": [],

        "noCascadeFSI": [], "QE_deex": [], "multipleNucleon_noCluster": [],
        "multipleNucleon_noCluster_src": [], "nuclearCluster": [], "nuclearCluster_src": [],
        "oneProton": [], "protonPion": [], "other": [],

        "noCascadeFSI_pre": [], "QE_deex_pre": [], "multipleNucleon_noCluster_pre": [],
        "nuclearCluster_pre": [], "oneProton_pre": [], "protonPion_pre": [],
        "mu_only": [], "no_protons": [], "other_pre": [],

        "src_noCascadeFSI_pt": [], "src_qeDeEX_pt": [], "src_elasticProton_pt": [],
        "src_multipleNucleons_pt": [], "src_nuclearClusters_pt": [],
        "src_protonPion_pt": [], "src_muOnly_pt": [], "src_other_pt": [],
        "src_neutronPion_pt": [],

        "src_noCascadeFSI_at": [], "src_qeDeEX_at": [], "src_elasticProton_at": [],
        "src_multipleNucleons_at": [], "src_nuclearClusters_at": [],
        "src_protonPion_at": [], "src_muOnly_at": [], "src_other_at": [],
        "src_neutronPion_at": [],

        "src_noCascadeFSI": [], "src_qeDeEX": [], "src_elasticProton": [],
        "src_multipleNucleons": [], "src_nuclearClusters": [],
        "src_protonPion": [], "src_muOnly": [], "src_other": [],
        "src_neutronPion": [], 

        "src_elasticProton_T": [],"src_multipleNucleons_T": [], 
        "src_nuclearClusters_T": [],  "src_protonPion_T": [],

        "2p2h_noCascadeFSI": [], "2p2h_qeDeEX": [], "2p2h_elasticProton": [],
        "2p2h_multipleNucleons": [], "2p2h_nuclearClusters": [],
        "2p2h_protonPion": [], "2p2h_muOnly": [], "2p2h_other": [],
        "2p2h_neutronPion": [],

        "ex_mf": [], "ex_noCascadeFSI": [],"ex_QE_deex": [], "ex_multipleNucleon_noCluster": [],"ex_nuclearCluster": [],
        "ex_oneProton": [], "ex_protonPion": [],"ex_no_protons": [], "ex_src":[], "ex_2p2h": [],
        "ex_other": []
    }

    h_src = ROOT.TH2D("h_src", "SRC Events;Missing Energy (MeV);Neutron Momentum (MeV/c)", 
                100, 0, 500, 100, 0, 300)

    h_nonsrc = ROOT.TH2D("h_nonsrc", "Non-SRC Events;Missing Energy (MeV);Neutron Momentum (MeV/c)", 
                        100, 0, 500, 100, 0, 300)
    

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
            output_dic["ex_mf"].append(excitation_E)
            if (intChannel == intChannel_CCQE.muOnly) or (intChannel == intChannel_CCQE.neutronPion):
                output_dic["no_protons"].append(prefsi_proton_mom)
                output_dic["ex_no_protons"].append(excitation_E) 
                continue
                
            HMP_proton = nvect_class.HMPMom
            DPT = nvect_class.DpT
            dat = nvect_class.DaT

            if intChannel == intChannel_CCQE.noCascadeFSI:
                output_dic["noCascadeFSI"].append(HMP_proton)
                output_dic["noCascadeFSI_at"].append(dat)
                output_dic["noCascadeFSI_pt"].append(DPT)
                output_dic["ex_noCascadeFSI"].append(excitation_E) # Added
                if(prefsi_proton_mom != 0):
                    output_dic["noCascadeFSI_pre"].append(prefsi_proton_mom)
                else:
                    output_dic["noCascadeFSI_pre"].append(HMP_proton)

            elif intChannel == intChannel_CCQE.qeDeEX:
                output_dic["QE_deex"].append(HMP_proton)
                output_dic["QE_deex_at"].append(dat)
                output_dic["QE_deex_pt"].append(DPT)
                output_dic["QE_deex_pre"].append(prefsi_proton_mom)
                output_dic["ex_QE_deex"].append(excitation_E) # Added
            
            elif intChannel == intChannel_CCQE.oneProton:

                output_dic["oneProton"].append(HMP_proton)
                output_dic["oneProton_at"].append(dat)
                output_dic["oneProton_pt"].append(DPT)
                output_dic["oneProton_pre"].append(prefsi_proton_mom)
                output_dic["ex_oneProton"].append(excitation_E) # Added
            
            elif intChannel == intChannel_CCQE.multipleNucleon:
                output_dic["multipleNucleon_noCluster"].append(HMP_proton)
                output_dic["multipleNucleon_noCluster_at"].append(dat)
                output_dic["multipleNucleon_noCluster_pt"].append(DPT)
                output_dic["multipleNucleon_noCluster_pre"].append(prefsi_proton_mom)
                output_dic["ex_multipleNucleon_noCluster"].append(excitation_E) # Added
            
            elif intChannel == intChannel_CCQE.nuclearCluster:
                output_dic["nuclearCluster"].append(HMP_proton)
                output_dic["nuclearCluster_at"].append(dat)
                output_dic["nuclearCluster_pt"].append(DPT)
                output_dic["nuclearCluster_pre"].append(prefsi_proton_mom)
                output_dic["ex_nuclearCluster"].append(excitation_E) # Added

            elif intChannel == intChannel_CCQE.protonPion: 
                output_dic["protonPion"].append(HMP_proton)
                output_dic["protonPion_at"].append(dat)
                output_dic["protonPion_pt"].append(DPT)
                output_dic["protonPion_pre"].append(prefsi_proton_mom)  
                output_dic["ex_protonPion"].append(excitation_E) # Added

            elif intChannel == intChannel_CCQE.other:
                other_counter += 1 
                nvect_class.Print()
                output_dic["other"].append(HMP_proton)
                output_dic["other_at"].append(dat)
                output_dic["other_pt"].append(DPT)
                output_dic["other_pre"].append(prefsi_proton_mom)
                output_dic["ex_other"].append(excitation_E) # Added

        elif eventType == EventType.SRC:
            excitation_E = nvect_class.excitation_E_SRC()
            output_dic["ex_src"].append(excitation_E)
            if (intChannel == intChannel_CCQE.muOnly) or (intChannel == intChannel_CCQE.neutronPion):
                continue

            if nvect_class.HMPMom is None:    
                continue

            HMP_proton = nvect_class.HMPMom 
            DPT = 0 
            dat = 0 

            if intChannel == intChannel_CC0pi.noCascadeFSI:
                output_dic["src_noCascadeFSI"].append(HMP_proton)
                output_dic["src_noCascadeFSI_at"].append(dat)
                output_dic["src_noCascadeFSI_pt"].append(DPT)

            elif intChannel == intChannel_CC0pi.qeDeEX:
                output_dic["src_qeDeEX"].append(HMP_proton)
                output_dic["src_qeDeEX_at"].append(dat)
                output_dic["src_qeDeEX_pt"].append(DPT)
            
            elif intChannel == intChannel_CC0pi.elasticProton:
                if nvect_class.istransparent == True:
                    print("here")
                    output_dic["src_elasticProton_T"].append(HMP_proton)
                else:
                    output_dic["src_elasticProton"].append(HMP_proton)
                    output_dic["src_elasticProton_at"].append(dat)
                    output_dic["src_elasticProton_pt"].append(DPT)


       
            elif intChannel == intChannel_CC0pi.multipleNucleons:
                if nvect_class.istransparent == True:
                    print("here")
                    output_dic["src_multipleNucleons_T"].append(HMP_proton)
                else:
                    output_dic["src_multipleNucleons"].append(HMP_proton)
                    output_dic["src_multipleNucleons_at"].append(dat)
                    output_dic["src_multipleNucleons_pt"].append(DPT)
                    
            elif intChannel == intChannel_CC0pi.nuclearClusters: 
                if nvect_class.istransparent == True:
                    print("here")
                    output_dic["src_nuclearClusters_T"].append(HMP_proton)
                else:
                    output_dic["src_nuclearClusters"].append(HMP_proton)
                    output_dic["src_nuclearClusters_at"].append(dat)
                    output_dic["src_nuclearClusters_pt"].append(DPT)
                
            elif intChannel == intChannel_CC0pi.protonPion: 
                output_dic["src_protonPion"].append(HMP_proton)
                output_dic["src_protonPion_at"].append(dat)
                output_dic["src_protonPion_pt"].append(DPT)  

            elif intChannel == intChannel_CCQE.other:
                other_counter += 1 
                nvect_class.Print()
                output_dic["src_other"].append(HMP_proton)
                output_dic["src_other_at"].append(dat)
                output_dic["src_other_pt"].append(DPT)

        elif eventType == EventType.twop2h:
            excitation_E = nvect_class.excitation_E_2p2h()
            output_dic["ex_2p2h"].append(excitation_E)
            if (intChannel == intChannel_CCQE.muOnly) or (intChannel == intChannel_CCQE.neutronPion):
                continue

            if nvect_class.HMPMom is None:    
                #nvect_class.Print()
                continue

            HMP_proton = nvect_class.HMPMom 
            DPT = 0 
            dat = 0 

            if intChannel == intChannel_CC0pi.noCascadeFSI:
                output_dic["2p2h_noCascadeFSI"].append(HMP_proton)

            elif intChannel == intChannel_CC0pi.qeDeEX:
                output_dic["2p2h_qeDeEX"].append(HMP_proton)
            
            elif intChannel == intChannel_CC0pi.elasticProton:
                output_dic["2p2h_elasticProton"].append(HMP_proton)
       
            elif intChannel == intChannel_CC0pi.multipleNucleons:
                output_dic["2p2h_multipleNucleons"].append(HMP_proton)
                
            elif intChannel == intChannel_CC0pi.nuclearClusters: 
                output_dic["2p2h_nuclearClusters"].append(HMP_proton)
                
            elif intChannel == intChannel_CC0pi.protonPion: 
                output_dic["2p2h_protonPion"].append(HMP_proton)

            elif intChannel == intChannel_CCQE.other:
                other_counter += 1 
                #nvect_class.Print()
                output_dic["2p2h_other"].append(HMP_proton)
        else:
            print("how are we here")
            #nvect_class.Print()

        if eventType == EventType.SRC:
            h_src.Fill(miss_p, miss_E)
        elif eventType == EventType.MF:
            h_nonsrc.Fill(miss_p, miss_E)

    return output_dic, h_src,h_nonsrc

def INCL_processing(filename):
    
    filename_chunk = filename.split(".")[0].split("out_")[1]

    f = ROOT.TFile(filename)
    t = f.Get("neuttree")
    i = 0
    
    other_counter = 0
    evt_counter = 0

    output_dic,h_src,h_nonsrc = process_events(t)

    print("other counter", other_counter)  

    file = ROOT.TFile("incl_analysis_test_{}.root".format(filename_chunk),"RECREATE")

    c0 = ROOT.TCanvas("c0", "SRC Comparison", 1200, 600)
    c0.Divide(2, 1)
    c0.cd(1)
    h_src.Draw("COLZ") 
    c0.cd(2)
    h_nonsrc.Draw("COLZ")


    c1 = ROOT.TCanvas("c1", "Proton Momentum Distributions", 800, 600)
    c1.cd() 
    pad1 = ROOT.TPad("pad1", "pad1", 0, 0.3, 1, 1.0)
    pad1.SetBottomMargin(0.02) 
    pad1.SetLeftMargin(0.12)
    pad1.SetTicks(1, 1)        
    pad1.Draw()

    pad2 = ROOT.TPad("pad2", "pad2", 0, 0.0, 1, 0.3)
    pad2.SetTopMargin(0.01)  
    pad2.SetBottomMargin(0.25) 
    pad2.SetLeftMargin(0.12)
    pad2.SetTicks(1, 1)
    pad2.Draw()

    n_bins = 40
    x_min = 0
    x_max = 1500   

    h_nucCluster = create_histo("h_nuc", "Nuclear Clusters",       ROOT.kViolet-4, 3001, output_dic["nuclearCluster"],n_bins,x_min,x_max)
    #h_nucCluster_src = create_histo("h_nuc_src", "Nuclear Clusters, SRC",       ROOT.kViolet-1, 3001, output_dic["nuclearCluster_src"],n_bins,x_min,x_max)
    h_noCascade  = create_histo("h_noC", "No Cascade FSI",         ROOT.kOrange-4, 3001, output_dic["noCascadeFSI"],n_bins,x_min,x_max)
    h_QE         = create_histo("h_qe",  "QE Proton + De-ex",      ROOT.kOrange+7, 3001, output_dic["QE_deex"],n_bins,x_min,x_max)
    h_multiNuc   = create_histo("h_mul", "Multiple Nucleons",      ROOT.kRed+2,    3001, output_dic["multipleNucleon_noCluster"],n_bins,x_min,x_max)
    #h_multiNuc_src = create_histo("h_mul_src", "Multiple Nucleons, SRC",      ROOT.kRed,    3001, output_dic["multipleNucleon_noCluster_src"],n_bins,x_min,x_max)
    h_oneProton  = create_histo("h_one", "One Proton",             ROOT.kGreen,  3001, output_dic["oneProton"],n_bins,x_min,x_max)
    h_protonPion = create_histo("h_pion","Proton + Pion",          ROOT.kYellow-9, 3001, output_dic["protonPion"],n_bins,x_min,x_max)
    h_other      = create_histo("h_oth", "Other",                  ROOT.kGray,     3001, output_dic["other"],n_bins,x_min,x_max)

    pad1.cd() 

    hs = ROOT.THStack("hs", "Proton Momentum by Channel;Momentum [MeV/c];Events")
    hs.Add(h_protonPion)
    hs.Add(h_other)
    hs.Add(h_oneProton)
    hs.Add(h_multiNuc)
    #hs.Add(h_multiNuc_src)
    hs.Add(h_nucCluster)
    #hs.Add(h_nucCluster_src)
    hs.Add(h_noCascade)
    hs.Add(h_QE)

    hs.Draw("hist") 
    hs.GetXaxis().SetLabelSize(0)
    hs.SetTitle("; ;Number of Events") 


    hs.GetYaxis().SetLabelSize(0.04)
    hs.GetYaxis().SetTitleSize(0.05)
    hs.GetYaxis().SetTitleOffset(1.1)

    max_h = hs.GetMaximum()
    hs.SetMaximum(max_h * 1.2) 

    total_integral = sum(hist.Integral() for hist in hs.GetHists())

    print(f"Total visible events in the stack: {total_integral}")


    legend = ROOT.TLegend(0.55, 0.55, 0.88, 0.88)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0) # Transparent
    legend.AddEntry(h_noCascade, "no cascade FSI", "f")
    legend.AddEntry(h_QE, "QE proton + de-excitation", "f")
    legend.AddEntry(h_multiNuc, "multiple nucleons", "f")
    #legend.AddEntry(h_multiNuc_src, "multiple nucleons, SRC", "f")
    legend.AddEntry(h_nucCluster, "nuclear clusters", "f")
    #legend.AddEntry(h_nucCluster_src, "nuclear clusters", "f")
    legend.AddEntry(h_oneProton, "one proton", "f")
    legend.AddEntry(h_protonPion, "proton + pion", "f")
    legend.Draw()


    pad2.cd() 

    h_total = h_nucCluster.Clone("h_total")
    #h_total.Add(h_nucCluster_src)
    h_total.Add(h_noCascade)
    h_total.Add(h_QE)
    h_total.Add(h_multiNuc)
   # h_total.Add(h_multiNuc_src)
    h_total.Add(h_oneProton)
    h_total.Add(h_protonPion)
    h_total.Add(h_other)


    def create_ratio(h_in, h_tot):
        h_ratio = h_in.Clone(h_in.GetName() + "_ratio")
        h_ratio.Divide(h_tot) 


        h_ratio.SetFillStyle(0)
        h_ratio.SetLineWidth(2)
        h_ratio.SetLineColor(h_in.GetLineColor()) 
        return h_ratio

    r_nuc   = create_ratio(h_nucCluster, h_total)
    
    r_noCas = create_ratio(h_noCascade, h_total)
    r_qe    = create_ratio(h_QE, h_total)
    r_mul   = create_ratio(h_multiNuc, h_total) 
    #r_mul_src   = create_ratio(h_multiNuc_src, h_total)
    r_one   = create_ratio(h_oneProton, h_total)
    r_pion  = create_ratio(h_protonPion, h_total)


    r_nuc.Draw("HIST") 
    r_nuc.SetTitle(";P_{aft} (GeV/c);Fractional Contributions")


    y_axis = r_nuc.GetYaxis()
    y_axis.SetRangeUser(0, 1.1)
    y_axis.SetNdivisions(505)
    y_axis.SetLabelSize(0.08)
    y_axis.SetTitleSize(0.09)
    y_axis.SetTitleOffset(0.5)

    x_axis = r_nuc.GetXaxis()
    x_axis.SetLabelSize(0.08)
    x_axis.SetTitleSize(0.10)
    x_axis.SetTitleOffset(1.0)

    r_noCas.Draw("HIST SAME")
    r_qe.Draw("HIST SAME")
    r_mul.Draw("HIST SAME")
    #r_mul_src.Draw("hist SAME")
    r_one.Draw("HIST SAME")
    r_pion.Draw("HIST SAME")

    c1.Update()
    #c1.SaveAs("proton_fractional_plot_{}_newdefs.png".format(filename_chunk))

    c2 = ROOT.TCanvas("c2", "deltaPT", 800, 600)
    c2.cd()

    n_bins = 50    
    x_min = 0
    x_max = 1000   

    h_nucCluster = create_histo("h_nuc_pt", "Nuclear Clusters",       ROOT.kViolet-4, 3001, output_dic["nuclearCluster_pt"],n_bins,x_min,x_max)
    h_noCascade  = create_histo("h_noC_pt", "No Cascade FSI",         ROOT.kOrange-4, 3001, output_dic["noCascadeFSI_pt"],n_bins,x_min,x_max)
    h_QE         = create_histo("h_qe_pt",  "QE Proton + De-ex",      ROOT.kOrange+7, 3001, output_dic["QE_deex_pt"],n_bins,x_min,x_max)
    h_multiNuc   = create_histo("h_mul_pt", "Multiple Nucleons",      ROOT.kRed+2,    3001, output_dic["multipleNucleon_noCluster_pt"],n_bins,x_min,x_max)
    h_oneProton  = create_histo("h_one_pt", "One Proton",             ROOT.kGreen,  3001, output_dic["oneProton_pt"],n_bins,x_min,x_max)
    h_protonPion = create_histo("h_pion_pt","Proton + Pion",          ROOT.kYellow-9, 3001, output_dic["protonPion_pt"],n_bins,x_min,x_max)
    h_other      = create_histo("h_oth_pt", "Other",                  ROOT.kGray,     3001, output_dic["other_pt"],n_bins,x_min,x_max)

    hs1 = ROOT.THStack("hs_pt", "DPT by Channel;Momentum [MeV/c];Events")
    hs1.Add(h_protonPion)
    hs1.Add(h_other)
    hs1.Add(h_oneProton)
    hs1.Add(h_multiNuc)
    hs1.Add(h_nucCluster)
    hs1.Add(h_noCascade)
    hs1.Add(h_QE)


    hs1.Draw("hist") 
    hs1.SetTitle("; DPT;Number of Events") 

    hs1.GetYaxis().SetLabelSize(0.03)
    hs1.GetYaxis().SetTitleSize(0.04)
    hs1.GetYaxis().SetTitleOffset(1.1)

    hs1.GetXaxis().SetTitle("DPT")

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

    legend1.Draw()




    c3 = ROOT.TCanvas("c3", "deltaPT", 800, 600)
    c3.cd()

    n_bins = 50    
    x_min = 0
    x_max = 180  

    h_nucCluster = create_histo("h_nuc_pt", "Nuclear Clusters",       ROOT.kViolet-4, 3001, output_dic["nuclearCluster_at"],n_bins,x_min,x_max)
    h_noCascade  = create_histo("h_noC_pt", "No Cascade FSI",         ROOT.kOrange-4, 3001, output_dic["noCascadeFSI_at"],n_bins,x_min,x_max)
    h_QE         = create_histo("h_qe_pt",  "QE Proton + De-ex",      ROOT.kOrange+7, 3001, output_dic["QE_deex_at"],n_bins,x_min,x_max)
    h_multiNuc   = create_histo("h_mul_pt", "Multiple Nucleons",      ROOT.kRed+2,    3001, output_dic["multipleNucleon_noCluster_at"],n_bins,x_min,x_max)
    h_oneProton  = create_histo("h_one_pt", "One Proton",             ROOT.kGreen,  3001, output_dic["oneProton_at"],n_bins,x_min,x_max)
    h_protonPion = create_histo("h_pion_pt","Proton + Pion",          ROOT.kYellow-9, 3001, output_dic["protonPion_at"],n_bins,x_min,x_max)
    h_other      = create_histo("h_oth_pt", "Other",                  ROOT.kGray,     3001, output_dic["other_at"],n_bins,x_min,x_max)

    hs2 = ROOT.THStack("hs_at", "DaT by Channel;Momentum [MeV/c];Events")
    hs2.Add(h_noCascade)
    hs2.Add(h_QE)
    hs2.Add(h_oneProton)
    hs2.Add(h_protonPion)
    hs2.Add(h_multiNuc)
    hs2.Add(h_other)
    hs2.Add(h_nucCluster)


    hs2.Draw("hist") 
    hs2.SetTitle("; DaT;Number of Events") 

    hs2.GetYaxis().SetLabelSize(0.03)
    hs2.GetYaxis().SetTitleSize(0.04)
    hs2.GetYaxis().SetTitleOffset(1.1)

    hs2.GetXaxis().SetTitle("DaT")

    max_h = hs2.GetMaximum()
    hs2.SetMaximum(max_h * 1.2) 

    legend2 = ROOT.TLegend(0.55, 0.55, 0.88, 0.88)
    legend2.SetBorderSize(0)
    legend2.SetFillStyle(0) # Transparent
    legend2.AddEntry(h_noCascade, "no cascade FSI", "f")
    legend2.AddEntry(h_QE, "QE proton + de-excitation", "f")
    legend2.AddEntry(h_multiNuc, "multiple nucleons", "f")
    legend2.AddEntry(h_nucCluster, "nuclear clusters", "f")
    legend2.AddEntry(h_oneProton, "one proton", "f")
    legend2.AddEntry(h_protonPion, "proton + pion", "f")

    legend2.Draw()

    c4 = ROOT.TCanvas("c4", "Proton Momentum Distributions", 800, 600)
    c4.cd() 
    pad3 = ROOT.TPad("pad3", "pad3", 0, 0.3, 1, 1.0)
    pad3.SetBottomMargin(0.02) 
    pad3.SetLeftMargin(0.12)
    pad3.SetTicks(1, 1)        
    pad3.Draw()

    pad4 = ROOT.TPad("pad4", "pad4", 0, 0.0, 1, 0.3)
    pad4.SetTopMargin(0.01)  
    pad4.SetBottomMargin(0.25) 
    pad4.SetLeftMargin(0.12)
    pad4.SetTicks(1, 1)
    pad4.Draw()

    n_bins = 90    
    x_min = 0
    x_max = 1500   

    h_nucCluster_pre = create_histo("h_nuc_prefsi", "Nuclear Clusters",       ROOT.kViolet-4, 3001, output_dic["nuclearCluster_pre"],n_bins,x_min,x_max)
    h_noCascade_pre  = create_histo("h_noC_prefsi", "No Cascade FSI",         ROOT.kOrange-4, 3001, output_dic["noCascadeFSI_pre"],n_bins,x_min,x_max)
    h_QE_pre         = create_histo("h_qe_prefsi",  "QE Proton + De-ex",      ROOT.kOrange+7, 3001, output_dic["QE_deex_pre"],n_bins,x_min,x_max)
    h_multiNuc_pre   = create_histo("h_mul_prefsi", "Multiple Nucleons",      ROOT.kRed+2,    3001, output_dic["multipleNucleon_noCluster_pre"],n_bins,x_min,x_max)
    h_oneProton_pre  = create_histo("h_one_prefsi", "One Proton",             ROOT.kGreen,    3001, output_dic["oneProton_pre"],n_bins,x_min,x_max)
    h_protonPion_pre = create_histo("h_pion_prefsi","Proton + Pion",          ROOT.kYellow-9, 3001, output_dic["protonPion_pre"],n_bins,x_min,x_max)
    h_muonly_pre     = create_histo("h_muonly_prefsi","muonly",               ROOT.kViolet-9,  3001, output_dic["mu_only"],n_bins,x_min,x_max)
    h_noproton_pre   = create_histo("h_noproton_prefsi","no protons",         ROOT.kGreen-9,  3001, output_dic["no_protons"],n_bins,x_min,x_max)
    h_other_pre      = create_histo("h_oth_prefsi", "Other",                  ROOT.kGray,     3001, output_dic["other_pre"],n_bins,x_min,x_max)
    pad3.cd() 

    hs3 = ROOT.THStack("hs3", "Proton Momentum by Channel;Momentum [MeV/c];Events")
    hs3.Add(h_protonPion_pre)
    hs3.Add(h_other_pre)
    hs3.Add(h_oneProton_pre)
    hs3.Add(h_multiNuc_pre)
    hs3.Add(h_nucCluster_pre)
    hs3.Add(h_noCascade_pre)
    hs3.Add(h_QE_pre)
    hs3.Add(h_muonly_pre)
    hs3.Add(h_noproton_pre)

    hs3.Draw("hist") 
    hs3.GetXaxis().SetLabelSize(0)
    hs3.SetTitle("; ;Number of Events") 


    hs3.GetYaxis().SetLabelSize(0.04)
    hs3.GetYaxis().SetTitleSize(0.05)
    hs3.GetYaxis().SetTitleOffset(1.1)

    max_h = hs3.GetMaximum()
    hs3.SetMaximum(max_h * 1.2) 

    total_integral = sum(hist.Integral() for hist in hs2.GetHists())

    print(f"Total visible events in the stack: {total_integral}")


    legend4 = ROOT.TLegend(0.55, 0.55, 0.88, 0.88)
    legend4.SetBorderSize(0)
    legend4.SetFillStyle(0) # Transparent
    legend4.AddEntry(h_noCascade_pre, "no cascade FSI", "f")
    legend4.AddEntry(h_QE_pre, "QE proton + de-excitation", "f")
    legend4.AddEntry(h_multiNuc_pre, "multiple nucleons", "f")
    legend4.AddEntry(h_nucCluster_pre, "nuclear clusters", "f")
    legend4.AddEntry(h_oneProton_pre, "one proton", "f")
    legend4.AddEntry(h_protonPion_pre, "proton + pion", "f")
    legend4.AddEntry(h_muonly_pre, "mu only", "f")
    legend4.AddEntry(h_noproton_pre, "no protons", "f")
    legend4.Draw()


    pad4.cd() 

    h_total2 = h_nucCluster_pre.Clone("h_total")
    h_total2.Add(h_noCascade_pre)
    h_total2.Add(h_QE_pre)
    h_total2.Add(h_multiNuc_pre)
    h_total2.Add(h_oneProton_pre)
    h_total2.Add(h_protonPion_pre)
    h_total2.Add(h_other_pre)

    r_nuc1   = create_ratio(h_nucCluster_pre, h_total2)
    r_noCas1 = create_ratio(h_noCascade_pre, h_total2)
    r_qe1    = create_ratio(h_QE_pre, h_total2)
    r_mul1   = create_ratio(h_multiNuc_pre, h_total2)
    r_one1   = create_ratio(h_oneProton_pre, h_total2)
    r_pion1  = create_ratio(h_protonPion_pre, h_total2)
    r_muonly1= create_ratio(h_muonly_pre, h_total2)
    r_noproton =create_ratio(h_noproton_pre, h_total2)


    r_nuc1.Draw("HIST") 
    r_nuc1.SetTitle(";P_{aft} (GeV/c);Fractional Contributions")


    y_axis = r_nuc1.GetYaxis()
    y_axis.SetRangeUser(0, 1.1)
    y_axis.SetNdivisions(505)
    y_axis.SetLabelSize(0.08)
    y_axis.SetTitleSize(0.09)
    y_axis.SetTitleOffset(0.5)

    x_axis = r_nuc1.GetXaxis()
    x_axis.SetLabelSize(0.08)
    x_axis.SetTitleSize(0.10)
    x_axis.SetTitleOffset(1.0)


    r_noCas1.Draw("HIST SAME")
    r_qe1.Draw("HIST SAME")
    r_mul1.Draw("HIST SAME")
    r_one1.Draw("HIST SAME")
    r_pion1.Draw("HIST SAME")
    r_muonly1.Draw("HIST SAME")
    r_noproton.Draw("HIST SAME")

    c4.Update()


    c6 = ROOT.TCanvas("c6", "SRC Proton Momentum Distributions", 800, 600)
    c6.cd() 
    
    # Appended _c6 to pad names to avoid overwriting your c1 pads if run in the same script
    pad1_c6 = ROOT.TPad("pad1_c6", "pad1", 0, 0.3, 1, 1.0)
    pad1_c6.SetBottomMargin(0.02) 
    pad1_c6.SetLeftMargin(0.12)
    pad1_c6.SetTicks(1, 1)        
    pad1_c6.Draw()

    pad2_c6 = ROOT.TPad("pad2_c6", "pad2", 0, 0.0, 1, 0.3)
    pad2_c6.SetTopMargin(0.01)  
    pad2_c6.SetBottomMargin(0.25) 
    pad2_c6.SetLeftMargin(0.12)
    pad2_c6.SetTicks(1, 1)
    pad2_c6.Draw()

    n_bins = 20
    x_min = 0
    x_max = 1500   

    # Creating histograms using the updated src_ lists
    h_src_noCascade    = create_histo("h_src_noC",       "No Cascade FSI",      ROOT.kOrange-4, 3001, output_dic["src_noCascadeFSI"],    n_bins, x_min, x_max)
    h_src_qe           = create_histo("h_src_qe",        "QE Proton + De-ex",   ROOT.kOrange+7, 3001, output_dic["src_qeDeEX"],          n_bins, x_min, x_max)
    h_src_elasticProton= create_histo("h_src_elasticP",  "Elastic Proton",      ROOT.kGreen+2,  3001, output_dic["src_elasticProton"],   n_bins, x_min, x_max)
    h_src_multiNuc     = create_histo("h_src_mul",       "Multiple Nucleons",   ROOT.kRed+2,    3001, output_dic["src_multipleNucleons"],n_bins, x_min, x_max)
    h_src_nucCluster   = create_histo("h_src_nuc",       "Nuclear Clusters",    ROOT.kViolet-4, 3001, output_dic["src_nuclearClusters"], n_bins, x_min, x_max)
    h_src_protonPion   = create_histo("h_src_pion",      "Proton + Pion",       ROOT.kYellow-9, 3001, output_dic["src_protonPion"],      n_bins, x_min, x_max)
    h_src_other        = create_histo("h_src_oth",       "Other",               ROOT.kGray,     3001, output_dic["src_other"],           n_bins, x_min, x_max)
    
    h_src_elasticProton_T = create_histo("h_src_elasticP_t",  "Elastic Proton, transparent",      ROOT.kGreen-4,  3001, output_dic["src_elasticProton_T"],   n_bins, x_min, x_max)
    h_src_multiNuc_T     = create_histo("h_src_mul_T",       "Multiple Nucleons, transparent",   ROOT.kRed-4,    3001, output_dic["src_multipleNucleons_T"],n_bins, x_min, x_max)
    h_src_nucCluster_T   = create_histo("h_src_nuc_T",       "Nuclear Clusters, transparent",    ROOT.kViolet+4, 3001, output_dic["src_nuclearClusters_T"], n_bins, x_min, x_max)
   

    pad1_c6.cd() 

    hs_src = ROOT.THStack("hs_src", "SRC Proton Momentum by Channel;Momentum [MeV/c];Events")
    
    # Adding to stack (bottom to top)
    hs_src.Add(h_src_other)
    hs_src.Add(h_src_protonPion)
    hs_src.Add(h_src_nucCluster)
    hs_src.Add(h_src_nucCluster_T)
    hs_src.Add(h_src_multiNuc)
    hs_src.Add(h_src_multiNuc_T)   
    hs_src.Add(h_src_elasticProton)
    hs_src.Add(h_src_elasticProton_T)
    hs_src.Add(h_src_noCascade)
    hs_src.Add(h_src_qe)

    hs_src.Draw("hist") 
    hs_src.GetXaxis().SetLabelSize(0)
    hs_src.SetTitle("; ;Number of Events") 

    hs_src.GetYaxis().SetLabelSize(0.04)
    hs_src.GetYaxis().SetTitleSize(0.05)
    hs_src.GetYaxis().SetTitleOffset(1.1)

    max_h_src = hs_src.GetMaximum()
    hs_src.SetMaximum(max_h_src * 1.2) 

    total_integral_src = sum(hist.Integral() for hist in hs_src.GetHists())
    print(f"Total visible events in the SRC stack: {total_integral_src}")

    # Adjusted legend size
    legend_src = ROOT.TLegend(0.45, 0.45, 0.88, 0.88)
    legend_src.SetBorderSize(0)
    legend_src.SetFillStyle(0) # Transparent
    legend_src.AddEntry(h_src_qe,            "QE proton + de-excitation", "f")
    legend_src.AddEntry(h_src_noCascade,     "No cascade FSI",            "f")
    legend_src.AddEntry(h_src_elasticProton, "Elastic proton",            "f")
    legend_src.AddEntry(h_src_multiNuc,      "Multiple nucleons",         "f")
    legend_src.AddEntry(h_src_nucCluster,    "Nuclear clusters",          "f")
    legend_src.AddEntry(h_src_elasticProton_T, "Elastic proton + transparent","f")
    legend_src.AddEntry(h_src_multiNuc_T,      "Multiple nucleons + transparent", "f")
    legend_src.AddEntry(h_src_nucCluster_T,    "Nuclear clusters + transparent", "f")
    legend_src.AddEntry(h_src_protonPion,    "Proton + pion",             "f")
    legend_src.AddEntry(h_src_other,         "Other",                     "f")
    legend_src.Draw()

    pad2_c6.cd() 

    # Create total histogram for ratio calculation
    h_total_src = h_src_nucCluster.Clone("h_total_src")
    h_total_src.Add(h_src_noCascade)
    h_total_src.Add(h_src_qe)
    h_total_src.Add(h_src_elasticProton)
    h_total_src.Add(h_src_multiNuc)
    h_total_src.Add(h_src_protonPion)
    h_total_src.Add(h_src_other)
    h_total_src.Add(h_src_elasticProton_T)
    h_total_src.Add(h_src_multiNuc_T)
    h_total_src.Add(h_src_nucCluster_T)
    

    # Create ratios
    r_src_noCascade     = create_ratio(h_src_noCascade, h_total_src)
    r_src_qe            = create_ratio(h_src_qe, h_total_src)
    r_src_elasticProton = create_ratio(h_src_elasticProton, h_total_src)
    r_src_multiNuc      = create_ratio(h_src_multiNuc, h_total_src)
    r_src_nucCluster    = create_ratio(h_src_nucCluster, h_total_src)
    r_src_protonPion    = create_ratio(h_src_protonPion, h_total_src)
    r_src_other         = create_ratio(h_src_other, h_total_src)
    r_src_elasticProton_T = create_ratio(h_src_elasticProton_T, h_total_src)
    r_src_multiNuc_T      = create_ratio(h_src_multiNuc_T, h_total_src)
    r_src_nucCluster_T    = create_ratio(h_src_nucCluster_T, h_total_src)

    # Draw ratios
    r_src_noCascade.Draw("HIST") 
    r_src_noCascade.SetTitle(";P_{aft} (MeV/c);Fractional Contributions") 

    y_axis_src = r_src_noCascade.GetYaxis()
    y_axis_src.SetRangeUser(0, 1.1)
    y_axis_src.SetNdivisions(505)
    y_axis_src.SetLabelSize(0.08)
    y_axis_src.SetTitleSize(0.09)
    y_axis_src.SetTitleOffset(0.5)

    x_axis_src = r_src_noCascade.GetXaxis()
    x_axis_src.SetLabelSize(0.08)
    x_axis_src.SetTitleSize(0.10)
    x_axis_src.SetTitleOffset(1.0)

    # Draw the remaining ratios on top
    r_src_qe.Draw("HIST SAME")
    r_src_elasticProton.Draw("HIST SAME")
    r_src_multiNuc.Draw("HIST SAME")
    r_src_nucCluster.Draw("HIST SAME")
    r_src_protonPion.Draw("HIST SAME")
    r_src_other.Draw("HIST SAME")
    r_src_elasticProton_T.Draw("HIST SAME")
    r_src_multiNuc_T.Draw("HIST SAME")
    r_src_nucCluster_T.Draw("HIST SAME")
    
    c6.Update()

    c7 = ROOT.TCanvas("c7", "Excitation_energy_distribution_ex", 800, 600)
    c7.cd() 
    pad7 = ROOT.TPad("pad7", "pad7", 0, 0.3, 1, 1.0)
    pad7.SetBottomMargin(0.02) 
    pad7.SetLeftMargin(0.12)
    pad7.SetTicks(1, 1)        
    pad7.Draw()

    n_bins = 300
    x_min = -5
    x_max = 60

    # Extracting the lists directly from the dictionary
    h_ex_nucCluster = create_histo("h_ex_nuc", "Nuclear Clusters",       ROOT.kViolet-4, 3001, output_dic["ex_nuclearCluster"], n_bins, x_min, x_max)
    h_ex_noCascade  = create_histo("h_ex_noC", "No Cascade FSI",         ROOT.kOrange-4, 3001, output_dic["ex_noCascadeFSI"], n_bins, x_min, x_max)
    h_ex_QE         = create_histo("h_ex_qe",  "QE Proton + De-ex",      ROOT.kOrange+7, 3001, output_dic["ex_QE_deex"], n_bins, x_min, x_max)
    h_ex_multiNuc   = create_histo("h_ex_mul", "Multiple Nucleons",      ROOT.kRed+2,    3001, output_dic["ex_multipleNucleon_noCluster"], n_bins, x_min, x_max)
    h_ex_oneProton  = create_histo("h_ex_one", "One Proton",             ROOT.kGreen,    3001, output_dic["ex_oneProton"], n_bins, x_min, x_max)
    h_ex_protonPion = create_histo("h_ex_pion","Proton + Pion",          ROOT.kYellow-9, 3001, output_dic["ex_protonPion"], n_bins, x_min, x_max)
    h_ex_other      = create_histo("h_ex_oth", "No protons",             ROOT.kGray,     3001, output_dic["ex_no_protons"], n_bins, x_min, x_max)

    hs7 = ROOT.THStack("hs7_ExcitationE", "excitationE_distribution_ex")
    hs7.Add(h_ex_protonPion)
    hs7.Add(h_ex_other)
    hs7.Add(h_ex_oneProton)
    hs7.Add(h_ex_multiNuc)
    hs7.Add(h_ex_nucCluster)
    hs7.Add(h_ex_noCascade)
    hs7.Add(h_ex_QE)

    hs7.Draw("hist") 
    hs7.SetTitle("; excitation E [MeV];Number of Events") 

    hs7.GetYaxis().SetLabelSize(0.03)
    hs7.GetYaxis().SetTitleSize(0.04)
    hs7.GetYaxis().SetTitleOffset(1.1)

    hs7.GetXaxis().SetTitle("Excitation Energy [MeV]")

    max_h = hs7.GetMaximum()
    hs7.SetMaximum(max_h * 1.2) 

    legend7 = ROOT.TLegend(0.55, 0.55, 0.88, 0.88)
    legend7.SetBorderSize(0)
    legend7.SetFillStyle(0) # Transparent
    legend7.AddEntry(h_ex_noCascade, "no cascade FSI", "f")
    legend7.AddEntry(h_ex_QE, "QE proton + de-excitation", "f")
    legend7.AddEntry(h_ex_multiNuc, "multiple nucleons", "f")
    legend7.AddEntry(h_ex_nucCluster, "nuclear clusters", "f")
    legend7.AddEntry(h_ex_oneProton, "one proton", "f")
    legend7.AddEntry(h_ex_protonPion, "proton + pion", "f")
    legend7.AddEntry(h_ex_other, "no protons", "f")

    legend7.Draw()

    c9 = ROOT.TCanvas("c9", "Excitation_energy_distribution_ex", 800, 600)
    c9.cd() 
    pad9 = ROOT.TPad("pad7", "pad7", 0, 0.3, 1, 1.0)
    pad9.SetBottomMargin(0.02) 
    pad9.SetLeftMargin(0.12)
    pad9.SetTicks(1, 1)        
    pad9.Draw()

    n_bins = 150
    x_min = -50
    x_max = 60

    # Extracting the lists directly from the dictionary
    h_ex_MF         = create_histo("h_ex_mf", "Mean Field",             ROOT.kOrange+7, 3001, output_dic["ex_mf"], n_bins, x_min, x_max)
    h_ex_src        = create_histo("h_ex_src", "SRC",                   ROOT.kBlue,     3001, output_dic["ex_src"], n_bins, x_min, x_max)
    h_ex_2p2h       = create_histo("h_ex_2p2h", "2p2h",                 ROOT.kTeal+3,    3001, output_dic["ex_2p2h"], n_bins, x_min, x_max)



    hs9 = ROOT.THStack("hs9_ExcitationE", "excitationE_distribution_ex")
    hs9.Add(h_ex_MF)
    hs9.Add(h_ex_src)
    hs9.Add(h_ex_2p2h)

    hs9.Draw("hist") 
    hs9.SetTitle("; excitation E [MeV];Number of Events") 

    hs9.GetYaxis().SetLabelSize(0.03)
    hs9.GetYaxis().SetTitleSize(0.04)
    hs9.GetYaxis().SetTitleOffset(1.1)

    hs9.GetXaxis().SetTitle("Excitation Energy [MeV]")

    max_h = hs9.GetMaximum()
    hs9.SetMaximum(max_h * 1.2) 

    legend9 = ROOT.TLegend(0.55, 0.55, 0.88, 0.88)
    legend9.SetBorderSize(0)
    legend9.SetFillStyle(0) # Transparent
    legend9.AddEntry(h_ex_MF, "Mean Field", "f")
    legend9.AddEntry(h_ex_src, "src", "f")
    legend9.AddEntry(h_ex_2p2h, "2p2h", "f")

    legend9.Draw()


    c8 = ROOT.TCanvas("c8", "2p2h Proton Momentum Distributions", 800, 600)
    c8.cd() 
    
    # Appended _c8 to pad names to avoid overwriting your c1 pads if run in the same script
    pad1_c8 = ROOT.TPad("pad1_c8", "pad1", 0, 0.3, 1, 1.0)
    pad1_c8.SetBottomMargin(0.02) 
    pad1_c8.SetLeftMargin(0.12)
    pad1_c8.SetTicks(1, 1)        
    pad1_c8.Draw()

    pad2_c8 = ROOT.TPad("pad2_c8", "pad2", 0, 0.0, 1, 0.3)
    pad2_c8.SetTopMargin(0.01)  
    pad2_c8.SetBottomMargin(0.25) 
    pad2_c8.SetLeftMargin(0.12)
    pad2_c8.SetTicks(1, 1)
    pad2_c8.Draw()

    n_bins = 20
    x_min = 0
    x_max = 1500   

    # Creating histograms using the updated 2p2h_ lists
    h_2p2h_noCascade    = create_histo("h_2p2h_noC",       "No Cascade FSI",      ROOT.kOrange-4, 3001, output_dic["2p2h_noCascadeFSI"],    n_bins, x_min, x_max)
    h_2p2h_qe           = create_histo("h_2p2h_qe",        "QE Proton + De-ex",   ROOT.kOrange+7, 3001, output_dic["2p2h_qeDeEX"],          n_bins, x_min, x_max)
    h_2p2h_elasticProton= create_histo("h_2p2h_elasticP",  "Elastic Proton",      ROOT.kGreen+2,  3001, output_dic["2p2h_elasticProton"],   n_bins, x_min, x_max)
    h_2p2h_multiNuc     = create_histo("h_2p2h_mul",       "Multiple Nucleons",   ROOT.kRed+2,    3001, output_dic["2p2h_multipleNucleons"],n_bins, x_min, x_max)
    h_2p2h_nucCluster   = create_histo("h_2p2h_nuc",       "Nuclear Clusters",    ROOT.kViolet-4, 3001, output_dic["2p2h_nuclearClusters"], n_bins, x_min, x_max)
    h_2p2h_protonPion   = create_histo("h_2p2h_pion",      "Proton + Pion",       ROOT.kYellow-9, 3001, output_dic["2p2h_protonPion"],      n_bins, x_min, x_max)
    h_2p2h_other        = create_histo("h_2p2h_oth",       "Other",               ROOT.kGray,     3001, output_dic["2p2h_other"],           n_bins, x_min, x_max)
    
    pad1_c8.cd() 

    hs_2p2h = ROOT.THStack("hs_2p2h", "2p2h Proton Momentum by Channel;Momentum [MeV/c];Events")
    
    # Adding to stack (bottom to top)
    hs_2p2h.Add(h_2p2h_other)
    hs_2p2h.Add(h_2p2h_protonPion)
    hs_2p2h.Add(h_2p2h_nucCluster)
    hs_2p2h.Add(h_2p2h_multiNuc)
    hs_2p2h.Add(h_2p2h_elasticProton)
    hs_2p2h.Add(h_2p2h_noCascade)
    hs_2p2h.Add(h_2p2h_qe)

    hs_2p2h.Draw("hist") 
    hs_2p2h.GetXaxis().SetLabelSize(0)
    hs_2p2h.SetTitle("; ;Number of Events") 

    hs_2p2h.GetYaxis().SetLabelSize(0.04)
    hs_2p2h.GetYaxis().SetTitleSize(0.05)
    hs_2p2h.GetYaxis().SetTitleOffset(1.1)

    max_h_2p2h = hs_2p2h.GetMaximum()
    hs_2p2h.SetMaximum(max_h_2p2h * 1.2) 

    total_integral_2p2h = sum(hist.Integral() for hist in hs_2p2h.GetHists())
    print(f"Total visible events in the 2p2h stack: {total_integral_2p2h}")

    # Adjusted legend size
    legend_2p2h = ROOT.TLegend(0.45, 0.45, 0.88, 0.88)
    legend_2p2h.SetBorderSize(0)
    legend_2p2h.SetFillStyle(0) # Transparent
    legend_2p2h.AddEntry(h_2p2h_qe,            "QE proton + de-excitation", "f")
    legend_2p2h.AddEntry(h_2p2h_noCascade,     "No cascade FSI",            "f")
    legend_2p2h.AddEntry(h_2p2h_elasticProton, "Elastic proton",            "f")
    legend_2p2h.AddEntry(h_2p2h_multiNuc,      "Multiple nucleons",         "f")
    legend_2p2h.AddEntry(h_2p2h_nucCluster,    "Nuclear clusters",          "f")
    legend_2p2h.AddEntry(h_2p2h_protonPion,    "Proton + pion",             "f")
    legend_2p2h.AddEntry(h_2p2h_other,         "Other",                     "f")
    legend_2p2h.Draw()

    pad2_c8.cd() 

    # Create total histogram for ratio calculation
    h_total_2p2h = h_2p2h_nucCluster.Clone("h_total_2p2h")
    h_total_2p2h.Add(h_2p2h_noCascade)
    h_total_2p2h.Add(h_2p2h_qe)
    h_total_2p2h.Add(h_2p2h_elasticProton)
    h_total_2p2h.Add(h_2p2h_multiNuc)
    h_total_2p2h.Add(h_2p2h_protonPion)
    h_total_2p2h.Add(h_2p2h_other)

    # Create ratios
    r_2p2h_noCascade     = create_ratio(h_2p2h_noCascade, h_total_2p2h)
    r_2p2h_qe            = create_ratio(h_2p2h_qe, h_total_2p2h)
    r_2p2h_elasticProton = create_ratio(h_2p2h_elasticProton, h_total_2p2h)
    r_2p2h_multiNuc      = create_ratio(h_2p2h_multiNuc, h_total_2p2h)
    r_2p2h_nucCluster    = create_ratio(h_2p2h_nucCluster, h_total_2p2h)
    r_2p2h_protonPion    = create_ratio(h_2p2h_protonPion, h_total_2p2h)
    r_2p2h_other         = create_ratio(h_2p2h_other, h_total_2p2h)

    # Draw ratios
    r_2p2h_noCascade.Draw("HIST") 
    r_2p2h_noCascade.SetTitle(";P_{aft} (MeV/c);Fractional Contributions") 

    y_axis_2p2h = r_2p2h_noCascade.GetYaxis()
    y_axis_2p2h.SetRangeUser(0, 1.1)
    y_axis_2p2h.SetNdivisions(505)
    y_axis_2p2h.SetLabelSize(0.08)
    y_axis_2p2h.SetTitleSize(0.09)
    y_axis_2p2h.SetTitleOffset(0.5)

    x_axis_2p2h = r_2p2h_noCascade.GetXaxis()
    x_axis_2p2h.SetLabelSize(0.08)
    x_axis_2p2h.SetTitleSize(0.10)
    x_axis_2p2h.SetTitleOffset(1.0)

    # Draw the remaining ratios on top
    r_2p2h_qe.Draw("HIST SAME")
    r_2p2h_elasticProton.Draw("HIST SAME")
    r_2p2h_multiNuc.Draw("HIST SAME")
    r_2p2h_nucCluster.Draw("HIST SAME")
    r_2p2h_protonPion.Draw("HIST SAME")
    r_2p2h_other.Draw("HIST SAME")
    
    c8.Update()


    file.cd()
    c0.Write()
    c1.Write()
    c2.Write()
    c3.Write()
    c4.Write()
    c6.Write()
    c7.Write()
    c9.Write()
    c8.Write()
    file.Close()
    #c1.SaveAs("proton_channels_{}.root".format(filename_chunk)) 

    input("Saved as {}".format("incl_analysis_test_{}.root".format(filename_chunk)))
