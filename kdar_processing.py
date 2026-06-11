import ROOT
from incl_analysis_functions import *
import sys 
import numpy as np
import math
ROOT.gSystem.Load("libNEUTROOTClass.so")
ROOT.gSystem.Load("libNEUTOutput.so")
ROOT.gSystem.Load("libNEUTReWeight.so")
from enum import Enum
            
class kdar_measurements(nvect_reader):
    def __init__(self, nvect_,flavour = 2212):
        super().__init__(nvect_,flavour)


    def cluster_mom(self,flavour):
        cluster_mom = []
        for i in range(self.nopart):
            pinfo = self.nvect.PartInfo(i)
            if ((pinfo.fPID == abs(flavour)) and  (pinfo.fIsAlive == 1)):
                cluster_mom.append(np.linalg.norm(np.array([pinfo.fP.X(),pinfo.fP.Y(), pinfo.fP.Z()])))
        
        return cluster_mom
            
            
    def remnant_flavour(self):
        for i in range(self.nopart):
            pinfo = self.nvect.PartInfo(i)
            if ((pinfo.fPID > 1000) and (pinfo.fPID != abs(2212)) and  (pinfo.fPID != abs(2112)) and (pinfo.fPID <10000)):
                return pinfo.fPID 


    def Missing_energy_kdar(self, D2= False):
        E_nu = 0.0
        E_lep = 0.0 
        T_had = []
        proton_counter = 0
        for i in range(self.nopart):
            pinfo = self.nvect.PartInfo(i)

            if pinfo.fPID == abs(2212) and (pinfo.fIsAlive == 1):
                T_had.append(np.sqrt(np.linalg.norm(np.array([pinfo.fP.X(),pinfo.fP.Y(), pinfo.fP.Z()]))**2 + 938.272**2) - 938.272)
                proton_counter += 1
            if pinfo.fPID == abs(14):
                E_nu = np.linalg.norm(np.array([pinfo.fP.X(),pinfo.fP.Y(), pinfo.fP.Z()])) 
            if pinfo.fPID == abs(13):
                T_lep = np.sqrt(np.linalg.norm(np.array([pinfo.fP.X(),pinfo.fP.Y(), pinfo.fP.Z()]))**2 + 105.7**2) - 105.7
            if (D2 == True and pinfo.fPID == abs(1002) and (pinfo.fIsAlive == 1)):
                T_had.append(np.sqrt(np.linalg.norm(np.array([pinfo.fP.X(),pinfo.fP.Y(), pinfo.fP.Z()]))**2 + (2.014*931.5)**2) - (2*931.5))
            
        E_vis = T_lep + sum(T_had)

        E_miss = 131.2 - E_vis

        if E_vis > 20 and E_vis < 150 and proton_counter != 0:
            return E_vis, E_miss
        
        return None , None

def kdar_event_loop(t):


    output_dic = {
    "noCascadeFSI": [],
    "QE_deex": [],
    "multipleNucleon_noCluster": [],
    "nuclearCluster": [],
    "oneProton": [],
    "protonPion": [],
    "no_protons": [],
    "other": [],
    
    "deuteronMom": [],
    "He3Mom": [],
    "alphaMom": [],
    "protMom": [],
    
    "vis_E_D2": [],
    "vis_E": [],
    "missE_kdar": [],
    "missE_D2": []
}

    for event in t:
        nvect = event.vectorbranch
        nvect_class = kdar_measurements(nvect)
        
        prefsi_proton_mom = nvect_class.PreFSIProtMom
        eventType = nvect_class.eventType
        intChannel = nvect_class.intChannel
        miss_E = nvect_class.E_miss
        miss_p = nvect_class.P_miss

        if (intChannel == intChannel_CCQE.muOnly) or (intChannel == intChannel_CCQE.neutronPion):
                continue
        
        excitation_E, miss_E_kdar = nvect_class.Missing_energy_kdar()
        missing_E_D2, miss_E_d2 = nvect_class.Missing_energy_kdar(D2=True)

        if excitation_E is None or missing_E_D2 is None:
            continue
        
        # Already correct in your original code
        output_dic["vis_E_D2"].append(missing_E_D2)
        output_dic["vis_E"].append(excitation_E)

        # Updated to use dictionary keys
        output_dic["missE_kdar"].append(miss_E_kdar)
        output_dic["missE_D2"].append(miss_E_d2)

        deuteron_mom = nvect_class.cluster_mom(1002)
        if deuteron_mom != False:
            output_dic["deuteronMom"].extend(deuteron_mom)

        alpha_mom = nvect_class.cluster_mom(2004)
        if alpha_mom != False:
            output_dic["alphaMom"].extend(alpha_mom)

        He3_mom = nvect_class.cluster_mom(2003)
        if He3_mom != False:
            output_dic["He3Mom"].extend(He3_mom)

        proton_mom = nvect_class.cluster_mom(2212)
        if proton_mom != False:
            output_dic["protMom"].extend(proton_mom)
        

        if eventType == EventType.MF:
            
            if intChannel == intChannel_CCQE.noCascadeFSI:
                output_dic["noCascadeFSI"].append(excitation_E)

            elif intChannel == intChannel_CCQE.qeDeEX:
                output_dic["QE_deex"].append(excitation_E)
            
            elif intChannel == intChannel_CCQE.oneProton:
                output_dic["oneProton"].append(excitation_E)
            
            elif intChannel == intChannel_CCQE.multipleNucleon:
                output_dic["multipleNucleon_noCluster"].append(excitation_E)

            elif intChannel == intChannel_CCQE.nuclearCluster:
                output_dic["nuclearCluster"].append(excitation_E)

            elif intChannel == intChannel_CCQE.protonPion: 
                output_dic["protonPion"].append(excitation_E)

            elif intChannel == intChannel_CCQE.other:
                other_counter += 1 
                output_dic["other"].append(excitation_E)

    return output_dic 

def missing_energy_kdar(filename1,filename_list=None):
    filename_chunk = filename1.split(".")[0].split("out_")[1]

    f = ROOT.TFile(filename1)
    t = f.Get("neuttree")
    output_dic_dic = {}
    filename_chunk_list = []
    hist_list = []
    output_dic = kdar_event_loop(t)

    if filename_list:
        for file in filename_list:
            filename_chunk_temp = (file.split("out_KDAR_")[1].split("_CCQE")[0]).lower()
            filename_chunk_list.append(filename_chunk_temp)
            f2 = ROOT.TFile(file)
            t2 = f2.Get("neuttree")
            output_dic_dic[filename_chunk_temp] = kdar_event_loop(t2)

        file = ROOT.TFile("kdar_file_comparisons_{}.root".format(filename_chunk),"RECREATE")
    else:
        file = ROOT.TFile("kdar_file_{}.root".format(filename_chunk),"RECREATE")



    c0 = ROOT.TCanvas("c0", "deuteron_momentum", 800, 600)
    c0.cd() 
    h = ROOT.TH1F("deuteron_momentum", "deuteron_momentum", 20, 0, 400)
    for val in output_dic["deuteronMom"]:
        h.Fill(val)

    h2 = ROOT.TH1F("He3mom", "he3mom", 20, 0, 400)
    for val in output_dic["He3Mom"]:
        h2.Fill(val)

    h3 = ROOT.TH1F("Alphamom", "Alphamom", 20, 0, 400)
    for val in output_dic["alphaMom"]:
        h3.Fill(val)

    h4 = ROOT.TH1F("protmom", "protmom", 20, 0, 400)
    for val in output_dic["protMom"]:
        h4.Fill(val)

    # 2. Normalize the histograms (protect against empty histograms)
    if h.Integral() > 0: h.Scale(1.0 / h.Integral())
    if h2.Integral() > 0: h2.Scale(1.0 / h2.Integral())
    if h3.Integral() > 0: h3.Scale(1.0 / h3.Integral())
    if h4.Integral() > 0: h4.Scale(1.0 / h4.Integral())

    # 3. Apply different colors, line styles, and thicknesses
    # Deuteron (Blue, Solid)
    h.SetLineColor(ROOT.kBlue)
    h.SetLineStyle(1) 
    h.SetLineWidth(2)

    # Helium-3 (Red, Dashed)
    h2.SetLineColor(ROOT.kRed)
    h2.SetLineStyle(2) 
    h2.SetLineWidth(2)

    # Alpha (Dark Green, Dotted)
    h3.SetLineColor(ROOT.kGreen + 2)
    h3.SetLineStyle(3) 
    h3.SetLineWidth(2)

    h4.SetLineColor(ROOT.kViolet + 2)
    h4.SetLineStyle(3) 
    h4.SetLineWidth(2)

    # Update titles and axis formatting
    h.SetTitle("; Cluster Momentum [MeV/c]; Normalized Events / 20 MeV") 

    h.GetYaxis().SetLabelSize(0.03)
    h.GetYaxis().SetTitleSize(0.04)
    h.GetYaxis().SetTitleOffset(1.1)
    h.GetXaxis().SetTitle("Cluster Momentum [MeV/c]")

    # Calculate the global maximum across all three histograms to prevent clipping
    max_h = max(h.GetMaximum(), h2.GetMaximum(), h3.GetMaximum())
    h.SetMaximum(max_h * 1.2) 

    # Draw histograms
    h.Draw("hist") 
    h2.Draw("hist same")
    h3.Draw("hist same")
    h4.Draw("hist same")

    # 4. Add a legend
    leg = ROOT.TLegend(0.65, 0.70, 0.88, 0.88) # (x1, y1, x2, y2) in NDC coordinates
    leg.SetBorderSize(0) # Removes the box around the legend for a cleaner look
    leg.AddEntry(h, "Deuteron", "l")
    leg.AddEntry(h2, "^{3}He", "l")
    leg.AddEntry(h3, "Alpha", "l")
    leg.AddEntry(h4, "Proton", "l")
    leg.Draw()


    c1 = ROOT.TCanvas("c1", "missing_E_distribution", 800, 600)
    c1.cd() 
    pad1 = ROOT.TPad("pad1", "pad1", 0, 0.3, 1, 1.0)
    pad1.SetBottomMargin(0.02) 
    pad1.SetLeftMargin(0.12)
    pad1.SetTicks(1, 1)        
    pad1.Draw()

    x_min = 0
    x_max = 300

    n_bins = int((x_max - x_min)/4)
    print(n_bins)

    h_nucCluster = create_histo("h_nuc", "Nuclear Clusters",       ROOT.kViolet-4, 3001, output_dic["nuclearCluster"],n_bins,x_min,x_max)
    h_noCascade  = create_histo("h_noC", "No Cascade FSI",         ROOT.kOrange-4, 3001, output_dic["noCascadeFSI"],n_bins,x_min,x_max)
    h_QE         = create_histo("h_qe",  "QE Proton + De-ex",      ROOT.kOrange+7, 3001, output_dic["QE_deex"],n_bins,x_min,x_max)
    h_multiNuc   = create_histo("h_mul", "Multiple Nucleons",      ROOT.kRed+2,    3001, output_dic["multipleNucleon_noCluster"],n_bins,x_min,x_max)
    h_oneProton  = create_histo("h_one", "One Proton",             ROOT.kGreen,  3001, output_dic["oneProton"],n_bins,x_min,x_max)
    h_protonPion = create_histo("h_pion","Proton + Pion",          ROOT.kYellow-9, 3001, output_dic["protonPion"],n_bins,x_min,x_max)
    h_other      = create_histo("h_oth", "No protons",             ROOT.kGray,     3001, output_dic["no_protons"],n_bins,x_min,x_max)


    hs1 = ROOT.THStack("hs_missing_E_kdar", "missing_E_kdar_distribution")
    hs1.Add(h_protonPion)
    hs1.Add(h_other)
    hs1.Add(h_oneProton)
    hs1.Add(h_multiNuc)
    hs1.Add(h_nucCluster)
    hs1.Add(h_noCascade)
    hs1.Add(h_QE)

    hs1.Draw("hist") 
    hs1.SetTitle("; Missing E [MeV];Number of Events") 

    hs1.GetYaxis().SetLabelSize(0.03)
    hs1.GetYaxis().SetTitleSize(0.04)
    hs1.GetYaxis().SetTitleOffset(1.1)

    hs1.GetXaxis().SetTitle("missing_E_kdar [MeV]")

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

    legend1.Draw()

    c2 = ROOT.TCanvas("c2", "normalized_visible_energy", 800, 600)
    c2.cd() 

    h_Evis = create_histo("h_NoD2", "Without D2", ROOT.kViolet-4, 0, output_dic["vis_E"], n_bins, x_min, x_max)
    h_E_vis_D2 = create_histo("h_withD2", "with D2", ROOT.kOrange-4, 0, output_dic["vis_E_D2"], n_bins, x_min, x_max)


    if h_Evis.Integral() > 0:
        h_Evis.Scale(1.0 / h_Evis.Integral())

    if h_E_vis_D2.Integral() > 0:
        h_E_vis_D2.Scale(1.0 / h_E_vis_D2.Integral())

    h_Evis.SetLineWidth(2)
    h_E_vis_D2.SetLineWidth(2)

    hs2 = ROOT.THStack("hs_missing_E_kdar", "Normalized Visible Energy")
    hs2.Add(h_Evis)
    hs2.Add(h_E_vis_D2)

    hs2.Draw("nostack hist") 

    hs2.GetXaxis().SetTitle("Visible Energy [MeV]")
    hs2.GetYaxis().SetTitle("Arbitrary Units (Normalized)")

    max_h = hs2.GetMaximum("nostack")
    hs2.SetMaximum(max_h * 1.2) 

    legend2 = ROOT.TLegend(0.55, 0.70, 0.88, 0.88)
    legend2.SetBorderSize(0)
    legend2.SetFillStyle(0) 
    legend2.AddEntry(h_Evis, "No D2 (Normalized)", "l")
    legend2.AddEntry(h_E_vis_D2, "with D2 (Normalized)", "l")

    legend2.Draw()
    c2.Update()


    x_min = -15.0
    x_max = 110.0
    n_bins = int((x_max - x_min) / 5.0) 

    c3 = ROOT.TCanvas("c3", "normalized_missing_energy", 800, 600)
    c3.cd() 

    h_missE    = create_histo("h_missE_noD2", "Without D2", ROOT.kBlue-7,  0, output_dic["missE_kdar"], n_bins, x_min, x_max, line=True)
    h_missE_D2 = create_histo("h_missE_withD2", "With D2",    ROOT.kRed-7,   1, output_dic["missE_D2"],   n_bins, x_min, x_max,line=True)
    if filename_chunk_list:
        counter = 0
        for file_name in filename_chunk_list:
            counter += 1
            hist_list.append( create_histo("h_missE_{}".format(counter), "With D2", ( 2 + ((counter) % 48) ),counter+2, (output_dic_dic[file_name])["missE_D2"], n_bins, x_min, x_max,line=True)) 
   

    if h_missE.Integral() > 0:
        h_missE.Scale(1.0 / h_missE.Integral())

    if h_missE_D2.Integral() > 0:
        h_missE_D2.Scale(1.0 / h_missE_D2.Integral())

    if filename_chunk_list:
        counter = 0
        for hist in hist_list:
                if hist.Integral() > 0:
                    hist.Scale(1.0 / hist.Integral())


    h_missE.SetLineWidth(3)
    h_missE_D2.SetLineWidth(3)
    if filename_chunk_list:
        counter = 0
        for hist in hist_list:
            hist.SetLineWidth(3)


    hs3 = ROOT.THStack("hs3", "kDar Missing Energy Comparison")
    hs3.Add(h_missE)
    hs3.Add(h_missE_D2)
    if filename_chunk_list:
        counter = 0
        for hist in hist_list:
            hs3.Add(hist)


    hs3.Draw("nostack hist") 
    hs3.GetXaxis().SetTitle("Missing Energy [MeV]")
    hs3.GetYaxis().SetTitle("Probability Density")
    hs3.GetYaxis().SetTitleOffset(1.2)


    max_h3 = hs3.GetMaximum("nostack")
    hs3.SetMaximum(max_h3 * 1.2) 


    data_bin_centers = [
        -7.5, -2.5, 2.5, 7.5, 12.5, 17.5, 22.5, 27.5, 32.5, 37.5, 
        42.5, 47.5, 52.5, 57.5, 62.5, 67.5, 72.5, 77.5, 82.5, 87.5, 92.5
    ]
    data_values = [
        0.005, 0.018, 0.015, 0.015, 0.050, 0.330, 0.160, 0.057, 0.060, 0.065, 
        0.055, 0.045, 0.030, 0.020, 0.020, 0.010, 0.010, 0.008, 0.008, 0.005, 0.005
    ]
    data_errors = [
        0.005, 0.015, 0.010, 0.010, 0.010, 0.035, 0.020, 0.015, 0.015, 0.015, 
        0.015, 0.015, 0.010, 0.010, 0.010, 0.008, 0.008, 0.006, 0.006, 0.005, 0.005
    ]

    h_data = ROOT.TH1F("h_data", "Data (Stat. Error)", n_bins, x_min, x_max)

    for center, val, err in zip(data_bin_centers, data_values, data_errors):
        bin_idx = h_data.FindBin(center)
        h_data.SetBinContent(bin_idx, val)
        h_data.SetBinError(bin_idx, err)

    h_data.SetLineColor(ROOT.kBlack)
    h_data.SetLineWidth(2)
    h_data.SetMarkerColor(ROOT.kBlack)
    h_data.SetMarkerStyle(1)

    # Draw the data on top of the existing THStack
    h_data.Draw("SAME HIST P")
    h_data.Draw("SAME E")


    legend3 = ROOT.TLegend(0.45, 0.65, 0.88, 0.88)
    legend3.SetBorderSize(0)
    legend3.SetFillStyle(0) 

    legend3.AddEntry(h_data, "Data (Stat. Error)", "le")
    legend3.AddEntry(h_missE, "No D2 (Normalised)", "l")
    legend3.AddEntry(h_missE_D2, "With D2 (Normalised)", "l")
    if filename_chunk_list:
        counter = 0
        for hist, filename_chunk in zip(hist_list, filename_chunk_list):
            legend3.AddEntry(hist, filename_chunk, "l")
    legend3.Draw()
    c3.Update()

    file.cd()
    c0.Write()
    c1.Write()
    c2.Write()
    c3.Write()
    file.Close()