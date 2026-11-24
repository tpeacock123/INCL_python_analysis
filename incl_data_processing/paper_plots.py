import numpy as np
import array
import sys 
import argparse
import ROOT
import os
from incl_analysis_functions import *


ROOT.gROOT.GetColor(ROOT.kBlue).SetRGB(87 / 255.0, 144 / 255.0, 252 / 255.0)  # #5790FC
ROOT.gROOT.GetColor(ROOT.kOrange + 7).SetRGB(
    248 / 255.0, 156 / 255.0, 32 / 255.0
)  # #F89C20
ROOT.gROOT.GetColor(ROOT.kRed).SetRGB(228 / 255.0, 37 / 255.0, 54 / 255.0)  # #E42536
ROOT.gROOT.GetColor(ROOT.kMagenta).SetRGB(
    150 / 255.0, 74 / 255.0, 139 / 255.0
)  # #964A8B

ROOT.gROOT.GetColor(ROOT.kAzure).SetRGB(
    156 / 255.0, 156 / 255.0, 161 / 255.0
)  # #9c9ca1

font_id = 132
ROOT.gStyle.SetTextFont(font_id)
ROOT.gStyle.SetLegendFont(font_id)
ROOT.gStyle.SetLabelFont(font_id, "XYZ")  # Axis tick labels
ROOT.gStyle.SetTitleFont(font_id, "XYZ")  # Axis titles
ROOT.gStyle.SetTitleFont(font_id, "")     # Histogram/Pad titles

ROOT.TGaxis.SetMaxDigits(3)
ROOT.gStyle.SetLabelSize(5, "XYZ")        # Axis tick labels
ROOT.gStyle.SetTitleSize(5, "XYZ")

ROOT.gStyle.SetStatFont(font_id)
ROOT.gStyle.SetPadGridX(True)
ROOT.gStyle.SetPadGridY(True)
ROOT.gStyle.SetNdivisions(304, "XY")


ROOT.gSystem.Load("libNEUTROOTClass.so")
ROOT.gSystem.Load("libNEUTOutput.so")
ROOT.gSystem.Load("libNEUTReWeight.so")
ROOT.TH1.AddDirectory(False)
ROOT.TColor.SetGrayscale(False)
ROOT.gROOT.SetBatch(True)



def spectral_function2d(filename, additional_files, name, n_events=500000):
    # Prevent ROOT from spawning interactive GUI windows (Fixes hanging)
    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(0)

    additional_filename = additional_files[0]
    third_filename = additional_files[1]

    filename_chunk = filename.split(".")[0].split("out_")[1]

    # --- File 1 ---
    f = ROOT.TFile(filename)
    t = f.Get("neuttree")
    # Speed up PyROOT loop by disabling unused branches
    t.SetBranchStatus("*", 0)
    t.SetBranchStatus("vectorbranch", 1)

    h_miss = ROOT.TH2D(
        "h_miss", "; P_{N} [MeV/c]; E_{RMV} [MeV/c]", 150, 0, 450, 80, 0, 80
    )
    h1_emiss = ROOT.TH1D(
        "h1_emiss", "; E_{RMV} [MeV/c]; [1/N_{tot} Events/MeV]", 150, 0, 75
    )
    for i, event in enumerate(t):
        if i == n_events:
            break
        nvect_class = nvect_reader(event.vectorbranch)
        h_miss.Fill(nvect_class.P_miss, nvect_class.E_miss)
        h1_emiss.Fill(nvect_class.E_miss + 1.29)

    # --- File 2 ---
    f2 = ROOT.TFile(additional_filename)
    t2 = f2.Get("neuttree")
    t2.SetBranchStatus("*", 0)
    t2.SetBranchStatus("vectorbranch", 1)

    h2_miss = ROOT.TH2D(
        "h2_miss", "; P_{N} [MeV/c]; E_{RMV} [MeV/c]", 150, 0, 450, 80, 0, 80
    )
    h2_emiss = ROOT.TH1D(
        "h2_emiss", "; E_{RMV} [MeV/c]; [1/N_{tot} Events/MeV]", 150, 0, 75
    )
    for i, event in enumerate(t2):
        if i == n_events:
            break
        nvect_class = nvect_reader(event.vectorbranch)
        h2_miss.Fill(nvect_class.P_miss, nvect_class.E_miss)
        h2_emiss.Fill(nvect_class.E_miss)

    # --- File 3 ---
    f3 = ROOT.TFile(third_filename)
    t3 = f3.Get("neuttree")
    t3.SetBranchStatus("*", 0)
    t3.SetBranchStatus("vectorbranch", 1)

    h3_miss = ROOT.TH2D(
        "h3_miss", "; P_{N} [MeV/c]; E_{RMV} [MeV/c]", 150, 0, 450, 80, 0, 80
    )
    h3_emiss = ROOT.TH1D(
        "h3_emiss", "; E_{RMV} [MeV/c]; [1/N_{tot} Events/MeV]",  150, 0, 75
    )
    for i, event in enumerate(t3):
        if i == n_events:
            break
        nvect_class = nvect_reader(event.vectorbranch)
        h3_miss.Fill(nvect_class.P_miss, nvect_class.E_miss)
        h3_emiss.Fill(nvect_class.E_miss)

    # Scale max and min across all three 2D histograms
    global_max = max(
        h_miss.GetMaximum(), h2_miss.GetMaximum(), h3_miss.GetMaximum()
    )

    for h in (h_miss, h2_miss, h3_miss):
        h.SetMaximum(global_max)
        h.SetMinimum(0.00001)

    ROOT.gStyle.SetPalette(ROOT.kBlueRedYellow)
    ROOT.gStyle.SetNumberContours(99)

    def format_axes(obj):
        """Formats axis labels, titles, and line colors for histograms and THStack."""
        xaxis = obj.GetXaxis()
        yaxis = obj.GetYaxis()
        for axis in [xaxis, yaxis]:
            axis.SetLabelColor(ROOT.kBlack)
            axis.SetTitleColor(ROOT.kBlack)
            axis.SetAxisColor(ROOT.kBlack)
            axis.SetLabelSize(0.04)
            axis.SetTitleSize(0.045)
            axis.CenterTitle(True)
        xaxis.SetTitleOffset(1.1)
        yaxis.SetTitleOffset(1.2)

    # --- PAGE 1: 2D Spectral Functions ---
    c1 = ROOT.TCanvas("c1", "c1", 2400, 700)
    c1.Divide(3, 1)

    for i, h in enumerate([h_miss, h2_miss, h3_miss], start=1):
        pad = c1.cd(i)
        pad.SetFillColor(ROOT.kWhite)
        pad.SetFrameFillColor(ROOT.kWhite)
        pad.SetLeftMargin(0.14)
        pad.SetRightMargin(0.15)
        pad.SetBottomMargin(0.14)
        pad.SetTopMargin(0.08)

        format_axes(h)
        h.Draw("COLZ")

    c1.Update()

    # --- PAGE 2: 1D Overlapping Line Histograms of E_miss ---
    c2 = ROOT.TCanvas("c2", "c2", 600, 600)
    c2.SetLeftMargin(0.14)
    c2.SetRightMargin(0.05)
    c2.SetBottomMargin(0.14)
    c2.SetTopMargin(0.08)

    emiss_hists = [h1_emiss, h2_emiss, h3_emiss]
    labels = ["SF", "ED-RMF", "Nieves LFG"]
    colors = [
        ROOT.kBlue,
        ROOT.kRed,
        ROOT.kOrange + 7,
    ]

    hs_emiss = ROOT.THStack("hs_emiss", "; E_{RMV} [MeV/c]; [1/N_{tot} Events/MeV]")
    legend = ROOT.TLegend(0.58, 0.65, 0.90, 0.88)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.SetTextSize(0.04)

    for h, color, label in zip(emiss_hists, colors, labels):
        if h.Integral() > 0:
            h.Scale(1.0 / h.Integral())

        h.SetFillStyle(0)
        h.SetLineColor(color)
        h.SetLineWidth(2)
        hs_emiss.Add(h)
        legend.AddEntry(h, label, "l")


    c2.cd()
    hs_emiss.Draw("NOSTACK HIST")
    c2.Update()

    format_axes(hs_emiss)
    legend.Draw()
    c2.Update()

    # --- SAVE TO MULTI-PAGE PDF ---
    base_filename = f"spectral_function_{name[0]}"

    pdf_filename = f"{base_filename}.pdf"
    c1.Print(f"{pdf_filename}(") 
    c2.Print(f"{pdf_filename})")  

    root_filename = f"{base_filename}.root"
    root_file = ROOT.TFile(root_filename, "RECREATE")
    c1.Write("c1")  
    c2.Write("c2")  
    root_file.Close()



def SRC_plot(filename, additional_filename):
    ROOT.gStyle.SetOptStat(0)

    def format_axes(hist):
        for axis in [hist.GetXaxis(), hist.GetYaxis()]:
            axis.SetLabelColor(ROOT.kBlack)
            axis.SetTitleColor(ROOT.kBlack)
            axis.SetAxisColor(ROOT.kWhite)
            # Increased label and title size
            axis.SetTitleSize(0.05)
            axis.SetLabelSize(0.045)
            axis.SetTitleOffset(1.1)

    filename_chunk = filename.split(".")[0].split("out_")[1] if "out_" in filename else "output"
    f = ROOT.TFile(filename)
    t = f.Get("neuttree")   

    p_mf_1, e_mf_1 = array.array('d'), array.array('d')
    p_src_1, e_src_1 = array.array('d'), array.array('d')

    for i, event in enumerate(t):
        if i == 50000:
            break
        nvect = event.vectorbranch
        nvect_class = nvect_reader(nvect)
        if nvect_class.eventType == EventType.MF:
            print(nvect_class.P_miss, nvect_class.E_miss)
            p_mf_1.append(nvect_class.P_miss)
            e_mf_1.append(nvect_class.E_miss)
        elif nvect_class.eventType == EventType.SRC:
            print("SRC ", nvect_class.P_miss, nvect_class.E_miss)
            p_src_1.append(nvect_class.P_miss)
            e_src_1.append(nvect_class.E_miss)

    print("NEUT", len(p_src_1) / (len(p_mf_1) + len(p_src_1)))

    g1_mf = ROOT.TGraph(len(p_mf_1), p_mf_1, e_mf_1)
    g1_src = ROOT.TGraph(len(p_src_1), p_src_1, e_src_1)

    f2 = ROOT.TFile(additional_filename)
    t2 = f2.Get("neuttree")  

    p_mf_2, e_mf_2 = array.array('d'), array.array('d')
    p_src_2, e_src_2 = array.array('d'), array.array('d')

    for i, event in enumerate(t):
        if i == 50000:
            break
        nvect = event.vectorbranch
        nvect_class = nvect_reader(nvect)
        if nvect_class.eventType == EventType.MF:
            p_mf_2.append(nvect_class.P_miss)
            e_mf_2.append(nvect_class.E_miss)
        elif nvect_class.eventType == EventType.SRC:
            p_src_2.append(nvect_class.P_miss)
            e_src_2.append(nvect_class.E_miss)

    print("NuWro", len(p_src_2) / (len(p_mf_2) + len(p_src_2)))

    g2_mf = ROOT.TGraph(len(p_mf_2), p_mf_2, e_mf_2)
    g2_src = ROOT.TGraph(len(p_src_2), p_src_2, e_src_2)

    color_mf = ROOT.kAzure - 1
    color_src = ROOT.kRed - 4

    # Increased marker size from 0.4 to 0.6
    for g in [g1_mf, g2_mf]:
        g.SetMarkerStyle(20) 
        g.SetMarkerSize(0.6)
        g.SetMarkerColorAlpha(color_mf, 0.5) 

    for g in [g1_src, g2_src]:
        g.SetMarkerStyle(20)
        g.SetMarkerSize(0.6)
        g.SetMarkerColorAlpha(color_src, 0.5)

    h_dummy1 = ROOT.TH2D("hd1", "; P_{N} [MeV/c]; E_{RMV} [MeV/c]", 150, 0, 600, 150, 0, 250)
    h_dummy2 = ROOT.TH2D("hd2", "; P_{N} [MeV/c]; E_{RMV} [MeV/c]", 150, 0, 600, 150, 0, 250)

    leg_dummy_mf = ROOT.TGraph()
    leg_dummy_mf.SetMarkerStyle(20)
    leg_dummy_mf.SetMarkerColor(color_mf)
    leg_dummy_mf.SetMarkerSize(1.5) 

    leg_dummy_src = ROOT.TGraph()
    leg_dummy_src.SetMarkerStyle(20)
    leg_dummy_src.SetMarkerColor(color_src)
    leg_dummy_src.SetMarkerSize(1.5)

    def create_legend():
        leg = ROOT.TLegend(0.60, 0.75, 0.88, 0.88)
        leg.AddEntry(leg_dummy_mf, "Mean Field", "p")
        leg.AddEntry(leg_dummy_src, "SRC", "p")
        leg.SetTextColor(ROOT.kWhite)
        leg.SetFillColorAlpha(ROOT.kBlack, 0.7) 
        leg.SetLineColor(ROOT.kWhite)
        leg.SetBorderSize(1)
        return leg

    c1 = ROOT.TCanvas("c1", "c1", 1600, 700) 
    c1.Divide(2, 1)

    pad1 = c1.cd(1)
    pad1.SetFillColor(ROOT.kWhite) 
    pad1.SetFrameFillColor(ROOT.kBlack)
    format_axes(h_dummy1)
    
    h_dummy1.Draw()           
    g1_mf.Draw("P SAME")      
    g1_src.Draw("P SAME")

    pad2 = c1.cd(2)
    pad2.SetFillColor(ROOT.kWhite) 
    pad2.SetFrameFillColor(ROOT.kBlack)
    format_axes(h_dummy2)
    
    h_dummy2.Draw()
    g2_mf.Draw("P SAME")
    g2_src.Draw("P SAME")

    leg2 = create_legend()
    leg2.Draw()

    c1.Update()
    c1.SaveAs("spectral_scatter2_{}.pdf".format(filename_chunk))

def transparency(filename, filename_list=[]):

    hist_list = []
    filename_list.append(filename)

    styles = [
        (ROOT.kBlack+3, 3),    # Dotted Black
        (ROOT.kRed, 7),      # Dashed Red
        (ROOT.kBlue, 9),     # Dash-Dot Blue
        (ROOT.kGreen+2, 5),  # Long-Dash Green
        (ROOT.kMagenta, 2),  # Short-Dash Magenta
        (ROOT.kCyan, 1)      # Solid Cyan (fallback)
    ]

    for filename in filename_list:
        filename_chunk = filename.split(".")[0].split("out_")[-1]

        f = ROOT.TFile(filename)
        t = f.Get("neuttree")   

        # Changed X-axis range to 0.0 -> 1.0 (MeV)
        h_trans = ROOT.TH1D("h_trans_{}".format(filename_chunk), "Numerator", 30, 0.0, 500.0)
        h_total = ROOT.TH1D("h_total_{}".format(filename_chunk), "Denominator", 30, 0.0, 500.0)
        h_trans.Sumw2()
        h_total.Sumw2()
        counter = 0
        for event in t:
            counter += 1
            if counter == 10000:
                break

            nvect = event.vectorbranch
            nvect_class = nvect_reader(nvect) 
            if nvect_class.eventType == EventType.MF: 

                fsiProtonKE_MeV = np.sqrt(nvect_class.fsiProton**2 + 938.27**2) - 938.27
                if fsiProtonKE_MeV == 0.0:
                    continue
                
                if nvect_class.intChannel == intChannel_CCQE.qeDeEX or nvect_class.intChannel ==  intChannel_CCQE.noCascadeFSIPhoton or nvect_class.intChannel ==  intChannel_CCQE.noCascadeFSI:
                    h_trans.Fill(fsiProtonKE_MeV)
                    h_total.Fill(fsiProtonKE_MeV)

                else:
                    h_total.Fill(fsiProtonKE_MeV)


        h_ratio = h_trans.Clone("h_ratio_{}".format(filename_chunk))

        h_ratio.SetTitle(";T_{p} [GeV];Transparency p C")
        h_ratio.Divide(h_trans, h_total, 1.0, 1.0, "B")
        

        h_ratio.SetDirectory(0) 
        
        hist_list.append((filename_chunk, h_ratio))
        f.Close()

    ROOT.gStyle.SetOptTitle(0)

    out_file = ROOT.TFile("paper_transparency_combined_output.root", "RECREATE")


    canvas = ROOT.TCanvas("c1", "Transparency Plot", 800, 600)
    canvas.SetTickx(1)
    canvas.SetTicky(1)

    canvas.SetBottomMargin(0.15) 
    canvas.SetLeftMargin(0.12)


    legend = ROOT.TLegend(0.15, 0.18, 0.85, 0.35) 
    legend.SetNColumns(2) 
    legend.SetBorderSize(1)
    legend.SetTextSize(0.04)

    for i, (name, h) in enumerate(hist_list):
        color, style = styles[i % len(styles)]
        
        h.SetLineColor(color)
        h.SetLineStyle(style)
        h.SetLineWidth(2)
        
        h.SetStats(0) 
        h.SetMinimum(0.0) 
        h.SetMaximum(1.0) 

        h.GetXaxis().SetTitleSize(0.06)
        h.GetXaxis().SetTitleOffset(1.0)
        h.GetYaxis().SetTitleSize(0.06)
        h.GetYaxis().SetTitleOffset(0.9)


        if i == 0:
            h.Draw("HIST") 
        else:
            h.Draw("HIST SAME")

        legend.AddEntry(h, name, "l")

    legend.Draw()
    out_file.cd() 
    canvas.Write("transparency_canvas") 

    out_file.Close()


def ccqe_combined_plots(filename, name, no_events = 500000, filename_list=None):
    if filename_list is None:
        filename_list = []
        
    files_to_process = [filename] + filename_list

    print("no events", no_events)
    
    # Original lists
    hist_list_proton = []
    hist_list_mu = []
    hist_list_ex = []
    hist_list_expostfsi = []
    hist_list_avail = []
    hist_list_nop_avail = []
    hist_list_non_incl_avail = [] 
    
    # KDAR Missing Energy
    #hist_list_missing_kdar = []

    # Transverse Kinematics
    hist_list_dpt = []
    hist_list_dat = []

    # Invisible & Total Energy Lists
    hist_list_invisible = []
    hist_list_total_e = []  # ADDED: Total energy list
    
    # Momentum lists
    hist_list_momp_tot  = []
    hist_list_momp_casc = [] 
    hist_list_momp_deex = [] 
    hist_list_momgamma  = []
    hist_list_momn_tot  = []
    hist_list_momn_casc = [] 
    hist_list_momn_deex = [] 
    hist_list_momd = []
    hist_list_momt = []
    hist_list_momhe3 = []
    hist_list_moma_casc = []
    hist_list_moma_deex = []
    hist_list_mompi0 = []
    hist_list_mompipm = []

    # Multiplicity lists
    hist_list_mult_p = []
    hist_list_mult_gamma = []
    hist_list_mult_n = []
    hist_list_mult_d = []
    hist_list_mult_t = []
    hist_list_mult_he3 = []
    hist_list_mult_a = []
    hist_list_mult_pi0 = []
    hist_list_mult_pipm = []

    """
        styles = [
            (ROOT.kBlack, 1),      # Dashed Red
            (ROOT.kBlue, 9),     # Dash-Dot Blue
            (ROOT.kRed, 9),     # Dash-Dot Blue
            (ROOT.kGreen+2, 5),
            (ROOT.kOrange-3,9), # Nieves Color 
            (ROOT.kBlue, 9),     #
    # Long-Dash Green
            (ROOT.kBlue, 9),     # Dash-Dot Blue
            (ROOT.kRed, 9),     # Dash-Dot Blue
            (ROOT.kMagenta, 2),  # Short-Dash Magenta
            (ROOT.kCyan, 1)      # Solid Cyan (fallback)
        ]
    """
    #legend_name = ["No Subprimary Nucleon", "NEUT SRC Definition", "NuWro SRC Definition"]#
    #legend_name = ["No Kinematic Threshold", "14 MeV Threshold", "28 MeV Threshold", "NEUT SRC"]
    #legend_name = ["sf*", "LFG", "sf"]
    #legend_name = ["sf*", "LFG"]
    #legend_name = ["SF","Nieves LFG"] 

    """
    legend_name = ["SF", 
                   "SF*", 
                   ]
    styles = [
                (ROOT.kMagenta, 1),  # Uses index 2 (modified to Petroff Blue)
                (ROOT.kBlue,1), 
            ]
    """

    legend_name = ["PB Hybrid = 1", "PB Hybrid = 0","PB on", "PB off"]
    styles = [ (ROOT.kRed, 1), (ROOT.kBlue, 1), (ROOT.kRed, 2), (ROOT.kBlue, 2)]


    """
    legend_name = ["NEUT Cascade - NuWro SRC", 
                   "NEUT Cascade - Box SRC", 
                   "INCL + ABLA  - NuWro SRC",
                   "INCL + ABLA  - Box SRC",
                   ]
    styles = [
                (ROOT.kRed, 4),  # Uses index 2 (modified to Petroff Blue)
                (ROOT.kBlue, 4),
                (ROOT.kRed, 1),  
                (ROOT.kBlue, 1),
            ]
    """


    """
    legend_name = ["NEUT Cascade", "INCL + ABLA"
                
                ]
    styles = [(ROOT.kAzure,1),
                (ROOT.kBlue, 1),  # Uses index 2 (modified to Petroff Blue)

                #(ROOT.kBlue, 4),  # Modified to Petroff Grape
                #(ROOT.kOrange + 7, 1),  # Modified to Petroff Red
                #(ROOT.kOrange + 7, 4),  # Modified to Petroff Yellow
            ]
    """
    """
    legend_name = ["No SRC Contribution", 
                   "NEUT SRC",
                   "NuWro SRC, no Kin. Cut",
                   "NuWro SRC, 14 MeV Kin. Cut", 
                   "NuWro SRC, iso Kin. Cut", 
                ]
    styles = [
                
                (ROOT.kAzure,1),
                (ROOT.kBlue, 1), 
                (ROOT.kOrange+7,1),
                (ROOT.kRed,1), 
                (ROOT.kMagenta,1), 
            ]
    """

    """
    legend_name = ["INCL + ABLA",
                   "NEUT Cascade",
                   "No Cascade"
                ]

    styles = [
                
                (ROOT.kAzure,1),
                (ROOT.kBlue, 1), 
                (ROOT.kOrange+7,1),
            ]
    """
    """
    legend_name = ["SF*", 
                   "ED-RMF",
                   "Nieves LFG", 
                   ]
    styles = [
                (ROOT.kBlue, 1),  # Uses index 2 (modified to Petroff Blue)
                (ROOT.kRed, 1),  # Modified to Petroff Grape
                (ROOT.kOrange + 7, 1)  # Modified to Petroff Red
    ]
    """
    """
    legend_name = ["tom incl interface", "main"]
    styles = [
                    (ROOT.kBlue, 1),  # Uses index 2 (modified to Petroff Blue)
                    (ROOT.kRed, 1), 
     ] # Modified to Petroff Grape
    """

    count = 0 
    
    for fname in files_to_process:
        print("processing file {}".format(fname))
        filename_chunk = legend_name[count]
        count += 1
        
        f = ROOT.TFile(fname)
        t = f.Get("neuttree")   

# --- Leading Proton & Excitation ---
        h_leadin = ROOT.TH1D(f"h_leading_{filename_chunk}", "Leading Proton Distribution", 50, 0.0, 1500)
        h_leadin.Sumw2()
        h_leadin.SetTitle(";P_{N} [MeV/c]; [Protons / Evt / 30 MeV/c]")

        h_mu = ROOT.TH1D(f"h_leading_{filename_chunk}", "Muon Distribution", 50, 0.0, 1500)
        h_mu.Sumw2()
        h_mu.SetTitle(";P_{#mu} [MeV/c]; [Muons / Evt / 30 MeV/c]")

        h_ex = ROOT.TH1D(f"h_excitation_E_{filename_chunk}", "Excitation Energy Distribution", 200, -20, 80)
        h_ex.Sumw2()
        h_ex.SetTitle(";E^{*}_{init} [MeV]; Probability density [MeV^{-1}]")

        hpostfsi_ex = ROOT.TH1D(f"h_excitation_Epostfsi_{filename_chunk}", "Excitation Energy post fsi Distribution", 200, -20, 80)
        hpostfsi_ex.Sumw2()
        hpostfsi_ex.SetTitle(";E^{*}_{postfsi} [MeV]; Probability density [MeV^{-1}]")
    

        # --- Energy Histograms ---
        h_eavail = ROOT.TH1D(f"h_avail_E_{filename_chunk}", "Hadronic Available Energy Distribution", 40, 0.0, 400)
        h_eavail.Sumw2()
        h_eavail.SetTitle(";E_{Had} [MeV]; [Events / Evt / 10 MeV]")

        h_notpeavail = ROOT.TH1D(f"h_avail_E_noP_{filename_chunk}", "Hadronic Available Energy Distribution (no p)", 100, 0.0, 1000)
        h_notpeavail.Sumw2()
        h_notpeavail.SetTitle(";(E_{Had} - T_{p}) [MeV]; [Events / Evt / 10 MeV]")

        h_non_incl_avail = ROOT.TH1D(f"h_avail_E_nonINCL_{filename_chunk}", "Non-INCL Available Energy Distribution", 100, 0.0, 1000)
        h_non_incl_avail.Sumw2()
        h_non_incl_avail.SetTitle(";E_{Non-INCL} [MeV]; [Events / Evt / 10 MeV]")

        h_invisible = ROOT.TH1D(f"h_invisible_E_{filename_chunk}", "Invisible Energy Distribution", 40, 0.0, 400)
        h_invisible.Sumw2()
        h_invisible.SetTitle(";E_{Inv} [MeV]; [Events / Evt / 10 MeV]")

        h_total_e = ROOT.TH1D(f"h_tot_E_{filename_chunk}", "Total Energy Distribution", 50, 0.0, 1000)
        h_total_e.Sumw2()
        h_total_e.SetTitle(";E_{Total} [MeV]; [Events / Evt / 20 MeV]")

        # --- Transverse Kinematics ---
        h_dpt = ROOT.TH1D(f"h_dpt_{filename_chunk}", "#delta p_{T} Distribution", 50, 0.0, 1000.0)
        h_dpt.Sumw2()
        h_dpt.SetTitle(";#delta p_{T} [MeV/c]; [#delta p_{T} / Evt / 20 MeV/c]")

        h_dat = ROOT.TH1D(f"h_dat_{filename_chunk}", "#delta #alpha_{T} Distribution", 36, 0.0, 180.0)
        h_dat.Sumw2()
        h_dat.SetTitle(";#delta#alpha_{T} [deg]; [#delta#alpha_{T} / Evt / 5 deg]")

        # --- Particle Momentum Distributions ---
        h_momp_tot = ROOT.TH1D(f"h_momp_tot_{filename_chunk}", "Proton Momentum Distribution0", 50, 0.0, 1000)
        h_momp_tot.Sumw2()
        h_momp_tot.SetTitle(";p_{p_{all}} [MeV/c]; [Protons / Evt / 20 MeV/c]")
        
        h_momp_casc = ROOT.TH1D(f"h_momp_casc_{filename_chunk}", "Proton Momentum Distribution1", 50, 0.0, 1000)
        h_momp_casc.Sumw2()
        h_momp_casc.SetTitle(";p_{p_{all}} [MeV/c]; [Protons / Evt / 20 MeV/c]")

        h_momp_deex = ROOT.TH1D(f"h_momp_deex_{filename_chunk}", "Proton Momentum Distribution2", 50, 0.0, 1000)
        h_momp_deex.Sumw2()
        h_momp_deex.SetTitle(";p_{p_{all}} [MeV/c]; [Protons / Evt / 20 MeV/c]")

        h_momgamma = ROOT.TH1D(f"h_momgamma_{filename_chunk}", "Total Photon Momentum Distribution", 100, 0.0, 30)
        h_momgamma.Sumw2()
        h_momgamma.SetTitle(";P_{#gamma}^{Total} [MeV/c]; [Photons / Evt / 0.3 MeV/c]")

        h_momn_tot = ROOT.TH1D(f"h_momn_tot_{filename_chunk}", "neutron Momentum Distribution0", 50, 0.0, 1000)
        h_momn_tot.Sumw2()
        h_momn_tot.SetTitle(";p_{n_{all}} [MeV/c]; [Neutron / Evt / 20 MeV/c]")
        
        h_momn_casc = ROOT.TH1D(f"h_momn_casc_{filename_chunk}", "Neutron Momentum Distribution1", 50, 0.0, 1000)
        h_momn_casc.Sumw2()
        h_momn_casc.SetTitle(";p_{n_{all}} [MeV/c]; [Neutrons / Evt / 20 MeV/c]")

        h_momn_deex = ROOT.TH1D(f"h_momn_deex_{filename_chunk}", "Neutron Momentum Distribution2", 50, 0.0, 1000)
        h_momn_deex.Sumw2()
        h_momn_deex.SetTitle(";p_{n_{all}} [MeV/c]; [Neutrons / Evt / 20 MeV/c]")

        h_momd = ROOT.TH1D(f"h_momd_{filename_chunk}", "Total Deuteron Momentum Distribution", 50, 0.0, 800)
        h_momd.Sumw2()
        h_momd.SetTitle(";P_{d}^{Total} [MeV/c]; [Deuterons / Evt / 8 MeV/c]")

        h_momt = ROOT.TH1D(f"h_momt_{filename_chunk}", "Total Triton Momentum Distribution", 50, 0.0, 800)
        h_momt.Sumw2()
        h_momt.SetTitle(";P_{t}^{Total} [MeV/c]; [Tritons / Evt / 8 MeV/c]")

        h_momhe3 = ROOT.TH1D(f"h_momhe3_{filename_chunk}", "Total Helium-3 Momentum Distribution", 50, 0.0, 800)
        h_momhe3.Sumw2()
        h_momhe3.SetTitle(";P_{^{3}He}^{Total} [MeV/c]; [Helium-3 / Evt / 8 MeV/c]")

        h_moma_casc = ROOT.TH1D(f"h_momacasc_{filename_chunk}", "Total Alpha Momentum Distribution", 50, 0.0, 1000)
        h_moma_casc.Sumw2()
        h_moma_casc.SetTitle(";p_{#alpha} [MeV/c]; [Alphas / Evt / 10 MeV/c]")

        h_moma_deex = ROOT.TH1D(f"h_momadeex_{filename_chunk}", "Total Alpha Momentum Distribution", 50, 0.0, 1000)
        h_moma_deex.Sumw2()
        h_moma_deex.SetTitle(";p_{#alpha} [MeV/c]; [Alphas / Evt / 10 MeV/c]")

        h_mompi0 = ROOT.TH1D(f"h_mompi0_{filename_chunk}", "Total #pi^{0} Momentum Distribution", 100, 0.0, 500)
        h_mompi0.Sumw2()
        h_mompi0.SetTitle(";P_{#pi^{0}}^{Total} [MeV/c]; [#pi^{0} / Evt / 5 MeV/c]")

        h_mompipm = ROOT.TH1D(f"h_mompipm_{filename_chunk}", "Total #pi^{#pm} Momentum Distribution", 100, 0.0, 500)
        h_mompipm.Sumw2()
        h_mompipm.SetTitle(";P_{#pi^{#pm}}^{Total} [MeV/c]; [#pi^{#pm} / Evt / 5 MeV/c]")


        # --- Multiplicity histograms ---
        h_mult_p = ROOT.TH1D(f"h_mult_p_{filename_chunk}", "Proton Multiplicity", 15, 0.0, 15)
        h_mult_p.Sumw2()
        h_mult_p.SetTitle(";N_{p};#frac{N_{events}}{N_{total}}")

        h_mult_gamma = ROOT.TH1D(f"h_mult_gamma_{filename_chunk}", "Photon Multiplicity", 20, 0.0, 20)
        h_mult_gamma.Sumw2()
        h_mult_gamma.SetTitle(";N_{#gamma};#frac{N_{events}}{N_{total}}")

        h_mult_n = ROOT.TH1D(f"h_mult_n_{filename_chunk}", "Neutron Multiplicity", 15, 0.0, 15)
        h_mult_n.Sumw2()
        h_mult_n.SetTitle(";N_{n};#frac{N_{events}}{N_{total}}")

        h_mult_d = ROOT.TH1D(f"h_mult_d_{filename_chunk}", "Deuteron Multiplicity", 10, 0.0, 10)
        h_mult_d.Sumw2()
        h_mult_d.SetTitle(";N_{d};#frac{N_{events}}{N_{total}}")

        h_mult_t = ROOT.TH1D(f"h_mult_t_{filename_chunk}", "Triton Multiplicity", 10, 0.0, 10)
        h_mult_t.Sumw2()
        h_mult_t.SetTitle(";N_{t};#frac{N_{events}}{N_{total}}")

        h_mult_he3 = ROOT.TH1D(f"h_mult_he3_{filename_chunk}", "Helium-3 Multiplicity", 10, 0.0, 10)
        h_mult_he3.Sumw2()
        h_mult_he3.SetTitle(";N_{^{3}He};#frac{N_{events}}{N_{total}}")

        h_mult_a = ROOT.TH1D(f"h_mult_a_{filename_chunk}", "Alpha Multiplicity", 10, 0.0, 10)
        h_mult_a.Sumw2()
        h_mult_a.SetTitle(";N_{#alpha};#frac{N_{events}}{N_{total}}")

        h_mult_pi0 = ROOT.TH1D(f"h_mult_pi0_{filename_chunk}", "#pi^{0} Multiplicity", 10, 0.0, 10)
        h_mult_pi0.Sumw2()
        h_mult_pi0.SetTitle(";N_{#pi^{0}};#frac{N_{events}}{N_{total}}")

        h_mult_pipm = ROOT.TH1D(f"h_mult_pipm_{filename_chunk}", "#pi^{#pm} Multiplicity", 10, 0.0, 10)
        h_mult_pipm.Sumw2()
        h_mult_pipm.SetTitle(";N_{#pi^{#pm}};#frac{N_{events}}{N_{total}}")

        event_no = 0
        EVENT_COUNT = 0
        SRC_COUNT = 0
        for event in t:
            nvect = event.vectorbranch
            nvect_class = nvect_reader(nvect) 
            event_no +=1
            if event_no % 100000 == 0:
                print("event number ", event_no)
            if event_no == no_events:
                break
            
            if nvect_class.eventType == EventType.SRC or nvect_class.eventType == EventType.MF: 
               
                dpt, dat = nvect_class.get_deltaPT()
                
                h_dpt.Fill(dpt)
                h_dat.Fill(dat)
                postFSI_ex = nvect_class.post_casc_E_star()
                if not (postFSI_ex ==  None): 
                
                    hpostfsi_ex.Fill(nvect_class.post_casc_E_star())

                if nvect_class.intChannel == intChannel_CCQE.other:
                    print(nvect_class.eventType.name)

                if nvect_class.eventType == EventType.MF:
                    EVENT_COUNT += 1
                #    if nvect_class.excitation_E_CCQE(carbon=True) < 0.0:
                #        h_ex.Fill(0.0)
                #    else:
                    h_ex.Fill(nvect_class.excitation_E_CCQE(carbon=True))
                elif nvect_class.eventType == EventType.SRC:
                    EVENT_COUNT += 1
                    SRC_COUNT += 1
                    #if nvect_class.excitation_E_SRC(carbon=True) < 0.0:
                    #    h_ex.Fill(0.0)
                    #else:
                    h_ex.Fill(nvect_class.excitation_E_SRC(carbon=True))
               

               # Missing_energy_kdar, E_Miss = nvect_class.Missing_energy_kdar2()
               # if Missing_energy_kdar:
               #     h_missing_kdar.Fill(E_Miss)

                particles, energies, flags = nvect_class.particles(flag=True,lepton=True)

                eavail = 0
                for particle, energy in zip(particles, energies):
                    if particle == 2112 or particle == 13:
                        continue
                    eavail += energy

                h_eavail.Fill(eavail) 

                # Counters for per-event totals
                non_proton_evail = 0
                muon_mom = 0
                non_INCL_eavail = 0 
                invisible_energy = 0
                total_E = 0  # ADDED: Initialized per event
                total_proton_tot_mom = []
                total_proton_cascade_mom = []
                total_proton_deex_mom = []
                total_gamma_mom = 0
                total_neutron_tot_mom  = []
                total_neutron_casc_mom = []
                total_neutron_deex_mom = []
                total_deuteron_mom = 0
                total_triton_mom = 0
                total_he3_mom = 0
                total_alpha_cascade_mom = []
                total_alpha_deex_mom = []
                total_pi0_mom = 0
                total_pipm_mom = 0

                mult_p = 0
                mult_gamma = 0
                mult_n = 0
                mult_d = 0
                mult_t = 0
                mult_he3 = 0
                mult_a = 0
                mult_pi0 = 0
                mult_pipm = 0
                
                for particle, energy,flag in zip(particles, energies,flags):
                    if  particle == 13:
                        muon_mom = np.sqrt((energy+139.27)**2 - 139.27)
                    total_E += energy
                    if particle != 2212 and particle != 2112:
                        non_proton_evail += energy

                    if (particle < 1000000 and particle > 20) and particle != 2112:
                        #print(particle, energy)
                        non_INCL_eavail += energy

                    if particle != 2212 and abs(particle) != 211:
                        invisible_energy += energy
                    
                    # Convert Kinetic Energy to Momentum (MeV/c)
                    if particle == 2212: # Proton
                        total_proton_tot_mom.append(np.sqrt((energy+938.27)**2 - 938.27**2))
                        if flag == 8:
                            total_proton_cascade_mom.append(np.sqrt((energy+938.27)**2 - 938.27**2))
                        elif flag == 10:
                            total_proton_deex_mom.append(np.sqrt((energy+938.27)**2 - 938.27**2))
                        #total_proton_mom += energy
                        mult_p += 1
                    elif particle == 22: # Gamma
                        total_gamma_mom += energy # p = E for photon
                        mult_gamma += 1

                    elif particle == 2112: # Neutron
                        total_neutron_tot_mom.append(np.sqrt((energy+939.57)**2 - 939.57**2))
                        if flag == 8 or flag == 0:
                            total_neutron_casc_mom.append(np.sqrt((energy+939.57)**2 - 939.57**2))
                        elif flag == 10:
                            total_neutron_deex_mom.append(np.sqrt((energy+939.57)**2 - 939.57**2))
                        mult_n += 1
                    elif particle == 1000010020: # Deuteron
                        if flag == 10:
                            total_deuteron_mom += np.sqrt((energy+1875.61)**2 - 1875.61**2)
                        mult_d += 1
                    elif particle == 1000010030: # Triton
                        if flag == 10:
                            total_triton_mom += np.sqrt((energy+2808.92)**2 - 2808.92**2)
                        mult_t += 1
                    elif particle == 1000020030: # Helium-3
                        if flag == 10:
                            total_he3_mom += np.sqrt((energy+2808.39)**2 - 2808.39**2)
                        mult_he3 += 1
                    elif particle == 1000020040: # Alpha
                        if flag == 8:
                            total_alpha_cascade_mom.append(np.sqrt((energy+3727.38)**2 - 3727.38**2))
                        if flag == 10:
                            total_alpha_deex_mom.append( np.sqrt((energy+3727.38)**2 - 3727.38**2))
                        mult_a += 1
                    elif particle == 111: # Pi0
                        total_pi0_mom += np.sqrt((energy+134.98)**2 - 134.98**2)
                        mult_pi0 += 1
                    elif abs(particle) == 211: # Pi+-
                        total_pipm_mom += np.sqrt((energy+139.57)**2 - 139.57**2)
                        mult_pipm += 1
                
                # Fill energy histograms (skipping zeros)
                if non_proton_evail > 0: h_notpeavail.Fill(non_proton_evail)
                if non_INCL_eavail > 0: h_non_incl_avail.Fill(non_INCL_eavail)
                if invisible_energy > 0: h_invisible.Fill(invisible_energy)
                if total_E > 0: h_total_e.Fill(total_E)  # ADDED: Fill total energy

                # Fill momentum histograms
                if len(total_proton_tot_mom) > 0: 
                    for proton in total_proton_tot_mom:
                        h_momp_tot.Fill(proton)

                if len(total_proton_cascade_mom) > 0: 
                    for proton in total_proton_cascade_mom:
                        h_momp_casc.Fill(proton)

                if len(total_proton_deex_mom) > 0: 
                    for proton in total_proton_deex_mom:
                        h_momp_deex.Fill(proton)

                if len(total_neutron_tot_mom) > 0: 
                    for neutron in total_neutron_tot_mom:
                        h_momn_tot.Fill(neutron)

                if len(total_neutron_casc_mom) > 0: 
                    for neutron in total_neutron_casc_mom:
                        h_momn_casc.Fill(neutron)

                if len(total_neutron_deex_mom) > 0: 
                    for neutron in total_neutron_deex_mom:
                        h_momn_deex.Fill(neutron)

                if muon_mom > 0: h_mu.Fill(muon_mom)

                if total_gamma_mom > 0: h_momgamma.Fill(total_gamma_mom)

                if total_deuteron_mom > 0: h_momd.Fill(total_deuteron_mom)
                if total_triton_mom > 0: h_momt.Fill(total_triton_mom)
                if total_he3_mom > 0: h_momhe3.Fill(total_he3_mom)

                if len(total_alpha_cascade_mom) > 0:
                    for alpha in total_alpha_cascade_mom:
                        h_moma_casc.Fill(alpha)
                if len(total_alpha_deex_mom) > 0:
                    for alpha in total_alpha_deex_mom:
                        h_moma_deex.Fill(alpha)
                if total_pi0_mom > 0: h_mompi0.Fill(total_pi0_mom)
                if total_pipm_mom > 0: h_mompipm.Fill(total_pipm_mom)

                # Fill multiplicity histograms
                h_mult_p.Fill(mult_p)
                h_mult_gamma.Fill(mult_gamma)
                h_mult_n.Fill(mult_n)
                h_mult_d.Fill(mult_d)
                h_mult_t.Fill(mult_t)
                h_mult_he3.Fill(mult_he3)
                h_mult_a.Fill(mult_a)
                h_mult_pi0.Fill(mult_pi0)
                h_mult_pipm.Fill(mult_pipm)

                if nvect_class.HMPMom:
                    h_leadin.Fill(nvect_class.HMPMom)

        print(SRC_COUNT/EVENT_COUNT)
        # Set directories
        h_leadin.SetDirectory(0) 
        h_mu.SetDirectory(0)
        h_ex.SetDirectory(0) 
        hpostfsi_ex.SetDirectory(0)
        h_eavail.SetDirectory(0) 
        h_notpeavail.SetDirectory(0) 
        h_non_incl_avail.SetDirectory(0)
        h_invisible.SetDirectory(0)
        h_total_e.SetDirectory(0)  # ADDED: Detach from file directory
       # h_missing_kdar.SetDirectory(0)

        h_dpt.SetDirectory(0)
        h_dat.SetDirectory(0)
        h_momp_tot.SetDirectory(0)
        h_momp_casc.SetDirectory(0)
        h_momp_deex.SetDirectory(0)
        h_momgamma.SetDirectory(0)
        h_momn_tot.SetDirectory(0)
        h_momn_casc.SetDirectory(0)
        h_momn_deex.SetDirectory(0)
        h_momd.SetDirectory(0)
        h_momt.SetDirectory(0)
        h_momhe3.SetDirectory(0)
        h_moma_casc.SetDirectory(0)
        h_moma_deex.SetDirectory(0)
        h_mompi0.SetDirectory(0)
        h_mompipm.SetDirectory(0)

        h_mult_p.SetDirectory(0)
        h_mult_gamma.SetDirectory(0)
        h_mult_n.SetDirectory(0)
        h_mult_d.SetDirectory(0)
        h_mult_t.SetDirectory(0)
        h_mult_he3.SetDirectory(0)
        h_mult_a.SetDirectory(0)
        h_mult_pi0.SetDirectory(0)
        h_mult_pipm.SetDirectory(0)
        
        # Append to lists
        hist_list_proton.append((filename_chunk, h_leadin))
        hist_list_mu.append((filename_chunk,h_mu))
        hist_list_ex.append((filename_chunk, h_ex))
        hist_list_expostfsi.append((filename_chunk,hpostfsi_ex))
        hist_list_avail.append((filename_chunk, h_eavail))
        hist_list_nop_avail.append((filename_chunk, h_notpeavail))
        hist_list_non_incl_avail.append((filename_chunk, h_non_incl_avail))
        hist_list_invisible.append((filename_chunk, h_invisible))
        hist_list_total_e.append((filename_chunk, h_total_e))  # ADDED: Append to list
       # hist_list_missing_kdar.append((filename_chunk, h_missing_kdar))

        hist_list_dpt.append((filename_chunk, h_dpt))
        hist_list_dat.append((filename_chunk, h_dat))

        hist_list_momp_tot.append((filename_chunk,   h_momp_tot))        
        hist_list_momp_casc.append((filename_chunk, h_momp_casc))
        hist_list_momp_deex.append((filename_chunk, h_momp_deex))

        hist_list_momgamma.append((filename_chunk, h_momgamma))

        hist_list_momn_tot.append((filename_chunk, h_momn_tot))
        hist_list_momn_casc.append((filename_chunk, h_momn_casc))
        hist_list_momn_deex.append((filename_chunk, h_momn_deex))
        hist_list_momd.append((filename_chunk, h_momd))
        hist_list_momt.append((filename_chunk, h_momt))
        hist_list_momhe3.append((filename_chunk, h_momhe3))
        hist_list_moma_casc.append((filename_chunk, h_moma_casc))
        hist_list_moma_deex.append((filename_chunk,h_moma_deex))
        hist_list_mompi0.append((filename_chunk, h_mompi0))
        hist_list_mompipm.append((filename_chunk, h_mompipm))

        hist_list_mult_p.append((filename_chunk, h_mult_p))
        hist_list_mult_gamma.append((filename_chunk, h_mult_gamma))
        hist_list_mult_n.append((filename_chunk, h_mult_n))
        hist_list_mult_d.append((filename_chunk, h_mult_d))
        hist_list_mult_t.append((filename_chunk, h_mult_t))
        hist_list_mult_he3.append((filename_chunk, h_mult_he3))
        hist_list_mult_a.append((filename_chunk, h_mult_a))
        hist_list_mult_pi0.append((filename_chunk, h_mult_pi0))
        hist_list_mult_pipm.append((filename_chunk, h_mult_pipm))

        f.Close()

    ROOT.gStyle.SetOptTitle(0)
    out_file = ROOT.TFile("paper_combined_comparisons_{}.root".format(name), "RECREATE")
    latex_font = 132

    # Canvas helper definitions remain standard...
    def draw_canvas(canvas_name, hist_list, ymax=None, ymin=None, log = False):
        c = ROOT.TCanvas(canvas_name, canvas_name, 650, 600)
        c.SetGrayscale(False)
        c.SetLeftMargin(0.17)   # Default is ~0.12
        c.SetRightMargin(0.08)  # Default is ~0.10
        c.SetTopMargin(0.06)    # Default is ~0.10
        c.SetBottomMargin(0.12) # Default is ~0.12
        c.SetLogy(log)

        ROOT.gPad.SetTickx(1)
        ROOT.gPad.SetTicky(1)
        #ROOT.gPad.SetBottomMargin(0.15) 
        #ROOT.gPad.SetLeftMargin(0.12)
        

        legend = ROOT.TLegend(0.52, 0.67, 0.86, 0.90)
        legend.SetNColumns(1) 
        legend.SetTextFont(latex_font) 
        legend.SetBorderSize(0) 
        legend.SetFillStyle(1001)  # Solid fill
        legend.SetFillColor(ROOT.kWhite)
        legend.SetTextSize(0.03) 
        

        valid_hists = []
        for name, h in hist_list:
            #print(name,h)
            #print(event_no)
            #print(h.Integral())
            h.Scale(1/ event_no)
            pass
            #print(h.Integral())

            valid_hists.append((name, h))

        if ymax is None and valid_hists:
            global_max = max(h.GetMaximum() for _, h in valid_hists)
            ymax = global_max * 1.2

        for i, (name, h) in enumerate(valid_hists):
            h_name = h.GetName()
            if ("dat" in h_name or "dpt" in h_name or "excitation_E" in h_name) and "no cascade" in name:
                continue

            
            color, style = styles[i]
            
            h.SetLineColor(color)
            h.SetLineStyle(style)
            h.SetLineWidth(4) 
            h.SetStats(0) 

            h.GetXaxis().CenterTitle(True)
            h.GetYaxis().CenterTitle(True)

            h.GetYaxis().SetMaxDigits(4)

            # Force X-axis to NEVER use scientific exponent notation
            h.GetXaxis().SetNoExponent(True)

            h.GetXaxis().SetTitleFont(latex_font)
            h.GetXaxis().SetLabelFont(latex_font)
            h.GetYaxis().SetTitleFont(latex_font)
            h.GetYaxis().SetLabelFont(latex_font)
            h.GetXaxis().SetTitleSize(0.05) 
            h.GetXaxis().SetLabelSize(0.05) 
            h.GetXaxis().SetTitleOffset(1.0)

            h.GetXaxis().SetNdivisions(505, True)

            h.GetYaxis().SetTitleSize(0.05) 
            h.GetYaxis().SetLabelSize(0.05) 
            h.GetYaxis().SetTitleOffset(1.6)

            if i == 0:
                h.SetMaximum(ymax)
                if ymin is not None:
                    h.SetMinimum(ymin)
                h.Draw("HIST") 
            else:
                h.Draw("HIST SAME")
                
            legend.AddEntry(h, name, "l")

        legend.Draw()
        return c, legend 

    def draw_canvas_casc_deex(canvas_name, hist_list_casc, hist_list_deex, ymax=None, ymin=None):


        c = ROOT.TCanvas(canvas_name, canvas_name, 650, 600)
        c.SetGrayscale(False)
        c.SetLeftMargin(0.17)   # Default is ~0.12
        c.SetRightMargin(0.08)  # Default is ~0.10
        c.SetTopMargin(0.06)    # Default is ~0.10
        c.SetBottomMargin(0.12) # Default is ~0.12

        ROOT.gPad.SetTickx(1)
        ROOT.gPad.SetTicky(1)

        legend = ROOT.TLegend(0.22, 0.67, 0.56, 0.90)
        legend.SetNColumns(1) 
        legend.SetTextFont(latex_font) 
        legend.SetBorderSize(0) 
        legend.SetFillStyle(1001)  # Solid fill
        legend.SetFillColor(ROOT.kWhite)
        legend.SetTextSize(0.03) 

        # Process and scale histograms for both lists
        def process_hists(hist_list):
            valid = []
            for name, h in hist_list:
                #print(name, h)
                #print("Entries:", h.GetEntries())  # Prints raw entry count
                h.Scale(1 / event_no)
                pass
                #print(h.Integral())
                valid.append((name, h))
            return valid
        
        valid_casc = process_hists(hist_list_casc)
        print("valid casc ",valid_casc)
        valid_deex = process_hists(hist_list_deex)

        # Determine global maximum across both lists if ymax is not explicitly provided
        all_valid = valid_casc + valid_deex
  
        if ymax is None and all_valid:
            global_max = max(h.GetMaximum() for _, h in all_valid)
            ymax = global_max * 1.2

        def apply_hist_style(h, color, style):
            h.SetLineColor(color)
            h.SetLineStyle(style)
            h.SetLineWidth(4) 
            h.SetStats(0) 

            h.GetXaxis().CenterTitle(True)
            h.GetYaxis().CenterTitle(True)
            h.GetYaxis().SetMaxDigits(4)
            h.GetXaxis().SetNoExponent(True)

            h.GetXaxis().SetTitleFont(latex_font)
            h.GetXaxis().SetLabelFont(latex_font)
            h.GetYaxis().SetTitleFont(latex_font)
            h.GetYaxis().SetLabelFont(latex_font)

            h.GetXaxis().SetTitleSize(0.05) 
            h.GetXaxis().SetLabelSize(0.05) 
            h.GetXaxis().SetTitleOffset(1.0)

            h.GetYaxis().SetTitleSize(0.05) 
            h.GetYaxis().SetLabelSize(0.05) 
            h.GetYaxis().SetTitleOffset(1.6)

        first_draw = True
        num_samples = max(len(valid_deex), len(valid_casc))

        for i in range(num_samples):
            color, _ = styles[i]

            # 1. Plot De-excitation histogram (Solid Line)
            if i < len(valid_deex):
                name_deex, h_deex = valid_deex[i]
                h_name = h_deex.GetName()

                if not (("dat" in h_name or "dpt" in h_name or "excitation_E" in h_name) and "no cascade" in name_deex):
                    apply_hist_style(h_deex, color=color, style=1)  # 1 = Solid

                    draw_opt = "HIST" if first_draw else "HIST SAME"
                    if first_draw:
                        h_deex.SetMaximum(ymax)
                        if ymin is not None:
                            h_deex.SetMinimum(ymin)
                        first_draw = False

                    h_deex.Draw(draw_opt)
                    label = name_deex if "deex" in name_deex.lower() else f"{name_deex} (de-excitation)"
                    legend.AddEntry(h_deex, label, "l")

            # 2. Plot Cascade histogram (Dashed Line, Same Color)
            if i < len(valid_casc):
                name_casc, h_casc = valid_casc[i]
                h_name = h_casc.GetName()

                if not (("dat" in h_name or "dpt" in h_name or "excitation_E" in h_name) and "no cascade" in name_casc):
                    apply_hist_style(h_casc, color=color, style=2)  # 2 = Dashed

                    draw_opt = "HIST" if first_draw else "HIST SAME"
                    if first_draw:
                        h_casc.SetMaximum(ymax)
                        if ymin is not None:
                            h_casc.SetMinimum(ymin)
                        first_draw = False

                    h_casc.Draw(draw_opt)
                    label = name_casc if "casc" in name_casc.lower() else f"{name_casc} (cascade)"
                    legend.AddEntry(h_casc, label, "l")

        legend.Draw()
        return c, legend

    def draw_canvas_with_ratio(canvas_name, hist_list, ymax=None, ymin=None, ratio_min = 0.5, ratio_max = 1.5, log=False):
        c = ROOT.TCanvas(canvas_name, canvas_name, 600, 750)
        c.SetGrayscale(False)
        c.SetLeftMargin(0.12)   # Default is ~0.12
        c.SetRightMargin(0.08)  # Default is ~0.10
        c.SetTopMargin(0.06)    # Default is ~0.10
        c.SetBottomMargin(0.0) # Default is ~0.12

        pad1 = ROOT.TPad("pad1", "pad1", 0.0, 0.4, 1.0, 0.98)
        pad1.SetBottomMargin(0.02)
        pad1.SetLeftMargin(0.12)   # Default is ~0.12
        pad1.SetRightMargin(0.08)  # Default is ~0.10
        pad1.SetTopMargin(0.06)    # Default is ~0.10
        pad1.SetBottomMargin(0.02) # Default is ~0.12

        pad1.SetTickx(1)
        pad1.SetTicky(1)
        pad1.Draw()
        pad1.cd()
        if log == True:
            pad1.SetLogy(True)

        legend = ROOT.TLegend(0.47, 0.67, 0.8, 0.9)
        #legend = ROOT.TLegend(0.45, 0.50, 0.8, 0.80)
        legend.SetNColumns(1) 
        legend.SetTextFont(latex_font) 
        legend.SetBorderSize(0) 
        legend.SetFillStyle(1001)  # Solid fill
        legend.SetFillColor(ROOT.kWhite)
        legend.SetTextSize(0.04) 

        valid_hists = []
        for name, h in hist_list:
            if ymax is None:
                integral = h.Integral()
                if integral > 0:
                    h.Scale(1/ no_events)
                    pass

            valid_hists.append((name, h))

        if ymax is None and valid_hists:
            global_max = max(h.GetMaximum() for _, h in valid_hists)
            ymax = global_max * 1.2


        xaxis_title = hist_list[0][1].GetXaxis().GetTitle()

        for i, (name, h) in enumerate(valid_hists):

            
            color, style = styles[i % len(styles)]
            #print(color)
            h.SetLineColor(color)
            h.SetLineStyle(style)
            h.SetLineWidth(3) 
            h.SetStats(0) 

            h.GetXaxis().SetLabelSize(0)
            h.GetXaxis().SetTitleSize(0)
            
            h.GetYaxis().SetTitleFont(latex_font)
            h.GetYaxis().SetLabelFont(latex_font)
            h.GetYaxis().SetTitleOffset(1.2)

            h.GetYaxis().SetTitleSize(0.05) 
            h.GetYaxis().SetLabelSize(0.05) 

            h.GetXaxis().SetNdivisions(505, True)

            h.GetXaxis().CenterTitle(True)
            h.GetYaxis().CenterTitle(True)

            h.GetYaxis().SetMaxDigits(3)
            h.GetXaxis().SetNoExponent(True)
            

            print(
                f"Hist: {h.GetName()} | Set Color: {color} | Real LineColor: {h.GetLineColor()}"
            )

            if i == 0:
                h.SetMaximum(ymax)
                if ymin is not None:
                    h.SetMinimum(ymin)
                h.Draw("HIST") 
            else:
                h.Draw("HIST SAME")
                
            legend.AddEntry(h, name, "l")

        legend.Draw()

        c.cd()
        pad2 = ROOT.TPad("pad2", "pad2", 0.0, 0.15, 1.0, 0.4)
        pad2.SetLeftMargin(0.12)   # Default is ~0.12
        pad2.SetRightMargin(0.08)  # Default is ~0.10
        pad2.SetTopMargin(0.06)    # Default is ~0.10
        pad2.SetBottomMargin(0.40) # Default is ~0.12
        pad2.SetTickx(1)
        pad2.SetTicky(1)
        pad2.Draw()
        pad2.cd()

        h_base = hist_list[0][1]
        ratios = []
        for i in range(1, len(hist_list)):
            name, h = hist_list[i]
            h_ratio = h.Clone(f"{h.GetName()}_ratio")
            h_ratio.Divide(h_base)

            color, style = styles[i % len(styles)]
            h_ratio.SetLineColor(color)
            h_ratio.SetLineStyle(style)
            h_ratio.SetLineWidth(3)
            h_ratio.SetStats(0)

            h_ratio.SetTitle("")
            h_ratio.GetYaxis().SetTitle("NEUT/INCL")
            h_ratio.GetYaxis().CenterTitle()
            h_ratio.GetYaxis().SetNdivisions(505)
            h_ratio.GetYaxis().SetTitleSize(0.11)
            h_ratio.GetYaxis().SetTitleFont(latex_font)
            h_ratio.GetYaxis().SetTitleOffset(0.5)
            h_ratio.GetYaxis().SetLabelFont(latex_font)
            h_ratio.GetYaxis().SetLabelSize(0.1)
            h_ratio.GetXaxis().SetNdivisions(505, True)
            
            h_ratio.GetXaxis().SetTitle(xaxis_title)
            h_ratio.GetXaxis().SetTitleSize(0.13)
            h_ratio.GetXaxis().SetTitleFont(latex_font)
            h_ratio.GetXaxis().SetTitleOffset(1.0)
            h_ratio.GetXaxis().SetLabelFont(latex_font)
            h_ratio.GetXaxis().SetLabelSize(0.12)



            h_ratio.GetXaxis().CenterTitle(True)
            h_ratio.GetYaxis().CenterTitle(True)

            if i == 1:
                h_ratio.SetMinimum(ratio_min)
                h_ratio.SetMaximum(ratio_max)
                h_ratio.Draw("HIST")
            else:
                h_ratio.Draw("HIST SAME")
            ratios.append(h_ratio)
            
        line = ROOT.TLine(h_base.GetXaxis().GetXmin(), 1, h_base.GetXaxis().GetXmax(), 1)
        line.SetLineColor(ROOT.kBlack)
        line.SetLineStyle(2)
        line.Draw("SAME")

        c.pad1 = pad1
        c.pad2 = pad2
        c.ratios = ratios
        c.line = line

        pad1.Modified()
        pad1.Update()
        pad2.Modified()
        pad2.Update()
        
        c.cd()
        c.Modified()
        c.Update()
        
        return c, legend

    def draw_canvas_with_ratio2(canvas_name, hist_list, ymax=None, ymin=None, ratio_min=0.5, ratio_max=1.5):
        c = ROOT.TCanvas(canvas_name, canvas_name, 600, 750)
        c.SetGrayscale(False)
        c.SetLeftMargin(0.12)
        c.SetRightMargin(0.08)
        c.SetTopMargin(0.06)
        c.SetBottomMargin(0.0)

        # Upper Pad
        pad1 = ROOT.TPad("pad1", "pad1", 0.0, 0.4, 1.0, 0.98)
        pad1.SetLeftMargin(0.12)
        pad1.SetRightMargin(0.08)
        pad1.SetTopMargin(0.06)
        pad1.SetBottomMargin(0.02)
        pad1.SetTickx(1)
        pad1.SetTicky(1)
        pad1.Draw()
        pad1.cd()

        legend = ROOT.TLegend(0.45, 0.50, 0.8, 0.80)
        legend.SetNColumns(1) 
        legend.SetTextFont(latex_font) 
        legend.SetBorderSize(0) 
        legend.SetFillStyle(1001)
        legend.SetFillColor(ROOT.kWhite)
        legend.SetTextSize(0.04) 

        if ymax is None:
            global_max = max(h.GetMaximum() for _, h in hist_list)
            ymax = global_max * 1.2

        xaxis_title = hist_list[0][1].GetXaxis().GetTitle()

        # Draw all histograms on top pad
        for i, (name, h) in enumerate(hist_list):
            color, style = styles[i % len(styles)]

            integral = h.Integral()
            if integral > 0:
                h.Scale(1.0 / integral)
            h.SetLineColor(color)
            h.SetLineStyle(style)
            h.SetLineWidth(3) 
            h.SetStats(0) 

            h.GetXaxis().SetLabelSize(0)
            h.GetXaxis().SetTitleSize(0)
            
            h.GetYaxis().SetTitleFont(latex_font)
            h.GetYaxis().SetLabelFont(latex_font)
            h.GetYaxis().SetTitleOffset(1.2)
            h.GetYaxis().SetTitleSize(0.05) 
            h.GetYaxis().SetLabelSize(0.05) 

            h.GetXaxis().CenterTitle(True)
            h.GetYaxis().CenterTitle(True)
            h.GetYaxis().SetMaxDigits(3)
            h.GetXaxis().SetNoExponent(True)

            if i == 0:
                h.SetMaximum(ymax)
                if ymin is not None:
                    h.SetMinimum(ymin)
                h.Draw("HIST") 
            else:
                h.Draw("HIST SAME")
                
            legend.AddEntry(h, name, "l")

        legend.Draw()

        # Lower Pad
        c.cd()
        pad2 = ROOT.TPad("pad2", "pad2", 0.0, 0.15, 1.0, 0.4)
        pad2.SetLeftMargin(0.12)
        pad2.SetRightMargin(0.08)
        pad2.SetTopMargin(0.06)
        pad2.SetBottomMargin(0.40)
        pad2.SetTickx(1)
        pad2.SetTicky(1)
        pad2.Draw()
        pad2.cd()

        ratio_pairs = [(1, 0), (3, 2)]


        ratios = []
        for idx, (num_idx, den_idx) in enumerate(ratio_pairs):
            name_num, h_num = hist_list[num_idx]
            name_den, h_den = hist_list[den_idx]

            h_ratio = h_num.Clone(f"{h_num.GetName()}_ratio_{num_idx}_{den_idx}")
            h_ratio.Divide(h_den)

            color, style = styles[num_idx % len(styles)]
            h_ratio.SetLineColor(color)
            h_ratio.SetLineStyle(style)
            h_ratio.SetLineWidth(3)
            h_ratio.SetStats(0)

            h_ratio.SetTitle("")
            h_ratio.GetYaxis().SetTitle("#frac{NEUT SRC}{NuWro SRC}")
            h_ratio.GetYaxis().CenterTitle()
            h_ratio.GetYaxis().SetNdivisions(505)
            h_ratio.GetYaxis().SetTitleSize(0.11)
            h_ratio.GetYaxis().SetTitleFont(latex_font)
            h_ratio.GetYaxis().SetTitleOffset(0.5)
            h_ratio.GetYaxis().SetLabelFont(latex_font)
            h_ratio.GetYaxis().SetLabelSize(0.1)
            
            h_ratio.GetXaxis().SetTitle(xaxis_title)
            h_ratio.GetXaxis().SetTitleSize(0.13)
            h_ratio.GetXaxis().SetTitleFont(latex_font)
            h_ratio.GetXaxis().SetTitleOffset(1.0)
            h_ratio.GetXaxis().SetLabelFont(latex_font)
            h_ratio.GetXaxis().SetLabelSize(0.12)

            h_ratio.GetXaxis().CenterTitle(True)
            h_ratio.GetYaxis().CenterTitle(True)

            if idx == 0:
                h_ratio.SetMinimum(ratio_min)
                h_ratio.SetMaximum(ratio_max)
                h_ratio.Draw("HIST")
            else:
                h_ratio.Draw("HIST SAME")

            ratios.append(h_ratio)


        h_base = hist_list[0][1]
        line = ROOT.TLine(h_base.GetXaxis().GetXmin(), 1, h_base.GetXaxis().GetXmax(), 1)
        line.SetLineColor(ROOT.kBlack)
        line.SetLineStyle(2)
        line.Draw("SAME")

        c.pad1 = pad1
        c.pad2 = pad2
        c.ratios = ratios
        c.line = line

        pad1.Modified()
        pad1.Update()
        pad2.Modified()
        pad2.Update()
        
        c.cd()
        c.Modified()
        c.Update()
        
        return c, legend


    def draw_side_by_side_canvas(canvas_name, hist_list_energy, hist_list_mult, ymax_left=None, ymax_right=None):
        c = ROOT.TCanvas(canvas_name, canvas_name, 1400, 750)
        c.Divide(2, 1)

        def configure_and_draw(pad_idx, hist_list, ymax=None):
            c.cd(pad_idx)
            ROOT.gPad.SetTickx(1)
            ROOT.gPad.SetTicky(1)
            ROOT.gPad.SetBottomMargin(0.15) 
            ROOT.gPad.SetLeftMargin(0.12)

            legend = ROOT.TLegend(0.35, 0.7, 0.60, 0.80) 
            legend.SetNColumns(1) 
            legend.SetTextFont(latex_font) 
            legend.SetBorderSize(0) 
            legend.SetFillStyle(0)  
            legend.SetTextSize(0.04) 

            if ymax is None:
                global_max = max(h.GetMaximum() for _, h in hist_list)
                ymax = global_max * 1.2

            for i, (name, h) in enumerate(hist_list):
                color, style = styles[i % len(styles)]
                
                h.SetLineColor(color)
                h.SetLineStyle(style)
                h.SetLineWidth(2) 
                h.SetStats(0) 

                h.GetXaxis().SetTitleFont(latex_font)
                h.GetXaxis().SetLabelFont(latex_font)
                h.GetYaxis().SetTitleFont(latex_font)
                h.GetYaxis().SetLabelFont(latex_font)
                h.GetXaxis().SetTitleSize(0.05) 
                h.GetXaxis().SetLabelSize(0.05) 
                h.GetXaxis().SetTitleOffset(1.0)
                h.GetYaxis().SetTitleOffset(1.1)
                h.GetYaxis().SetTitleSize(0.05) 
                h.GetYaxis().SetLabelSize(0.05) 

                h.GetXaxis().CenterTitle(True)
                h.GetYaxis().CenterTitle(True)

                h.GetYaxis().SetMaxDigits(3)
                h.GetXaxis().SetNoExponent(True)

                if i == 0:
                    h.SetMaximum(ymax)
                    h.Draw("HIST") 
                else:
                    h.Draw("HIST SAME")
                    
                legend.AddEntry(h, name, "l")
            legend.Draw()
            return legend

        leg_e = configure_and_draw(1, hist_list_energy, ymax_left)
        leg_m = configure_and_draw(2, hist_list_mult, ymax_right)
        
        c.legends = (leg_e, leg_m)
        return c, c.legends
    
    # Aggregated property canvases
    c1, leg1 = draw_canvas("c_leading_proton", hist_list_proton)#, ymax=31000)
    #c1.SaveAs("c_leading_proton.pdf")
    c20,leg20 = draw_canvas("c_leading_mu", hist_list_mu)#, ymax=31000)
    c2, leg2 = draw_canvas("c_excitation_energy", hist_list_ex,log=True)
    c21, leg21 = draw_canvas("c_excitation_energy_postfsi", hist_list_expostfsi ,log=True)
    c3, leg3 = draw_canvas_with_ratio("c_avail_energy", hist_list_avail)#, ymax=34000)
    c4, leg4 = draw_canvas_with_ratio("c_nop_avail_energy", hist_list_nop_avail)#, ymax=30000)
    c13, leg13 = draw_canvas_with_ratio("c_non_incl_avail_energy", hist_list_non_incl_avail)#, ymax=35000)
    
    # Invisible & Total Energy Canvases
    c18, leg18 = draw_canvas_with_ratio("c_invisible_energy", hist_list_invisible)
    c19, leg19 = draw_canvas_with_ratio("c_total_energy", hist_list_total_e)  # ADDED: Total Energy canvas

    # KDAR canvas
    #c14, leg14 = draw_kdar_canvas("c_missing_energy_kdar", hist_list_missing_kdar) 

    # Transverse kinematics canvases
    c16, leg16 = draw_canvas_with_ratio("c_dpt", hist_list_dpt)#, ymax=65000)
    c17, leg17 = draw_canvas_with_ratio("c_dat", hist_list_dat)#, ymax=26000)


    c22_mom, leg22_mom = draw_canvas("c_tot_proton_all", hist_list_momp_tot)
    c5_mom, leg5_mom = draw_canvas_casc_deex("c_tot_proton_mom", hist_list_momp_casc, hist_list_momp_deex)
    c5_mult, leg5_mult = draw_canvas("c_tot_proton_mult", hist_list_mult_p)
    
    # Neutron plots (Momentum with ratio, Multiplicity without ratio)
    c23_mom, leg23_mom = draw_canvas("c_tot_neutron_all", hist_list_momn_tot)
    c7_mom, leg7_mom = draw_canvas_casc_deex("c_tot_neutron_mom", hist_list_momn_casc, hist_list_momn_deex)
    c7_mult, leg7_mult = draw_canvas("c_tot_neutron_mult", hist_list_mult_n)

    # Other Particle Specific Split Canvases (Momentum + Multiplicity)
    c6, legs6 = draw_side_by_side_canvas("c_tot_gamma_mom_and_mult", hist_list_momgamma, hist_list_mult_gamma)
    c8, legs8 = draw_side_by_side_canvas("c_tot_deuteron_mom_and_mult", hist_list_momd, hist_list_mult_d)#, ymax_left=4500, ymax_right=550000)
    c9, legs9 = draw_side_by_side_canvas("c_tot_triton_mom_and_mult", hist_list_momt, hist_list_mult_t)
    c15, legs15 = draw_side_by_side_canvas("c_tot_he3_mom_and_mult", hist_list_momhe3, hist_list_mult_he3)#, ymax_left=1600, ymax_right=550000)
    #c10, legs10 = draw_side_by_side_canvas("c_tot_alpha_mom_and_mult", hist_list_moma, hist_list_mult_a)#, ymax_left=7000, ymax_right=350000)
    c10, legs10 = draw_canvas("c_tot_alpha_mom", hist_list_moma_deex)
    c11, legs11 = draw_side_by_side_canvas("c_tot_pi0_mom_and_mult", hist_list_mompi0, hist_list_mult_pi0)
    c12, legs12 = draw_side_by_side_canvas("c_tot_pipm_mom_and_mult", hist_list_mompipm, hist_list_mult_pipm)


    canvases = [
        (c1, "c_leading_proton"),
        (c2, "c_excitation_energy"),
        (c21, "c_excitation_energy_postfsi"),
        (c3, "c_avail_energy"),
        (c4, "c_nop_avail_energy"),
        (c13, "c_non_incl_avail_energy"),
        (c18, "c_invisible_energy"),
        (c19, "c_total_energy"),
        (c16, "c_dpt"),
        (c17, "c_dat"),
        (c22_mom, "c_tot_proton_all"),
        (c5_mom, "c_tot_proton_mom"),
        (c5_mult, "c_tot_proton_mult"),
        (c23_mom, "c_tot_neutron_all"),
        (c7_mom, "c_tot_neutron_mom"),
        (c7_mult, "c_tot_neutron_mult"),
        (c6, "c_tot_gamma_mom_and_mult"),
        (c8, "c_tot_deuteron_mom_and_mult"),
        (c9, "c_tot_triton_mom_and_mult"),
        (c15, "c_tot_he3_mom_and_mult"),
        (c10, "c_tot_alpha_mom_and_mult"),
        (c11, "c_tot_pi0_mom_and_mult"),
        (c12, "c_tot_pipm_mom_and_mult"),
        (c20, "c_leading_mu")
    ]





    output_dir = "CCQE_combined_plots"
    os.makedirs(output_dir, exist_ok=True)

    output_pdf = os.path.join(output_dir, f"paper_plots_{name[0]}.pdf")

    canvases = [c1, c2, c21, c3, c4, c13, c18, c19, c16, c17,  c22_mom, c5_mom, c5_mult, c6, c23_mom,c7_mom, c7_mult, c8, c9, c15, c10, c11, c12,c20]

    for i, canvas in enumerate(canvases):
        canvas.cd()
        canvas.Update()
        
        if i == 0:
            canvas.Print(f"{output_pdf}(")
        elif i == len(canvases) - 1:
            canvas.Print(f"{output_pdf})")
        else:
            canvas.Print(output_pdf)

  
    out_file.cd() 

    for color_idx in [1179, 1180, 1181, 1182]:
        tcolor_obj = ROOT.gROOT.GetColor(color_idx)
        if tcolor_obj:
            tcolor_obj.Write()

    c1.Write()
    c2.Write()
    c21.Write()
    c3.Write()
    c4.Write()
    c5_mult.Write()
    c5_mom.Write()
    c6.Write()
    c23_mom.Write()
    c7_mult.Write()
    c7_mom.Write()
    c8.Write()
    c9.Write()
    c15.Write()
    c10.Write()
    c11.Write()
    c12.Write()
    c13.Write()
    #c14.Write()
    c16.Write()
    c17.Write()
    c18.Write()
    c19.Write()  # ADDED: Write Total Energy canvas
    c20.Write()
    
    out_file.Close()


def pion_combined_plots(filename, name, filename_list=None):
    if filename_list is None:
        filename_list = []
        
    files_to_process = [filename] + filename_list
    
    # Original lists
    hist_list_proton = []
    hist_list_ex = []
    hist_list_avail = []
    hist_list_nop_avail = []
    hist_list_non_incl_avail = [] 
    
    # KDAR Missing Energy
    #hist_list_missing_kdar = []

    # Transverse Kinematics
    hist_list_dpt = []
    hist_list_dat = []

    # Invisible & Total Energy Lists
    hist_list_invisible = []
    hist_list_total_e = []  # ADDED: Total energy list
    
    # Momentum lists
    hist_list_momp = [] 
    hist_list_momgamma = []
    hist_list_momn = []
    hist_list_momd = []
    hist_list_momt = []
    hist_list_momhe3 = []
    hist_list_moma = []
    hist_list_mompi0 = []
    hist_list_mompip = []
    hist_list_mompim = []

    # Multiplicity lists
    hist_list_mult_p = []
    hist_list_mult_gamma = []
    hist_list_mult_n = []
    hist_list_mult_d = []
    hist_list_mult_t = []
    hist_list_mult_he3 = []
    hist_list_mult_a = []
    hist_list_mult_pi0 = []
    hist_list_mult_pip = []
    hist_list_mult_pim = []

    """
        styles = [
            (ROOT.kBlack, 1),      # Dashed Red
            (ROOT.kBlue, 9),     # Dash-Dot Blue
            (ROOT.kRed, 9),     # Dash-Dot Blue
            (ROOT.kGreen+2, 5),
            (ROOT.kOrange-3,9), # Nieves Color 
            (ROOT.kBlue, 9),     #
    # Long-Dash Green
            (ROOT.kBlue, 9),     # Dash-Dot Blue
            (ROOT.kRed, 9),     # Dash-Dot Blue
            (ROOT.kMagenta, 2),  # Short-Dash Magenta
            (ROOT.kCyan, 1)      # Solid Cyan (fallback)
        ]
    """
    #legend_name = ["No Subprimary Nucleon", "NEUT SRC Definition", "NuWro SRC Definition"]#
    #legend_name = ["No Kinematic Threshold", "14 MeV Threshold", "28 MeV Threshold", "NEUT SRC"]
    #legend_name = ["sf*", "LFG", "sf"]
    #legend_name = ["sf*", "LFG"]
    #legend_name = ["SF","Nieves LFG"] 

    """
    legend_name = ["SF", 
                   #"SF no cascade",
                   "Nieves LFG", 
                   #"Nieves LFG no cascade"
                   ]
    styles = [
                (ROOT.kBlue, 1),  # Uses index 2 (modified to Petroff Blue)
                #(ROOT.kBlue, 4),  # Modified to Petroff Grape
                (ROOT.kOrange + 7, 1),  # Modified to Petroff Red
                #(ROOT.kOrange + 7, 4),  # Modified to Petroff Yellow
            ]
    """

    """
    legend_name = ["INCL + ABLA", 
                "NEUT Cascade", 
                ]
    styles = [
                (ROOT.kBlue, 1),  # Uses index 2 (modified to Petroff Blue)
                (ROOT.kAzure,1)
                #(ROOT.kBlue, 4),  # Modified to Petroff Grape
                #(ROOT.kOrange + 7, 1),  # Modified to Petroff Red
                #(ROOT.kOrange + 7, 4),  # Modified to Petroff Yellow
            ]
    """
    """
    legend_name = ["No SRC Contribution", 
                   "NEUT SRC",
                   "NuWro SRC, no Kin. Cut",
                   "NuWro SRC, 14 MeV Kin. Cut", 
                   "NuWro SRC, iso Kin. Cut", 
                ]
    styles = [
                
                (ROOT.kAzure,1),
                (ROOT.kBlue, 1), 
                (ROOT.kOrange+7,1),
                (ROOT.kRed,1), 
                (ROOT.kMagenta,1), 
            ]
    """


    legend_name = [
                    "NEUT Cascade",
                    "INCL + ABLA",        
                ]

    styles = [
                
                (ROOT.kAzure,1),
                (ROOT.kBlue, 1), 
                (ROOT.kOrange+7,1),
            ]


                #(ROOT.kBlue, 4),  # Modified to Petroff Grape
                #(ROOT.kOrange + 7, 1),  # Modified to Petroff Red
                #(ROOT.kOrange + 7, 4),  # Modified to Petroff Yellow


    #styles = [
    #    (ROOT.kBlue, 1),  # Uses index 2 (modified to Petroff Blue)
    #    #(ROOT.kMagenta, 8),  # Modified to Petroff Grape
    #    (ROOT.kRed, 9),  # Modified to Petroff Red
    #    #(ROOT.kOrange + 7, 7),  # Modified to Petroff Yellow
    #]

   
    
    #legend_name = ["INCL + ABLA ","NEUT Cascade"]
    

    count = 0 
    
    for fname in files_to_process:
        print("processing file {}".format(fname))
        filename_chunk = legend_name[count]
        count += 1
        
        f = ROOT.TFile(fname)
        t = f.Get("neuttree")   

        h_leadin = ROOT.TH1D(f"h_leading_{filename_chunk}", "Leading Proton Distribution", 50, 0.0, 1500)
        h_leadin.Sumw2()
        h_leadin.SetTitle(";P_{N} [MeV/c];Counts")

        h_ex = ROOT.TH1D(f"h_excitation_E_{filename_chunk}", "Excitation Energy Distribution", 52, -4.0, 100)
        h_ex.Sumw2()
        h_ex.SetTitle(";E_{x} [MeV];Counts")

        h_eavail = ROOT.TH1D(f"h_avail_E_{filename_chunk}", "Hadronic Available Energy Distribution", 40, 0.0, 400)
        h_eavail.Sumw2()
        h_eavail.SetTitle(";E_{Had} [MeV];Counts")

        h_notpeavail = ROOT.TH1D(f"h_avail_E_noP_{filename_chunk}", "Hadronic Available Energy Distribution (no p)", 100, 0.0, 1000)
        h_notpeavail.Sumw2()
        h_notpeavail.SetTitle(";(E_{Had} - T_{p}) [MeV];Counts")

        h_non_incl_avail = ROOT.TH1D(f"h_avail_E_nonINCL_{filename_chunk}", "Non-INCL Available Energy Distribution", 100, 0.0, 1000)
        h_non_incl_avail.Sumw2()
        h_non_incl_avail.SetTitle(";E_{Non-INCL} [MeV];Counts")

        # Invisible Energy Histogram
        h_invisible = ROOT.TH1D(f"h_invisible_E_{filename_chunk}", "Invisible Energy Distribution", 40, 0.0, 400)
        h_invisible.Sumw2()
        h_invisible.SetTitle(";E_{Inv} [MeV];Counts")

        # ADDED: Total Energy Histogram
        h_total_e = ROOT.TH1D(f"h_tot_E_{filename_chunk}", "Total Energy Distribution", 50, 0.0, 1000)
        h_total_e.Sumw2()
        h_total_e.SetTitle(";E_{Total} [MeV];Counts")

        h_dpt = ROOT.TH1D(f"h_dpt_{filename_chunk}", "#delta p_{T} Distribution", 50, 0.0, 1000.0)
        h_dpt.Sumw2()
        h_dpt.SetTitle(";#delta p _{T} [MeV/c];Counts")

        h_dat = ROOT.TH1D(f"h_dat_{filename_chunk}", "#delta #alpha_{T} Distribution", 36, 0.0, 180.0)
        h_dat.Sumw2()
        h_dat.SetTitle(";#delta #alpha _{T} [deg];Counts")

        # --- Momentum histograms ---
        h_momp = ROOT.TH1D(f"h_momp_{filename_chunk}", "Total Proton Momentum Distribution", 100, 0.0, 2000)
        h_momp.Sumw2()
        h_momp.SetTitle(";P_{p}^{Total} [MeV/c];Counts")

        h_momgamma = ROOT.TH1D(f"h_momgamma_{filename_chunk}", "Total Photon Momentum Distribution", 100, 0.0, 30)
        h_momgamma.Sumw2()
        h_momgamma.SetTitle(";P_{#gamma}^{Total} [MeV/c];Counts")

        h_momn = ROOT.TH1D(f"h_momn_{filename_chunk}", "Total Neutron Momentum Distribution", 50, 0.0, 1000)
        h_momn.Sumw2()
        h_momn.SetTitle(";P_{n}^{Total} [MeV/c];Counts")

        h_momd = ROOT.TH1D(f"h_momd_{filename_chunk}", "Total Deuteron Momentum Distribution", 100, 0.0, 800)
        h_momd.Sumw2()
        h_momd.SetTitle(";P_{d}^{Total} [MeV/c];Counts")

        h_momt = ROOT.TH1D(f"h_momt_{filename_chunk}", "Total Triton Momentum Distribution", 100, 0.0, 800)
        h_momt.Sumw2()
        h_momt.SetTitle(";P_{t}^{Total} [MeV/c];Counts")

        h_momhe3 = ROOT.TH1D(f"h_momhe3_{filename_chunk}", "Total Helium-3 Momentum Distribution", 100, 0.0, 800)
        h_momhe3.Sumw2()
        h_momhe3.SetTitle(";P_{^{3}He}^{Total} [MeV/c];Counts")

        h_moma = ROOT.TH1D(f"h_moma_{filename_chunk}", "Total Alpha Momentum Distribution", 100, 0.0, 1000)
        h_moma.Sumw2()
        h_moma.SetTitle(";P_{#alpha}^{Total} [MeV/c];Counts")
        
        h_mompi0 = ROOT.TH1D(f"h_mompi0_{filename_chunk}", "Total #pi^{0} Momentum Distribution", 100, 0.0, 1000)
        h_mompi0.Sumw2()
        h_mompi0.SetTitle(";P_{#pi^{0}}^{Total} [MeV/c];Counts")
        
        h_mompip = ROOT.TH1D(f"h_mompipm_{filename_chunk}", "Total #pi^{+} Momentum Distribution", 100, 0.0, 1000)
        h_mompip.Sumw2()
        h_mompip.SetTitle(";P_{#pi^{+}}^{Total} [MeV/c];Counts")

        h_mompim = ROOT.TH1D(f"h_mompipm_{filename_chunk}", "Total #pi^{-} Momentum Distribution", 100, 0.0, 500)
        h_mompim.Sumw2()
        h_mompim.SetTitle(";P_{#pi^{-}}^{Total} [MeV/c];Counts")

        # --- Multiplicity histograms ---
        h_mult_p = ROOT.TH1D(f"h_mult_p_{filename_chunk}", "Proton Multiplicity", 15, 0.0, 15)
        h_mult_p.Sumw2()
        h_mult_p.SetTitle(";N_{p};Counts")

        h_mult_gamma = ROOT.TH1D(f"h_mult_gamma_{filename_chunk}", "Photon Multiplicity", 20, 0.0, 20)
        h_mult_gamma.Sumw2()
        h_mult_gamma.SetTitle(";N_{#gamma};Counts")

        h_mult_n = ROOT.TH1D(f"h_mult_n_{filename_chunk}", "Neutron Multiplicity", 15, 0.0, 15)
        h_mult_n.Sumw2()
        h_mult_n.SetTitle(";N_{n};Counts")

        h_mult_d = ROOT.TH1D(f"h_mult_d_{filename_chunk}", "Deuteron Multiplicity", 10, 0.0, 10)
        h_mult_d.Sumw2()
        h_mult_d.SetTitle(";N_{d};Counts")

        h_mult_t = ROOT.TH1D(f"h_mult_t_{filename_chunk}", "Triton Multiplicity", 10, 0.0, 10)
        h_mult_t.Sumw2()
        h_mult_t.SetTitle(";N_{t};Counts")

        h_mult_he3 = ROOT.TH1D(f"h_mult_he3_{filename_chunk}", "Helium-3 Multiplicity", 10, 0.0, 10)
        h_mult_he3.Sumw2()
        h_mult_he3.SetTitle(";N_{^{3}He};Counts")

        h_mult_a = ROOT.TH1D(f"h_mult_a_{filename_chunk}", "Alpha Multiplicity", 10, 0.0, 10)
        h_mult_a.Sumw2()
        h_mult_a.SetTitle(";N_{#alpha};Counts")

        h_mult_pi0 = ROOT.TH1D(f"h_mult_pi0_{filename_chunk}", "#pi^{0} Multiplicity", 10, 0.0, 10)
        h_mult_pi0.Sumw2()
        h_mult_pi0.SetTitle(";N_{#pi^{0}};Counts")

        h_mult_pip = ROOT.TH1D(f"h_mult_pip_{filename_chunk}", "#pi^{+} Multiplicity", 10, 0.0, 10)
        h_mult_pip.Sumw2()
        h_mult_pip.SetTitle(";N_{#pi^{+}};Counts")

        h_mult_pim = ROOT.TH1D(f"h_mult_pip_{filename_chunk}", "#pi^{-} Multiplicity", 10, 0.0, 10)
        h_mult_pim.Sumw2()
        h_mult_pim.SetTitle(";N_{#pi^{-}};Counts")
        event_no = 0
        for event in t:
            nvect = event.vectorbranch
            nvect_class = nvect_reader(nvect) 
            event_no +=1
            if event_no % 100000 == 0:
                print("event number ", event_no)
            if event_no == 5000:
                break
            
            if nvect_class.eventType == EventType.piProd: 
                dpt, dat = nvect_class.get_deltaPT()
                
                h_dpt.Fill(dpt)
                h_dat.Fill(dat)

                if nvect_class.eventType == EventType.piProd:
                    if nvect_class.excitation_E_pion() < 0.0:
                        h_ex.Fill(0.0)
                    else:
                        h_ex.Fill(nvect_class.excitation_E_pion())

            
               # Missing_energy_kdar, E_Miss = nvect_class.Missing_energy_kdar2()
               # if Missing_energy_kdar:
               #     h_missing_kdar.Fill(E_Miss)

                particles, energies = nvect_class.particles()
                eavail = 0
                for particle, energy in zip(particles, energies):
                    if particle == 2112 and particle == 13:
                        continue
                    eavail += energy

                h_eavail.Fill(eavail) 

                # Counters for per-event totals
                non_proton_evail = 0
                non_INCL_eavail = 0 
                invisible_energy = 0
                total_E = 0  # ADDED: Initialized per event
                
                total_proton_mom = 0
                total_gamma_mom = 0
                total_neutron_mom = 0
                total_deuteron_mom = 0
                total_triton_mom = 0
                total_he3_mom = 0
                total_alpha_mom = 0
                total_pi0_mom = 0
                total_pip_mom = 0
                total_pim_mom = 0

                mult_p = 0
                mult_gamma = 0
                mult_n = 0
                mult_d = 0
                mult_t = 0
                mult_he3 = 0
                mult_a = 0
                mult_pi0 = 0
                mult_pip = 0
                mult_pim = 0
                
                for particle, energy in zip(particles, energies):
                    total_E += energy
                    if particle != 2212 and particle != 2112:
                        non_proton_evail += energy

                    if (particle < 1000000 and particle > 20) and particle != 2112:
                        #print(particle, energy)
                        non_INCL_eavail += energy

                    if particle != 2212 and abs(particle) != 211:
                        invisible_energy += energy
                    
                    # Convert Kinetic Energy to Momentum (MeV/c)
                    if particle == 2212: # Proton
                        total_proton_mom += np.sqrt((energy+938.27)**2 - 938.27**2)
                        #total_proton_mom += energy
                        mult_p += 1
                    elif particle == 22: # Gamma
                        total_gamma_mom += energy # p = E for photon
                        mult_gamma += 1
                    elif particle == 2112: # Neutron
                        total_neutron_mom += np.sqrt((energy+939.57)**2 - 939.57**2)
                        mult_n += 1
                    elif particle == 1000010020: # Deuteron
                        total_deuteron_mom += np.sqrt((energy+1875.61)**2 - 1875.61**2)
                        mult_d += 1
                    elif particle == 1000010030: # Triton
                        total_triton_mom += np.sqrt((energy+2808.92)**2 - 2808.92**2)
                        mult_t += 1
                    elif particle == 1000020030: # Helium-3
                        total_he3_mom += np.sqrt((energy+2808.39)**2 - 2808.39**2)
                        mult_he3 += 1
                    elif particle == 1000020040: # Alpha
                        total_alpha_mom += np.sqrt((energy+3727.38)**2 - 3727.38**2) 
                        mult_a += 1
                    elif particle == 111: # Pi0
                        total_pi0_mom += np.sqrt((energy+134.98)**2 - 134.98**2)
                        mult_pi0 += 1
                    elif particle == 211: # Pi+-
                        total_pip_mom += np.sqrt((energy+139.57)**2 - 139.57**2)
                        mult_pip += 1
                    elif particle == -211: # Pi+-
                        total_pim_mom += np.sqrt((energy+139.57)**2 - 139.57**2)
                        mult_pim += 1
                
                # Fill energy histograms (skipping zeros)
                if non_proton_evail > 0: h_notpeavail.Fill(non_proton_evail)
                if non_INCL_eavail > 0: h_non_incl_avail.Fill(non_INCL_eavail)
                if invisible_energy > 0: h_invisible.Fill(invisible_energy)
                if total_E > 0: h_total_e.Fill(total_E)  # ADDED: Fill total energy

                # Fill momentum histograms
                if total_proton_mom > 0: h_momp.Fill(total_proton_mom)
                if total_gamma_mom > 0: h_momgamma.Fill(total_gamma_mom)
                if total_neutron_mom > 0: h_momn.Fill(total_neutron_mom)
                if total_deuteron_mom > 0: h_momd.Fill(total_deuteron_mom)
                if total_triton_mom > 0: h_momt.Fill(total_triton_mom)
                if total_he3_mom > 0: h_momhe3.Fill(total_he3_mom)
                if total_alpha_mom > 0: h_moma.Fill(total_alpha_mom)
                if total_pi0_mom > 0: h_mompi0.Fill(total_pi0_mom)
                if total_pip_mom > 0: h_mompip.Fill(total_pip_mom)
                if total_pim_mom > 0: h_mompim.Fill(total_pim_mom)

                # Fill multiplicity histograms
                h_mult_p.Fill(mult_p)
                h_mult_gamma.Fill(mult_gamma)
                h_mult_n.Fill(mult_n)
                h_mult_d.Fill(mult_d)
                h_mult_t.Fill(mult_t)
                h_mult_he3.Fill(mult_he3)
                h_mult_a.Fill(mult_a)
                h_mult_pi0.Fill(mult_pi0)
                h_mult_pip.Fill(mult_pip)
                h_mult_pim.Fill(mult_pim)

                if nvect_class.HMPMom:
                    h_leadin.Fill(nvect_class.HMPMom)

        # Set directories
        h_leadin.SetDirectory(0) 
        h_ex.SetDirectory(0) 
        h_eavail.SetDirectory(0) 
        h_notpeavail.SetDirectory(0) 
        h_non_incl_avail.SetDirectory(0)
        h_invisible.SetDirectory(0)
        h_total_e.SetDirectory(0)  # ADDED: Detach from file directory
       # h_missing_kdar.SetDirectory(0)

        h_dpt.SetDirectory(0)
        h_dat.SetDirectory(0)
        
        h_momp.SetDirectory(0)
        h_momgamma.SetDirectory(0)
        h_momn.SetDirectory(0)
        h_momd.SetDirectory(0)
        h_momt.SetDirectory(0)
        h_momhe3.SetDirectory(0)
        h_moma.SetDirectory(0)
        h_mompi0.SetDirectory(0)
        h_mompip.SetDirectory(0)
        h_mompim.SetDirectory(0)

        h_mult_p.SetDirectory(0)
        h_mult_gamma.SetDirectory(0)
        h_mult_n.SetDirectory(0)
        h_mult_d.SetDirectory(0)
        h_mult_t.SetDirectory(0)
        h_mult_he3.SetDirectory(0)
        h_mult_a.SetDirectory(0)
        h_mult_pi0.SetDirectory(0)
        h_mult_pip.SetDirectory(0)
        h_mult_pim.SetDirectory(0)
        
        # Append to lists
        hist_list_proton.append((filename_chunk, h_leadin))
        hist_list_ex.append((filename_chunk, h_ex))
        hist_list_avail.append((filename_chunk, h_eavail))
        hist_list_nop_avail.append((filename_chunk, h_notpeavail))
        hist_list_non_incl_avail.append((filename_chunk, h_non_incl_avail))
        hist_list_invisible.append((filename_chunk, h_invisible))
        hist_list_total_e.append((filename_chunk, h_total_e))  # ADDED: Append to list
       # hist_list_missing_kdar.append((filename_chunk, h_missing_kdar))

        hist_list_dpt.append((filename_chunk, h_dpt))
        hist_list_dat.append((filename_chunk, h_dat))
        
        hist_list_momp.append((filename_chunk, h_momp))
        hist_list_momgamma.append((filename_chunk, h_momgamma))
        hist_list_momn.append((filename_chunk, h_momn))
        hist_list_momd.append((filename_chunk, h_momd))
        hist_list_momt.append((filename_chunk, h_momt))
        hist_list_momhe3.append((filename_chunk, h_momhe3))
        hist_list_moma.append((filename_chunk, h_moma))
        hist_list_mompi0.append((filename_chunk, h_mompi0))
        hist_list_mompip.append((filename_chunk, h_mompip))
        hist_list_mompim.append((filename_chunk, h_mompim))

        hist_list_mult_p.append((filename_chunk, h_mult_p))
        hist_list_mult_gamma.append((filename_chunk, h_mult_gamma))
        hist_list_mult_n.append((filename_chunk, h_mult_n))
        hist_list_mult_d.append((filename_chunk, h_mult_d))
        hist_list_mult_t.append((filename_chunk, h_mult_t))
        hist_list_mult_he3.append((filename_chunk, h_mult_he3))
        hist_list_mult_a.append((filename_chunk, h_mult_a))
        hist_list_mult_pi0.append((filename_chunk, h_mult_pi0))
        hist_list_mult_pip.append((filename_chunk, h_mult_pip))
        hist_list_mult_pim.append((filename_chunk, h_mult_pim))

        f.Close()

    ROOT.gStyle.SetOptTitle(0)
    out_file = ROOT.TFile("paper_combined_comparisons_{}.root".format(name), "RECREATE")
    latex_font = 132

    # Canvas helper definitions remain standard...
    def draw_canvas(canvas_name, hist_list, ymax=None, ymin=None):
        c = ROOT.TCanvas(canvas_name, canvas_name, 650, 600)
        c.SetGrayscale(False)
        c.SetLeftMargin(0.12)   # Default is ~0.12
        c.SetRightMargin(0.08)  # Default is ~0.10
        c.SetTopMargin(0.06)    # Default is ~0.10
        c.SetBottomMargin(0.12) # Default is ~0.12


        ROOT.gPad.SetTickx(1)
        ROOT.gPad.SetTicky(1)
        #ROOT.gPad.SetBottomMargin(0.15) 
        #ROOT.gPad.SetLeftMargin(0.12)
        
        
        legend = ROOT.TLegend(0.55, 0.75, 0.82, 0.85)
        legend.SetNColumns(1) 
        legend.SetTextFont(latex_font) 
        legend.SetBorderSize(0) 
        legend.SetFillStyle(1001)  # Solid fill
        legend.SetFillColor(ROOT.kWhite)
        legend.SetTextSize(0.03) 

        if ymax is None:
            global_max = max(h.GetMaximum() for _, h in hist_list)
            ymax = global_max * 1.2

        for i, (name, h) in enumerate(hist_list):
            color, style = styles[i]
            
            h.SetLineColor(color)
            h.SetLineStyle(style)
            h.SetLineWidth(4) 
            h.SetStats(0) 

            h.GetXaxis().CenterTitle(True)
            h.GetYaxis().CenterTitle(True)

            h.GetYaxis().SetMaxDigits(3)

            # Force X-axis to NEVER use scientific exponent notation
            h.GetXaxis().SetNoExponent(True)

            h.GetXaxis().SetTitleFont(latex_font)
            h.GetXaxis().SetLabelFont(latex_font)
            h.GetYaxis().SetTitleFont(latex_font)
            h.GetYaxis().SetLabelFont(latex_font)
            h.GetXaxis().SetTitleSize(0.05) 
            h.GetXaxis().SetLabelSize(0.05) 
            h.GetXaxis().SetTitleOffset(1.0)

            h.GetYaxis().SetTitleSize(0.05) 
            h.GetYaxis().SetLabelSize(0.05) 
            h.GetYaxis().SetTitleOffset(1.2)

            if i == 0:
                h.SetMaximum(ymax)
                if ymin is not None:
                    h.SetMinimum(ymin)
                h.Draw("HIST") 
            else:
                h.Draw("HIST SAME")
                
            legend.AddEntry(h, name, "l")

        legend.Draw()
        return c, legend 
        
    def draw_canvas_with_ratio(canvas_name, hist_list, ymax=None, ymin=None, ratio_min = 0.5, ratio_max = 1.5):
        c = ROOT.TCanvas(canvas_name, canvas_name, 600, 750)
        c.SetGrayscale(False)
        c.SetLeftMargin(0.12)   # Default is ~0.12
        c.SetRightMargin(0.08)  # Default is ~0.10
        c.SetTopMargin(0.06)    # Default is ~0.10
        c.SetBottomMargin(0.0) # Default is ~0.12

        pad1 = ROOT.TPad("pad1", "pad1", 0.0, 0.4, 1.0, 0.98)
        pad1.SetBottomMargin(0.02)
        pad1.SetLeftMargin(0.12)   # Default is ~0.12
        pad1.SetRightMargin(0.08)  # Default is ~0.10
        pad1.SetTopMargin(0.06)    # Default is ~0.10
        pad1.SetBottomMargin(0.02) # Default is ~0.12

        pad1.SetTickx(1)
        pad1.SetTicky(1)
        pad1.Draw()
        pad1.cd()

        legend = ROOT.TLegend(0.42, 0.17, 0.8, 0.9)
        legend.SetNColumns(1) 
        legend.SetTextFont(latex_font) 
        legend.SetBorderSize(0) 
        legend.SetFillStyle(1001)  # Solid fill
        legend.SetFillColor(ROOT.kWhite)
        legend.SetTextSize(0.04) 

        if ymax is None:
            global_max = max(h.GetMaximum() for _, h in hist_list)
            ymax = global_max * 1.2

        xaxis_title = hist_list[0][1].GetXaxis().GetTitle()

        for i, (name, h) in enumerate(hist_list):
            color, style = styles[i % len(styles)]
            #print(color)
            h.SetLineColor(color)
            h.SetLineStyle(style)
            h.SetLineWidth(3) 
            h.SetStats(0) 

            h.GetXaxis().SetLabelSize(0)
            h.GetXaxis().SetTitleSize(0)
            
            h.GetYaxis().SetTitleFont(latex_font)
            h.GetYaxis().SetLabelFont(latex_font)
            h.GetYaxis().SetTitleOffset(1.2)

            h.GetYaxis().SetTitleSize(0.05) 
            h.GetYaxis().SetLabelSize(0.05) 

            h.GetXaxis().CenterTitle(True)
            h.GetYaxis().CenterTitle(True)

            h.GetYaxis().SetMaxDigits(3)
            h.GetXaxis().SetNoExponent(True)
            

            print(
                f"Hist: {h.GetName()} | Set Color: {color} | Real LineColor: {h.GetLineColor()}"
            )

            if i == 0:
                h.SetMaximum(ymax)
                if ymin is not None:
                    h.SetMinimum(ymin)
                h.Draw("HIST") 
            else:
                h.Draw("HIST SAME")
                
            legend.AddEntry(h, name, "l")

        legend.Draw()

        c.cd()
        pad2 = ROOT.TPad("pad2", "pad2", 0.0, 0.15, 1.0, 0.4)
        pad2.SetLeftMargin(0.12)   # Default is ~0.12
        pad2.SetRightMargin(0.08)  # Default is ~0.10
        pad2.SetTopMargin(0.06)    # Default is ~0.10
        pad2.SetBottomMargin(0.40) # Default is ~0.12
        pad2.SetTickx(1)
        pad2.SetTicky(1)
        pad2.Draw()
        pad2.cd()

        h_base = hist_list[0][1]
        ratios = []
        for i in range(1, len(hist_list)):
            name, h = hist_list[i]
            h_ratio = h.Clone(f"{h.GetName()}_ratio")
            h_ratio.Divide(h_base)

            color, style = styles[i % len(styles)]
            h_ratio.SetLineColor(color)
            h_ratio.SetLineStyle(style)
            h_ratio.SetLineWidth(3)
            h_ratio.SetStats(0)

            h_ratio.SetTitle("")
            h_ratio.GetYaxis().SetTitle("Ratio to {}".format("No SRC"))
            h_ratio.GetYaxis().CenterTitle()
            h_ratio.GetYaxis().SetNdivisions(505)
            h_ratio.GetYaxis().SetTitleSize(0.11)
            h_ratio.GetYaxis().SetTitleFont(latex_font)
            h_ratio.GetYaxis().SetTitleOffset(0.5)
            h_ratio.GetYaxis().SetLabelFont(latex_font)
            h_ratio.GetYaxis().SetLabelSize(0.1)
            
            h_ratio.GetXaxis().SetTitle(xaxis_title)
            h_ratio.GetXaxis().SetTitleSize(0.13)
            h_ratio.GetXaxis().SetTitleFont(latex_font)
            h_ratio.GetXaxis().SetTitleOffset(1.0)
            h_ratio.GetXaxis().SetLabelFont(latex_font)
            h_ratio.GetXaxis().SetLabelSize(0.12)

            h_ratio.GetXaxis().CenterTitle(True)
            h_ratio.GetYaxis().CenterTitle(True)

            if i == 1:
                h_ratio.SetMinimum(ratio_min)
                h_ratio.SetMaximum(ratio_max)
                h_ratio.Draw("HIST")
            else:
                h_ratio.Draw("HIST SAME")
            ratios.append(h_ratio)
            
        line = ROOT.TLine(h_base.GetXaxis().GetXmin(), 1, h_base.GetXaxis().GetXmax(), 1)
        line.SetLineColor(ROOT.kBlack)
        line.SetLineStyle(2)
        line.Draw("SAME")

        c.pad1 = pad1
        c.pad2 = pad2
        c.ratios = ratios
        c.line = line

        pad1.Modified()
        pad1.Update()
        pad2.Modified()
        pad2.Update()
        
        c.cd()
        c.Modified()
        c.Update()
        
        return c, legend

    def draw_side_by_side_canvas(canvas_name, hist_list_energy, hist_list_mult, ymax_left=None, ymax_right=None):
        c = ROOT.TCanvas(canvas_name, canvas_name, 1400, 750)
        c.Divide(2, 1)

        def configure_and_draw(pad_idx, hist_list, ymax=None):
            c.cd(pad_idx)
            ROOT.gPad.SetTickx(1)
            ROOT.gPad.SetTicky(1)
            ROOT.gPad.SetBottomMargin(0.15) 
            ROOT.gPad.SetLeftMargin(0.12)

            legend = ROOT.TLegend(0.35, 0.7, 0.60, 0.85) 
            legend.SetNColumns(1) 
            legend.SetTextFont(latex_font) 
            legend.SetBorderSize(0) 
            legend.SetFillStyle(0)  
            legend.SetTextSize(0.04) 

            if ymax is None:
                global_max = max(h.GetMaximum() for _, h in hist_list)
                ymax = global_max * 1.2

            for i, (name, h) in enumerate(hist_list):
                color, style = styles[i % len(styles)]
                
                h.SetLineColor(color)
                h.SetLineStyle(style)
                h.SetLineWidth(2) 
                h.SetStats(0) 

                h.GetXaxis().SetTitleFont(latex_font)
                h.GetXaxis().SetLabelFont(latex_font)
                h.GetYaxis().SetTitleFont(latex_font)
                h.GetYaxis().SetLabelFont(latex_font)
                h.GetXaxis().SetTitleSize(0.05) 
                h.GetXaxis().SetLabelSize(0.05) 
                h.GetXaxis().SetTitleOffset(1.0)
                h.GetYaxis().SetTitleOffset(1.1)
                h.GetYaxis().SetTitleSize(0.05) 
                h.GetYaxis().SetLabelSize(0.05) 

                h.GetXaxis().CenterTitle(True)
                h.GetYaxis().CenterTitle(True)

                h.GetYaxis().SetMaxDigits(3)
                h.GetXaxis().SetNoExponent(True)

                if i == 0:
                    h.SetMaximum(ymax)
                    h.Draw("HIST") 
                else:
                    h.Draw("HIST SAME")
                    
                legend.AddEntry(h, name, "l")
            legend.Draw()
            return legend

        leg_e = configure_and_draw(1, hist_list_energy, ymax_left)
        leg_m = configure_and_draw(2, hist_list_mult, ymax_right)
        
        c.legends = (leg_e, leg_m)
        return c, c.legends
    
    # Aggregated property canvases
    c1, leg1 = draw_canvas("c_leading_proton", hist_list_proton)#, ymax=31000)
    c1.SaveAs("c_leading_proton.pdf")
    c2, leg2 = draw_canvas_with_ratio("c_excitation_energy", hist_list_ex,ratio_min=0.0,ratio_max=2.0)
    c3, leg3 = draw_canvas_with_ratio("c_avail_energy", hist_list_avail)#, ymax=34000)
    c4, leg4 = draw_canvas_with_ratio("c_nop_avail_energy", hist_list_nop_avail)#, ymax=30000)
    c13, leg13 = draw_canvas_with_ratio("c_non_incl_avail_energy", hist_list_non_incl_avail)#, ymax=35000)
    
    # Invisible & Total Energy Canvases
    c18, leg18 = draw_canvas_with_ratio("c_invisible_energy", hist_list_invisible)
    c19, leg19 = draw_canvas_with_ratio("c_total_energy", hist_list_total_e)  # ADDED: Total Energy canvas

    # KDAR canvas
    #c14, leg14 = draw_kdar_canvas("c_missing_energy_kdar", hist_list_missing_kdar) 

    # Transverse kinematics canvases
    c16, leg16 = draw_canvas_with_ratio("c_dpt", hist_list_dpt)#, ymax=65000)
    c17, leg17 = draw_canvas_with_ratio("c_dat", hist_list_dat)#, ymax=26000)

    c5_mom, leg5_mom = draw_canvas_with_ratio("c_tot_proton_mom", hist_list_momp)
    c5_mult, leg5_mult = draw_canvas("c_tot_proton_mult", hist_list_mult_p)
    
    # Neutron plots (Momentum with ratio, Multiplicity without ratio)
    c7_mom, leg7_mom = draw_canvas_with_ratio("c_tot_neutron_mom", hist_list_momn)
    c7_mult, leg7_mult = draw_canvas("c_tot_neutron_mult", hist_list_mult_n)

    # Other Particle Specific Split Canvases (Momentum + Multiplicity)
    c6, legs6 = draw_side_by_side_canvas("c_tot_gamma_mom_and_mult", hist_list_momgamma, hist_list_mult_gamma)
    c8, legs8 = draw_side_by_side_canvas("c_tot_deuteron_mom_and_mult", hist_list_momd, hist_list_mult_d)#, ymax_left=4500, ymax_right=550000)
    c9, legs9 = draw_side_by_side_canvas("c_tot_triton_mom_and_mult", hist_list_momt, hist_list_mult_t)
    c15, legs15 = draw_side_by_side_canvas("c_tot_he3_mom_and_mult", hist_list_momhe3, hist_list_mult_he3)#, ymax_left=1600, ymax_right=550000)
    c10, legs10 = draw_side_by_side_canvas("c_tot_alpha_mom_and_mult", hist_list_moma, hist_list_mult_a)#, ymax_left=7000, ymax_right=350000)
    c11, legs11 = draw_side_by_side_canvas("c_tot_pi0_mom_and_mult", hist_list_mompi0, hist_list_mult_pi0)
    c12, legs12 = draw_side_by_side_canvas("c_tot_pip_mom_and_mult", hist_list_mompip, hist_list_mult_pip)
    c13, legs13 = draw_side_by_side_canvas("c_tot_pim_mom_and_mult", hist_list_mompim, hist_list_mult_pim)


    c1, leg1 = draw_canvas("c_leading_proton", hist_list_proton)#, ymax=31000)

    canvases = [
        (c1, "c_leading_proton"),
        (c2, "c_excitation_energy"),
        (c3, "c_avail_energy"),
        (c4, "c_nop_avail_energy"),
        (c13, "c_non_incl_avail_energy"),
        (c18, "c_invisible_energy"),
        (c19, "c_total_energy"),
        (c16, "c_dpt"),
        (c17, "c_dat"),
        (c5_mom, "c_tot_proton_mom"),
        (c5_mult, "c_tot_proton_mult"),
        (c7_mom, "c_tot_neutron_mom"),
        (c7_mult, "c_tot_neutron_mult"),
        (c6, "c_tot_gamma_mom_and_mult"),
        (c8, "c_tot_deuteron_mom_and_mult"),
        (c9, "c_tot_triton_mom_and_mult"),
        (c15, "c_tot_he3_mom_and_mult"),
        (c10, "c_tot_alpha_mom_and_mult"),
        (c11, "c_tot_pi0_mom_and_mult"),
        (c12, "c_tot_pip_mom_and_mult"),
        (c13, "c_tot_pip_mom_and_mult"),
    ]



    output_pdf = "fruuuuuuuu_PI.pdf"
    canvases = [c1, c2, c3, c4, c13, c18, c19, c16, c17, c5_mom, c5_mult, c6, c7_mom, c7_mult, c8, c9, c15, c10, c11, c12,c13]

    for i, canvas in enumerate(canvases):
        canvas.cd()
        canvas.Update()
        
        if i == 0:
            canvas.Print(f"{output_pdf}(")
        elif i == len(canvases) - 1:
            canvas.Print(f"{output_pdf})")
        else:
            canvas.Print(output_pdf)

  
    out_file.cd() 

    for color_idx in [1179, 1180, 1181, 1182]:
        tcolor_obj = ROOT.gROOT.GetColor(color_idx)
        if tcolor_obj:
            tcolor_obj.Write()

    c1.Write()
    c2.Write()
    c3.Write()
    c4.Write()
    c5_mult.Write()
    c5_mom.Write()
    c6.Write()
    c7_mult.Write()
    c7_mom.Write()
    c8.Write()
    c9.Write()
    c15.Write()
    c10.Write()
    c11.Write()
    c12.Write()
    c13.Write()
    #c14.Write()
    c16.Write()
    c17.Write()
    c18.Write()
    c19.Write()  # ADDED: Write Total Energy canvas
    
    out_file.Close()



def ccqe_momentum_plots_side_by_side(filename, out_name, filename_list=None):
    label1 = "NEUT Cascade"
    label2 = "INCL + ABLA" 
    
    if filename_list is None:
        filename_list = []
        
    files_to_process = [filename] + filename_list

    # Petroff 10 Colorblind Palette via TColor.GetColor (older ROOT compatible)
    particle_specs = {
        "Proton":   {"pdg": [2212],       "mass": 938.27,  "color": ROOT.TColor.GetColor("#3f90da"), "style": 1},
        "Neutron":  {"pdg": [2112],       "mass": 939.57,  "color": ROOT.TColor.GetColor("#ffa90e"), "style": 1},
        "Deuteron": {"pdg": [1000010020], "mass": 1875.61, "color": ROOT.TColor.GetColor("#bd1f01"), "style": 1},
        "Triton":   {"pdg": [1000010030], "mass": 2808.92, "color": ROOT.TColor.GetColor("#832db6"), "style": 1},
        "He-3":     {"pdg": [1000020030], "mass": 2808.39, "color": ROOT.TColor.GetColor("#92dadd"), "style": 1},
        "Alpha":    {"pdg": [1000020040], "mass": 3727.38, "color": ROOT.TColor.GetColor("#a96b59"), "style": 1},
        "Pi0":      {"pdg": [111],        "mass": 134.98,  "color": ROOT.TColor.GetColor("#b8ab6f"), "style": 1},
        "Pi+":      {"pdg": [211],        "mass": 139.57,  "color": ROOT.TColor.GetColor("#e76300"), "style": 1},
        "Pi-":      {"pdg": [-211],       "mass": 139.57,  "color": ROOT.TColor.GetColor("#94a4a2"), "style": 1},
    }

    charged_specs = [2212, 1000010020, 1000010030, 1000020030, 1000020040, 211, -211]

    pdg_map = {}
    for name, spec in particle_specs.items():
        for pdg in spec["pdg"]:
            pdg_map[pdg] = (name, spec["mass"])

    charged_names = [name for name, spec in particle_specs.items() if any(p in charged_specs for p in spec["pdg"])]

    # -------------------------------------------------------------------------
    # 1. FILE PROCESSING
    # -------------------------------------------------------------------------
    def process_file(fname, file_tag):
        print(f"Processing file: {fname}")
        f = ROOT.TFile.Open(fname)
        if not f or f.IsZombie():
            raise RuntimeError(f"Could not open file: {fname}")
            
        t = f.Get("neuttree")

        # 1a. Total scalar momentum 1D histograms
        hists_mom = {}
        for name, spec in particle_specs.items():
            h = ROOT.TH1D(f"h_mom_{name}_{file_tag}", f"; p_{{tot}} [MeV/c];Counts", 50, 0.0, 1500.0)
            h.Sumw2()
            h.SetLineColor(spec["color"])
            h.SetLineStyle(spec["style"])
            h.SetLineWidth(2)
            h.SetStats(0)
            hists_mom[name] = h

        # 1b. Leading particle momentum 1D histograms
        hists_lead_mom = {}
        for name, spec in particle_specs.items():
            h = ROOT.TH1D(f"h_lead_mom_{name}_{file_tag}", f"; p_{{lead}} [MeV/c];Counts", 50, 0.0, 1500.0)
            h.Sumw2()
            h.SetLineColor(spec["color"])
            h.SetLineStyle(spec["style"])
            h.SetLineWidth(2)
            h.SetStats(0)
            hists_lead_mom[name] = h

        # 1c. TProfile: Fraction vs E_avail
        profs_frac = {}
        for name in charged_names:
            spec = particle_specs[name]
            p = ROOT.TProfile(f"p_frac_{name}_{file_tag}", "; E_{avail} [MeV]; Mean Energy Fraction", 150, 0.0, 1500.0)
            p.SetLineColor(spec["color"])
            p.SetLineStyle(spec["style"])
            p.SetLineWidth(2)
            p.SetStats(0)
            profs_frac[name] = p

        # 1d. Weighted E_avail 1D histograms for stacking
        hists_eavail_frac = {}
        for name in charged_names:
            spec = particle_specs[name]
            h = ROOT.TH1D(f"h_eavail_frac_{name}_{file_tag}", "; E_{avail} [MeV]; Counts", 150, 0.0, 1500.0)
            h.Sumw2()
            h.SetFillColor(spec["color"])
            h.SetLineColor(ROOT.kBlack)
            h.SetLineWidth(1)
            h.SetStats(0)
            hists_eavail_frac[name] = h

        event_no = 0
        for event in t:
            nvect = event.vectorbranch
            nvect_class = nvect_reader(nvect) 
            event_no += 1
            
            if event_no % 100000 == 0:
                print(f"[{file_tag}] Event number {event_no}")
            if event_no == 500000:
                break
            print(nvect_class.eventType)
            if nvect_class.eventType in ( EventType.NC, EventType.MF, EventType.SRC,EventType.piProd,EventType.twop2h):
                particles, energies = nvect_class.particles()
                
                event_species_mom = {p_name: 0.0 for p_name in particle_specs}
                event_species_lead_mom = {p_name: 0.0 for p_name in particle_specs}
                event_species_energy = {p_name: 0.0 for p_name in particle_specs}
                E_avail = 0.0

                for particle, energy in zip(particles, energies):
                    if particle in pdg_map:
                        p_name, mass = pdg_map[particle]
                        mom = np.sqrt((energy + mass)**2 - mass**2) if mass > 0 else energy
                        
                        # Total species momentum
                        event_species_mom[p_name] += mom
                        
                        # Leading particle momentum per species
                        if mom > event_species_lead_mom[p_name]:
                            event_species_lead_mom[p_name] = mom

                        if particle in charged_specs:
                            # Total energy for created pions (kinetic + rest mass 139.57 MeV), kinetic for nucleons
                            particle_energy = (energy + mass) if particle in [211, -211] else energy
                            E_avail += particle_energy
                            event_species_energy[p_name] += particle_energy

                # Fill total momentum TH1D histograms
                for p_name, tot_mom in event_species_mom.items():
                    if tot_mom > 0:
                        hists_mom[p_name].Fill(tot_mom)

                # Fill leading momentum TH1D histograms
                for p_name, lead_mom in event_species_lead_mom.items():
                    if lead_mom > 0:
                        hists_lead_mom[p_name].Fill(lead_mom)

                # Fill energy fraction TProfiles and weighted E_avail histograms
                if E_avail > 0:
                    for p_name in charged_names:
                        frac = event_species_energy[p_name] / E_avail
                        profs_frac[p_name].Fill(E_avail, frac)
                        if frac > 0:
                            hists_eavail_frac[p_name].Fill(E_avail, frac)

        for h in hists_mom.values():
            h.SetDirectory(0)
        for h in hists_lead_mom.values():
            h.SetDirectory(0)
        for p in profs_frac.values():
            p.SetDirectory(0)
        for h in hists_eavail_frac.values():
            h.SetDirectory(0)

        f.Close()
        return hists_mom, hists_lead_mom, profs_frac, hists_eavail_frac

    hists_m1, hists_lm1, profs_f1, hists_e1 = process_file(files_to_process[0], "f1")
    hists_m2, hists_lm2, profs_f2, hists_e2 = process_file(files_to_process[1], "f2")

    ROOT.gStyle.SetOptTitle(0)
    latex_font = 132

    # Canvas helper
    def create_top_bottom_canvas(canvas_name, title):
        c = ROOT.TCanvas(canvas_name, title, 1000, 1000)
        pad1 = ROOT.TPad(f"{canvas_name}_pad1", "Pad 1 (Top)",    0.0, 0.535, 1.0, 1.0)
        pad2 = ROOT.TPad(f"{canvas_name}_pad2", "Pad 2 (Bottom)", 0.0, 0.0,   1.0, 0.535)

        pad1.SetLeftMargin(0.14)   
        pad1.SetRightMargin(0.05)  
        pad1.SetTopMargin(0.08)
        pad1.SetBottomMargin(0.0)

        pad2.SetLeftMargin(0.14)   
        pad2.SetRightMargin(0.05) 
        pad2.SetTopMargin(0.0)     
        pad2.SetBottomMargin(0.15) 

        pad1.Draw()
        pad2.Draw()
        return c, pad1, pad2

    def draw_mom_pad(pad, hists_dict, pad_title, ymax=None, is_bottom_pad=False):
        pad.cd()
        ROOT.gPad.SetTickx(1)
        ROOT.gPad.SetTicky(1)

        if ymax is None:
            ymax = shared_ymax

        is_first = True
        for name, h in hists_dict.items():
            h.GetXaxis().SetNoExponent(True)
            h.GetYaxis().SetTitleFont(latex_font)
            h.GetYaxis().SetLabelFont(latex_font)
            h.GetYaxis().SetTitleSize(0.06)
            h.GetYaxis().SetLabelSize(0.05)
            h.GetYaxis().SetTitleOffset(1.3)
            h.GetYaxis().CenterTitle(True)

            h.SetMaximum(ymax)
            h.SetMinimum(0)

            if not is_bottom_pad:
                h.GetXaxis().SetLabelSize(0)
                h.GetXaxis().SetTitleSize(0)
                h.GetXaxis().SetTitle("")
                h.GetYaxis().ChangeLabel(1, -1, 0)
            else:
                h.GetXaxis().SetTitleFont(latex_font)
                h.GetXaxis().SetLabelFont(latex_font)
                h.GetXaxis().SetTitleSize(0.06)
                h.GetXaxis().SetLabelSize(0.05)
                h.GetXaxis().SetTitleOffset(1.1)
                h.GetXaxis().CenterTitle(True)
                h.GetYaxis().ChangeLabel(-1, -1, 0)

            if is_first:
                h.Draw("HIST")
                is_first = False
            else:
                h.Draw("HIST SAME")

        legend = None
        if not is_bottom_pad:
            legend = ROOT.TLegend(0.70, 0.40, 0.90, 0.90)
            legend.SetTextFont(latex_font)
            legend.SetBorderSize(0)
            legend.SetFillStyle(1001)
            legend.SetFillColor(ROOT.kWhite)
            legend.SetTextSize(0.04)
            for name, h in hists_dict.items():
                legend.AddEntry(h, name, "l")
            legend.Draw()

        latex = ROOT.TLatex()
        latex.SetNDC()
        latex.SetTextFont(133)
        latex.SetTextSize(26)
        x_pos = 0.18
        y_pos = 0.92 if is_bottom_pad else 0.84
        latex.DrawLatex(x_pos, y_pos, pad_title)

        return legend

    # -------------------------------------------------------------------------
    # 2a. CANVAS 1: TOTAL SCALAR MOMENTUM DISTRIBUTIONS
    # -------------------------------------------------------------------------
    max_f1 = max(h.GetMaximum() for h in hists_m1.values())
    max_f2 = max(h.GetMaximum() for h in hists_m2.values())
    shared_ymax = max(max_f1, max_f2) * 1.25
    if shared_ymax <= 0:
        shared_ymax = 1.0

    c_mom, pad1_m, pad2_m = create_top_bottom_canvas("c_momentum", "Momentum Distributions")
    c_mom.leg1 = draw_mom_pad(pad1_m, hists_m1, label1, ymax=shared_ymax, is_bottom_pad=False)
    c_mom.leg2 = draw_mom_pad(pad2_m, hists_m2, label2, ymax=shared_ymax, is_bottom_pad=True)
    c_mom.Update()
    c_mom.SaveAs(f"{out_name}_momentum.pdf")

    # -------------------------------------------------------------------------
    # 2b. CANVAS 2: LEADING PARTICLE MOMENTUM DISTRIBUTIONS
    # -------------------------------------------------------------------------
    max_lm1 = max(h.GetMaximum() for h in hists_lm1.values())
    max_lm2 = max(h.GetMaximum() for h in hists_lm2.values())
    shared_ymax_lead = max(max_lm1, max_lm2) * 1.25
    if shared_ymax_lead <= 0:
        shared_ymax_lead = 1.0

    c_lead_mom, pad1_lm, pad2_lm = create_top_bottom_canvas("c_leading_momentum", "Leading Momentum Distributions")
    c_lead_mom.leg1 = draw_mom_pad(pad1_lm, hists_lm1, label1, ymax=shared_ymax_lead, is_bottom_pad=False)
    c_lead_mom.leg2 = draw_mom_pad(pad2_lm, hists_lm2, label2, ymax=shared_ymax_lead, is_bottom_pad=True)
    c_lead_mom.Update()
    c_lead_mom.SaveAs(f"{out_name}_leading_momentum.pdf")

    # -------------------------------------------------------------------------
    # 3. CANVAS 3: CHARGED ENERGY FRACTION VS E_AVAIL PROFILES
    # -------------------------------------------------------------------------
    c_frac, pad1_f, pad2_f = create_top_bottom_canvas("c_energy_fractions", "Charged Energy Fractions")

    def draw_frac_pad(pad, profs_dict, pad_title, is_bottom_pad=False):
        pad.cd()
        ROOT.gPad.SetTickx(1)
        ROOT.gPad.SetTicky(1)

        is_first = True
        for name, p in profs_dict.items():
            p.GetXaxis().SetNoExponent(True)
            p.GetYaxis().SetTitleFont(latex_font)
            p.GetYaxis().SetLabelFont(latex_font)
            p.GetYaxis().SetTitleSize(0.045)
            p.GetYaxis().SetLabelSize(0.04)
            p.GetYaxis().SetTitleOffset(1.3)
            p.GetYaxis().CenterTitle(True)

            p.SetMaximum(1.0)
            p.SetMinimum(0.0)

            if not is_bottom_pad:
                p.GetXaxis().SetLabelSize(0)
                p.GetXaxis().SetTitleSize(0)
                p.GetXaxis().SetTitle("")
                p.GetYaxis().ChangeLabel(1, -1, 0)
            else:
                p.GetXaxis().SetTitleFont(latex_font)
                p.GetXaxis().SetLabelFont(latex_font)
                p.GetXaxis().SetTitleSize(0.045)
                p.GetXaxis().SetLabelSize(0.04)
                p.GetXaxis().SetTitleOffset(1.1)
                p.GetXaxis().CenterTitle(True)
                p.GetYaxis().ChangeLabel(-1, -1, 0)

            if is_first:
                p.Draw("HIST")
                is_first = False
            else:
                p.Draw("HIST SAME")

        legend = None
        if not is_bottom_pad:
            legend = ROOT.TLegend(0.70, 0.40, 0.90, 0.90)
            legend.SetTextFont(latex_font)
            legend.SetBorderSize(0)
            legend.SetFillStyle(1001)
            legend.SetFillColor(ROOT.kWhite)
            legend.SetTextSize(0.04)
            for name, p in profs_dict.items():
                legend.AddEntry(p, name, "l")
            legend.Draw()

        latex = ROOT.TLatex()
        latex.SetNDC()
        latex.SetTextFont(133)
        latex.SetTextSize(26)
        x_pos = 0.18
        y_pos = 0.92 if is_bottom_pad else 0.84
        latex.DrawLatex(x_pos, y_pos, pad_title)

        return legend

    c_frac.leg1 = draw_frac_pad(pad1_f, profs_f1, label1, is_bottom_pad=False)
    c_frac.leg2 = draw_frac_pad(pad2_f, profs_f2, label2, is_bottom_pad=True)
    c_frac.Update()
    c_frac.SaveAs(f"{out_name}_energy_fractions.pdf")

    # -------------------------------------------------------------------------
    # 4. CANVAS 4: STACKED E_AVAIL DISTRIBUTION (WEIGHTED BY SPECIES FRACTION)
    # -------------------------------------------------------------------------
    c_eavail, pad1_e, pad2_e = create_top_bottom_canvas("c_eavail_stacked", "Stacked E_avail Breakdown")

    sum_e1 = hists_e1[charged_names[0]].Clone("sum_e1")
    sum_e1.Reset()
    for h in hists_e1.values():
        sum_e1.Add(h)

    sum_e2 = hists_e2[charged_names[0]].Clone("sum_e2")
    sum_e2.Reset()
    for h in hists_e2.values():
        sum_e2.Add(h)

    stack_ymax = max(sum_e1.GetMaximum(), sum_e2.GetMaximum()) * 1.25

    def draw_stacked_eavail_pad(pad, hists_eavail_dict, pad_title, file_tag, is_bottom_pad=False):
        pad.cd()
        ROOT.gPad.SetTickx(1)
        ROOT.gPad.SetTicky(1)

        hstack = ROOT.THStack(f"hs_eavail_{file_tag}", "; E_{avail} [MeV]; Counts")
        for name in charged_names:
            hstack.Add(hists_eavail_dict[name])

        hstack.Draw("HIST")
        hstack.SetMaximum(stack_ymax)
        hstack.SetMinimum(0.0)

        hstack.GetYaxis().SetTitleFont(latex_font)
        hstack.GetYaxis().SetLabelFont(latex_font)
        hstack.GetYaxis().SetTitleSize(0.06)
        hstack.GetYaxis().SetLabelSize(0.05)
        hstack.GetYaxis().SetTitleOffset(1.3)
        hstack.GetYaxis().CenterTitle(True)

        if not is_bottom_pad:
            hstack.GetXaxis().SetLabelSize(0)
            hstack.GetXaxis().SetTitleSize(0)
            hstack.GetXaxis().SetTitle("")
            hstack.GetYaxis().ChangeLabel(1, -1, 0)
        else:
            hstack.GetXaxis().SetNoExponent(True)
            hstack.GetXaxis().SetTitleFont(latex_font)
            hstack.GetXaxis().SetLabelFont(latex_font)
            hstack.GetXaxis().SetTitleSize(0.06)
            hstack.GetXaxis().SetLabelSize(0.05)
            hstack.GetXaxis().SetTitleOffset(1.1)
            hstack.GetXaxis().CenterTitle(True)
            hstack.GetYaxis().ChangeLabel(-1, -1, 0)

        legend = None
        if not is_bottom_pad:
            legend = ROOT.TLegend(0.70, 0.40, 0.90, 0.90)
            legend.SetTextFont(latex_font)
            legend.SetBorderSize(0)
            legend.SetFillStyle(1001)
            legend.SetFillColor(ROOT.kWhite)
            legend.SetTextSize(0.05)
            for name in reversed(charged_names):
                legend.AddEntry(hists_eavail_dict[name], name, "f")
            legend.Draw()

        latex = ROOT.TLatex()
        latex.SetNDC()
        latex.SetTextFont(133)
        latex.SetTextSize(26)
        x_pos = 0.18
        y_pos = 0.92 if is_bottom_pad else 0.84
        latex.DrawLatex(x_pos, y_pos, pad_title)

        return hstack, legend

    hs1, leg1_e = draw_stacked_eavail_pad(pad1_e, hists_e1, label1, "f1", is_bottom_pad=False)
    hs2, leg2_e = draw_stacked_eavail_pad(pad2_e, hists_e2, label2, "f2", is_bottom_pad=True)

    c_eavail.hs1 = hs1
    c_eavail.hs2 = hs2
    c_eavail.leg1 = leg1_e
    c_eavail.Update()
    c_eavail.SaveAs(f"{out_name}_eavail_stacked.pdf")

    # -------------------------------------------------------------------------
    # 5. WRITE ALL CANVASES & HISTOGRAMS TO ROOT FILE
    # -------------------------------------------------------------------------
    out_file = ROOT.TFile(f"{out_name}.root", "RECREATE")
    
    c_mom.Write("c_momentum")
    c_lead_mom.Write("c_leading_momentum")
    c_frac.Write("c_energy_fractions")
    c_eavail.Write("c_eavail_stacked")

    for h in hists_m1.values():
        h.Write()
    for h in hists_m2.values():
        h.Write()

    for h in hists_lm1.values():
        h.Write()
    for h in hists_lm2.values():
        h.Write()

    for p in profs_f1.values():
        p.Write()
    for p in profs_f2.values():
        p.Write()

    for h in hists_e1.values():
        h.Write()
    for h in hists_e2.values():
        h.Write()

    out_file.Close()

    print(f"Saved total momentum plot to {out_name}_momentum.pdf")
    print(f"Saved leading momentum plot to {out_name}_leading_momentum.pdf")
    print(f"Saved fraction profiles plot to {out_name}_energy_fractions.pdf")
    print(f"Saved stacked E_avail plot to {out_name[0]}_eavail_stacked.pdf")
    print(f"Saved all canvases & ROOT objects to {out_name}.root")

    return c_mom, c_lead_mom, c_frac, c_eavail



def ccqe_individual_particle_momentum_plots(filename, out_name, filename_list=None):
    label1 = "NEUT Cascade"
    label2 = "INCL + ABLA" 

    #label1 = "FB"
    #label2 = "abla++"
    
    if filename_list is None:
        filename_list = []
        
    files_to_process = [filename] + filename_list

    # Particle specifications and color palette
    particle_specs = {
        "Proton":   {"pdg": [2212],       "mass": 938.27,  "color": ROOT.TColor.GetColor("#3f90da"), "style": 1},
        "Neutron":  {"pdg": [2112],       "mass": 939.57,  "color": ROOT.TColor.GetColor("#ffa90e"), "style": 1},
        "Deuteron": {"pdg": [1000010020], "mass": 1875.61, "color": ROOT.TColor.GetColor("#bd1f01"), "style": 1},
        "Triton":   {"pdg": [1000010030], "mass": 2808.92, "color": ROOT.TColor.GetColor("#832db6"), "style": 1},
        "He-3":     {"pdg": [1000020030], "mass": 2808.39, "color": ROOT.TColor.GetColor("#92dadd"), "style": 1},
        "Alpha":    {"pdg": [1000020040], "mass": 3727.38, "color": ROOT.TColor.GetColor("#a96b59"), "style": 1},
    }

    # Stack rank definitions (1st highest momentum down to 5th+ highest)
    rank_specs = [
        {"label": "1st (Leading)", "color": ROOT.TColor.GetColor("#1f77b4")},
        {"label": "2nd",           "color": ROOT.TColor.GetColor("#ff7f0e")},
        {"label": "3rd",           "color": ROOT.TColor.GetColor("#2ca02c")},
        {"label": "4th",           "color": ROOT.TColor.GetColor("#d62728")},
        {"label": "5th+",          "color": ROOT.TColor.GetColor("#9467bd")},
    ]
    MAX_RANKS = len(rank_specs)

    pdg_map = {}
    for name, spec in particle_specs.items():
        for pdg in spec["pdg"]:
            pdg_map[pdg] = (name, spec["mass"])

    # -------------------------------------------------------------------------
    # 1. FILE PROCESSING
    # -------------------------------------------------------------------------
    def process_file(fname, file_tag):
        print(f"Processing file: {fname}")
        f = ROOT.TFile.Open(fname)
        if not f or f.IsZombie():
            raise RuntimeError(f"Could not open file: {fname}")
            
        t = f.Get("neuttree")

        # 1a. Histogram for ALL individual particles of each species
        hists_all_mom = {}
        for name, spec in particle_specs.items():
            h = ROOT.TH1D(f"h_all_mom_{name}_{file_tag}", f"; p [MeV/c]; Particle Count", 50, 0.0, 1500.0)
            h.Sumw2()
            h.SetLineColor(spec["color"])
            h.SetLineStyle(spec["style"])
            h.SetLineWidth(2)
            h.SetStats(0)
            hists_all_mom[name] = h

        # 1b. Histograms for 1st, 2nd, 3rd, 4th, 5th+ momentum rank per species
        hists_rank_mom = {name: [] for name in particle_specs}
        for name in particle_specs:
            for r_idx, r_spec in enumerate(rank_specs):
                h = ROOT.TH1D(f"h_rank_{r_idx}_{name}_{file_tag}", f"; p [MeV/c]; Particle Count", 50, 0.0, 1500.0)
                h.Sumw2()
                h.SetFillColor(r_spec["color"])
                h.SetLineColor(ROOT.kBlack)
                h.SetLineWidth(1)
                h.SetStats(0)
                hists_rank_mom[name].append(h)

        event_no = 0
        for event in t:
            nvect = event.vectorbranch
            nvect_class = nvect_reader(nvect) 
            event_no += 1
            
            if event_no % 100000 == 0:
                print(f"[{file_tag}] Event number {event_no}")
            if event_no == 500000:
                break

            if nvect_class.eventType in ( EventType.SRC, EventType.MF, EventType.twop2h ,EventType.piProd, EventType.NC_mf,EventType.NC_SRC):
                particles, energies = nvect_class.particles()
                
                # Store momentum list per particle species for this event
                event_species_moms = {p_name: [] for p_name in particle_specs}

                for particle, energy in zip(particles, energies):
                    if particle in pdg_map:
                        p_name, mass = pdg_map[particle]
                        mom = np.sqrt((energy + mass)**2 - mass**2) if mass > 0 else energy
                        event_species_moms[p_name].append(mom)

                # Fill individual particle histograms and rank histograms
                for p_name, moms in event_species_moms.items():
                    # Fill individual momentum histogram for every particle present
                    for m in moms:
                        hists_all_mom[p_name].Fill(m)

                    # Sort particle momenta descending to determine rank contributions
                    sorted_moms = sorted(moms, reverse=True)
                    for idx, mom_val in enumerate(sorted_moms):
                        rank_idx = min(idx, MAX_RANKS - 1)
                        hists_rank_mom[p_name][rank_idx].Fill(mom_val)

        for h in hists_all_mom.values():
            h.SetDirectory(0)
        for h_list in hists_rank_mom.values():
            for h in h_list:
                h.SetDirectory(0)

        f.Close()
        return hists_all_mom, hists_rank_mom

    hists_all1, hists_rank1 = process_file(files_to_process[0], "f1")
    hists_all2, hists_rank2 = process_file(files_to_process[1], "f2")

    ROOT.gStyle.SetOptTitle(0)
    latex_font = 132

    # Canvas helper
    def create_top_bottom_canvas(canvas_name, title):
        c = ROOT.TCanvas(canvas_name, title, 1000, 1000)
        pad1 = ROOT.TPad(f"{canvas_name}_pad1", "Pad 1 (Top)",    0.0, 0.535, 1.0, 1.0)
        pad2 = ROOT.TPad(f"{canvas_name}_pad2", "Pad 2 (Bottom)", 0.0, 0.0,   1.0, 0.535)

        pad1.SetLeftMargin(0.14)   
        pad1.SetRightMargin(0.05)  
        pad1.SetTopMargin(0.08)
        pad1.SetBottomMargin(0.0)

        pad2.SetLeftMargin(0.14)   
        pad2.SetRightMargin(0.05) 
        pad2.SetTopMargin(0.0)     
        pad2.SetBottomMargin(0.15) 

        pad1.Draw()
        pad2.Draw()
        return c, pad1, pad2

    # -------------------------------------------------------------------------
    # 2. CANVAS 1: ALL INDIVIDUAL PARTICLE MOMENTA (ALL SPECIES OVERLAID)
    # -------------------------------------------------------------------------
    max_a1 = max(h.GetMaximum() for h in hists_all1.values())
    max_a2 = max(h.GetMaximum() for h in hists_all2.values())
    shared_ymax_all = max(max_a1, max_a2) * 1.25
    if shared_ymax_all <= 0:
        shared_ymax_all = 1.0

    c_all_mom, pad1_a, pad2_a = create_top_bottom_canvas("c_all_particle_momentum", "All Particle Momentum")

    def draw_all_mom_pad(pad, hists_dict, pad_title, is_bottom_pad=False):
        pad.cd()
        ROOT.gPad.SetTickx(1)
        ROOT.gPad.SetTicky(1)

        is_first = True
        for name, h in hists_dict.items():
            h.GetXaxis().SetNoExponent(True)
            h.GetYaxis().SetTitleFont(latex_font)
            h.GetYaxis().SetLabelFont(latex_font)
            h.GetYaxis().SetTitleSize(0.08)
            h.GetYaxis().SetLabelSize(0.07)
            h.GetYaxis().SetTitleOffset(0.9)
            h.GetYaxis().CenterTitle(True)

            h.SetMaximum(shared_ymax_all)
            h.SetMinimum(0)

            if not is_bottom_pad:
                h.GetXaxis().SetLabelSize(0)
                h.GetXaxis().SetTitleSize(0)
                h.GetXaxis().SetTitle("")
                h.GetYaxis().ChangeLabel(1, -1, 0)
            else:
                h.GetXaxis().SetTitleFont(latex_font)
                h.GetXaxis().SetLabelFont(latex_font)
                h.GetXaxis().SetTitleSize(0.08)
                h.GetXaxis().SetLabelSize(0.07)
                h.GetXaxis().SetTitleOffset(0.8)
                h.GetXaxis().CenterTitle(True)
                #h.GetYaxis().ChangeLabel(-1, -1, 0)

            if is_first:
                h.Draw("HIST")
                is_first = False
            else:
                h.Draw("HIST SAME")

        legend = None
        if not is_bottom_pad:
            legend = ROOT.TLegend(0.70, 0.40, 0.90, 0.90)
            legend.SetTextFont(latex_font)
            legend.SetBorderSize(0)
            legend.SetFillStyle(1001)
            legend.SetFillColor(ROOT.kWhite)
            legend.SetTextSize(0.06)
            for name, h in hists_dict.items():
                legend.AddEntry(h, name, "l")
            legend.Draw()

        latex = ROOT.TLatex()
        latex.SetNDC()
        latex.SetTextFont(133)
        latex.SetTextSize(26)
        x_pos = 0.18
        y_pos = 0.92 if is_bottom_pad else 0.84
        latex.DrawLatex(x_pos, y_pos, pad_title)

        return legend

    c_all_mom.leg1 = draw_all_mom_pad(pad1_a, hists_all1, label1, is_bottom_pad=False)
    c_all_mom.leg2 = draw_all_mom_pad(pad2_a, hists_all2, label2, is_bottom_pad=True)
    c_all_mom.Update()

    # -------------------------------------------------------------------------
    # 3. CANVAS SET 2: STACKED RANK CONTRIBUTIONS PER SPECIES
    # -------------------------------------------------------------------------
    rank_canvases = {}

    for p_name in particle_specs:
        sum_h1 = hists_all1[p_name]
        sum_h2 = hists_all2[p_name]
        species_ymax = max(sum_h1.GetMaximum(), sum_h2.GetMaximum()) * 1
        if species_ymax <= 0:
            species_ymax = 1.0

        c_spec, pad1_s, pad2_s = create_top_bottom_canvas(f"c_rank_stack_{p_name}", f"Rank Breakdown: {p_name}")

        def draw_rank_stack_pad(pad, rank_hists_list, pad_title, file_tag, is_bottom_pad=False):
            pad.cd()
            ROOT.gPad.SetTickx(1)
            ROOT.gPad.SetTicky(1)

            hstack = ROOT.THStack(f"hs_rank_{p_name}_{file_tag}", f"; p [MeV/c]; Particle Count")
            for h_rank in rank_hists_list:
                hstack.Add(h_rank)

            hstack.Draw("HIST")
            hstack.SetMaximum(species_ymax)
            hstack.SetMinimum(0.0)

            hstack.GetYaxis().SetTitleFont(latex_font)
            hstack.GetYaxis().SetLabelFont(latex_font)
            hstack.GetYaxis().SetTitleSize(0.08)
            hstack.GetYaxis().SetLabelSize(0.07)
            hstack.GetYaxis().SetTitleOffset(0.9)
            hstack.GetYaxis().CenterTitle(True)

            if not is_bottom_pad:
                hstack.GetXaxis().SetLabelSize(0)
                hstack.GetXaxis().SetTitleSize(0)
                hstack.GetXaxis().SetTitle("")
                hstack.GetYaxis().ChangeLabel(1, -1, 0)
            else:
                hstack.GetXaxis().SetNoExponent(True)
                hstack.GetXaxis().SetTitleFont(latex_font)
                hstack.GetXaxis().SetLabelFont(latex_font)
                hstack.GetXaxis().SetTitleSize(0.08)
                hstack.GetXaxis().SetLabelSize(0.07)
                hstack.GetXaxis().SetTitleOffset(0.8)
                hstack.GetXaxis().CenterTitle(True)
                hstack.GetYaxis().ChangeLabel(-1, -1, 0)

            legend = None
            if not is_bottom_pad:
                legend = ROOT.TLegend(0.68, 0.45, 0.90, 0.90)
                legend.SetTextFont(latex_font)
                legend.SetBorderSize(0)
                legend.SetFillStyle(1001)
                legend.SetFillColor(ROOT.kWhite)
                legend.SetTextSize(0.07)
                for r_idx in reversed(range(MAX_RANKS)):
                    legend.AddEntry(rank_hists_list[r_idx], rank_specs[r_idx]["label"], "f")
                legend.Draw()

            latex = ROOT.TLatex()
            latex.SetNDC()
            latex.SetTextFont(133)
            latex.SetTextSize(26)
            x_pos = 0.18
            y_pos = 0.92 if is_bottom_pad else 0.84
            latex.DrawLatex(x_pos, y_pos, f"{p_name} - {pad_title}")

            return hstack, legend

        hs1, leg1_s = draw_rank_stack_pad(pad1_s, hists_rank1[p_name], label1, "f1", is_bottom_pad=False)
        hs2, leg2_s = draw_rank_stack_pad(pad2_s, hists_rank2[p_name], label2, "f2", is_bottom_pad=True)

        c_spec.hs1 = hs1
        c_spec.hs2 = hs2
        c_spec.leg1 = leg1_s
        c_spec.Update()
        rank_canvases[p_name] = c_spec

    # -------------------------------------------------------------------------
    # 4. SAVE ALL CANVASES TO A SINGLE MULTI-PAGE PDF
    # -------------------------------------------------------------------------

    output_dir = "CCQE_stacked_momentum_dist"
    os.makedirs(output_dir, exist_ok=True)
    
    pdf_filename  = os.path.join(output_dir, f"all_particles_{out_name[0]}_all_plots.pdf")

    all_canvases = [c_all_mom] + list(rank_canvases.values())

    # Open multi-page PDF container
    all_canvases[0].Print(f"{pdf_filename}[")

    # Print each canvas as a page into the combined PDF
    for c in all_canvases:
        c.Print(pdf_filename)

    # Close multi-page PDF container
    all_canvases[-1].Print(f"{pdf_filename}]")

    # -------------------------------------------------------------------------
    # 5. WRITE ALL CANVASES & HISTOGRAMS TO ROOT FILE
    # -------------------------------------------------------------------------
    out_file = ROOT.TFile(f"{out_name}_all_particles.root", "RECREATE")
    
    c_all_mom.Write("c_all_particle_momentum")
    for p_name, c_spec in rank_canvases.items():
        c_spec.Write(f"c_rank_stack_{p_name}")

    for h in hists_all1.values():
        h.Write()
    for h in hists_all2.values():
        h.Write()

    for h_list in hists_rank1.values():
        for h in h_list:
            h.Write()
    for h_list in hists_rank2.values():
        for h in h_list:
            h.Write()

    out_file.Close()


    print(f"Saved all {len(all_canvases)} canvas pages to single PDF: {pdf_filename}")
    print(f"Saved all canvases & ROOT objects to ROOT file: {out_name}_all_particles.root")

    return c_all_mom, rank_canvases
def proton_momentum_plots(filename, name, filename_list=None):
    if filename_list is None:
        filename_list = []

    files_to_process = [filename] + filename_list
    
    legend_name = [
                    "NEUT Cascade",
                    #"INCL (NO ABLA)",
                    "INCL + ABLA",        
                ]

    styles = [
                
                (ROOT.kAzure,1),
                (ROOT.kBlue, 1), 
               # (ROOT.kOrange+7,1),
            ]

    """
    legend_name = ["NEUT Cascade - NuWro SRC", 
                   "NEUT Cascade - Box SRC", 
                   "INCL + ABLA  - NuWro SRC",
                   "INCL + ABLA  - Box SRC",
                   ]
    
    styles = [
                (ROOT.kRed, 4),  # Uses index 2 (modified to Petroff Blue)
                (ROOT.kBlue, 4),
                (ROOT.kRed, 1),  
                (ROOT.kBlue, 1),
            ]
    """
    # Lists to store (label, histogram) for each rank
    hist_list_p1 = []
    hist_list_p2 = []
    hist_list_p3 = []
    hist_list_p4plus = []

    m_p = 938.27  # Proton mass in MeV/c^2

    for count, fname in enumerate(files_to_process):
        print(f"Processing file: {fname}")
        chunk_label = (
            legend_name[count]
            if count < len(legend_name)
            else f"File {count + 1}"
        )

        f = ROOT.TFile.Open(fname)
        t = f.Get("neuttree")

        # Rank-based histograms for this file
        h_p1 = ROOT.TH1D(
            f"h_p1_{count}",
            "1^{st} Highest Proton;p_{p} [MeV/c];Counts",
            150,
            0.0,
            1500,
        )
        h_p2 = ROOT.TH1D(
            f"h_p2_{count}",
            "2^{nd} Highest Proton;p_{p} [MeV/c];Counts",
            80,
            0.0,
            800,
        )
        h_p3 = ROOT.TH1D(
            f"h_p3_{count}",
            "3^{rd} Highest Proton;p_{p} [MeV/c];Counts",
            60,
            0.0,
            600,
        )
        h_p4plus = ROOT.TH1D(
            f"h_p4plus_{count}",
            "4^{+} Highest Protons;p_{p} [MeV/c];Counts",
            40,
            0.0,
            400,
        )

        for h in [h_p1, h_p2, h_p3, h_p4plus]:
            h.Sumw2()
            h.SetDirectory(0)

        event_no = 0
        for event in t:
            event_no += 1
            if event_no % 100000 == 0:
                print(f"  Event number {event_no}")
            if event_no == 500000:
                break

            nvect = event.vectorbranch
            nvect_class = nvect_reader(nvect)

            if nvect_class.eventType in (EventType.MF,EventType.twop2h,EventType.piProd):
                particles, energies = nvect_class.particles()

                # Collect and calculate momentum for all protons in the event
                proton_momenta = []
                for particle, energy in zip(particles, energies):
                    if particle == 2212:  # Proton PDG
                        p_mom = np.sqrt((energy + m_p)**2 - m_p**2)
                        proton_momenta.append(p_mom)

                proton_momenta.sort(reverse=True)

                # Fill rank histograms
                if len(proton_momenta) > 0:
                    h_p1.Fill(proton_momenta[0])
                if len(proton_momenta) > 1:
                    h_p2.Fill(proton_momenta[1])
                if len(proton_momenta) > 2:
                    h_p3.Fill(proton_momenta[2])
                if len(proton_momenta) > 3:
                    for p in proton_momenta[3:]:
                        h_p4plus.Fill(p)

        f.Close()

        hist_list_p1.append((chunk_label, h_p1))
        hist_list_p2.append((chunk_label, h_p2))
        hist_list_p3.append((chunk_label, h_p3))
        hist_list_p4plus.append((chunk_label, h_p4plus))

    # --- Plotting on a single canvas with 4 pads side-by-side ---
    ROOT.gStyle.SetOptTitle(0)
    latex_font = 132

    # Aspect ratio set for 4 horizontal pads
    canvas = ROOT.TCanvas(
        "c_proton_rank_momenta", "Proton Momentum Distributions", 2000, 600
    )
    canvas.Divide(4, 1)

    pad_data = [
        ("1^{st} Highest Proton", hist_list_p1),
        ("2^{nd} Highest Proton", hist_list_p2),
        ("3^{rd} Highest Proton", hist_list_p3),
        ("4^{+} Highest Protons", hist_list_p4plus),
    ]

    legends = []  # Keep references to prevent garbage collection

    for pad_idx, (pad_title, h_list) in enumerate(pad_data, start=1):
        pad = canvas.cd(pad_idx)
        pad.SetLeftMargin(0.14)
        pad.SetRightMargin(0.04)
        pad.SetTopMargin(0.10)
        pad.SetBottomMargin(0.14)
        pad.SetTickx(1)
        pad.SetTicky(1)

        # Scale Y axis based on pad maximum
        global_max = max(h.GetMaximum() for _, h in h_list) if h_list else 1.0

        legend = ROOT.TLegend(0.38, 0.65, 0.92, 0.88)
        legend.SetNColumns(1)
        legend.SetTextFont(latex_font)
        legend.SetBorderSize(0)
        legend.SetFillStyle(0)
        legend.SetTextSize(0.035)

        for i, (label, h) in enumerate(h_list):
            color, style = styles[i % len(styles)]
            h.SetLineColor(color)
            h.SetLineStyle(style)
            h.SetLineWidth(3)
            h.SetStats(0)

            # Axis formatting
            h.GetXaxis().SetTitle("p_{p} [MeV/c]")
            h.GetYaxis().SetTitle("Counts")
            h.GetXaxis().CenterTitle(True)
            h.GetYaxis().CenterTitle(True)
            h.GetXaxis().SetNoExponent(True)

            h.GetXaxis().SetTitleFont(latex_font)
            h.GetXaxis().SetLabelFont(latex_font)
            h.GetYaxis().SetTitleFont(latex_font)
            h.GetYaxis().SetLabelFont(latex_font)

            h.GetXaxis().SetTitleSize(0.05)
            h.GetXaxis().SetLabelSize(0.045)
            h.GetYaxis().SetTitleSize(0.05)
            h.GetYaxis().SetLabelSize(0.045)
            h.GetXaxis().SetTitleOffset(1.1)
            h.GetYaxis().SetTitleOffset(1.3)

            if i == 0:
                h.SetMaximum(global_max * 1.25)
                h.SetMinimum(0)
                h.Draw("HIST")
            else:
                h.Draw("HIST SAME")

            legend.AddEntry(h, label, "l")

        legend.Draw()
        legends.append(legend)

        # Pad Header Title
        title_latex = ROOT.TLatex()
        title_latex.SetNDC()
        title_latex.SetTextFont(latex_font)
        title_latex.SetTextSize(0.05)
        title_latex.SetTextAlign(22)  # 2 = Horizontal Center, 2 = Vertical Center
        title_latex.DrawLatex(0.50, 0.94, pad_title)

    canvas.Update()
    canvas.SaveAs(f"proton_momenta_ranked_{name}.pdf")

    out_file = ROOT.TFile(f"proton_momenta_ranked_{name}.root", "RECREATE")
    canvas.Write()
    out_file.Close()

    return canvas
        

def monoenergetic_CCQE_combined_plots():

    """
    Processes monoenergetic ROOT files for energies 250, 500, 750, 1000, and 1250 MeV
    across different model implementations (excluding KDAR files).
    
    Generates a canvas for each variable containing 5 subplots (one per energy),
    overlaying the model predictions in each subplot.
    """
    energies = ["250", "500", "750","1000", "1250"]
    
    # Models and their corresponding filename templates
    models = [
        ("NEUT SRC cut", "out_PAPER_INCL_BOX_MONO{energy}.root"),
        ("NuWro SRC Cut", "out_PAPER_INCL_NUWRO_MONO{energy}.root"),
        ("NEUT Default Cascade", "out_PAPER_NEUTCASC_BOX_MONO{energy}.root")
    ]
    
    styles = [
        (ROOT.kRed, 7),      # Dashed Red
        (ROOT.kBlue, 9),     # Dash-Dot Blue
        (ROOT.kGreen+2, 5),  # Long-Dash Green
        (ROOT.kMagenta, 2),  # Short-Dash Magenta
        (ROOT.kCyan, 1)      # Solid Cyan
    ]

    # Specification for all variables: (Title, Bins, X-min, X-max, X-title)
    var_specs = {
        "leading_proton": ("Leading Proton Momentum", 50, 0.0, 1500, "P_{N} [MeV]"),
        "excitation_E": ("Excitation Energy Distribution", 40, 0.0, 100, "E_{ex} [MeV]"),
        "avail_E": ("Hadronic Available Energy Distribution", 40, 0.0, 400, "E_{Had} [MeV]"),
        "avail_E_noP": ("Hadronic Available Energy (no p)", 40, 0.0, 400, "(E_{Had} - T_{p}) [MeV]"),
        "avail_E_nonINCL": ("Non-INCL Available Energy Distribution", 40, 0.0, 400, "E_{Non-INCL} [MeV]"),
        
        # Energy distributions
        "totp_E": ("Total Proton Energy", 40, 0.0, 400, "E_{p}^{Total} [MeV]"),
        "totgamma_E": ("Total Photon Energy", 100, 0.0, 30, "E_{#gamma}^{Total} [MeV]"),
        "totn_E": ("Total Neutron Energy", 400, 0.0, 400, "E_{n}^{Total} [MeV]"),
        "totd_E": ("Total Deuteron Energy", 100, 0.0, 100, "E_{d}^{Total} [MeV]"),
        "tott_E": ("Total Triton Energy", 100, 0.0, 100, "E_{t}^{Total} [MeV]"),
        "tota_E": ("Total Alpha Energy", 100, 0.0, 100, "E_{#alpha}^{Total} [MeV]"),
        "totpi0_E": ("Total #pi^{0} Energy", 100, 0.0, 100, "E_{#pi^{0}}^{Total} [MeV]"),
        "totpipm_E": ("Total #pi^{#pm} Energy", 100, 0.0, 100, "E_{#pi^{#pm}}^{Total} [MeV]"),
        
        # Multiplicity distributions
        "mult_p": ("Proton Multiplicity", 15, 0.0, 15, "N_{p}"),
        "mult_gamma": ("Photon Multiplicity", 20, 0.0, 20, "N_{#gamma}"),
        "mult_n": ("Neutron Multiplicity", 15, 0.0, 15, "N_{n}"),
        "mult_d": ("Deuteron Multiplicity", 10, 0.0, 10, "N_{d}"),
        "mult_t": ("Triton Multiplicity", 10, 0.0, 10, "N_{t}"),
        "mult_a": ("Alpha Multiplicity", 10, 0.0, 10, "N_{#alpha}"),
        "mult_pi0": ("#pi^{0} Multiplicity", 10, 0.0, 10, "N_{#pi^{0}}"),
        "mult_pipm": ("#pi^{#pm} Multiplicity", 10, 0.0, 10, "N_{#pi^{#pm}}")
    }

    # Nested dictionary to hold TH1D objects: hists[var_key][energy][model_label]
    hists = {vk: {e: {} for e in energies} for vk in var_specs}

    # --- Step 1: Read Files and Fill Histograms ---
    for energy in energies:
        for model_label, fname_pattern in models:
            fname = fname_pattern.format(energy=energy)
            
            f = ROOT.TFile.Open(fname)
            if not f or f.IsZombie():
                print(f"Warning: Could not open {fname}, skipping...")
                continue

            t = f.Get("neuttree")
            if not t:
                f.Close()
                continue

            # Instantiate TH1D histograms for this file
            file_hists = {}
            for vk, (title, nbins, xlow, xhigh, xtitle) in var_specs.items():
                hname = f"h_{vk}_{energy}_{model_label.replace(' ', '_')}"
                h = ROOT.TH1D(hname, title, nbins, xlow, xhigh)
                h.Sumw2()
                h.SetTitle(f";{xtitle};Counts")
                h.SetDirectory(0)
                file_hists[vk] = h

            for event in t:
                nvect = event.vectorbranch
                nvect_class = nvect_reader(nvect) 

                if nvect_class.eventType == EventType.MF or nvect_class.eventType == EventType.SRC:
                    file_hists["excitation_E"].Fill(nvect_class.excitation_E_CCQE())

                    particles, energies_list = nvect_class.particles()
                    file_hists["avail_E"].Fill(sum(energies_list))

                    non_proton_evail = 0
                    non_INCL_eavail = 0 
                    
                    total_proton_energy = 0
                    total_gamma_energy = 0
                    total_neutron_energy = 0
                    total_deuteron_energy = 0
                    total_triton_energy = 0
                    total_alpha_energy = 0
                    total_pi0_energy = 0
                    total_pipm_energy = 0

                    mult_p = 0
                    mult_gamma = 0
                    mult_n = 0
                    mult_d = 0
                    mult_t = 0
                    mult_a = 0
                    mult_pi0 = 0
                    mult_pipm = 0

                    for particle, energy_val in zip(particles, energies_list):
                        if particle != 2212:
                            non_proton_evail += energy_val

                        if particle < 1000000 and particle > 10:
                            non_INCL_eavail += energy_val

                        if particle == 2212:
                            total_proton_energy += energy_val
                            mult_p += 1
                        elif particle == 22:
                            total_gamma_energy += energy_val
                            mult_gamma += 1
                        elif particle == 2112:
                            total_neutron_energy += energy_val
                            mult_n += 1
                        elif particle == 1000010020:
                            total_deuteron_energy += energy_val
                            mult_d += 1
                        elif particle == 1000010030:
                            total_triton_energy += energy_val
                            mult_t += 1
                        elif particle == 1000020040:
                            total_alpha_energy += energy_val
                            mult_a += 1
                        elif particle == 111:
                            total_pi0_energy += energy_val
                            mult_pi0 += 1
                        elif abs(particle) == 211:
                            total_pipm_energy += energy_val
                            mult_pipm += 1

                    # Fill energy histograms
                    if non_proton_evail > 0: file_hists["avail_E_noP"].Fill(non_proton_evail)
                    if non_INCL_eavail > 0: file_hists["avail_E_nonINCL"].Fill(non_INCL_eavail)
                    if total_proton_energy > 0: file_hists["totp_E"].Fill(total_proton_energy)
                    if total_gamma_energy > 0: file_hists["totgamma_E"].Fill(total_gamma_energy)
                    if total_neutron_energy > 0: file_hists["totn_E"].Fill(total_neutron_energy)
                    if total_deuteron_energy > 0: file_hists["totd_E"].Fill(total_deuteron_energy)
                    if total_triton_energy > 0: file_hists["tott_E"].Fill(total_triton_energy)
                    if total_alpha_energy > 0: file_hists["tota_E"].Fill(total_alpha_energy)
                    if total_pi0_energy > 0: file_hists["totpi0_E"].Fill(total_pi0_energy)
                    if total_pipm_energy > 0: file_hists["totpipm_E"].Fill(total_pipm_energy)

                    # Fill multiplicity histograms
                    file_hists["mult_p"].Fill(mult_p)
                    file_hists["mult_gamma"].Fill(mult_gamma)
                    file_hists["mult_n"].Fill(mult_n)
                    file_hists["mult_d"].Fill(mult_d)
                    file_hists["mult_t"].Fill(mult_t)
                    file_hists["mult_a"].Fill(mult_a)
                    file_hists["mult_pi0"].Fill(mult_pi0)
                    file_hists["mult_pipm"].Fill(mult_pipm)

                    if nvect_class.HMPMom:
                        file_hists["leading_proton"].Fill(nvect_class.HMPMom)

            f.Close()

            # Store in master dictionary
            for vk in var_specs:
                hists[vk][energy][model_label] = file_hists[vk]

    # --- Step 2: Draw 5 Subplots per Canvas ---
    ROOT.gStyle.SetOptTitle(0)
    latex_font = 132
    out_file = ROOT.TFile("monoenergetic_combined_comparisons.root", "RECREATE")

    canvases = []

    for vk, (var_title, nbins, xlow, xhigh, xtitle) in var_specs.items():
        c_name = f"c_mono_{vk}"
        c = ROOT.TCanvas(c_name, f"Monoenergetic - {var_title}", 1500, 1000)
        c.Divide(3, 2)  # 6-pad grid: 5 for energies, pad 6 for canvas info

        legends = []
        latexs = []

        for idx, energy in enumerate(energies):
            c.cd(idx + 1)
            ROOT.gPad.SetTickx(1)
            ROOT.gPad.SetTicky(1)
            ROOT.gPad.SetBottomMargin(0.15)
            ROOT.gPad.SetLeftMargin(0.14)

            # Subplot Title (Monoenergetic Energy Label)
            lat = ROOT.TLatex()
            lat.SetNDC()
            lat.SetTextFont(latex_font)
            lat.SetTextSize(0.055)

            legend = ROOT.TLegend(0.50, 0.65, 0.88, 0.85)
            legend.SetTextFont(latex_font)
            legend.SetBorderSize(0)
            legend.SetFillStyle(0)
            legend.SetTextSize(0.038)

            # Find global max across models for the current energy pad
            global_max = 0.0
            energy_hists = hists[vk][energy]
            for model_label, h in energy_hists.items():
                if h.GetMaximum() > global_max:
                    global_max = h.GetMaximum()

            for i, (model_label, h) in enumerate(energy_hists.items()):
                color, style = styles[i % len(styles)]
                h.SetLineColor(color)
                h.SetLineStyle(style)
                h.SetLineWidth(3)
                h.SetStats(0)

                h.GetXaxis().SetTitleFont(latex_font)
                h.GetXaxis().SetLabelFont(latex_font)
                h.GetYaxis().SetTitleFont(latex_font)
                h.GetYaxis().SetLabelFont(latex_font)

                h.GetXaxis().SetTitleSize(0.05)
                h.GetXaxis().SetLabelSize(0.04)
                h.GetXaxis().SetTitleOffset(1.0)
                h.GetYaxis().SetTitleSize(0.05)
                h.GetYaxis().SetLabelSize(0.04)
                h.GetYaxis().SetTitleOffset(1.1)

                if i == 0:
                    h.SetMaximum(global_max * 1.25)
                    h.Draw("HIST")
                else:
                    h.Draw("HIST SAME")

                legend.AddEntry(h, model_label, "l")

            legend.Draw()
            lat.DrawLatex(0.18, 0.82, f"E_{{#nu}} = {energy} MeV")
            
            legends.append(legend)
            latexs.append(lat)

        # 6th pad reserved for variable title header
        c.cd(6)
        lat_title = ROOT.TLatex()
        lat_title.SetNDC()
        lat_title.SetTextFont(latex_font)
        lat_title.SetTextSize(0.065)
        lat_title.DrawLatex(0.1, 0.5, var_title)
        latexs.append(lat_title)

        # Retain references to avoid PyROOT garbage collection
        c.legends = legends
        c.latexs = latexs
        canvases.append(c)

        out_file.cd()
        c.Write()

    out_file.Close()
    print("Done! Results saved to monoenergetic_combined_comparisons.root")


def ccqe_stacked_plots(filename, filename_list=None):
    if filename_list is None:
        filename_list = []
        
    files_to_process = [filename] + filename_list
    n_files = len(files_to_process)
    file_chunks = []

    # Dictionary to hold pair dictionaries per file: {var_name: [{'src': h, 'mf': h}, ...]}
    var_names = [
        'leadin', 'ex', 'eavail', 'notpeavail', 'non_incl_avail', 'missing_kdar',
        'momp', 'momgamma', 'momn', 'momd', 'momt', 'moma', 'mompi0', 'mompipm',
        'mult_p', 'mult_gamma', 'mult_n', 'mult_d', 'mult_t', 'mult_a', 'mult_pi0', 'mult_pipm'
    ]
    file_hists = {v: [] for v in var_names}

    def style_histo(h, color, fill_style=3001):
        h.SetLineColor(color - 2)  
        h.SetLineWidth(3)
        h.SetFillColor(color)
        h.SetFillStyle(fill_style) 
        h.SetStats(0)
        h.Sumw2()
        h.SetDirectory(0)

    # 1. Fill Histograms (SRC vs MF split for each file)
    for fname in files_to_process:
        filename_chunk = fname.split(".")[0].split("out_")[1] if "out_" in fname else fname.split(".")[0]
        file_chunks.append(filename_chunk)

        f = ROOT.TFile(fname)
        t = f.Get("neuttree")   

        # Create SRC and MF histogram pairs for all variables
        h_dict = {}

        # Definition parameters: (nbins, xmin, xmax, title)
        hist_defs = {
            'leadin':          (50, 0.0, 1500, "Leading Proton Distribution;P_{N} [MeV/c];Counts"),
            'ex':              (52, -4.0, 100, "Excitation Energy Distribution;E_{ex} [MeV];Counts"),
            'eavail':          (40, 0.0, 400,  "Hadronic Available Energy Distribution;E_{Had} [MeV];Counts"),
            'notpeavail':      (40, 0.0, 400,  "Hadronic Available Energy (no p);(E_{Had} - T_{p}) [MeV];Counts"),
            'non_incl_avail':  (40, 0.0, 400,  "Non-INCL Available Energy;E_{Non-INCL} [MeV];Counts"),
            'missing_kdar':    (21, -10.0, 95.0, "KDAR Missing Energy;E_{miss}^{KDAR} [MeV];Counts"),
            'momp':            (40, 0.0, 1000, "Total Proton Momentum;P_{p}^{Total} [MeV/c];Counts"),
            'momgamma':        (100, 0.0, 30,  "Total Photon Momentum;P_{#gamma}^{Total} [MeV/c];Counts"),
            'momn':            (100, 0.0, 1000, "Total Neutron Momentum;P_{n}^{Total} [MeV/c];Counts"),
            'momd':            (100, 0.0, 800, "Total Deuteron Momentum;P_{d}^{Total} [MeV/c];Counts"),
            'momt':            (100, 0.0, 800, "Total Triton Momentum;P_{t}^{Total} [MeV/c];Counts"),
            'moma':            (100, 0.0, 1000, "Total Alpha Momentum;P_{#alpha}^{Total} [MeV/c];Counts"),
            'mompi0':          (100, 0.0, 500, "Total #pi^{0} Momentum;P_{#pi^{0}}^{Total} [MeV/c];Counts"),
            'mompipm':         (100, 0.0, 500, "Total #pi^{#pm} Momentum;P_{#pi^{#pm}}^{Total} [MeV/c];Counts"),
            'mult_p':          (15, 0.0, 15,   "Proton Multiplicity;N_{p};Counts"),
            'mult_gamma':      (20, 0.0, 20,   "Photon Multiplicity;N_{#gamma};Counts"),
            'mult_n':          (15, 0.0, 15,   "Neutron Multiplicity;N_{n};Counts"),
            'mult_d':          (10, 0.0, 10,   "Deuteron Multiplicity;N_{d};Counts"),
            'mult_t':          (10, 0.0, 10,   "Triton Multiplicity;N_{t};Counts"),
            'mult_a':          (10, 0.0, 10,   "Alpha Multiplicity;N_{#alpha};Counts"),
            'mult_pi0':        (10, 0.0, 10,   "#pi^{0} Multiplicity;N_{#pi^{0}};Counts"),
            'mult_pipm':       (10, 0.0, 10,   "#pi^{#pm} Multiplicity;N_{#pi^{#pm}};Counts"),
        }

        for v, (nbins, xmin, xmax, title) in hist_defs.items():
            h_src = ROOT.TH1D(f"h_{v}_src_{filename_chunk}", f"SRC {title}", nbins, xmin, xmax)
            h_mf  = ROOT.TH1D(f"h_{v}_mf_{filename_chunk}",  f"MF {title}",  nbins, xmin, xmax)
            
            style_histo(h_src, ROOT.kBlue, 3001)   # Red fill with red-2 line
            style_histo(h_mf, ROOT.kOrange+7, 3001)   # Blue fill with blue-2 line
            
            h_dict[v] = {'src': h_src, 'mf': h_mf, 'title': title}

        # Event Loop
        for event in t:
            nvect = event.vectorbranch
            nvect_class = nvect_reader(nvect) 
            
            evt_type = nvect_class.eventType
            if evt_type not in [EventType.SRC, EventType.MF]:
                continue

            tag = 'src' if evt_type == EventType.SRC else 'mf'

            # Excitation Energy
            if evt_type == EventType.MF:
                e_ex = nvect_class.excitation_E_CCQE()
                h_dict['ex'][tag].Fill(0.0 if e_ex < 0.0 else e_ex)
            elif evt_type == EventType.SRC:
                e_ex = nvect_class.excitation_E_SRC()
                h_dict['ex'][tag].Fill(0.0 if e_ex < 0.0 else e_ex)

            # Missing Energy KDAR
            Missing_energy_kdar, E_Miss = nvect_class.Missing_energy_kdar2()
            if Missing_energy_kdar:
                h_dict['missing_kdar'][tag].Fill(E_Miss)

            # Particles analysis
            particles, energies = nvect_class.particles()
            eavail = sum([e for p, e in zip(particles, energies) if p != 2112])
            h_dict['eavail'][tag].Fill(eavail)

            non_proton_eavail = 0
            non_INCL_eavail = 0 
            
            moms = {'p': 0.0, 'gamma': 0.0, 'n': 0.0, 'd': 0.0, 't': 0.0, 'a': 0.0, 'pi0': 0.0, 'pipm': 0.0}
            mults = {'p': 0, 'gamma': 0, 'n': 0, 'd': 0, 't': 0, 'a': 0, 'pi0': 0, 'pipm': 0}

            for particle, energy in zip(particles, energies):
                if particle not in [2212, 2112]:
                    non_proton_eavail += energy

                if (10 < particle < 1000000) and particle != 2112:
                    non_INCL_eavail += energy
                
                # Particle kinematics
                if particle == 2212: # Proton
                    moms['p'] += np.sqrt((energy+938.27)**2 - 938.27**2)
                    mults['p'] += 1
                elif particle == 22: # Gamma
                    moms['gamma'] += energy
                    mults['gamma'] += 1
                elif particle == 2112: # Neutron
                    moms['n'] += np.sqrt((energy+939.57)**2 - 939.57**2)
                    mults['n'] += 1
                elif particle == 1000010020: # Deuteron
                    moms['d'] += np.sqrt((energy+1875.61)**2 - 1875.61**2)
                    mults['d'] += 1
                elif particle == 1000010030: # Triton
                    moms['t'] += np.sqrt((energy+2808.92)**2 - 2808.92**2)
                    mults['t'] += 1
                elif particle == 1000020040: # Alpha
                    moms['a'] += np.sqrt((energy+3727.38)**2 - 3727.38**2) 
                    mults['a'] += 1
                elif particle == 111: # Pi0
                    moms['pi0'] += np.sqrt((energy+134.98)**2 - 134.98**2)
                    mults['pi0'] += 1
                elif abs(particle) == 211: # Pi+-
                    moms['pipm'] += np.sqrt((energy+139.57)**2 - 139.57**2)
                    mults['pipm'] += 1

            if non_proton_eavail > 0: h_dict['notpeavail'][tag].Fill(non_proton_eavail)
            if non_INCL_eavail > 0:   h_dict['non_incl_avail'][tag].Fill(non_INCL_eavail)

            for key in moms:
                if moms[key] > 0:
                    h_dict[f'mom{key}'][tag].Fill(moms[key])
                h_dict[f'mult_{key}'][tag].Fill(mults[key])

            if getattr(nvect_class, 'HMPMom', None):
                h_dict['leadin'][tag].Fill(nvect_class.HMPMom)

        for v in var_names:
            file_hists[v].append(h_dict[v])

        f.Close()

    # 2. Plotting Phase
    ROOT.gStyle.SetOptTitle(1)
    out_file = ROOT.TFile("paper_combined_comparisons_stacked.root", "RECREATE")
    latex_font = 132

    # Helper function to render standard 1-row stacked canvases across files
    def draw_single_row_stacked(canvas_name, var_key):
        c = ROOT.TCanvas(canvas_name, canvas_name, 600 * n_files, 500)
        c.Divide(n_files, 1)
        c.stacks = []
        c.legends = []

        stacks = []
        for i, chunk in enumerate(file_chunks):
            pair = file_hists[var_key][i]
            hs = ROOT.THStack(f"hs_{var_key}_{chunk}", f"{chunk} {pair['title']}")
            hs.Add(pair['mf'])
            hs.Add(pair['src'])
            stacks.append(hs)

        max_y = max([hs.GetMaximum() for hs in stacks]) * 1.25 if stacks else 1.0

        for i, chunk in enumerate(file_chunks):
            pad = c.cd(i + 1)
            pad.SetTickx(1)
            pad.SetTicky(1)
            pad.SetLeftMargin(0.14)
            pad.SetBottomMargin(0.15)

            hs = stacks[i]
            hs.SetMaximum(max_y)
            hs.Draw("HIST")
            c.stacks.append(hs)

            if i == 0:
                leg = ROOT.TLegend(0.60, 0.72, 0.88, 0.88)
                leg.SetBorderSize(0)
                leg.SetFillStyle(0)
                leg.SetTextFont(latex_font)
                leg.AddEntry(file_hists[var_key][i]['src'], "SRC", "f")
                leg.AddEntry(file_hists[var_key][i]['mf'], "Mean Field", "f")
                leg.Draw()
                c.legends.append(leg)

        return c

    # Helper function to render side-by-side Momentum (Top) + Multiplicity (Bottom)
    def draw_side_by_side_stacked(canvas_name, p_code):
        c = ROOT.TCanvas(canvas_name, canvas_name, 600 * n_files, 800)
        c.Divide(n_files, 2)
        c.stacks = []
        c.legends = []

        mom_key = f"mom{p_code}"
        mult_key = f"mult_{p_code}"

        mom_stacks = []
        mult_stacks = []

        for i, chunk in enumerate(file_chunks):
            pair_mom = file_hists[mom_key][i]
            pair_mult = file_hists[mult_key][i]

            hs_mom = ROOT.THStack(f"hs_{mom_key}_{chunk}", f"{chunk} {pair_mom['title']}")
            hs_mom.Add(pair_mom['mf'])
            hs_mom.Add(pair_mom['src'])
            mom_stacks.append(hs_mom)

            hs_mult = ROOT.THStack(f"hs_{mult_key}_{chunk}", f"{chunk} {pair_mult['title']}")
            hs_mult.Add(pair_mult['mf'])
            hs_mult.Add(pair_mult['src'])
            mult_stacks.append(hs_mult)

        max_mom = max([hs.GetMaximum() for hs in mom_stacks]) * 1.25 if mom_stacks else 1.0
        max_mult = max([hs.GetMaximum() for hs in mult_stacks]) * 1.25 if mult_stacks else 1.0

        for i, chunk in enumerate(file_chunks):
            # Top row: Momentum
            pad_mom = c.cd(i + 1)
            pad_mom.SetTickx(1)
            pad_mom.SetTicky(1)
            pad_mom.SetLeftMargin(0.14)
            pad_mom.SetBottomMargin(0.15)

            hs_mom = mom_stacks[i]
            hs_mom.SetMaximum(max_mom)
            hs_mom.Draw("HIST")
            c.stacks.append(hs_mom)

            if i == 0:
                leg1 = ROOT.TLegend(0.60, 0.72, 0.88, 0.88)
                leg1.SetBorderSize(0)
                leg1.SetFillStyle(0)
                leg1.SetTextFont(latex_font)
                leg1.AddEntry(file_hists[mom_key][i]['src'], "SRC", "f")
                leg1.AddEntry(file_hists[mom_key][i]['mf'], "Mean Field", "f")
                leg1.Draw()
                c.legends.append(leg1)

            # Bottom row: Multiplicity
            pad_mult = c.cd(i + 1 + n_files)
            pad_mult.SetTickx(1)
            pad_mult.SetTicky(1)
            pad_mult.SetLeftMargin(0.14)
            pad_mult.SetBottomMargin(0.15)

            hs_mult = mult_stacks[i]
            hs_mult.SetMaximum(max_mult)
            hs_mult.Draw("HIST")
            c.stacks.append(hs_mult)

        return c

    # Create & Save Canvases
    canvases = []

    # Single row event-level variables
    canvases.append(draw_single_row_stacked("c_leading_proton", "leadin"))
    canvases.append(draw_single_row_stacked("c_excitation_energy", "ex"))
    canvases.append(draw_single_row_stacked("c_avail_energy", "eavail"))
    canvases.append(draw_single_row_stacked("c_nop_avail_energy", "notpeavail"))
    canvases.append(draw_single_row_stacked("c_non_incl_avail_energy", "non_incl_avail"))
    canvases.append(draw_single_row_stacked("c_missing_energy_kdar", "missing_kdar"))

    # Particle-specific momentum + multiplicity stacked canvases
    particles = ['p', 'gamma', 'n', 'd', 't', 'a', 'pi0', 'pipm']
    for p in particles:
        canvases.append(draw_side_by_side_stacked(f"c_tot_{p}_mom_and_mult", p))

    # Write out to ROOT file
    out_file.cd()
    for c in canvases:
        c.Write()
    out_file.Close()

    print("All stacked canvases successfully written to 'paper_combined_comparisons_stacked.root'")
    return canvases



def e_reco_bias_combined_plots(filename, name, n_events, filename_list=None):
    if filename_list is None:
        filename_list = []

    files_to_process = [filename] + filename_list

    # Output histogram containers for the 3 cases
    hist_list_vis0_neu0 = []  # visible=False, neutron=False
    hist_list_vis1_neu0 = []  # visible=True,  neutron=False
    hist_list_vis1_neu1 = []  # visible=True,  neutron=True

    styles = [
        (ROOT.kRed, 7),  # Dashed Red
        (ROOT.kBlue, 9),  # Dash-Dot Blue
        (ROOT.kGreen + 2, 5),  # Long-Dash Green
        (ROOT.kMagenta, 2),  # Short-Dash Magenta
        (ROOT.kCyan, 1),  # Solid Cyan (fallback)
    ]
    legend_name = ["No Subprimary Nucleon", "NEUT SRC Cut", "NuWro SRC Cut"]
    count = 0

    for i, fname in enumerate(files_to_process):
        filename_chunk = legend_name[count]
        count += 1

        f = ROOT.TFile(fname)
        t = f.Get("neuttree")

        # Histograms for energy bias: (E_reco - E_nu) / E_nu
        h_bias_vis0_neu0 = ROOT.TH1D(
            f"h_bias_vis0_neu0_{filename_chunk}",
            "Energy Reco Bias (Invis, No Neutrons)",
            32,
            -0.6,
            0.2,
        )
        h_bias_vis0_neu0.Sumw2()
        h_bias_vis0_neu0.SetTitle(
            ";(E_{reco} - E_{#nu}) / E_{#nu};Counts"
        )

        h_bias_vis1_neu0 = ROOT.TH1D(
            f"h_bias_vis1_neu0_{filename_chunk}",
            "Energy Reco Bias (Vis, No Neutrons)",
            32,
            -0.6,
            0.2,
        )
        h_bias_vis1_neu0.Sumw2()
        h_bias_vis1_neu0.SetTitle(
            ";(E_{reco} - E_{#nu}) / E_{#nu};Counts"
        )

        h_bias_vis1_neu1 = ROOT.TH1D(
            f"h_bias_vis1_neu1_{filename_chunk}",
            "Energy Reco Bias (Vis, With Neutrons)",
            32,
            -0.6,
            0.2,
        )
        h_bias_vis1_neu1.Sumw2()
        h_bias_vis1_neu1.SetTitle(
            ";(E_{reco} - E_{#nu}) / E_{#nu};Counts"
        )

        for event in t:
            nvect = event.vectorbranch
            nvect_class = nvect_reader(nvect)

            if (
                nvect_class.eventType == EventType.SRC
                or nvect_class.eventType == EventType.MF
            ):
                # Calculate bias for the 3 conditions
                bias_vis0_neu0 = nvect_class.E_reco_bias(
                    visible=False, neutron=False
                )
                bias_vis1_neu0 = nvect_class.E_reco_bias(
                    visible=True, neutron=False
                )
                bias_vis1_neu1 = nvect_class.E_reco_bias(
                    visible=True, neutron=True
                )

                h_bias_vis0_neu0.Fill(bias_vis0_neu0)
                h_bias_vis1_neu0.Fill(bias_vis1_neu0)
                h_bias_vis1_neu1.Fill(bias_vis1_neu1)

        # Detach histograms from file
        h_bias_vis0_neu0.SetDirectory(0)
        h_bias_vis1_neu0.SetDirectory(0)
        h_bias_vis1_neu1.SetDirectory(0)

        hist_list_vis0_neu0.append((filename_chunk, h_bias_vis0_neu0))
        hist_list_vis1_neu0.append((filename_chunk, h_bias_vis1_neu0))
        hist_list_vis1_neu1.append((filename_chunk, h_bias_vis1_neu1))

        f.Close()

    ROOT.gStyle.SetOptTitle(0)
    out_file = ROOT.TFile(
        f"paper_e_reco_bias_{name}.root", "RECREATE"
    )
    latex_font = 132

    # Canvas helper function with ratio pad
    def draw_canvas_with_ratio(canvas_name, hist_list, ymax=None, ymin=None):
        c = ROOT.TCanvas(canvas_name, canvas_name, 800, 800)

        pad1 = ROOT.TPad("pad1", "pad1", 0, 0.3, 1, 1.0)
        pad1.SetBottomMargin(0.02)
        pad1.SetLeftMargin(0.12)
        pad1.SetTickx(1)
        pad1.SetTicky(1)
        pad1.Draw()
        pad1.cd()

        legend = ROOT.TLegend(0.60, 0.65, 0.85, 0.85)
        legend.SetNColumns(1)
        legend.SetTextFont(latex_font)
        legend.SetBorderSize(0)
        legend.SetFillStyle(0)
        legend.SetTextSize(0.04)

        if ymax is None:
            global_max = max(h.GetMaximum() for _, h in hist_list)
            ymax = global_max * 1.2

        xaxis_title = hist_list[0][1].GetXaxis().GetTitle()

        for i, (name, h) in enumerate(hist_list):
            color, style = styles[i % len(styles)]

            h.SetLineColor(color)
            h.SetLineStyle(style)
            h.SetLineWidth(3)
            h.SetStats(0)

            h.GetXaxis().SetLabelSize(0)
            h.GetXaxis().SetTitleSize(0)

            h.GetYaxis().SetTitleFont(latex_font)
            h.GetYaxis().SetLabelFont(latex_font)
            h.GetYaxis().SetTitleSize(0.05)
            h.GetYaxis().SetLabelSize(0.04)
            h.GetYaxis().SetTitleOffset(1.1)

            if i == 0:
                h.SetMaximum(ymax)
                if ymin is not None:
                    h.SetMinimum(ymin)
                h.Draw("HIST")
            else:
                h.Draw("HIST SAME")

            legend.AddEntry(h, name, "l")

        legend.Draw()

        c.cd()
        pad2 = ROOT.TPad("pad2", "pad2", 0, 0.0, 1, 0.3)
        pad2.SetTopMargin(0.02)
        pad2.SetBottomMargin(0.3)
        pad2.SetLeftMargin(0.12)
        pad2.SetTickx(1)
        pad2.SetTicky(1)
        pad2.Draw()
        pad2.cd()

        h_base = hist_list[0][1]
        ratios = []
        for i in range(1, len(hist_list)):
            name, h = hist_list[i]
            h_ratio = h.Clone(f"{h.GetName()}_ratio")
            h_ratio.Divide(h_base)

            color, style = styles[i % len(styles)]
            h_ratio.SetLineColor(color)
            h_ratio.SetLineStyle(style)
            h_ratio.SetLineWidth(3)
            h_ratio.SetStats(0)

            h_ratio.SetTitle("")
            h_ratio.GetYaxis().SetTitle("Ratio to ", legend_name[0])
            h_ratio.GetYaxis().CenterTitle()
            h_ratio.GetYaxis().SetNdivisions(505)
            h_ratio.GetYaxis().SetTitleSize(0.11)
            h_ratio.GetYaxis().SetTitleFont(latex_font)
            h_ratio.GetYaxis().SetTitleOffset(0.45)
            h_ratio.GetYaxis().SetLabelFont(latex_font)
            h_ratio.GetYaxis().SetLabelSize(0.1)

            h_ratio.GetXaxis().SetTitle(xaxis_title)
            h_ratio.GetXaxis().SetTitleSize(0.13)
            h_ratio.GetXaxis().SetTitleFont(latex_font)
            h_ratio.GetXaxis().SetTitleOffset(1.0)
            h_ratio.GetXaxis().SetLabelFont(latex_font)
            h_ratio.GetXaxis().SetLabelSize(0.12)

            if i == 1:
                h_ratio.SetMinimum(0.0)
                h_ratio.SetMaximum(2.0)
                h_ratio.Draw("HIST")
            else:
                h_ratio.Draw("HIST SAME")
            ratios.append(h_ratio)

        line = ROOT.TLine(
            h_base.GetXaxis().GetXmin(), 1, h_base.GetXaxis().GetXmax(), 1
        )
        line.SetLineColor(ROOT.kBlack)
        line.SetLineStyle(2)
        line.Draw("SAME")

        c.pad1 = pad1
        c.pad2 = pad2
        c.ratios = ratios
        c.line = line

        return c, legend

    # Draw 3 separate canvases for each configuration
    c_vis0_neu0, leg1 = draw_canvas_with_ratio(
        "c_e_reco_bias_vis0_neu0", hist_list_vis0_neu0, ymax=200000
    )
    c_vis1_neu0, leg2 = draw_canvas_with_ratio(
        "c_e_reco_bias_vis1_neu0", hist_list_vis1_neu0, ymax=200000
    )
    c_vis1_neu1, leg3 = draw_canvas_with_ratio(
        "c_e_reco_bias_vis1_neu1", hist_list_vis1_neu1, ymax=200000
    )

    out_file.cd()
    c_vis0_neu0.Write()
    c_vis1_neu0.Write()
    c_vis1_neu1.Write()

    out_file.Close()
def e_avail_2d_combined_plots(
    filename, name, filename_list=None, max_events=500000
):
    if filename_list is None:
        filename_list = []

    files_to_process = [filename] + filename_list
    num_files = len(files_to_process)
    legend_name = ["No Subprimary Nucleon", "NEUT SRC Cut", "NuWro SRC Cut"]

    settings = [
        (
            False,
            False,
            "vis0_tot0",
            "Available Energy",
        ),
        (
            True,
            False,
            "vis1_tot0",
            "Available Energy W/ Fragments",
        ),
        (
            True,
            True,
            "vis1_tot1",
            "Total Energy W/O Neutrinos",
        ),
    ]

    hists_2d = {s[2]: [] for s in settings}

    # Set ROOT global styles to remove titles/stats and apply palette
    ROOT.gStyle.SetOptStat(0)
    #ROOT.gStyle.SetOptTitle(0)
   # ROOT.gStyle.SetPalette(ROOT.kBlueRedYellow)
    ROOT.gStyle.SetNumberContours(99)

    count = 0
    for fname in files_to_process:
        file_label = (
            legend_name[count] if count < len(legend_name) else f"File_{count}"
        )
        count += 1

        f = ROOT.TFile(fname)
        t = f.Get("neuttree")
        total_events = t.GetEntries()

        print(
            f"--> Processing file: {fname} [{file_label}] ({total_events} total events)"
        )

        h_vis0_tot0 = ROOT.TH2D(
            f"h2d_vis0_tot0_{count}",
            ";T_{p} [MeV];E_{avail} [MeV]",
            100,
            0.0,
            1000.0,
            100,
            0.0,
            1000.0,
        )
        h_vis1_tot0 = ROOT.TH2D(
            f"h2d_vis1_tot0_{count}",
            ";T_{p} [MeV];E_{avail} [MeV]",
            100,
            0.0,
            1000.0,
            100,
            0.0,
            1000.0,
        )
        h_vis1_tot1 = ROOT.TH2D(
            f"h2d_vis1_tot1_{count}",
            ";T_{p} [MeV];E_{avail} [MeV]",
            100,
            0.0,
            1000.0,
            100,
            0.0,
            1000.0,
        )

        for i_evt, event in enumerate(t):
            if max_events is not None and i_evt >= max_events:
                print(
                    f"    [{file_label}] Reached {max_events} event limit. Stopping loop."
                )
                break

            if (i_evt + 1) % (max_events//10) == 0:
                print(
                    f"    [{file_label}] Processed {i_evt + 1}/{min(total_events, max_events)} events..."
                )

            nvect = event.vectorbranch
            nvect_class = nvect_reader(nvect)

            if (
                nvect_class.eventType == EventType.SRC
                or nvect_class.eventType == EventType.MF
            ):
                e_reco1, p_e1 = nvect_class.E_avail(visible=False, total=False)
                e_reco2, p_e2 = nvect_class.E_avail(visible=True, total=False)
                e_reco3, p_e3 = nvect_class.E_avail(visible=True, total=True)

                h_vis0_tot0.Fill(p_e1, e_reco1)
                h_vis1_tot0.Fill(p_e2, e_reco2)
                h_vis1_tot1.Fill(p_e3, e_reco3)

        h_vis0_tot0.SetDirectory(0)
        h_vis1_tot0.SetDirectory(0)
        h_vis1_tot1.SetDirectory(0)

        hists_2d["vis0_tot0"].append((file_label, h_vis0_tot0))
        hists_2d["vis1_tot0"].append((file_label, h_vis1_tot0))
        hists_2d["vis1_tot1"].append((file_label, h_vis1_tot1))

        f.Close()

    latex_font = 132

    def format_axes(hist):
        for axis in [hist.GetXaxis(), hist.GetYaxis()]:
            axis.SetLabelColor(ROOT.kBlack)
            axis.SetTitleColor(ROOT.kBlack)
            axis.SetAxisColor(ROOT.kWhite)

    out_file = ROOT.TFile(f"paper_e_avail_2d_{name}.root", "RECREATE")

    for vis, tot, key, setting_title in settings:
        canvas_name = f"c_2d_e_avail_{key}"
        #ROOT.gStyle.SetPalette(ROOT.kBlueRedYellow)
        c = ROOT.TCanvas(canvas_name, canvas_name, 550 * num_files, 500)
        c.Divide(num_files, 1)
       

        file_hists = hists_2d[key]
        for pad_idx, (file_label, h2d) in enumerate(file_hists, start=1):
            pad = c.cd(pad_idx)
            pad.SetFillColor(ROOT.kWhite)
            pad.SetFrameFillColor(ROOT.kBlack)
            pad.SetTickx(1)
            pad.SetTicky(1)
            pad.SetLeftMargin(0.14)
            pad.SetRightMargin(0.15)
            pad.SetBottomMargin(0.14)

            # Strip statistics boxes and titles
            h2d.SetStats(0)
            h2d.SetTitle(file_label)
            #h2d.SetMinimum(1.0)
            

            # Set a non-zero minimum threshold so ROOT scales the palette properly across all pads
            #h2d.SetMinimum(0.00001)

            h2d.GetXaxis().SetTitleFont(latex_font)
            h2d.GetXaxis().SetLabelFont(latex_font)
            h2d.GetYaxis().SetTitleFont(latex_font)
            h2d.GetYaxis().SetLabelFont(latex_font)
            h2d.GetZaxis().SetTitleFont(latex_font)
            h2d.GetZaxis().SetLabelFont(latex_font)

            h2d.GetXaxis().SetTitleSize(0.045)
            h2d.GetYaxis().SetTitleSize(0.045)
            h2d.GetXaxis().SetLabelSize(0.035)
            h2d.GetYaxis().SetLabelSize(0.035)


  

            format_axes(h2d)
            #ROOT.gStyle.SetPalette(ROOT.kBlueRedYellow)
            h2d.Draw("COLZ")
            pad.Modified()
            pad.Update()

        out_file.cd()
        c.Write()

    out_file.Close()


def neutrons_per_muon_mom(filename, name, filename_list=None, max_events=500000):

    if filename_list is None:
        filename_list = []

    files_to_process = [filename] + filename_list
    legend_names = ["No Subprimary Nucleon", "NEUT SRC Cut", "NuWro SRC Cut"]
    colors = [ROOT.kBlack, ROOT.kRed + 1, ROOT.kBlue + 1]
    latex_font = 132

    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(0)

    hists_p = []
    hists_costheta = []

    for count, fname in enumerate(files_to_process):
        file_label = (
            legend_names[count]
            if count < len(legend_names)
            else f"File_{count}"
        )

        f = ROOT.TFile(fname)
        t = f.Get("neuttree")
        total_events = t.GetEntries()

        print(
            f"--> Processing file: {fname} [{file_label}] ({total_events} total events)"
        )

        h_p = ROOT.TH1D(
            f"h_muon_p_{count}",
            ";p_{#mu} [MeV/c];Events",
            100,
            0.0,
            5000.0,
        )
        h_p.SetDirectory(0)

        h_costheta = ROOT.TH1D(
            f"h_cos_theta_{count}",
            ";cos(#theta_{#mu});Events",
            50,
            -1.0,
            1.0,
        )
        h_costheta.SetDirectory(0)

        for i_evt, event in enumerate(t):
            if max_events is not None and i_evt >= max_events:
                print(
                    f"    [{file_label}] Reached {max_events} event limit. Stopping loop."
                )
                break

            if max_events and (i_evt + 1) % (max_events // 10) == 0:
                print(
                    f"    [{file_label}] Processed {i_evt + 1}/{min(total_events, max_events)} events..."
                )

            nvect = event.vectorbranch
            nvect_class = nvect_reader(nvect)
            muon_p, cos_theta = nvect_class.Neutron_per_muon()
            
            if muon_p is not None:
                h_p.Fill(muon_p)
            if cos_theta is not None:
                h_costheta.Fill(cos_theta)

        f.Close()
        hists_p.append((file_label, h_p))
        hists_costheta.append((file_label, h_costheta))

    # Helper function to construct split-pad canvas with ratio plot
    def build_ratio_canvas(canvas_name, canvas_title, hists_list, x_axis_title):
        c = ROOT.TCanvas(canvas_name, canvas_title, 700, 800)

        pad1 = ROOT.TPad(f"pad1_{canvas_name}", "pad1", 0.0, 0.3, 1.0, 1.0)
        pad1.SetLeftMargin(0.14)
        pad1.SetRightMargin(0.06)
        pad1.SetTopMargin(0.08)
        pad1.SetBottomMargin(0.02)
        pad1.Draw()

        pad2 = ROOT.TPad(f"pad2_{canvas_name}", "pad2", 0.0, 0.0, 1.0, 0.3)
        pad2.SetLeftMargin(0.14)
        pad2.SetRightMargin(0.06)
        pad2.SetTopMargin(0.02)
        pad2.SetBottomMargin(0.35)
        pad2.SetGridy()
        pad2.Draw()

        leg = ROOT.TLegend(0.55, 0.68, 0.88, 0.88)
        leg.SetBorderSize(0)
        leg.SetFillStyle(0)
        leg.SetTextFont(latex_font)
        leg.SetTextSize(0.035)

        max_y = max(h.GetMaximum() for _, h in hists_list) if hists_list else 1.0

        # Draw Main Histograms (Top Pad)
        pad1.cd()
        for idx, (file_label, h) in enumerate(hists_list):
            color = colors[idx % len(colors)]
            h.SetLineColor(color)
            h.SetLineWidth(2)

            h.GetYaxis().SetTitleFont(latex_font)
            h.GetYaxis().SetLabelFont(latex_font)
            h.GetYaxis().SetTitleSize(0.045)
            h.GetYaxis().SetLabelSize(0.035)

            h.GetXaxis().SetLabelSize(0)
            h.GetXaxis().SetTitle("")

            h.SetMaximum(max_y * 1.15)

            draw_option = "HIST" if idx == 0 else "HIST SAME"
            h.Draw(draw_option)
            leg.AddEntry(h, file_label, "l")

        leg.Draw()

        # Draw Ratio Histograms (Bottom Pad)
        pad2.cd()
        ratio_hists = []
        ref_h = hists_list[0][1] if hists_list else None

        if ref_h:
            for idx, (file_label, h) in enumerate(hists_list):
                h_ratio = h.Clone(f"{h.GetName()}_ratio")
                h_ratio.Divide(ref_h)

                color = colors[idx % len(colors)]
                h_ratio.SetLineColor(color)
                h_ratio.SetLineWidth(2)

                h_ratio.GetXaxis().SetTitle(x_axis_title)
                h_ratio.GetYaxis().SetTitle("Ratio")

                h_ratio.GetXaxis().SetTitleFont(latex_font)
                h_ratio.GetYaxis().SetTitleFont(latex_font)
                h_ratio.GetXaxis().SetLabelFont(latex_font)
                h_ratio.GetYaxis().SetLabelFont(latex_font)

                h_ratio.GetXaxis().SetTitleSize(0.11)
                h_ratio.GetXaxis().SetLabelSize(0.09)
                h_ratio.GetXaxis().SetTitleOffset(1.2)
                h_ratio.GetYaxis().SetTitleSize(0.10)
                h_ratio.GetYaxis().SetLabelSize(0.08)
                h_ratio.GetYaxis().SetTitleOffset(0.55)
                h_ratio.GetYaxis().SetNdivisions(505)

                h_ratio.SetMinimum(0.5)
                h_ratio.SetMaximum(1.5)

                draw_option = "HIST" if idx == 0 else "HIST SAME"
                h_ratio.Draw(draw_option)
                ratio_hists.append(h_ratio)

        c.Update()
        return c, leg, ratio_hists, (pad1, pad2)

    # Build canvases for both variables
    c_p, leg_p, ratios_p, pads_p = build_ratio_canvas(
        f"c_muon_p_{name}", f"Muon Momentum - {name}", hists_p, "p_{#mu} [MeV/c]"
    )
    c_costheta, leg_costheta, ratios_costheta, pads_costheta = build_ratio_canvas(
        f"c_cos_theta_{name}", f"Cos(Theta) - {name}", hists_costheta, "cos(#theta_{#mu})"
    )

    # Save both canvases into output file
    out_file = ROOT.TFile(f"paper_neutron_per_muon_{name}.root", "RECREATE")
    c_p.Write()
    c_costheta.Write()
    out_file.Close()

    return {
        "momentum": (c_p, hists_p, leg_p, ratios_p),
        "costheta": (c_costheta, hists_costheta, leg_costheta, ratios_costheta),
    }



def neutrons_per_neutrino_E(filename, name, filename_list=None, max_events=500000):

    if filename_list is None:
        filename_list = []

    files_to_process = [filename] + filename_list
    legend_names = ["No Subprimary Nucleon", "NEUT SRC Cut", "NuWro SRC Cut"]
    colors = [ROOT.kBlack, ROOT.kRed + 1, ROOT.kBlue + 1]
    latex_font = 132

    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(0)

    hists_mult = []

    for count, fname in enumerate(files_to_process):
        file_label = (
            legend_names[count]
            if count < len(legend_names)
            else f"File_{count}"
        )

        f = ROOT.TFile(fname)
        t = f.Get("neuttree")
        total_events = t.GetEntries()

        print(
            f"--> Processing file: {fname} [{file_label}] ({total_events} total events)"
        )

        # Create a TProfile: 6 bins from 0 to 1200 MeV (200 MeV per bin)
        h_mult = ROOT.TProfile(
            f"h_avg_mult_E_{count}",
            ";E_{#nu} [MeV];Average Neutron Multiplicity",
            6,
            0.0,
            1200.0,
        )
        h_mult.SetDirectory(0)

        for i_evt, event in enumerate(t):
            if max_events is not None and i_evt >= max_events:
                print(
                    f"    [{file_label}] Reached {max_events} event limit. Stopping loop."
                )
                break

            if max_events and (i_evt + 1) % (max_events // 10) == 0:
                print(
                    f"    [{file_label}] Processed {i_evt + 1}/{min(total_events, max_events)} events..."
                )

            nvect = event.vectorbranch
            nvect_class = nvect_reader(nvect)
            
            neutrino_E, multiplicity = nvect_class.neutron_multiplicity()
       
            if neutrino_E is not None:
                h_mult.Fill(neutrino_E, multiplicity)

        f.Close()
        hists_mult.append((file_label, h_mult))


    # Helper function to construct split-pad canvas with ratio plot
    def build_ratio_canvas(canvas_name, canvas_title, hists_list, x_axis_title):
        c = ROOT.TCanvas(canvas_name, canvas_title, 700, 800)

        pad1 = ROOT.TPad(f"pad1_{canvas_name}", "pad1", 0.0, 0.3, 1.0, 1.0)
        pad1.SetLeftMargin(0.14)
        pad1.SetRightMargin(0.06)
        pad1.SetTopMargin(0.08)
        pad1.SetBottomMargin(0.02)
        pad1.Draw()

        pad2 = ROOT.TPad(f"pad2_{canvas_name}", "pad2", 0.0, 0.0, 1.0, 0.3)
        pad2.SetLeftMargin(0.14)
        pad2.SetRightMargin(0.06)
        pad2.SetTopMargin(0.02)
        pad2.SetBottomMargin(0.35)
        pad2.SetGridy()
        pad2.Draw()

        leg = ROOT.TLegend(0.55, 0.68, 0.88, 0.88)
        leg.SetBorderSize(0)
        leg.SetFillStyle(0)
        leg.SetTextFont(latex_font)
        leg.SetTextSize(0.035)

        #max_y = max(h.GetMaximum() for _, h in hists_list) if hists_list else 1.0

        pad1.cd()
        for idx, (file_label, h) in enumerate(hists_list):
            color = colors[idx % len(colors)]
            h.SetLineColor(color)
            h.SetLineWidth(2)
            h.SetMarkerColor(color)
            h.SetMarkerStyle(20)

            h.GetYaxis().SetTitleFont(latex_font)
            h.GetYaxis().SetLabelFont(latex_font)
            h.GetYaxis().SetTitleSize(0.045)
            h.GetYaxis().SetLabelSize(0.035)

            h.GetXaxis().SetLabelSize(0)
            h.GetXaxis().SetTitle("")

            h.SetMaximum(0.7) 
            h.SetMinimum(0)

            draw_option = "E1" if idx == 0 else "E1 SAME"
            h.Draw(draw_option)
            leg.AddEntry(h, file_label, "lep")

        leg.Draw()

        pad2.cd()
        ratio_hists = []
        
        if hists_list:
            # Safely project the reference TProfile to a standard TH1D for division
            ref_h_prof = hists_list[0][1]
            ref_h = ref_h_prof.ProjectionX(f"{ref_h_prof.GetName()}_px")

            for idx, (file_label, h_prof) in enumerate(hists_list):
                h_proj = h_prof.ProjectionX(f"{h_prof.GetName()}_px")
                
                h_ratio = h_proj.Clone(f"{h_prof.GetName()}_ratio")
                h_ratio.Divide(ref_h)

                color = colors[idx % len(colors)]
                h_ratio.SetLineColor(color)
                h_ratio.SetLineWidth(2)
                h_ratio.SetMarkerColor(color)
                h_ratio.SetMarkerStyle(20)

                h_ratio.GetXaxis().SetTitle(x_axis_title)
                h_ratio.GetYaxis().SetTitle("Ratio")

                h_ratio.GetXaxis().SetTitleFont(latex_font)
                h_ratio.GetYaxis().SetTitleFont(latex_font)
                h_ratio.GetXaxis().SetLabelFont(latex_font)
                h_ratio.GetYaxis().SetLabelFont(latex_font)

                h_ratio.GetXaxis().SetTitleSize(0.11)
                h_ratio.GetXaxis().SetLabelSize(0.09)
                h_ratio.GetXaxis().SetTitleOffset(1.2)
                h_ratio.GetYaxis().SetTitleSize(0.10)
                h_ratio.GetYaxis().SetLabelSize(0.08)
                h_ratio.GetYaxis().SetTitleOffset(0.55)
                h_ratio.GetYaxis().SetNdivisions(505)

                h_ratio.SetMinimum(0.5)
                h_ratio.SetMaximum(1.5)

                draw_option = "E1" if idx == 0 else "E1 SAME"
                h_ratio.Draw(draw_option)
                ratio_hists.append(h_ratio)

        c.Update()
        return c, leg, ratio_hists, (pad1, pad2)

    # 1. Build the ratio canvas (Projected for ratio calculation)
    c_mult, leg_mult, ratios_mult, pads_mult = build_ratio_canvas(
        f"c_avg_mult_E_{name}", f"Avg Multiplicity vs E_nu - {name}", hists_mult, "E_{#nu} [MeV]"
    )

    # 2. Build the standard unprojected TProfile canvas
    c_unproj = ROOT.TCanvas(f"c_unproj_mult_E_{name}", f"Unprojected Avg Multiplicity vs E_nu - {name}", 800, 600)
    c_unproj.SetLeftMargin(0.12)
    c_unproj.SetBottomMargin(0.12)
    
    leg_unproj = ROOT.TLegend(0.55, 0.70, 0.88, 0.88)
    leg_unproj.SetBorderSize(0)
    leg_unproj.SetFillStyle(0)
    leg_unproj.SetTextFont(latex_font)
    leg_unproj.SetTextSize(0.035)

    max_y_unproj = max(h.GetMaximum() for _, h in hists_mult) if hists_mult else 1.0

    c_unproj.cd()
    for idx, (file_label, h) in enumerate(hists_mult):
        color = colors[idx % len(colors)]
        h.SetLineColor(color)
        h.SetLineWidth(2)
        h.SetMarkerColor(color)
        h.SetMarkerStyle(20)

        # Formatting for a single-pad canvas
        h.GetYaxis().SetTitleFont(latex_font)
        h.GetYaxis().SetLabelFont(latex_font)
        h.GetYaxis().SetTitleSize(0.045)
        h.GetYaxis().SetLabelSize(0.035)
        h.GetXaxis().SetTitleFont(latex_font)
        h.GetXaxis().SetLabelFont(latex_font)
        h.GetXaxis().SetTitleSize(0.045)
        h.GetXaxis().SetLabelSize(0.035)
        
        # Ensure the axes have titles drawn properly
        h.GetXaxis().SetTitle("E_{#nu} [MeV]")
        
        h.SetMaximum(0.7)
        h.SetMinimum(0)

        draw_option = "E1" if idx == 0 else "E1 SAME"
        h.Draw(draw_option)
        leg_unproj.AddEntry(h, file_label, "lep")

    leg_unproj.Draw()
    c_unproj.Update()

    # Save both canvases to output file
    out_file = ROOT.TFile(f"paper_avg_multiplicity_per_E_{name}.root", "RECREATE")
    c_mult.Write()     # Split-pad canvas
    c_unproj.Write()   # Simple unprojected TProfile canvas
    out_file.Close()

    return {
        "multiplicity_ratio": (c_mult, hists_mult, leg_mult, ratios_mult),
        "multiplicity_unproj": (c_unproj, leg_unproj)
    }