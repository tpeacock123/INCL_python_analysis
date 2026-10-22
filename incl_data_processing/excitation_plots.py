from incl_analysis_functions import *
import ROOT
import numpy as np
import array
import sys 
import argparse
import os


font_id = 132
ROOT.gStyle.SetTextFont(font_id)
ROOT.gStyle.SetLegendFont(font_id)
ROOT.gStyle.SetLabelFont(font_id, "XYZ")  # Axis tick labels
ROOT.gStyle.SetTitleFont(font_id, "XYZ")  # Axis titles
ROOT.gStyle.SetTitleFont(font_id, "")     # Histogram/Pad titles
ROOT.gStyle.SetStatFont(font_id)

# Create a dark gray/charcoal color (#333333 or #444444)
dark_gray = ROOT.TColor.GetColor("#333333")

# Apply color globally to axes, labels, titles, and frames
ROOT.gStyle.SetAxisColor(dark_gray, "XYZ")   # Ticks and axis line color
ROOT.gStyle.SetLabelColor(dark_gray, "XYZ")  # Axis numbers color
ROOT.gStyle.SetTitleColor(dark_gray, "XYZ")  # Axis titles color
ROOT.gStyle.SetFrameLineColor(dark_gray)     # Plot border box color
ROOT.gStyle.SetTextColor(dark_gray)          # Default text color

ROOT.gStyle.SetTickLength(0.015, "XYZ")
ROOT.gStyle.SetNdivisions(306, "Y")
ROOT.gStyle.SetGridColor(ROOT.kGray)
ROOT.gStyle.SetOptTitle(0)

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

ROOT.gROOT.GetColor(ROOT.kYellow).SetRGB(
    122 / 255.0, 33 / 255.0, 221 / 255.0
)  # #7a21dd


ROOT.gSystem.Load("libNEUTROOTClass.so")
ROOT.gSystem.Load("libNEUTOutput.so")
ROOT.gSystem.Load("libNEUTReWeight.so")

ROOT.TH1.AddDirectory(False)



#ROOT.TGaxis.SetMaxDigits(4)

def excitation_plots(filename, name="excitation_combined_plot"):
    pre_FSI = False

    output_dic = {"ex_noCascadeFSI": [], "ex_noCascadeFSIPhoton": [], "ex_QE_deex": [], "ex_cascade_nodeex": [], "ex_cascade_nucleon": [], "ex_cascade_cluster": [], "ex_oneProton": [], "ex_protonPion": [], "ex_other": []}

    deex_pdg_map = {2112: "n", 2212: "p", 1000010020: "d", 1000020030: "3He", 1000020040: "#alpha"}

    colors = [ROOT.kYellow, ROOT.kBlue, ROOT.kRed, ROOT.kOrange + 7, ROOT.kMagenta]
    line_styles = [1, 1, 1, 1, 1]

    bin_edges = np.unique(np.concatenate([np.arange(-5.0, 0.1, 0.3), np.arange(0.0, 6.1, 0.3), np.arange(6.0, 60.1, 1.0)]))
    xbins = array.array("f", bin_edges)
    n_bins = len(xbins) - 1

    deex_nbins = 65
    deex_xmin = -5.5
    deex_xmax = 59.5

    y_axis_title = "Yield / Event"

    hists_flag8 = {}
    hists_flag10 = {}

    for idx, (pdg, label) in enumerate(deex_pdg_map.items()):
        color = colors[idx % len(colors)]
        style = line_styles[idx % len(line_styles)]

        h8 = ROOT.TH1D(f"h_deex_f8_{pdg}", f";;{y_axis_title}", deex_nbins, deex_xmin, deex_xmax)
        h8.SetLineColor(color)
        h8.SetLineStyle(style)
        h8.SetLineWidth(2)
        h8.SetStats(0)
        h8.GetYaxis().SetNoExponent(True)
        hists_flag8[pdg] = h8

        h10 = ROOT.TH1D(f"h_deex_f10_{pdg}", f";;{y_axis_title}", deex_nbins, deex_xmin, deex_xmax)
        h10.SetLineColor(color)
        h10.SetLineStyle(style)
        h10.SetLineWidth(2)
        h10.SetStats(0)
        h10.GetYaxis().SetNoExponent(True)
        hists_flag10[pdg] = h10

    f = ROOT.TFile.Open(filename)
    if not f or f.IsZombie():
        print(f"Error opening file {filename}")
        return

    t = f.Get("neuttree")
    total_mf_events = 0

    for i, event in enumerate(t):
        if i % 100000 == 0:
            print(f"Event number {i}")
        if i == 50000:
            break

        nvect = event.vectorbranch
        nvect_class = nvect_reader(nvect)

        eventType = nvect_class.eventType
        intChannel = nvect_class.intChannel

        if eventType == EventType.MF:
            total_mf_events += 1

            if pre_FSI:
                excitation_E = nvect_class.excitation_E_CCQE()
            else:
                excitation_E = nvect_class.post_casc_E_star()
                if excitation_E is None:
                    continue

            particles, energies, flags = nvect_class.particles(flag=True)
            for particle, flg in zip(particles, flags):
                if flg == 8 and particle in hists_flag8:
                    hists_flag8[particle].Fill(excitation_E)
                elif flg == 10 and particle in hists_flag10:
                    hists_flag10[particle].Fill(excitation_E)

            if (intChannel == intChannel_CCQE.muOnly) or (intChannel == intChannel_CCQE.neutronPion):
                continue

            if intChannel == intChannel_CCQE.noCascadeFSI:
                output_dic["ex_noCascadeFSI"].append(excitation_E)
            elif intChannel == intChannel_CCQE.noCascadeFSIPhoton:
                output_dic["ex_noCascadeFSIPhoton"].append(excitation_E)
            elif intChannel == intChannel_CCQE.qeDeEX:
                output_dic["ex_QE_deex"].append(excitation_E)
            elif intChannel == intChannel_CCQE.oneProton:
                output_dic["ex_oneProton"].append(excitation_E)
            elif intChannel == intChannel_CCQE.CascadeNoDeex:
                output_dic["ex_cascade_nodeex"].append(excitation_E)
            elif intChannel == intChannel_CCQE.CascadeNucleonDeex:
                output_dic["ex_cascade_nucleon"].append(excitation_E)
            elif intChannel == intChannel_CCQE.CascadeClusterDeex:
                output_dic["ex_cascade_cluster"].append(excitation_E)
            elif intChannel == intChannel_CCQE.protonPion:
                output_dic["ex_protonPion"].append(excitation_E)
            elif intChannel == intChannel_CCQE.other:
                output_dic["ex_other"].append(excitation_E)

    if total_mf_events > 0:
        for h in hists_flag8.values():
            h.Scale(1.0 / total_mf_events)
        for h in hists_flag10.values():
            h.Scale(1.0 / total_mf_events)

    max_y_pad1 = max([h.GetMaximum() for h in hists_flag8.values()] or [1.0])
    max_y_pad2 = max([h.GetMaximum() for h in hists_flag10.values()] or [1.0])

    c1 = ROOT.TCanvas("c1", "Excitation Energy & De-excitation Plots", 800, 1150)

    pad1 = ROOT.TPad("pad1", "pad1", 0.0, 0.7, 1.0, 1.0)
    pad1.SetTopMargin(0.08)
    pad1.SetBottomMargin(0.0)
    pad1.SetLeftMargin(0.12)
    pad1.SetRightMargin(0.05)
    pad1.SetTicks(1, 1)
    pad1.SetGrid(1, 1)
    pad1.Draw()

    pad2 = ROOT.TPad("pad2", "pad2", 0.0, 0.4, 1.0, 0.7)
    pad2.SetTopMargin(0.0)
    pad2.SetBottomMargin(0.0)
    pad2.SetLeftMargin(0.12)
    pad2.SetRightMargin(0.05)
    pad2.SetTicks(1, 1)
    pad2.SetGrid(1, 1)
    pad2.Draw()

    pad3 = ROOT.TPad("pad3", "pad3", 0.0, 0.0, 1.0, 0.4)
    pad3.SetTopMargin(0.0)
    pad3.SetBottomMargin(0.14)
    pad3.SetLeftMargin(0.12)
    pad3.SetRightMargin(0.05)
    pad3.SetTicks(1, 1)
    pad3.SetGrid(1, 1)
    pad3.Draw()

    # Histogram creation matching all updated output dictionary channels
    h_ex_noCascade = create_histo("h_ex_noC", "No Cascade", ROOT.kOrange - 4, 3001, output_dic["ex_noCascadeFSI"], n_bins, xbins)
    h_ex_noCascadePhoton = create_histo("h_ex_noC_p", "No Cascade FSI w/ photon", ROOT.kOrange - 1, 3001, output_dic["ex_noCascadeFSIPhoton"], n_bins, xbins)
    h_ex_QE = create_histo("h_ex_qe", "No Cascade w/ Deexcitation", ROOT.kOrange + 7, 3001, output_dic["ex_QE_deex"], n_bins, xbins)
    h_ex_cascNoDeex = create_histo("h_ex_casc_nodeex", "Cascade, No Deex", ROOT.kRed + 2, 3001, output_dic["ex_cascade_nodeex"], n_bins, xbins)
    h_ex_cascNucleon = create_histo("h_ex_casc_nuc", "Cascade, Nucleon Deex", ROOT.kViolet + 1, 3001, output_dic["ex_cascade_nucleon"], n_bins, xbins)
    h_ex_cascCluster = create_histo("h_ex_casc_clus", "Cascade with Clusters", ROOT.kViolet - 4, 3001, output_dic["ex_cascade_cluster"], n_bins, xbins)
    h_ex_oneProton = create_histo("h_ex_one", "Single Scattered Proton", ROOT.kGreen + 2, 3001, output_dic["ex_oneProton"], n_bins, xbins)
    h_ex_protonPion = create_histo("h_ex_pion", "Proton + Pion", ROOT.kYellow - 9, 3001, output_dic["ex_protonPion"], n_bins, xbins)
    h_ex_other = create_histo("h_ex_oth", "No Protons", ROOT.kGray, 3001, output_dic["ex_other"], n_bins, xbins)

    all_histograms = [h_ex_noCascade, h_ex_noCascadePhoton, h_ex_QE, h_ex_cascNoDeex, h_ex_cascNucleon, h_ex_cascCluster, h_ex_oneProton, h_ex_protonPion, h_ex_other]

    total_events = sum(h.Integral() for h in all_histograms)
    if total_events > 0:
        for h in all_histograms:
            h.Scale(1.0 / total_events, "width")

    hs1 = ROOT.THStack("hs1_ExcitationE", ";E_{x} [MeV];Normalised Counts")
    for h in [h_ex_other, h_ex_protonPion, h_ex_oneProton, h_ex_cascCluster, h_ex_cascNucleon, h_ex_cascNoDeex, h_ex_QE, h_ex_noCascadePhoton, h_ex_noCascade]:
        hs1.Add(h)

    # Pad 1 (Flag 8 De-excitation)
    pad1.cd()
    for idx, (pdg, h) in enumerate(hists_flag8.items()):
        h.GetYaxis().SetTitle(y_axis_title)
        h.GetYaxis().SetLabelSize(0.08)
        h.GetYaxis().SetTitleSize(0.10)
        h.GetYaxis().SetTitleOffset(0.6)
        h.GetYaxis().ChangeLabel(1, -1, -1, 31)
        h.GetYaxis().CenterTitle(True)
        h.GetYaxis().SetNoExponent(True)
        h.GetXaxis().SetLabelSize(0)

        if idx == 0:
            h.SetMaximum(max_y_pad1 * 1.25)
            h.SetMinimum(0)
            h.Draw("HIST")
        else:
            h.Draw("HIST SAME")

    label_f8 = ROOT.TLatex()
    label_f8.SetNDC()
    label_f8.SetTextSize(0.08)
    label_f8.SetTextColor(ROOT.kGray + 2)
    label_f8.DrawLatex(0.15, 0.80, "#it{Deexcitation}")

    # Pad 2 (Flag 10 Cascade)
    pad2.cd()
    leg_deex = ROOT.TLegend(0.78, 0.45, 0.88, 0.88)
    leg_deex.SetBorderSize(0)
    leg_deex.SetFillStyle(1001)
    leg_deex.SetTextSize(0.08)

    for idx, (pdg, h) in enumerate(hists_flag10.items()):
        h.GetYaxis().SetTitle(y_axis_title)
        h.GetYaxis().SetLabelSize(0.08)
        h.GetYaxis().SetTitleSize(0.10)
        h.GetYaxis().SetTitleOffset(0.6)
        h.GetYaxis().CenterTitle(True)
        h.GetYaxis().ChangeLabel(1, -1, -1, 31)
        h.GetYaxis().SetNoExponent(True)
        h.GetXaxis().SetLabelSize(0)

        if idx == 0:
            h.SetMaximum(max_y_pad2 * 1.25)
            h.SetMinimum(0)
            h.Draw("HIST")
        else:
            h.Draw("HIST SAME")

        leg_deex.AddEntry(h, deex_pdg_map[pdg], "l")

    leg_deex.Draw()

    label_f10 = ROOT.TLatex()
    label_f10.SetNDC()
    label_f10.SetTextSize(0.08)
    label_f10.SetTextColor(ROOT.kGray + 2)
    label_f10.DrawLatex(0.15, 0.80, "#it{Cascade}")

    # Pad 3 (Stacked Channel Breakdown)
    pad3.cd()
    pad3.SetLogy(1)
    hs1.Draw("hist")

    hs1.GetYaxis().SetLabelSize(0.06)
    hs1.GetYaxis().SetTitleSize(0.08)
    hs1.GetYaxis().SetTitleOffset(0.75)
    hs1.GetYaxis().CenterTitle(True)

    hs1.GetXaxis().SetLabelSize(0.07)
    hs1.GetXaxis().SetTitleSize(0.08)
    hs1.GetXaxis().SetTitleOffset(0.75)
    hs1.GetXaxis().CenterTitle(True)

    max_h = hs1.GetMaximum()
    hs1.SetMaximum(max_h * 10.0 if max_h > 0 else 1.1)

    legend1 = ROOT.TLegend(0.48, 0.55, 0.92, 0.93)
    legend1.SetBorderSize(0)
    legend1.SetFillColor(ROOT.kWhite)
    legend1.SetFillStyle(1001)
    legend1.SetTextSize(0.035)
    legend1.AddEntry(h_ex_noCascade, "No Cascade", "f")
    legend1.AddEntry(h_ex_noCascadePhoton, "No Cascade FSI w/ #gamma", "f")
    legend1.AddEntry(h_ex_QE, "No Cascade w/ Deexcitation", "f")
    legend1.AddEntry(h_ex_cascNoDeex, "Cascade, No Deex", "f")
    legend1.AddEntry(h_ex_cascNucleon, "Cascade, Nucleon Deex", "f")
    legend1.AddEntry(h_ex_cascCluster, "Cascade with Clusters", "f")
    legend1.AddEntry(h_ex_oneProton, "Single Scattered Proton", "f")
    legend1.AddEntry(h_ex_protonPion, "Proton + Pion", "f")
    legend1.AddEntry(h_ex_other, "No Protons", "f")
    legend1.Draw()

    for pad in [pad1, pad2, pad3]:
        pad.Modified()
        pad.Update()

    c1.Modified()
    c1.Update()
    output_dir = "triple_deex_plots"
    os.makedirs(output_dir, exist_ok=True)

    filepath = os.path.join(output_dir, f"triple_deex_{name}.pdf")
    c1.SaveAs(filepath)

    return c1



def _channel_migration(
    filename, output_filename="fsi_excitation_2d_plots"
):
    # Ensure string format without double extensions
    base_name = str(output_filename).replace(".root", "")
    root_filename = f"migration_{base_name}.root"

    channel_map = {
        intChannel_CCQE.noCascadeFSI: "noCascadeFSI",
        intChannel_CCQE.noCascadeFSIPhoton: "noCascadeFSIPhoton",
        intChannel_CCQE.qeDeEX: "QE_deex",
        intChannel_CCQE.oneProton: "oneProton",
        intChannel_CCQE.multipleNucleon: "multipleNucleon_noCluster",
        intChannel_CCQE.nuclearCluster: "nuclearCluster",
        intChannel_CCQE.protonPion: "protonPion",
        intChannel_CCQE.other: "other",
    }

    nbins_2d = 55
    xmin_2d = -5.0
    xmax_2d = 105.0

    # Initialize 2D histograms per channel
    hists_2d = {}
    for ch_enum, ch_name in channel_map.items():
        h2 = ROOT.TH2D(
            f"h2_pre_vs_post_{ch_name}",
            f"{ch_name};Post-Cascade E_{{ex}} [MeV];Pre-Cascade E_{{ex}} [MeV]",
            nbins_2d,
            xmin_2d,
            xmax_2d,
            nbins_2d,
            xmin_2d,
            xmax_2d,
        )
        h2.SetStats(0)
        hists_2d[ch_enum] = h2

    # Inclusive 2D histogram
    h2_inclusive = ROOT.TH2D(
        "h2_pre_vs_post_inclusive",
        f"Inclusive;Post-Cascade E_{{ex}} [MeV];Pre-Cascade E_{{ex}} [MeV]",
        nbins_2d,
        xmin_2d,
        xmax_2d,
        nbins_2d,
        xmin_2d,
        xmax_2d,
    )
    h2_inclusive.SetStats(0)

    f = ROOT.TFile.Open(str(filename))
    if not f or f.IsZombie():
        print(f"Error opening file {filename}")
        return None, None

    t = f.Get("neuttree")

    for i, event in enumerate(t):
        if i % 100000 == 0:
            print(f"Event number {i}")
        if i == 500000:
            break

        nvect = event.vectorbranch
        nvect_class = nvect_reader(nvect)

        eventType = nvect_class.eventType
        intChannel = nvect_class.intChannel

        if eventType == EventType.MF:
            if (intChannel == intChannel_CCQE.muOnly) or (
                intChannel == intChannel_CCQE.neutronPion
            ):
                continue

            pre_FSI_E = nvect_class.excitation_E_CCQE()
            post_FSI_E = nvect_class.post_casc_E_star()



            if pre_FSI_E is None or post_FSI_E is None:
                continue

            if post_FSI_E - pre_FSI_E < -5.0:
                nvect_class.Print()

            if intChannel in hists_2d:
                hists_2d[intChannel].Fill(post_FSI_E, pre_FSI_E)
                h2_inclusive.Fill(post_FSI_E, pre_FSI_E)

    # --- Canvas & Formatting Setup ---
    # Enable titles locally for the canvas grid
    ROOT.gStyle.SetOptTitle(1)

    all_plots = [h2_inclusive] + list(hists_2d.values())
    
    # Create large canvas (1800x1600 pixels) for a 3x3 grid
    c = ROOT.TCanvas("c_channel_migration", "Channel Migration 2D", 1800, 1600)
    c.Divide(3, 3)

    for idx, h2 in enumerate(all_plots, 1):
        pad = c.cd(idx)
        
        # Increase margins to avoid cutting off axis titles and colorbar
        pad.SetLeftMargin(0.14)
        pad.SetRightMargin(0.16)   # Space for colorbar palette
        pad.SetBottomMargin(0.14)
        pad.SetTopMargin(0.10)     # Space for pad title
        pad.SetLogz()

        # Adjust font sizes for divided pad scaling
        h2.GetXaxis().SetTitleSize(0.048)
        h2.GetXaxis().SetLabelSize(0.042)
        h2.GetXaxis().SetTitleOffset(1.1)

        h2.GetYaxis().SetTitleSize(0.048)
        h2.GetYaxis().SetLabelSize(0.042)
        h2.GetYaxis().SetTitleOffset(1.2)

        h2.GetZaxis().SetLabelSize(0.038)

        # Draw 2D histogram with color palette
        h2.Draw("COLZ")

    c.Update()

    # Save to image files and ROOT file

    c.SaveAs(f"migration_{base_name}.pdf")

    out_file = ROOT.TFile(root_filename, "RECREATE")
    h2_inclusive.Write()
    for h2 in hists_2d.values():
        h2.Write()
    c.Write("canvas_all_channels")
    out_file.Close()

    # Re-apply global opt title setting if needed elsewhere
    ROOT.gStyle.SetOptTitle(0)

    return hists_2d, h2_inclusive



def channel_migration(
    filename, output_filename="fsi_excitation_2d_plots"
):
    # Ensure string format without double extensions
    base_name = str(output_filename).replace(".root", "")
    root_filename = f"migration_{base_name}.root"

    channel_map = {
        intChannel_CCQE.noCascadeFSI: "No cascade, no de-excitation",
#        intChannel_CCQE.noCascadeFSIPhoton: "noCascadeFSIPhoton",
        intChannel_CCQE.qeDeEX: "No cascade, deexcitation",
        intChannel_CCQE.oneProton: "one proton",
        intChannel_CCQE.CascadeNoDeex: "Cascade, No Deexcitation",
        intChannel_CCQE.CascadeNucleonDeex: "Cascade, Deexcitation, no Clusters",
        intChannel_CCQE.CascadeClusterDeex: "Cascade, Deexcitation, Clusters",
#        intChannel_CCQE.protonPion: "protonPion",
#        intChannel_CCQE.other: "other",
    }

    nbins_2d = 55
    ymin_2d = -5.0
    ymax_2d = 105.0
    xmax = 2000
    xmin = 200
    

    # Initialize 2D histograms per channel
    hists_2d = {}
    for ch_enum, ch_name in channel_map.items():
        h2 = ROOT.TH2D(
            f"h2_pre_vs_post_{ch_name}",
            f"{ch_name}; Neutrino Energy [MeV];Post-Cascade E_{{ex}} [MeV]",
            nbins_2d,
            xmin,
            xmax,
            nbins_2d,
            ymin_2d,
            ymax_2d,
        )
        h2.SetStats(0)
        hists_2d[ch_enum] = h2

    # Inclusive 2D histogram
    h2_inclusive = ROOT.TH2D(
        "h2_pre_vs_post_inclusive",
        f"Inclusive; Neutrino Energy [MeV];Post-Cascade E_{{ex}} [MeV]",
        nbins_2d,
        xmin,
        xmax,
        nbins_2d,
        ymin_2d,
        ymax_2d,
    )
    h2_inclusive.SetStats(0)

    f = ROOT.TFile.Open(str(filename))
    if not f or f.IsZombie():
        print(f"Error opening file {filename}")
        return None, None

    t = f.Get("neuttree")

    for i, event in enumerate(t):
        if i % 100000 == 0:
            print(f"Event number {i}")
        if i == 50000:
            break

        nvect = event.vectorbranch
        nvect_class = nvect_reader(nvect)

        eventType = nvect_class.eventType
        intChannel = nvect_class.intChannel

        if eventType == EventType.MF:
            if (intChannel == intChannel_CCQE.muOnly) or (
                intChannel == intChannel_CCQE.neutronPion
            ):
                continue

            pre_FSI_E = nvect_class.neutrino_E
            post_FSI_E = nvect_class.post_casc_E_star()



            if pre_FSI_E is None or post_FSI_E is None:
                continue


            if intChannel in hists_2d:
                if (intChannel == intChannel_CCQE.CascadeNoDeex):
                    if post_FSI_E > 10.0:
                        
                        nvect_class.Print()

                hists_2d[intChannel].Fill(pre_FSI_E,post_FSI_E)
                h2_inclusive.Fill(pre_FSI_E,post_FSI_E)

    # --- Canvas & Formatting Setup ---
    # Enable titles locally for the canvas grid
    ROOT.gStyle.SetOptTitle(1)

    all_plots = [h2_inclusive] + list(hists_2d.values())
    
    # Create large canvas (1800x1600 pixels) for a 3x3 grid
    c = ROOT.TCanvas("c_channel_migration", "Channel Migration 2D", 1800, 1600)
    c.Divide(3, 3)

    for idx, h2 in enumerate(all_plots, 1):
        pad = c.cd(idx)
        
        # Increase margins to avoid cutting off axis titles and colorbar
        pad.SetLeftMargin(0.14)
        pad.SetRightMargin(0.16)   # Space for colorbar palette
        pad.SetBottomMargin(0.14)
        pad.SetTopMargin(0.10)     # Space for pad title
        pad.SetLogz()

        # Adjust font sizes for divided pad scaling
        h2.GetXaxis().SetTitleSize(0.048)
        h2.GetXaxis().SetLabelSize(0.042)
        h2.GetXaxis().SetTitleOffset(1.1)

        h2.GetYaxis().SetTitleSize(0.048)
        h2.GetYaxis().SetLabelSize(0.042)
        h2.GetYaxis().SetTitleOffset(1.2)

        h2.GetZaxis().SetLabelSize(0.038)

        # Draw 2D histogram with color palette
        h2.Draw("COLZ")

    c.Update()

    # Save to image files and ROOT file

    c.SaveAs(f"migration_{base_name}.pdf")

    out_file = ROOT.TFile(root_filename, "RECREATE")
    h2_inclusive.Write()
    for h2 in hists_2d.values():
        h2.Write()
    c.Write("canvas_all_channels")
    out_file.Close()

    # Re-apply global opt title setting if needed elsewhere
    ROOT.gStyle.SetOptTitle(0)

    return hists_2d, h2_inclusive