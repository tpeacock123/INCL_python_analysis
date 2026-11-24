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
#ROOT.TGaxis.SetMaxDigits(4)
def excitation_plots(filename, name="excitation_combined_plot"):
    pre_FSI    = False

    output_dic = {
        "ex_noCascadeFSI": [],
        "ex_noCascadeFSIPhoton": [],
        "ex_QE_deex": [],
        "ex_cascade_nucleon_nodeex": [],
        "ex_cascade_cluster_nodeex": [],
        "ex_cascade_nucleon": [],
        "ex_cascade_cluster": [],
        "ex_oneProton": [],
        "ex_protonAbsorption": [],
        "ex_protonPion": [],
        "ex_other": [],
    }

    deex_pdg_map = {
        2112: "n",
        2212: "p",
        1000010020: "d",
        1000020030: "3He",
        1000020040: "#alpha",
    }

    colors = [ROOT.kYellow, ROOT.kBlue, ROOT.kRed, ROOT.kOrange + 7, ROOT.kMagenta]
    line_styles = [1, 1, 1, 1, 1]

    deex_nbins = 82
    deex_xmin = -1
    deex_xmax = 81

    bin_edges = np.linspace(deex_xmin, deex_xmax, deex_nbins + 1)
    xbins = array.array("f", bin_edges)
    n_bins = deex_nbins

    y_axis_title = "[particle / Evt / 0.5 MeV]"

    hists_flag8 = {}
    hists_flag10 = {}

    for idx, (pdg, label) in enumerate(deex_pdg_map.items()):
        color = colors[idx % len(colors)]
        style = line_styles[idx % len(line_styles)]

        h8 = ROOT.TH1D(
            f"h_deex_f8_{pdg}", f";;{y_axis_title}", deex_nbins, deex_xmin, deex_xmax
        )
        h8.SetLineColor(color)
        h8.SetLineStyle(style)
        h8.SetLineWidth(2)
        h8.SetStats(0)
        h8.GetYaxis().SetNoExponent(True)
        h8.GetXaxis().SetNdivisions(505, True)
        hists_flag8[pdg] = h8

        h10 = ROOT.TH1D(
            f"h_deex_f10_{pdg}", f";;{y_axis_title}", deex_nbins, deex_xmin, deex_xmax
        )
        h10.SetLineColor(color)
        h10.SetLineStyle(style)
        h10.SetLineWidth(2)
        h10.SetStats(0)
        h10.GetYaxis().SetNoExponent(True)
        h10.GetXaxis().SetNdivisions(505, True)
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
        if i == 500000:
            break

        nvect = event.vectorbranch
        nvect_class = nvect_reader(nvect)

        eventType = nvect_class.eventType
        intChannel = nvect_class.intChannel

        if eventType == EventType.MF:
            total_mf_events += 1

            if pre_FSI:
                excitation_E = nvect_class.excitation_E_CCQE(carbon=False)
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

            if (intChannel == intChannel_CCQE.muOnly) or (
                intChannel == intChannel_CCQE.neutronPion
            ):
                continue

            if intChannel == intChannel_CCQE.noCascadeFSI:
                output_dic["ex_noCascadeFSI"].append(excitation_E)
            elif intChannel == intChannel_CCQE.noCascadeFSIPhoton:
                output_dic["ex_noCascadeFSIPhoton"].append(excitation_E)
            elif intChannel == intChannel_CCQE.qeDeEX:
                output_dic["ex_QE_deex"].append(excitation_E)
            elif intChannel == intChannel_CCQE.oneProton:
                output_dic["ex_oneProton"].append(excitation_E)
            elif intChannel == intChannel_CCQE.protonAbsorption:
                output_dic["ex_protonAbsorption"].append(excitation_E)
                if excitation_E < 10:
                    nvect_class.Print()
            elif intChannel == intChannel_CCQE.CascadeNucleonNoDeex:
                output_dic["ex_cascade_nucleon_nodeex"].append(excitation_E)
            elif intChannel == intChannel_CCQE.CascadeClusterNoDeex:
                output_dic["ex_cascade_cluster_nodeex"].append(excitation_E)
            elif intChannel == intChannel_CCQE.CascadeNucleonDeex:
                output_dic["ex_cascade_nucleon"].append(excitation_E)
            elif intChannel == intChannel_CCQE.CascadeClusterDeex:
                output_dic["ex_cascade_cluster"].append(excitation_E)
            elif intChannel == intChannel_CCQE.protonPion or intChannel == intChannel_CCQE.neutronPion:
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
    h_ex_cascNucleonNoDeex = create_histo("h_ex_casc_nuc_nodeex", "Cascade Nucleon, No Deex", ROOT.kRed + 2, 3001, output_dic["ex_cascade_nucleon_nodeex"], n_bins, xbins)
    h_ex_cascClusterNoDeex = create_histo("h_ex_casc_clus_nodeex", "Cascade Cluster, No Deex", ROOT.kRed - 4, 3001, output_dic["ex_cascade_cluster_nodeex"], n_bins, xbins)
    h_ex_cascNucleon = create_histo("h_ex_casc_nuc", "Cascade, Nucleon Deex", ROOT.kViolet + 4, 3001, output_dic["ex_cascade_nucleon"], n_bins, xbins)
    h_ex_cascCluster = create_histo("h_ex_casc_clus", "Cascade with Clusters", ROOT.kViolet - 4, 3001, output_dic["ex_cascade_cluster"], n_bins, xbins)
    h_ex_oneProton = create_histo("h_ex_one", "Single Scattered Proton", ROOT.kGreen + 2, 3001, output_dic["ex_oneProton"], n_bins, xbins)
    h_ex_protonAbsorption = create_histo("h_ex_p_abs", "Proton Absorption", ROOT.kGreen + 4, 3001, output_dic["ex_protonAbsorption"], n_bins, xbins)
    h_ex_protonPion = create_histo("h_ex_pion", "Proton + Pion", ROOT.kYellow - 9, 3001, output_dic["ex_protonPion"], n_bins, xbins)
    h_ex_other = create_histo("h_ex_oth", "No Protons", ROOT.kGray, 3001, output_dic["ex_other"], n_bins, xbins)

    all_histograms = [
        h_ex_noCascade,
        h_ex_noCascadePhoton,
        h_ex_QE,
        h_ex_cascNucleonNoDeex,
        h_ex_cascClusterNoDeex,
        h_ex_cascNucleon,
        h_ex_cascCluster,
        h_ex_oneProton,
        h_ex_protonAbsorption,
        h_ex_protonPion,
        h_ex_other,
    ]

    total_events = sum(h.Integral() for h in all_histograms)
    if total_events > 0:
        for h in all_histograms:
            h.Scale(1.0 / total_events, "width")

    hs1 = ROOT.THStack("hs1_ExcitationE", ";E*_{FSI} [MeV];Normalised Counts")
    for h in [
        h_ex_other,
        h_ex_protonPion,
        h_ex_cascCluster,
        h_ex_cascNucleon,
        h_ex_cascClusterNoDeex,
        h_ex_cascNucleonNoDeex,
        h_ex_oneProton,
        h_ex_protonAbsorption,
        h_ex_QE,
        h_ex_noCascadePhoton,
        h_ex_noCascade,
    ]:
        hs1.Add(h)

    # Pad 1 (Flag 8 De-excitation)
    pad1.cd()
    for idx, (pdg, h) in enumerate(hists_flag8.items()):
        h.GetYaxis().SetTitle(y_axis_title)
        h.GetYaxis().SetLabelSize(0.07)
        h.GetYaxis().SetTitleSize(0.08)
        h.GetYaxis().SetTitleOffset(0.7)
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
    label_f8.DrawLatex(0.15, 0.80, "#it{Cascade}")

    # Pad 2 (Flag 10 Cascade)
    pad2.cd()
    leg_deex = ROOT.TLegend(0.15, 0.45, 0.25, 0.88)
    leg_deex.SetBorderSize(0)
    leg_deex.SetFillStyle(1001)
    leg_deex.SetTextSize(0.08)

    for idx, (pdg, h) in enumerate(hists_flag10.items()):
        h.GetYaxis().SetTitle(y_axis_title)
        h.GetYaxis().SetLabelSize(0.07)
        h.GetYaxis().SetTitleSize(0.08)
        h.GetYaxis().SetTitleOffset(0.7)
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
    label_f10.DrawLatex(0.15, 0.90, "#it{De-excitation}")

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
    hs1.GetXaxis().SetNdivisions(505, True)

    max_h = hs1.GetMaximum()
    hs1.SetMaximum(1.5)
    hs1.SetMinimum(5e-4)

    legend1 = ROOT.TLegend(0.58, 0.55, 0.92, 0.93)
    legend1.SetBorderSize(0)
    legend1.SetFillColor(ROOT.kWhite)
    legend1.SetFillStyle(1001)
    legend1.SetTextSize(0.035)
    legend1.AddEntry(h_ex_noCascade, "No Cascade", "f")
#   legend1.AddEntry(h_ex_noCascadePhoton, "No Cascade FSI w/ #gamma", "f")
    legend1.AddEntry(h_ex_QE, "No Cascade w/ Deexcitation", "f")
    legend1.AddEntry(h_ex_cascNucleonNoDeex, "Cascade Nucleon, No Deex", "f")
#   legend1.AddEntry(h_ex_cascClusterNoDeex, "Cascade Cluster, No Deex", "f")
    legend1.AddEntry(h_ex_cascNucleon, "Cascade, Deex no Clusters", "f")
    legend1.AddEntry(h_ex_cascCluster, "Cascade, Deex Clusters", "f")
    legend1.AddEntry(h_ex_oneProton, "Single Scattered Proton", "f")
    legend1.AddEntry(h_ex_protonAbsorption, "Proton Absorption", "f")
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

    filepath = os.path.join(output_dir, f"triple_deex_{name}.root")
    c1.SaveAs(filepath)
    filepath = os.path.join(output_dir, f"triple_deex_{name[0]}.pdf")
    c1.SaveAs(filepath)

    return c1


def channel_migration(
    filename, output_filename="fsi_excitation_2d_plots"
):
    # Ensure string format without double extensions
    base_name = str(output_filename).replace(".root", "")
    root_filename = f"migration_{base_name}.root"

    channel_map = {
        intChannel_CCQE.noCascadeFSI: "No cascade, no de-excitation",
        intChannel_CCQE.noCascadeFSIPhoton: "noCascadeFSIPhoton",
        intChannel_CCQE.qeDeEX: "No cascade, deexcitation",
        intChannel_CCQE.oneProton: "one proton",
        intChannel_CCQE.CascadeNucleonNoDeex: "Cascade Nucleon, No Deexcitation",
        intChannel_CCQE.CascadeClusterNoDeex: "Cascade Cluster, No Deexcitation",
        intChannel_CCQE.CascadeNucleonDeex: "Cascade, Deexcitation, no Clusters",
        intChannel_CCQE.CascadeClusterDeex: "Cascade, Deexcitation, Clusters",
#        intChannel_CCQE.protonPion: "protonPion",
#        intChannel_CCQE.other: "other",
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

            pre_FSI_E = nvect_class.excitation_E_CCQE(carbon=False)
            post_FSI_E = nvect_class.post_casc_E_star()



            if pre_FSI_E is None or post_FSI_E is None:
                continue

            if intChannel == intChannel_CCQE.qeDeEX:
                if pre_FSI_E  < 1. :
                    print("pre fsi E", pre_FSI_E)
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


def _channel_migration(
    filename, output_filename="fsi_excitation_2d_plots"
):
    # Ensure string format without double extensions
    base_name = str(output_filename).replace(".root", "")
    root_filename = f"migration_{base_name}.root"

    # Split duplicate labels to distinguish Nucleon vs Cluster no-deex channels
    channel_map = {
        intChannel_CCQE.noCascadeFSI: "No cascade, no de-excitation",
#        intChannel_CCQE.noCascadeFSIPhoton: "noCascadeFSIPhoton",
        intChannel_CCQE.qeDeEX: "No cascade, deexcitation",
        intChannel_CCQE.oneProton: "one proton",
        intChannel_CCQE.CascadeNucleonNoDeex: "Cascade Nucleon, No Deexcitation",
        intChannel_CCQE.CascadeClusterNoDeex: "Cascade Cluster, No Deexcitation",
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
        "Inclusive; Neutrino Energy [MeV];Post-Cascade E_{{ex}} [MeV]",
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
                # Updated to inspect both split No-Deexcitation channels
                if intChannel == intChannel_CCQE.qeDeEX:
                    if pre_FSI_E  < 1. :
                        print("pre fsi E", pre_FSI_E)
                        nvect_class.Print()

                hists_2d[intChannel].Fill(pre_FSI_E, post_FSI_E)
                h2_inclusive.Fill(pre_FSI_E, post_FSI_E)

    # --- Canvas & Formatting Setup ---
    # Enable titles locally for the canvas grid
    ROOT.gStyle.SetOptTitle(1)

    all_plots = [h2_inclusive] + list(hists_2d.values())

    # Create large canvas (1800x1600 pixels) for a 3x3 grid (8 plots fit within 9 pads)
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

def src_mf_excitation_E(filename, name, no_events=500000, filename_list=None):
    if filename_list is None:
        filename_list = [filename]
    else:
        filename_list = list(filename_list)
        if filename not in filename_list:
            filename_list.insert(0, filename)
    
    colors = [ROOT.kRed, ROOT.kBlue]
    line_styles = [1, 1]

    deex_nbins = 100
    deex_xmin = -50.0
    deex_xmax = 50.0

    bin_edges = np.linspace(deex_xmin, deex_xmax, deex_nbins + 1)
    xbins = array.array("f", bin_edges)

    num_files = len(filename_list)
    pad_labels = ["NEUT SRC", "NuWo SRC"]

    file_data = []

    for idx, file in enumerate(filename_list):
        f = ROOT.TFile.Open(file)
        if not f or f.IsZombie():
            print(f"Error opening file {file}")
            continue

        t = f.Get("neuttree")
        if not t:
            print(f"Tree 'neuttree' not found in file {file}")
            f.Close()
            continue

        hist_mf = ROOT.TH1F(f"hist_mf_{idx}", "Mean Field", deex_nbins, xbins)
        hist_src = ROOT.TH1F(f"hist_src_{idx}", "SRC", deex_nbins, xbins)

        for i, event in enumerate(t):
            if i % 100000 == 0:
                print(f"File {idx+1} ({file}) - Event number {i}")
            if i >= no_events:
                break

            nvect = event.vectorbranch
            nvect_class = nvect_reader(nvect)
            eventType = nvect_class.eventType
            
            if eventType == EventType.MF:
                excitation_E = nvect_class.excitation_E_CCQE()
                hist_mf.Fill(excitation_E)
            elif eventType == EventType.SRC:
                excitation_E = nvect_class.excitation_E_SRC()
                hist_src.Fill(excitation_E)

        f.Close()

        # Area normalisation
        if hist_mf.Integral() > 0:
            hist_mf.Scale(1.0 / hist_mf.Integral())
        if hist_src.Integral() > 0:
            hist_src.Scale(1.0 / hist_src.Integral())

        # Line styling
        hist_mf.SetLineColor(colors[0])
        hist_mf.SetLineStyle(line_styles[0])
        hist_mf.SetLineWidth(2)

        hist_src.SetLineColor(colors[1])
        hist_src.SetLineStyle(line_styles[1])
        hist_src.SetLineWidth(2)

        hs = ROOT.THStack(f"hs_excitation_E_{idx}", ";Excitation Energy [MeV]; Probability Distribution")
        hs.Add(hist_mf)
        hs.Add(hist_src)

        file_data.append({
            "idx": idx,
            "hs": hs,
            "hist_mf": hist_mf,
            "hist_src": hist_src
        })

    if not file_data:
        return None

    # -------------------------------------------------------------
    # 2. Determine global max Y for synchronized axes
    # -------------------------------------------------------------
    global_max_y = max(data["hs"].GetMaximum() for data in file_data)

    # -------------------------------------------------------------
    # 3. Canvas setup with seamless adjacent pads
    # -------------------------------------------------------------
    c1 = ROOT.TCanvas("c1", "Excitation Energy & De-excitation Plots", 1000, 500)
    c1.Divide(num_files, 1, 0.0, 0.0)

    all_stacks = []
    all_hists = []
    all_legends = []
    all_labels = []

    for data in file_data:
        idx = data["idx"]
        hs = data["hs"]
        hist_mf = data["hist_mf"]
        hist_src = data["hist_src"]

        pad = c1.cd(idx + 1)
        ROOT.gPad.SetLogy(1)

        pad.SetBottomMargin(0.15)
        pad.SetTopMargin(0.00)

        # Set unified Y limits on log scale across both pads
        hs.SetMinimum(1e-2)
        hs.SetMaximum(global_max_y * 3.0)

        # Adjust side margins to seamlessly join Pad 1 right edge & Pad 2 left edge
        if idx == 0:
            pad.SetLeftMargin(0.14)
            pad.SetRightMargin(0.0)
        else:
            pad.SetLeftMargin(0.0)
            pad.SetRightMargin(0.05)

        hs.Draw("HIST")

        # Formatting X axis
        hs.GetXaxis().CenterTitle(True)
        hs.GetXaxis().SetNoExponent(True)
        hs.GetXaxis().SetTitleSize(0.055)
        hs.GetXaxis().SetLabelSize(0.045)
        hs.GetXaxis().SetTitleOffset(1.1)
        hs.GetXaxis().SetNdivisions(606, True)

        # Formatting Y axis
        if idx == 0:
            hs.GetYaxis().CenterTitle(True)
            hs.GetYaxis().SetMaxDigits(4)
            hs.GetYaxis().SetTitleSize(0.055)
            hs.GetYaxis().SetLabelSize(0.045)
            hs.GetYaxis().SetTitleOffset(1.2)
        else:
            hs.GetYaxis().SetTitle("")
            hs.GetYaxis().SetLabelSize(0)

        # Draw pad text label (shifted up along with the plot)
        label_text = pad_labels[idx] if idx < len(pad_labels) else f"Pad {idx+1}"
        latex = ROOT.TLatex()
        latex.SetNDC(True)
        latex.SetTextFont(52)
        latex.SetTextColor(ROOT.kGray + 2)
        latex.SetTextSize(0.045)
        x_pos = 0.18 if idx == 0 else 0.05
        latex.DrawLatex(x_pos, 0.88, label_text)

        # Draw legend (shifted up along with the plot)
        leg_x1 = 0.65 if idx == 0 else 0.60
        leg_x2 = 0.93 if idx == 0 else 0.88
        leg = ROOT.TLegend(leg_x1, 0.74, leg_x2, 0.92)
        leg.SetBorderSize(0)
        leg.SetFillStyle(0)
        leg.AddEntry(hist_mf, "MF", "l")
        leg.AddEntry(hist_src, "SRC", "l")
        leg.Draw()

        all_stacks.append(hs)
        all_hists.extend([hist_mf, hist_src])
        all_legends.append(leg)
        all_labels.append(latex)

    c1.Modified()
    c1.Update()

    # Output generation
    output_dir = "triple_deex_plots"
    os.makedirs(output_dir, exist_ok=True)
    
    name_str = name[0] if isinstance(name, (list, tuple)) else name
    filepath_root = os.path.join(output_dir, f"SRC_vs_MF_deex_{name_str}.root")
    c1.SaveAs(filepath_root)
    filepath_pdf = os.path.join(output_dir, f"SRC_vs_MF_deex_{name_str}.pdf")
    c1.SaveAs(filepath_pdf)

    c1._stacks = all_stacks
    c1._hists = all_hists
    c1._legends = all_legends
    c1._labels = all_labels

    return c1

def excitation_energy_differences(
    filename, output_filename="comps", n_events=500000, files= []
):
    # Define legend labels, line styles, and colors per file index
    styles = {
        0: {"label": "Post Fix", "linestyle": ROOT.kSolid, "color": ROOT.kBlue + 1},
        1: {"label": "pre Fix", "linestyle": ROOT.kDashed, "color": ROOT.kRed + 1},
        2: {"label": "Post Fix 2", "linestyle": ROOT.kDotted, "color": ROOT.kGreen + 2},
    }

    ROOT.gStyle.SetOptStat(0)

    # Combine primary filename with any additional files
    all_files = (
        list(filename) if isinstance(filename, (list, tuple)) else [filename]
    ) + list(files)

    def read_file_diffs(filepath):
        f = ROOT.TFile.Open(str(filepath))
        if not f or f.IsZombie():
            print(f"Error opening file {filepath}")
            return None

        t = f.Get("neuttree")

        deex_mf, nondeex_mf = [], []
        deex_src, nondeex_src = [], []

        for i, event in enumerate(t):
            if i % 100000 == 0:
                print(f"[{os.path.basename(str(filepath))}] Event number {i}")
            if i == n_events:
                break

            nvect = event.vectorbranch
            nvect_class = nvect_reader(nvect)
            eventType = nvect_class.eventType

            if eventType == EventType.MF:
                particles, energies, flags = nvect_class.particles(flag=True)
                pre_FSI_E = nvect_class.excitation_E_CCQE(carbon=True)
                post_FSI_E = nvect_class.post_casc_E_star()
                if 10 in flags:
                    deex_mf.append(post_FSI_E - pre_FSI_E)
                else:
                    nondeex_mf.append(post_FSI_E - pre_FSI_E)

            if eventType == EventType.SRC:
                particles, energies, flags = nvect_class.particles(flag=True)
                pre_FSI_E = nvect_class.excitation_E_SRC(carbon=True)
                post_FSI_E = nvect_class.post_casc_E_star()
                if 10 in flags:
                    deex_src.append(post_FSI_E - pre_FSI_E)
                else:
                    if post_FSI_E - pre_FSI_E:
                        nvect_class.Print()
                        input(" ")
                    nondeex_src.append(post_FSI_E - pre_FSI_E)

        f.Close()
        return {
            "nondeex_MF": nondeex_mf,
            "deex_MF": deex_mf,
            "nondeex_SRC": nondeex_src,
            "deex_SRC": deex_src,
        }

    def make_hist(name, data, style):
        x_min, x_max = -100.0, 50.0

        h = ROOT.TH1D(
            name,
            ";E_{post-FSI} - E_{pre-FSI};Events",
            50,
            x_min,
            x_max,
        )
        h.SetDirectory(0)
        h.SetStats(0)
        h.SetLineColor(style["color"])
        h.SetLineStyle(style["linestyle"])
        h.SetLineWidth(2)

        h.GetXaxis().CenterTitle(True)
        h.GetXaxis().SetNoExponent(True)
        h.GetXaxis().SetTitleSize(0.055)
        h.GetXaxis().SetLabelSize(0.045)
        h.GetXaxis().SetTitleOffset(1.05)
        h.GetXaxis().SetNdivisions(606, True)

        h.GetYaxis().CenterTitle(True)
        h.GetYaxis().SetMaxDigits(4)
        h.GetYaxis().SetTitleSize(0.055)
        h.GetYaxis().SetLabelSize(0.045)
        h.GetYaxis().SetTitleOffset(1.25)

        for val in data:
            h.Fill(val)
        return h

    # Ordered to match the 2x2 pad layout:
    # Pad 1: nondeex_MF | Pad 2: deex_MF
    # Pad 3: nondeex_SRC | Pad 4: deex_SRC
    selection_keys = ["nondeex_MF", "deex_MF", "nondeex_SRC", "deex_SRC"]
    pad_histograms = {key: [] for key in selection_keys}
    legend_entries = []

    for file_idx, filepath in enumerate(all_files):
        diffs = read_file_diffs(filepath)
        if diffs is None:
            continue

        style = styles.get(
            file_idx,
            {
                "label": os.path.basename(str(filepath)),
                "linestyle": (file_idx % 9) + 1,
                "color": file_idx + 1,
            },
        )

        for key in selection_keys:
            h = make_hist(f"h_{key}_file{file_idx}", diffs[key], style)
            pad_histograms[key].append(h)

        # Track one histogram per file for the shared legend
        legend_entries.append((pad_histograms["deex_MF"][-1], style["label"]))

    if not legend_entries:
        print("No valid files were processed.")
        return None, None

    # Enlarged canvas with space on left (6%) and top (6%) for row/column titles
    canvas = ROOT.TCanvas("c_fsi_diff", "FSI Excitation Energy Differences", 1350, 1350)

    left_offset = 0.06
    top_offset = 0.06
    grid_pad = ROOT.TPad("grid_pad", "", left_offset, 0.0, 1.0, 1.0 - top_offset)
    grid_pad.Draw()
    grid_pad.Divide(2, 2, 0.002, 0.002)

    for pad_idx, key in enumerate(selection_keys, start=1):
        pad = grid_pad.cd(pad_idx)
        pad.SetLeftMargin(0.15)
        pad.SetBottomMargin(0.13)
        pad.SetRightMargin(0.03)
        pad.SetTopMargin(0.04)
        pad.SetGrid()

        hists = pad_histograms[key]
        # Scale y-axis to the tallest histogram in the overlay so nothing clips
        max_y = max(h.GetMaximum() for h in hists)
        hists[0].SetMinimum(0.0)
        hists[0].SetMaximum(max_y * 1.15 if max_y > 0 else 1.0)

        for idx, hist in enumerate(hists):
            draw_opt = "HIST" if idx == 0 else "HIST SAME"
            hist.Draw(draw_opt)

        if pad_idx == 2:
            leg_height = 0.055 * len(legend_entries)
            leg = ROOT.TLegend(0.56, 0.93 - leg_height, 0.94, 0.93)
            leg.SetBorderSize(1)
            leg.SetFillColor(ROOT.kWhite)
            leg.SetFillStyle(1001)
            leg.SetTextFont(42)
            leg.SetTextSize(0.04)
            for hist, label in legend_entries:
                leg.AddEntry(hist, label, "l")
            leg.Draw()
            canvas._legend = leg

    canvas.cd()
    latex = ROOT.TLatex()
    latex.SetTextFont(62)
    latex.SetTextSize(0.032)
    latex.SetTextAlign(22)

    latex.DrawLatexNDC(0.323, 0.968, "(non deex)")
    latex.DrawLatexNDC(0.793, 0.968, "(deex)")

    latex.SetTextAngle(90)
    latex.DrawLatexNDC(0.03, 0.726, "MF")
    latex.DrawLatexNDC(0.03, 0.256, "SRC")

    canvas.Update()

    output_dir = "Excitation_e_diff"
    os.makedirs(output_dir, exist_ok=True)

    name_str = output_filename[0]
    filepath_root = os.path.join(output_dir, f"SexcitationE_diff_{name_str}.root")
    filepath_pdf = os.path.join(output_dir, f"excitationE_diff_{name_str}.pdf")

    canvas.SaveAs(filepath_pdf)

    out_root = ROOT.TFile(filepath_root, "RECREATE")
    canvas.Write()
    out_root.Close()

    return canvas, pad_histograms