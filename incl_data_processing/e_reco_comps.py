from incl_analysis_functions import *
from default_processing import *
from kdar_processing import *
from paper_plots import *
import numpy as np
from array import array
import sys 
import argparse

import ROOT
import numpy as np

ROOT.gSystem.Load("libNEUTROOTClass.so")
ROOT.gSystem.Load("libNEUTOutput.so")
ROOT.gSystem.Load("libNEUTReWeight.so")
ROOT.TH1.AddDirectory(False)


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

ROOT.gStyle.SetPadTickX(1)
ROOT.gStyle.SetPadTickY(1)

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
ROOT.gStyle.SetOptStat(0)


def e_reco_comps(filename, name, n_events, filename_list=None):
    if filename_list is None:
        filename_list = []

    files_to_process = [filename] + filename_list

    hist_list_bias = []
    hist_list_qe_bias = []

    styles = [
        (ROOT.kRed, 1),  # Solid Red
        (ROOT.kBlue, 1),  # Solid Blueout
        (ROOT.kRed, 7),  # Solid Red
        (ROOT.kBlue, 7),  # Solid Blue
    ]

    legend_name = ["Carbon, INCL+ABLA", "Oxygen, INCL+ABLA", "Carbon, NEUT Cascade", "Oxygen, NEUT Cascade", ]
    #legend_name = ["Carbon, INBC"]

    for count, fname in enumerate(files_to_process):
        filename_chunk = legend_name[count]

        f = ROOT.TFile(fname)
        t = f.Get("neuttree")

        # 1. Histogram for energy bias (Visible, No Neutrons)
        h_bias = ROOT.TH1D(
            f"h_bias_vis1_neu0_{count}",
            "Energy Reco Bias (Vis, No Neutrons)",
            31,
            -0.51,
            0.11,
        )
        h_bias.Sumw2()

        # 2. Histogram for QE energy reconstruction bias
        h_qe_bias = ROOT.TH1D(
            f"h_qe_bias_{count}",
            "QE Energy Reco Bias",
            22,
            -0.75,
            0.55,
        )
        h_qe_bias.Sumw2()

        for i, event in enumerate(t):
            if i % 100000 == 0:
                print(i, " of ", n_events)
            elif i == n_events:
                break

            nvect = event.vectorbranch
            nvect_class = nvect_reader(nvect)

            if nvect_class.eventType in (EventType.SRC, EventType.MF):
                bias = nvect_class.E_reco_bias(visible=False, energy_thresh=0.0)
                if count == 0 or count == 2:
                    qe_bias = nvect_class.E_reco_qe_bias(carbon = True)
                else:
                    qe_bias = nvect_class.E_reco_qe_bias()

                h_bias.Fill(bias)
                h_qe_bias.Fill(qe_bias)

        h_bias.SetDirectory(0)
        h_qe_bias.SetDirectory(0)

        hist_list_bias.append((filename_chunk, h_bias))
        hist_list_qe_bias.append((filename_chunk, h_qe_bias))
        f.Close()

    ROOT.gStyle.SetOptTitle(0)
    out_file = ROOT.TFile(f"paper_e_reco_bias_{name[0]}.root", "RECREATE")
    latex_font = 132

    def draw_canvas_with_ratio(
        canvas_name,
        hist_list,
        xtitle="(E^{rec,avail}_{#nu} - E^{True}_{#nu}) / E^{True}_{#nu}",
        x_range=None,
        ymax=None,
        ymin=None,
        ratio_ymin=0.5,
        ratio_ymax=1.5,
    ):
        c = ROOT.TCanvas(canvas_name, canvas_name, 800, 800)

        # Main Pad
        pad1 = ROOT.TPad("pad1", "pad1", 0, 0.3, 1, 1.0)
        pad1.SetBottomMargin(0.02)
        pad1.SetLeftMargin(0.12)
        pad1.SetTickx(1)
        pad1.SetTicky(1)
        pad1.Draw()
        pad1.cd()

        legend = ROOT.TLegend(0.20, 0.65, 0.45, 0.85)
        legend.SetNColumns(1)
        legend.SetTextFont(latex_font)
        legend.SetBorderSize(0)
        legend.SetFillStyle(0)
        legend.SetTextSize(0.04)

        valid_hists = []
        for h_name, h in hist_list:
            if ymax is None:
                integral = h.Integral()
                if integral > 0:
                    h.Scale(1.0 / n_events)
            valid_hists.append((h_name, h))

        if ymax is None and valid_hists:
            global_max = max(h.GetMaximum() for _, h in valid_hists)
            ymax = 0.3

        # Plot main histograms on pad1
        for i, (name_label, h) in enumerate(valid_hists):
            color, style = styles[i % len(styles)]

            h.SetLineColor(color)
            h.SetLineStyle(style)
            h.SetLineWidth(3)
            h.SetStats(0)

            if x_range is not None:
                h.GetXaxis().SetRangeUser(x_range[0], x_range[1])

            h.GetYaxis().CenterTitle()
            h.GetXaxis().SetLabelSize(0)
            h.GetXaxis().SetTitleSize(0)

            h.GetYaxis().SetTitleFont(latex_font)
            h.GetYaxis().SetLabelFont(latex_font)
            h.GetYaxis().SetTitle("[A.U]")
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

            legend.AddEntry(h, name_label, "l")

        legend.Draw()

        # Ratio Pad
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
            name_label, h = hist_list[i]
            h_ratio = h.Clone(f"{h.GetName()}_ratio")
            h_ratio.Divide(h_base)

            color, style = styles[i % len(styles)]
            h_ratio.SetLineColor(color)
            h_ratio.SetLineStyle(style)
            h_ratio.SetLineWidth(3)
            h_ratio.SetStats(0)

            h_ratio.SetTitle("")
            h_ratio.GetYaxis().SetTitle(f"Ratio to {hist_list[0][0]}")
            h_ratio.GetYaxis().CenterTitle()
            h_ratio.GetYaxis().SetNdivisions(505)
            h_ratio.GetYaxis().SetTitleSize(0.11)
            h_ratio.GetYaxis().SetTitleFont(latex_font)
            h_ratio.GetYaxis().SetTitleOffset(0.45)
            h_ratio.GetYaxis().SetLabelFont(latex_font)
            h_ratio.GetYaxis().SetLabelSize(0.1)

            if x_range is not None:
                h_ratio.GetXaxis().SetRangeUser(x_range[0], x_range[1])

            h_ratio.GetXaxis().SetTitle(xtitle)
            h_ratio.GetXaxis().SetTitleSize(0.11)
            h_ratio.GetXaxis().SetTitleFont(latex_font)
            h_ratio.GetXaxis().SetTitleOffset(1.0)
            h_ratio.GetXaxis().CenterTitle()
            h_ratio.GetXaxis().SetLabelFont(latex_font)
            h_ratio.GetXaxis().SetLabelSize(0.1)

            if i == 1:
                h_ratio.SetMinimum(ratio_ymin)
                h_ratio.SetMaximum(ratio_ymax)
                h_ratio.Draw("HIST")
            else:
                h_ratio.Draw("HIST SAME")
            ratios.append(h_ratio)

        # Unit ratio line
        xmin_line = x_range[0] if x_range else h_base.GetXaxis().GetXmin()
        xmax_line = x_range[1] if x_range else h_base.GetXaxis().GetXmax()

        line = ROOT.TLine(xmin_line, 1, xmax_line, 1)
        line.SetLineColor(ROOT.kRed)
        line.SetLineStyle(2)
        line.Draw("SAME")

        c.pad1 = pad1
        c.pad2 = pad2
        c.ratios = ratios
        c.line = line

        return c, legend

    # Page 1: Energy Reco Bias Canvas
    c_vis1_neu0, leg1 = draw_canvas_with_ratio(
        "c_e_reco_bias_vis1_neu0",
        hist_list_bias,
        xtitle="(E^{rec,avail}_{#nu} - E^{True}_{#nu}) / E^{True}_{#nu}",
    #    x_range=(-0.51, 0.11),
    )

    # Page 2: QE Energy Reco Bias Canvas
    c_qe_bias, leg2 = draw_canvas_with_ratio(
        "c_e_reco_qe_bias",
        hist_list_qe_bias,
        xtitle="(E^{rec, QE}_{#nu} - E^{True}_{#nu}) / E^{True}_{#nu}",
    #    x_range=(-0.51, 0.11),
    )

    # Save to ROOT file
    out_file.cd()
    c_vis1_neu0.Write()
    c_qe_bias.Write()
    out_file.Close()

    # Save both pages into a single PDF
    output_dir = "E_reco_comps"
    os.makedirs(output_dir, exist_ok=True)
    pdf_filename = f"{output_dir}/paper_e_reco_bias_{name[0]}.pdf"

    c_vis1_neu0.Print(f"{pdf_filename}(", "pdf")  # Open multi-page PDF with '('
    c_qe_bias.Print(f"{pdf_filename})", "pdf")  # Close multi-page PDF with ')'



def e_reco_threshold_comps(filename, name, n_events, filename_list=None):
    if filename_list is None:
        filename_list = []

    files_to_process = [filename] + filename_list
    # 8 energy threshold values (0 to 175 MeV in 25 MeV increments)
    thresholds = [0, 25, 50, 75, 100, 125, 150, 175]

    styles = [
        (ROOT.kRed, 1),  # Solid Red
        (ROOT.kBlue, 1),  # Solid Blue
    ]

    styles = [
        (ROOT.kRed, 7),  # Solid Red
        (ROOT.kBlue, 7),  # Solid Blue
    ]

    legend_name = ["Carbon, INCL+ABLA", "Oxygen, INCL+ABLA", "Carbon, NEUT Cascade", "Oxygen, NEUT Cascade", ]

    # Distinct color palette for the 8 thresholds on Page 2
    thresh_colors = [
        ROOT.kBlack,
        ROOT.kRed + 1,
        ROOT.kBlue + 1,
        ROOT.kGreen + 2,
        ROOT.kMagenta + 1,
        ROOT.kOrange + 7,
        ROOT.kCyan + 1,
        ROOT.kViolet + 1,
    ]

    # Dict to hold histograms: threshold -> list of (legend_label, TH1D)
    hist_dict = {t: [] for t in thresholds}

    for count, fname in enumerate(files_to_process):
        filename_chunk = (
            legend_name[count]
            if count < len(legend_name)
            else f"File_{count}"
        )

        f = ROOT.TFile(fname)
        t_tree = f.Get("neuttree")

        # Create histograms for each threshold for the current file
        h_bias_map = {}
        for t_val in thresholds:
            h = ROOT.TH1D(
                f"h_bias_thresh_{t_val}_file_{count}",
                f"Energy Reco Bias (Thresh {t_val} MeV)",
                31,
                -0.51,
                0.11,
            )
            h.Sumw2()
            h_bias_map[t_val] = h

        # Event processing loop
        for i, event in enumerate(t_tree):
            if i % 100000 == 0:
                print(i, " of ", n_events)
            elif i == n_events:
                break

            nvect = event.vectorbranch
            nvect_class = nvect_reader(nvect)

            if nvect_class.eventType in (EventType.SRC, EventType.MF):
                for t_val in thresholds:
                    bias = nvect_class.E_reco_bias(
                        visible=False, energy_thresh=t_val
                    )
                    h_bias_map[t_val].Fill(bias)

        for t_val in thresholds:
            h = h_bias_map[t_val]
            h.SetDirectory(0)
            hist_dict[t_val].append((filename_chunk, h))

        f.Close()

    ROOT.gStyle.SetOptTitle(0)
    latex_font = 132

    # Scale histograms and find global maximum Y across ALL 8 thresholds
    global_ymax = 0.0
    for t_val in thresholds:
        for _, h in hist_dict[t_val]:
            if n_events > 0 and h.Integral() > 0:
                h.Scale(1.0 / n_events)
            if h.GetMaximum() > global_ymax:
                global_ymax = h.GetMaximum()

    ymax_uniform = global_ymax * 1.25 if global_ymax > 0 else 1.0

    # =========================================================================
    # PAGE 1: 4x2 Grid of Target Comparisons (Carbon vs Oxygen per threshold)
    # =========================================================================
    c_page1 = ROOT.TCanvas(
        "c_e_reco_threshold_comps",
        "E Reco Bias Across Energy Thresholds",
        2000,
        1200,
    )
    c_page1.Divide(4, 2)

    subpads_p1 = []
    ratios_p1 = []
    lines_p1 = []
    legends_p1 = []
    latexs_p1 = []

    for pad_idx, t_val in enumerate(thresholds, start=1):
        cell_pad = c_page1.cd(pad_idx)

        pad1 = ROOT.TPad(
            f"pad1_p1_{pad_idx}", f"pad1_p1_{pad_idx}", 0.0, 0.3, 1.0, 1.0
        )
        pad2 = ROOT.TPad(
            f"pad2_p1_{pad_idx}", f"pad2_p1_{pad_idx}", 0.0, 0.0, 1.0, 0.3
        )

        pad1.SetTopMargin(0.08)
        pad1.SetBottomMargin(0.02)
        pad1.SetLeftMargin(0.16)
        pad1.SetRightMargin(0.05)
        pad1.SetTickx(1)
        pad1.SetTicky(1)

        pad2.SetTopMargin(0.02)
        pad2.SetBottomMargin(0.35)
        pad2.SetLeftMargin(0.16)
        pad2.SetRightMargin(0.05)
        pad2.SetTickx(1)
        pad2.SetTicky(1)

        cell_pad.cd()
        pad1.Draw()
        pad2.Draw()
        subpads_p1.extend([pad1, pad2])

        hist_list = hist_dict[t_val]

        # --- Pad 1: Main Histograms ---
        pad1.cd()
        legend = ROOT.TLegend(0.20, 0.52, 0.60, 0.75)
        legend.SetTextFont(latex_font)
        legend.SetBorderSize(0)
        legend.SetFillStyle(0)
        legend.SetTextSize(0.065)
        legends_p1.append(legend)

        for i, (name_label, h) in enumerate(hist_list):
            color, style = styles[i % len(styles)]
            h.SetLineColor(color)
            h.SetLineStyle(style)
            h.SetLineWidth(2)
            h.SetStats(0)

            h.GetXaxis().SetLabelSize(0)
            h.GetXaxis().SetTitleSize(0)

            h.GetYaxis().SetTitle("[A.U]")
            h.GetYaxis().SetTitleFont(latex_font)
            h.GetYaxis().SetLabelFont(latex_font)
            h.GetYaxis().SetTitleSize(0.07)
            h.GetYaxis().SetLabelSize(0.06)
            h.GetYaxis().SetTitleOffset(1.1)

            if i == 0:
                h.SetMinimum(0)
                h.SetMaximum(ymax_uniform)
                h.Draw("HIST")
            else:
                h.Draw("HIST SAME")

            legend.AddEntry(h, name_label, "l")

        legend.Draw()

        text = ROOT.TLatex()
        text.SetNDC()
        text.SetTextFont(latex_font)
        text.SetTextSize(0.07)
        text.DrawLatex(0.20, 0.80, f"E_{{thresh}} = {t_val} MeV")
        latexs_p1.append(text)

        # --- Pad 2: Target Ratio Plot ---
        pad2.cd()
        h_base = hist_list[0][1]

        for i in range(1, len(hist_list)):
            name_label, h = hist_list[i]
            h_ratio = h.Clone(f"{h.GetName()}_ratio_p1")
            h_ratio.Divide(h_base)

            color, style = styles[i % len(styles)]
            h_ratio.SetLineColor(color)
            h_ratio.SetLineStyle(style)
            h_ratio.SetLineWidth(2)
            h_ratio.SetStats(0)

            h_ratio.SetTitle("")
            h_ratio.GetYaxis().SetTitle(f"Ratio to {hist_list[0][0]}")
            h_ratio.GetYaxis().CenterTitle()
            h_ratio.GetYaxis().SetNdivisions(505)
            h_ratio.GetYaxis().SetTitleSize(0.12)
            h_ratio.GetYaxis().SetTitleFont(latex_font)
            h_ratio.GetYaxis().SetTitleOffset(0.55)
            h_ratio.GetYaxis().SetLabelFont(latex_font)
            h_ratio.GetYaxis().SetLabelSize(0.11)

            h_ratio.GetXaxis().SetTitle(
                "(E^{rec,avail}_{#nu} - E^{True}_{#nu}) / E^{True}_{#nu}"
            )
            h_ratio.GetXaxis().SetTitleSize(0.12)
            h_ratio.GetXaxis().SetTitleFont(latex_font)
            h_ratio.GetXaxis().SetTitleOffset(1.1)
            h_ratio.GetXaxis().CenterTitle()
            h_ratio.GetXaxis().SetLabelFont(latex_font)
            h_ratio.GetXaxis().SetLabelSize(0.11)

            if i == 1:
                h_ratio.SetMinimum(0.5)
                h_ratio.SetMaximum(1.5)
                h_ratio.Draw("HIST")
            else:
                h_ratio.Draw("HIST SAME")

            ratios_p1.append(h_ratio)

        xmin_line = h_base.GetXaxis().GetXmin()
        xmax_line = h_base.GetXaxis().GetXmax()
        line = ROOT.TLine(xmin_line, 1.0, xmax_line, 1.0)
        line.SetLineColor(ROOT.kRed)
        line.SetLineStyle(2)
        line.Draw("SAME")
        lines_p1.append(line)

    # =========================================================================
    # PAGE 2: Threshold Comparisons per Target (Left: Carbon, Right: Oxygen)
    # =========================================================================
    c_page2 = ROOT.TCanvas(
        "c_e_reco_threshold_targets",
        "Threshold Overlay per Target",
        1600,
        800,
    )
    c_page2.Divide(2, 1)

    subpads_p2 = []
    ratios_p2 = []
    lines_p2 = []
    legends_p2 = []
    latexs_p2 = []

    # Organize histograms by target index: 0 = Carbon, 1 = Oxygen
    target_names = [legend_name[0], legend_name[1]]

    for target_idx, target_label in enumerate(target_names):
        cell_pad = c_page2.cd(target_idx + 1)

        pad1 = ROOT.TPad(
            f"pad1_p2_{target_idx}", f"pad1_p2_{target_idx}", 0.0, 0.3, 1.0, 1.0
        )
        pad2 = ROOT.TPad(
            f"pad2_p2_{target_idx}", f"pad2_p2_{target_idx}", 0.0, 0.0, 1.0, 0.3
        )

        pad1.SetTopMargin(0.08)
        pad1.SetBottomMargin(0.02)
        pad1.SetLeftMargin(0.14)
        pad1.SetRightMargin(0.05)
        pad1.SetTickx(1)
        pad1.SetTicky(1)

        pad2.SetTopMargin(0.02)
        pad2.SetBottomMargin(0.35)
        pad2.SetLeftMargin(0.14)
        pad2.SetRightMargin(0.05)
        pad2.SetTickx(1)
        pad2.SetTicky(1)

        cell_pad.cd()
        pad1.Draw()
        pad2.Draw()
        subpads_p2.extend([pad1, pad2])

        # Retrieve histograms across all 8 thresholds for current target
        target_hists = [
            hist_dict[t_val][target_idx][1] for t_val in thresholds
        ]

        # --- Pad 1: Overlay All 8 Threshold Histograms ---
        pad1.cd()
        legend = ROOT.TLegend(0.18, 0.45, 0.58, 0.88)
        legend.SetNColumns(2)
        legend.SetTextFont(latex_font)
        legend.SetBorderSize(0)
        legend.SetFillStyle(0)
        legend.SetTextSize(0.045)
        legends_p2.append(legend)

        for i, (t_val, h) in enumerate(zip(thresholds, target_hists)):
            color = thresh_colors[i % len(thresh_colors)]
            h_clone = h.Clone(f"{h.GetName()}_p2_main")
            h_clone.SetLineColor(color)
            h_clone.SetLineWidth(2)
            h_clone.SetStats(0)

            h_clone.GetXaxis().SetLabelSize(0)
            h_clone.GetXaxis().SetTitleSize(0)

            h_clone.GetYaxis().SetTitle("[A.U]")
            h_clone.GetYaxis().SetTitleFont(latex_font)
            h_clone.GetYaxis().SetLabelFont(latex_font)
            h_clone.GetYaxis().SetTitleSize(0.06)
            h_clone.GetYaxis().SetLabelSize(0.05)
            h_clone.GetYaxis().SetTitleOffset(1.1)

            if i == 0:
                h_clone.SetMinimum(0)
                h_clone.SetMaximum(ymax_uniform)
                h_clone.Draw("HIST")
            else:
                h_clone.Draw("HIST SAME")

            legend.AddEntry(h_clone, f"E_{{thresh}} = {t_val} MeV", "l")
            ratios_p2.append(h_clone)  # Keep in memory

        legend.Draw()

        # Target label title
        text = ROOT.TLatex()
        text.SetNDC()
        text.SetTextFont(latex_font)
        text.SetTextSize(0.06)
        text.DrawLatex(0.20, 0.82, f"Target: {target_label}")
        latexs_p2.append(text)

        # --- Pad 2: Threshold Ratios (relative to baseline E_thresh = 0 MeV) ---
        pad2.cd()
        h_base_target = target_hists[0]

        for i, (t_val, h) in enumerate(zip(thresholds, target_hists)):
            h_ratio = h.Clone(f"{h.GetName()}_ratio_p2")
            h_ratio.Divide(h_base_target)

            color = thresh_colors[i % len(thresh_colors)]
            h_ratio.SetLineColor(color)
            h_ratio.SetLineWidth(2)
            h_ratio.SetStats(0)

            h_ratio.SetTitle("")
            h_ratio.GetYaxis().SetTitle("Ratio to 0 MeV")
            h_ratio.GetYaxis().CenterTitle()
            h_ratio.GetYaxis().SetNdivisions(505)
            h_ratio.GetYaxis().SetTitleSize(0.10)
            h_ratio.GetYaxis().SetTitleFont(latex_font)
            h_ratio.GetYaxis().SetTitleOffset(0.55)
            h_ratio.GetYaxis().SetLabelFont(latex_font)
            h_ratio.GetYaxis().SetLabelSize(0.09)

            h_ratio.GetXaxis().SetTitle(
                "(E^{rec,avail}_{#nu} - E^{True}_{#nu}) / E^{True}_{#nu}"
            )
            h_ratio.GetXaxis().SetTitleSize(0.10)
            h_ratio.GetXaxis().SetTitleFont(latex_font)
            h_ratio.GetXaxis().SetTitleOffset(1.1)
            h_ratio.GetXaxis().CenterTitle()
            h_ratio.GetXaxis().SetLabelFont(latex_font)
            h_ratio.GetXaxis().SetLabelSize(0.09)

            if i == 0:
                h_ratio.SetMinimum(0.5)
                h_ratio.SetMaximum(1.5)
                h_ratio.Draw("HIST")
            else:
                h_ratio.Draw("HIST SAME")

            ratios_p2.append(h_ratio)

        xmin_line = h_base_target.GetXaxis().GetXmin()
        xmax_line = h_base_target.GetXaxis().GetXmax()
        line = ROOT.TLine(xmin_line, 1.0, xmax_line, 1.0)
        line.SetLineColor(ROOT.kRed)
        line.SetLineStyle(2)
        line.Draw("SAME")
        lines_p2.append(line)

    # Save both canvases to ROOT file
    out_file = ROOT.TFile(
        f"paper_e_reco_threshold_comps_{name[0]}.root", "RECREATE"
    )
    c_page1.Write()
    c_page2.Write()
    out_file.Close()

    # Save to 2-page PDF
    output_dir = "E_reco_comps"
    os.makedirs(output_dir, exist_ok=True)
    pdf_filename = f"{output_dir}/paper_e_reco_threshold_comps_{name[0]}.pdf"

    c_page1.Print(f"{pdf_filename}(", "pdf")  # Open multi-page PDF with '('
    c_page2.Print(f"{pdf_filename})", "pdf")  # Close multi-page PDF with ')'

    # Store references on canvas objects to prevent C++ garbage collection
    c_page1.subpads = subpads_p1
    c_page1.legends = legends_p1
    c_page1.latexs = latexs_p1
    c_page1.ratios = ratios_p1
    c_page1.lines = lines_p1

    c_page2.subpads = subpads_p2
    c_page2.legends = legends_p2
    c_page2.latexs = latexs_p2
    c_page2.ratios = ratios_p2
    c_page2.lines = lines_p2

    return c_page1, c_page2