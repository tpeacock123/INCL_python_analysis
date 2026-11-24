from incl_analysis_functions import *
import ROOT
import numpy as np
import array
import sys 
import argparse
import os
import re 

from pathlib import Path


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




def transparency_neutcas():
    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(0)

    sample_names = ["NEUT_cascade", "INCL+ABLA"]
    colors = [ROOT.kAzure, ROOT.kBlue]
    markers = [1, 1]

    c1 = ROOT.TCanvas("c1", "Transparency", 800, 600)
    c1.SetLeftMargin(0.14)
    c1.SetBottomMargin(0.14)
    c1.SetRightMargin(0.06)
    c1.SetTopMargin(0.06)
    c1.SetGridy()

    y_min = max(0.50, 0.88 - 0.045 * 2)
    legend = ROOT.TLegend(0.58, y_min, 0.88, 0.88)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.SetTextFont(font_id)
    legend.SetTextSize(0.035)

    ratios = []
    tot_hists = []
    arb_hists = []

    


    for idx, (fpath, sname) in enumerate(zip(all_files, sample_names)):
        f = ROOT.TFile.Open(fpath)
        if not f or f.IsZombie():
            print(f"Error opening file: {fpath}")
            continue

        t = f.Get("neuttree")
        if not t:
            print(f"Error: 'neuttree' branch not found in {fpath}.")
            f.Close()
            continue

        t.SetBranchStatus("*", 0)
        t.SetBranchStatus("vectorbranch", 1)

        h1_tot = ROOT.TH1D(f"h1_tot_{idx}", f";T_{{p}} [MeV];Events", 80, 0, 400)
        h1_arb = ROOT.TH1D(f"h1_arb_{idx}", f";T_{{p}} [MeV];Events", 80, 0, 400)
        h1_tot.SetDirectory(0)
        h1_arb.SetDirectory(0)

        for i, event in enumerate(t):
            if i >= n_events:
                break
            nvect_class = nvect_reader(event.vectorbranch)
            HMPMom = nvect_class.HMPMom
            if HMPMom is None:
                continue

            T_p = math.sqrt(HMPMom**2 + 938.272**2) - 938.272  
            h1_tot.Fill(T_p)

            particle_list, energy_list, flag_list = nvect_class.particles(flag=True)

            if nvect_class.nopart == 4:
                h1_arb.Fill(T_p)
                continue
    
            for particle, energy, flag in zip(particle_list, energy_list, flag_list):
                if particle == 2212 and flag == 2:
                    h1_arb.Fill(T_p)

        f.Close()

        h_ratio = h1_arb.Clone(f"h_ratio_{idx}")
        h_ratio.SetDirectory(0)
        h_ratio.Divide(h1_arb, h1_tot, 1.0, 1.0, "B")

        color = colors[idx % len(colors)]
        marker = markers[idx % len(markers)]

        h_ratio.SetLineColor(color)
        h_ratio.SetMarkerColor(color)
        h_ratio.SetMarkerStyle(marker)
        h_ratio.SetMarkerSize(0.8)
        h_ratio.SetLineWidth(2)

        if idx == 0:
            h_ratio.SetTitle("")
            h_ratio.GetXaxis().SetTitle("T_{p} [MeV]")
            h_ratio.GetYaxis().SetTitle("Transparency")
            h_ratio.GetXaxis().SetTitleSize(0.045)
            h_ratio.GetYaxis().SetTitleSize(0.045)
            h_ratio.GetXaxis().SetLabelSize(0.04)
            h_ratio.GetYaxis().SetLabelSize(0.04)
            h_ratio.GetXaxis().SetTitleOffset(1.1)
            h_ratio.GetYaxis().SetTitleOffset(1.3)
            h_ratio.SetMinimum(-0.02)
            h_ratio.SetMaximum(1.1)
            h_ratio.Draw("EP")
        else:
            h_ratio.Draw("EP SAME")

        legend.AddEntry(h_ratio, f"{sname}", "lp")
        ratios.append(h_ratio)
        tot_hists.append((sname, h1_tot))
        arb_hists.append((sname, h1_arb))

    legend.Draw()

    out_dir = "neutcasc_transparency"
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, f"transparency_{name[0]}.pdf")

    # --- Page 1: Ratio plot (Open multi-page PDF with '(') ---
    c1.Print(f"{out_path}(", "pdf")

    # --- Pages 2+: Raw Histograms for each sample ---
    for i, ((sname, h_tot), (_, h_arb)) in enumerate(zip(tot_hists, arb_hists)):
        c1.Clear()
        c1.SetGridy(False)

        color = colors[i % len(colors)]

        # Auto-scale Y axis with padding
        max_y = max(h_tot.GetMaximum(), h_arb.GetMaximum()) * 1.15
        h_tot.SetMaximum(max_y)
        h_tot.SetMinimum(0)

        # Style Total (Denominator)
        h_tot.SetLineColor(dark_gray)
        h_tot.SetLineWidth(2)
        h_tot.GetXaxis().SetTitle("T_{p} [MeV]")
        h_tot.GetYaxis().SetTitle("Events")
        h_tot.GetXaxis().SetTitleSize(0.045)
        h_tot.GetYaxis().SetTitleSize(0.045)
        h_tot.GetXaxis().SetLabelSize(0.04)
        h_tot.GetYaxis().SetLabelSize(0.04)
        h_tot.GetXaxis().SetTitleOffset(1.1)
        h_tot.GetYaxis().SetTitleOffset(1.3)

        # Style Selected (Numerator)
        h_arb.SetLineColor(color)
        h_arb.SetLineWidth(2)

        h_tot.Draw("HIST")
        h_arb.Draw("HIST SAME")

        raw_legend = ROOT.TLegend(0.55, 0.72, 0.88, 0.88)
        raw_legend.SetBorderSize(0)
        raw_legend.SetFillStyle(0)
        raw_legend.SetTextFont(font_id)
        raw_legend.SetTextSize(0.035)
        raw_legend.SetHeader(f"Raw Counts ({sname})", "C")
        raw_legend.AddEntry(h_tot, "Total (h1_tot)", "l")
        raw_legend.AddEntry(h_arb, "Selected (h1_arb)", "l")
        raw_legend.Draw()

        # Close PDF stream on the last page with ')'
        if i == len(tot_hists) - 1:
            c1.Print(f"{out_path})", "pdf")
        else:
            c1.Print(out_path, "pdf")





import os
import re
from pathlib import Path
import ROOT



def transparency_calc(name):
    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(0)

    directory = Path("./Transparency")

    samples = [

        ("NEUT_carbon",   "NEUT cascade",   ROOT.kAzure , 23),  # 21 = Filled Square
        ("INCL_carbon",  "INCL+ABLA", ROOT.kBlue , 21),  # 21 = Filled Square

    ]

    mg = ROOT.TMultiGraph()
    legend = ROOT.TLegend(0.58, 0.75, 0.88, 0.88)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.SetTextSize(0.055)

    graphs = []

    for prefix, label, color, marker in samples:
        pattern = re.compile(rf"^out_neutcas_{prefix}_trans_(\d+\.\d+)\.root$")
        numbers = []
        trans_ratio = []

        for file in directory.glob(f"out_neutcas_{prefix}_trans_*.root"):
            match = pattern.match(file.name)
            if match:
                val = float(match.group(1))
                ratio = transparency_ratio(str(file))
                numbers.append(val*1000)
                trans_ratio.append(ratio)
                print(trans_ratio)

        if not numbers:
            print(f"Warning: No matching files found for pattern prefix '{prefix}'!")
            continue

        sorted_pairs = sorted(zip(numbers, trans_ratio))

        graph = ROOT.TGraph(len(sorted_pairs))
        for idx, (x, y) in enumerate(sorted_pairs):
            graph.SetPoint(idx, x, y)

        graph.SetLineColor(color)
        graph.SetLineWidth(2)
        graph.SetMarkerColor(color)
        graph.SetMarkerStyle(marker)
        graph.SetMarkerSize(1.0)

        mg.Add(graph)
        legend.AddEntry(graph, label, "lp")
        graphs.append(graph)

    if not graphs:
        print("Error: No data found for any dataset.")
        return

    c1 = ROOT.TCanvas("c1", "Transparency vs Value", 800, 600)
    c1.SetLeftMargin(0.14)
    c1.SetBottomMargin(0.14)
    c1.SetRightMargin(0.06)
    c1.SetTopMargin(0.06)
    c1.SetGridx()
    c1.SetGridy()


    mg.Draw("APL")
    mg.SetTitle("")
    mg.GetXaxis().SetTitle("T_{p} (Initial Proton) [MeV]")
    mg.GetXaxis().CenterTitle()
    mg.GetYaxis().SetTitle("Transparency")
    mg.GetYaxis().CenterTitle()
    mg.GetXaxis().SetTitleSize(0.06)
    mg.GetYaxis().SetTitleSize(0.06)
    mg.GetXaxis().SetLabelSize(0.06)
    mg.GetYaxis().SetLabelSize(0.06)
    mg.GetXaxis().SetTitleOffset(1)
    mg.GetYaxis().SetTitleOffset(1)
    mg.SetMinimum(0.4)
    mg.SetMaximum(1.05)

    legend.Draw()

    out_dir = "neutcasc_transparency"
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, f"transparency_{name[0]}.pdf")
    c1.SaveAs(out_path)
    print(f"Plot saved to {out_path}")


def transparency_ratio(file):
    f = ROOT.TFile.Open(file)
    if not f or f.IsZombie():
        print(f"Error opening file: {file}")
        return 0
    
    t = f.Get("neuttree")
    if not t:
        print(f"Error: 'neuttree' branch not found in {file}.")
        f.Close()
        return 0

    t.SetBranchStatus("*", 0)
    t.SetBranchStatus("vectorbranch", 1)
    no_FSI_count = 0
    total_count  = 0

    for i, event in enumerate(t):
        total_count += 1
        nvect_class = nvect_reader(event.vectorbranch)

        particle_list, energy_list, flag_list = nvect_class.particles(flag=True)
        if nvect_class.nopart == 4:
            no_FSI_count += 1
            continue
        for particle, energy, flag in zip(particle_list, energy_list, flag_list):
            if particle == 2212 and flag == 2:
                no_FSI_count+=1
                break

    return no_FSI_count/total_count

        



