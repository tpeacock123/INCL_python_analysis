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


import ROOT
def deex_multiplicity(filename, filename_list=None):
    if filename_list is None:
        filename_list = [filename]
    else:
        filename_list = list(filename_list)
        if filename not in filename_list:
            filename_list.append(filename)

    # Setup output ROOT file
    filename_chunk = filename.split(".")[0].split("out_")[-1]
    out_file = ROOT.TFile(f"multiplicity_{filename_chunk}_SRC.root", "RECREATE")

    pdg_map = {
        22: "#gamma",
        2212: "p",
        2112: "n",
        1000010020: "d",
        1000010030: "t",
        1000020040: "#alpha",
        211: "pi+",
        -211: "pi-",
        111: "pi0"
    }

    colors = [
        ROOT.kBlack, ROOT.kRed + 1, ROOT.kBlue + 1, 
        ROOT.kGreen + 2, ROOT.kMagenta + 1, ROOT.kOrange + 7, ROOT.kCyan + 2
    ]

    histograms_1d = []
    canvases_2d = []
    histograms_2d = [] 

    for idx, fname in enumerate(filename_list):
        fname_chunk = fname.split(".")[0].split("out_")[-1]

        f_in = ROOT.TFile(fname)
        t = f_in.Get("neuttree")   

        particle_counts = {label: 0 for label in pdg_map.values()}
        counter = 0
        nbins = len(pdg_map)


        h2_energy = ROOT.TH2D(
            f"h2_energy_{fname_chunk}", 
            f"Energy vs Particle Species ({fname_chunk});Particle Species;Energy (GeV)", 
            nbins, 0, nbins, 
            200, 0.0, 750 
        )
        h2_energy.SetStats(0)

        for i, label in enumerate(pdg_map.values(), start=1):
            h2_energy.GetXaxis().SetBinLabel(i, label)

        for event in t:
            nvect = event.vectorbranch
            nvect_class = nvect_reader(nvect) 

            if nvect_class.eventType == EventType.SRC or nvect_class.eventType == EventType.twop2h: 
                counter += 1
                particles, energies = nvect_class.particles()
                
                for particle, p_energy in zip(particles, energies):
                    if particle in pdg_map:
                        label = pdg_map[particle]
                        particle_counts[label] += 1
                        h2_energy.Fill(label, p_energy, 1.0)

        f_in.Close()

        if counter == 0:
            print(f"Warning: No events processed for {fname}")
            continue

        # Setup 1D Multiplicity Histogram
        multiplicities = {label: count / counter for label, count in particle_counts.items()}
        
        h_mult = ROOT.TH1D(
            f"h_mult_{fname_chunk}", 
            "Average Particle Multiplicity per Event;Particle Species;Average Multiplicity / Event", 
            nbins, 0, nbins
        )

        for i, (label, mult_val) in enumerate(multiplicities.items(), start=1):
            h_mult.GetXaxis().SetBinLabel(i, label)
            h_mult.SetBinContent(i, mult_val)

        color = colors[idx % len(colors)]
        h_mult.SetLineColor(color)
        h_mult.SetLineWidth(3)
        h_mult.SetStats(0)
        h_mult.GetXaxis().SetLabelSize(0.05)
        h_mult.GetYaxis().SetTitleOffset(1.2)

        histograms_1d.append((fname_chunk, h_mult))

        c2 = ROOT.TCanvas(f"c2_energy_{fname_chunk}", f"Energy - {fname_chunk}", 800, 600)
        c2.SetRightMargin(0.15) # Room for colorbar
        h2_energy.Draw("COLZ")
        c2.Update()
        
        canvases_2d.append(c2)
        histograms_2d.append(h2_energy)  # <--- ADDED THIS: Keeps reference alive during the loop

    if not histograms_1d:
        print("No valid ROOT objects created.")
        out_file.Close()
        return

    canvas_mult = ROOT.TCanvas("c_multiplicity_summary", "Particle Multiplicity Comparison", 800, 600)
    canvas_mult.SetGridy()

    legend = ROOT.TLegend(0.68, 0.70, 0.88, 0.88)
    legend.SetBorderSize(1)
    legend.SetFillColor(0)

    max_y = max(h.GetMaximum() for _, h in histograms_1d)

    canvas_mult.cd()
    for idx, (fname_chunk, h_mult) in enumerate(histograms_1d):
        if idx == 0:
            h_mult.SetMaximum(3.0)
            h_mult.Draw("HIST")
        else:
            h_mult.Draw("HIST SAME")

        legend.AddEntry(h_mult, fname_chunk, "l")

    legend.Draw()
    canvas_mult.Update()

    out_file.cd()
    canvas_mult.Write()  

    for c2 in canvases_2d:
        c2.Write()       

    #for h2 in histograms_2d:
    #    h2.Write()  # <--- ADDED THIS: Writes raw 2D histograms to file as well

    out_file.Close()

    return canvas_mult, canvases_2d


def nocasc_pi_ana(filename, filename_list=None):
    if filename_list is None:
        filename_list = [filename]
    else:
        filename_list = list(filename_list)
        if filename not in filename_list:
            filename_list.append(filename)

    # Setup output ROOT file
    filename_chunk = filename.split(".")[0].split("out_")[-1]
    out_file = ROOT.TFile(f"nocasc_pion_ana_{filename_chunk}.root", "RECREATE")

    colors = [
        ROOT.kBlack, ROOT.kRed + 1, ROOT.kBlue + 1, 
        ROOT.kGreen + 2, ROOT.kMagenta + 1, ROOT.kOrange + 7, ROOT.kCyan + 2
    ]

    # --- 1D HISTOGRAM CONFIGURATIONS ---
    var_configs = {
        # Original Kinematics
        "pion_count": {"title": "Pion Multiplicity;Number of Pions;Events", "bins": 8, "min": 0, "max": 8},
        "pion_momenta": {"title": "Pion Momenta;p_{#pi} (MeV/c);Entries", "bins": 50, "min": 0.0, "max": 1500},
        "initProton": {"title": "Initial Nucleon Momenta;p (MeV/c);Entries", "bins": 50, "min": 0.0, "max": 1500},
        "pion_angles": {"title": "Pion Scattering Angle;#theta_{#pi} (rad);Entries", "bins": 50, "min": 0.0, "max": 3.2},
        "nucleon_outgoing_momentum": {"title": "Outgoing Nucleon Momentum;p_{N} (MeV/c);Events", "bins": 50, "min": 0.0, "max": 1500},
        "lepton_momentum": {"title": "Outgoing Lepton Momentum;p_{l} (MeV/c);Events", "bins": 50, "min": 0.0, "max": 5000},
        "p_miss": {"title": "Missing Momentum;p_{miss} (MeV/c);Events", "bins": 50, "min": 0.0, "max": 1500},
        "e_miss": {"title": "Missing Energy;E_{miss} (MeV);Events", "bins": 110, "min": -100, "max": 100},

        # Angular Kinematics
        "cos_theta_l": {"title": "Lepton Cosine Angle;cos(#theta_{l});Events", "bins": 50, "min": -1.0, "max": 1.0},
        "theta_l_pi": {"title": "Lepton-Pion Opening Angle;#theta_{l#pi} (rad);Events", "bins": 50, "min": 0.0, "max": 3.2},
        "theta_N_pi": {"title": "Nucleon-Pion Opening Angle;#theta_{N#pi} (rad);Events", "bins": 50, "min": 0.0, "max": 3.2},

        # Global & Invariant Kinematics
        "Q2": {"title": "Four-Momentum Transfer;Q^{2} (MeV^{2});Events", "bins": 50, "min": 0.0, "max": 2.0e6},
        "W": {"title": "Hadronic Invariant Mass;W (MeV);Events", "bins": 50, "min": 900.0, "max": 2000.0},
        "nu": {"title": "Energy Transfer;#nu (MeV);Events", "bins": 50, "min": 0.0, "max": 1500.0},

        # Single-Transverse Kinematic Variables (STVs)
        "delta_pT": {"title": "Transverse Momentum Imbalance;#delta p_{T} (MeV/c);Events", "bins": 50, "min": 0.0, "max": 1000.0},
        "delta_alphaT": {"title": "Transverse Imbalance Angle;#delta#alpha_{T} (rad);Events", "bins": 50, "min": 0.0, "max": 3.2},
        "delta_phiT": {"title": "Transverse Azimuthal Imbalance;#delta#phi_{T} (rad);Events", "bins": 50, "min": 0.0, "max": 3.2}
    }

    # --- 2D HISTOGRAM CONFIGURATIONS ---
    var_2d_configs = {
        "emiss_vs_pmiss": {
            "x_var": "p_miss", "y_var": "e_miss",
            "title": "Missing Energy vs Missing Momentum;p_{miss} (MeV/c);E_{miss} (MeV)",
            "xbins": 50, "xmin": 0.0, "xmax": 800.0,
            "ybins": 50, "ymin": -20.0, "ymax": 120.0
        },
        "nu_vs_Q2": {
            "x_var": "Q2", "y_var": "nu",
            "title": "Energy Transfer vs Q^{2};Q^{2} (MeV^{2}); q_{0} (MeV)",
            "xbins": 50, "xmin": 0.0, "xmax": 2.0e6,
            "ybins": 50, "ymin": 0.0, "ymax": 1500.0
        },
        "pl_vs_costhetal": {
            "x_var": "cos_theta_l", "y_var": "lepton_momentum",
            "title": "Lepton Momentum vs Lepton Cosine Angle;cos(#theta_{l});p_{l} (MeV/c)",
            "xbins": 50, "xmin": -1.0, "xmax": 1.0,
            "ybins": 50, "ymin": 0.0, "ymax": 5000.0
        },
        "deltapT_vs_deltaalphaT": {
            "x_var": "delta_alphaT", "y_var": "delta_pT",
            "title": "Transverse Momentum Imbalance vs Imbalance Angle;#delta#alpha_{T} (rad);#delta p_{T} (MeV/c)",
            "xbins": 50, "xmin": 0.0, "xmax": 3.2,
            "ybins": 50, "ymin": 0.0, "ymax": 800.0
        },
        "W_vs_Q2": {
            "x_var": "Q2", "y_var": "W",
            "title": "Hadronic Invariant Mass vs Q^{2};Q^{2} (MeV^{2});W (MeV)",
            "xbins": 50, "xmin": 0.0, "xmax": 2.0e6,
            "ybins": 50, "ymin": 900.0, "ymax": 2000.0
        }
    }

    # Storage structures for histograms
    hists = {var: {} for var in var_configs}
    hists_2d = {var_2d: {} for var_2d in var_2d_configs}

    for idx, fname in enumerate(filename_list):
        fname_chunk = fname.split(".")[0].split("out_")[-1]
        
        # Initialize 1D histograms
        for var, config in var_configs.items():
            h_name = f"h_{var}_{fname_chunk}"
            h = ROOT.TH1D(h_name, config["title"], config["bins"], config["min"], config["max"])
            h.SetLineColor(colors[idx % len(colors)])
            h.SetLineWidth(3)
            h.SetStats(0)
            hists[var][fname_chunk] = h

        # Initialize 2D histograms
        for var_2d, config in var_2d_configs.items():
            h2_name = f"h2_{var_2d}_{fname_chunk}"
            h2 = ROOT.TH2D(
                h2_name, config["title"], 
                config["xbins"], config["xmin"], config["xmax"],
                config["ybins"], config["ymin"], config["ymax"]
            )
            h2.SetStats(0)
            hists_2d[var_2d][fname_chunk] = h2

        f_in = ROOT.TFile(fname)
        t = f_in.Get("neuttree")

        counter = 0

        for event in t:
            nvect = event.vectorbranch
            nvect_class = nvect_reader(nvect) 
            
            if nvect_class.nocasc_pi and nvect_class.nocasc_piprocessing is not None:
                counter += 1
                pi_data = nvect_class.nocasc_piprocessing
                
                # Fill 1D histograms
                for var in var_configs:
                    val = getattr(pi_data, var, None)
                    if val is not None:
                        hists[var][fname_chunk].Fill(val)

                # Fill 2D histograms
                for var_2d, config in var_2d_configs.items():
                    x_val = getattr(pi_data, config["x_var"], None)
                    y_val = getattr(pi_data, config["y_var"], None)
                    if x_val is not None and y_val is not None:
                        hists_2d[var_2d][fname_chunk].Fill(x_val, y_val)

        f_in.Close()

        if counter == 0:
            print(f"Warning: No valid nocasc_pi events processed for {fname}")

    canvases = {}
    out_file.cd()

    # --- DRAW & SAVE 1D CANVASES ---
    for var in var_configs:
        c = ROOT.TCanvas(f"c_{var}", f"Comparison of {var}", 800, 600)
        c.SetGridy()
        
        legend = ROOT.TLegend(0.68, 0.70, 0.88, 0.88)
        legend.SetBorderSize(1)
        legend.SetFillColor(0)

        max_y = max([h.GetMaximum() for h in hists[var].values()] or [0])
                
        is_first = True
        for fname_chunk, h in hists[var].items():
            if is_first:
                h.SetMaximum(max_y * 1.2)
                h.Draw("HIST")
                is_first = False
            else:
                h.Draw("HIST SAME")
            
            legend.AddEntry(h, fname_chunk, "l")
            
        legend.Draw()
        c.Update()
        c.Write()
        canvases[var] = c

    # --- DRAW & SAVE 2D CANVASES ---
    num_files = len(filename_list)
    for var_2d, config in var_2d_configs.items():
        c2 = ROOT.TCanvas(f"c_2d_{var_2d}", f"2D Comparison of {var_2d}", 600 * num_files, 500)
        
        if num_files > 1:
            c2.Divide(num_files, 1)

        for pad_idx, (fname_chunk, h2) in enumerate(hists_2d[var_2d].items()):
            if num_files > 1:
                c2.cd(pad_idx + 1)
            else:
                c2.cd()

            ROOT.gPad.SetRightMargin(0.15)  # Leave room for z-axis color bar
            ROOT.gPad.SetGridx()
            ROOT.gPad.SetGridy()

            # Set title to include the filename chunk (e.g., LFG or SF)
            raw_title = config["title"].split(";")
            h2.SetTitle(f"{raw_title[0]} ({fname_chunk});{raw_title[1]};{raw_title[2]}")
            h2.Draw("COLZ")

        c2.Update()
        c2.Write()
        canvases[f"2d_{var_2d}"] = c2

    out_file.Close()

    return canvases


def plot_deex_vs_excitation_energy(filename, filename_list=None, use_log_y=False):
    if filename_list is None:
        filename_list = [filename]
    else:
        filename_list = list(filename_list)
        if filename not in filename_list:
            filename_list.append(filename)

    filename_chunk = filename.split(".")[0].split("out_")[-1]
    out_file = ROOT.TFile(f"deex_counts_vs_Eexc_{filename_chunk}.root", "RECREATE")

    deex_pdg_map = {
        22: "#gamma",
        2212: "p",
        1000010020: "d",
        1000020040: "#alpha"
    }

    colors = [
        ROOT.kBlack, ROOT.kRed + 1, ROOT.kBlue + 1, 
        ROOT.kGreen + 2, ROOT.kMagenta + 1, ROOT.kOrange + 7
    ]
    
    # Distinct line styles: solid, dashed, dash-dot, long-dash, dotted, etc.
    line_styles = [1, 2, 7, 9, 3, 10]

    nbins, e_min, e_max = 100, 0.0, 100.0  # Energy range in MeV

    hists = {}
    for idx, (pdg, label) in enumerate(deex_pdg_map.items()):
        h = ROOT.TH1D(
            f"h_deex_count_{pdg}_{filename_chunk}",
            f";E_{{exc}} (MeV);Normalized Yield (A.U.)",
            nbins, e_min, e_max
        )
        color = colors[idx % len(colors)]
        h.SetLineColor(color)
        h.SetLineStyle(line_styles[idx % len(line_styles)])
        h.SetLineWidth(3)
        h.SetFillColorAlpha(color, 0.10)  # Subtle translucent fill
        h.SetStats(0)
        hists[pdg] = h

    # Process events
    for fname in filename_list:
        f_in = ROOT.TFile(fname)
        t = f_in.Get("neuttree")

        for event in t:
            nvect = event.vectorbranch
            nvect_class = nvect_reader(nvect)

            if nvect_class.eventType == EventType.MF and (nvect_class.intChannel == intChannel_CCQE.qeDeEX  or  nvect_class.intChannel == intChannel_CCQE.noCascadeFSIPhoton):
                E_exc = nvect_class.excitation_E_CCQE()
                nvect_class.Print()
                if E_exc is None:
                    continue

                particles, _ = nvect_class.particles()

                for particle in particles:
                    if particle in hists:
                        hists[particle].Fill(E_exc)

        f_in.Close()

    # Normalization
    ##for h in hists.values():
    #    integral = h.Integral()
    #    if integral > 0:
    #        h.Scale(1.0 / integral)

    # Setup Canvas
    canvas = ROOT.TCanvas(f"c_deex_counts_{filename_chunk}", "Deexcitation Counts vs Excitation Energy", 800, 600)
    canvas.SetGridx()
    canvas.SetGridy()
    canvas.SetLeftMargin(0.13)
    canvas.SetBottomMargin(0.12)

    if use_log_y:
        canvas.SetLogy()

    # Clean transparent legend
    legend = ROOT.TLegend(0.72, 0.58, 0.88, 0.88)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.SetTextSize(0.04)

    max_y = max([h.GetMaximum() for h in hists.values()] or [1.0])

    for idx, (pdg, h) in enumerate(hists.items()):
        # Set zoom display window to 0-30 MeV if desired
        h.GetXaxis().SetRangeUser(0, 30)
        
        if idx == 0:
            h.SetMaximum(max_y * (3.0 if use_log_y else 1.25))
            if use_log_y:
                h.SetMinimum(0.002)
            h.Draw("HIST")
        else:
            h.Draw("HIST SAME")

        legend.AddEntry(h, deex_pdg_map[pdg], "lf")

    legend.Draw()
    canvas.Update()

    out_file.cd()
    canvas.Write()
    out_file.Close()

    return canvas