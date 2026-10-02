from incl_analysis_functions import *
from default_processing import *
from kdar_processing import *
from paper_plots import *
from pion_processing import *
import numpy as np
from array import array
import sys 
import argparse
import ROOT

def nocasc_pi_ana(filename, name, filename_list=None, labels=None, legend_pos=(0.68, 0.70, 0.88, 0.88), legend_configs=None):
    """
    Parameters:
      - filename_list: List of ROOT files to analyze.
      - labels: Custom list of string labels for legend entries matching `filename_list`.
      - legend_pos: Default NDC coordinates (x1, y1, x2, y2) for legends.
      - legend_configs: Dict mapping variable names to custom NDC tuples, e.g.,
                        {"e_miss": (0.15, 0.70, 0.35, 0.88), "sigma_Enu": (0.20, 0.20, 0.40, 0.40)}
    """
    labels =  ["LFG (R/S)", "SF (R/S)", "LFG (DCC)", "SF (DCC)" ] 
    if filename_list is None:
        labels = []
        filename_list = [filename]
    else:
        filename_list = list(filename_list)
        if filename not in filename_list:
            filename_list.append(filename)

    # Handle custom or fallback sample labels
    if labels is None:
        labels = [f.split(".")[0].split("out_")[-1] for f in filename_list]
    elif len(labels) < len(filename_list):
        labels += [f.split(".")[0].split("out_")[-1] for f in filename_list[len(labels):]]

    if legend_configs is None:
        legend_configs = {}

    def build_legend(var_key):
        x1, y1, x2, y2 = legend_configs.get(var_key, legend_pos)
        leg = ROOT.TLegend(x1, y1, x2, y2)
        leg.SetBorderSize(1)
        leg.SetFillColor(0)
        return leg

    out_file = ROOT.TFile(f"nocasc_pion_ana_{name}.root", "RECREATE")

    colors = [
        ROOT.kBlack, ROOT.kRed + 1, ROOT.kBlue + 1, 
        ROOT.kGreen + 2, ROOT.kMagenta + 1, ROOT.kOrange + 7, ROOT.kCyan + 2
    ]

    # --- 1D HISTOGRAM CONFIGURATIONS ---
    var_configs = {
        "pion_count": {"title": "Pion Multiplicity;Number of Pions;Events", "bins": 8, "min": 0, "max": 8},
        "pion_momenta": {"title": "Pion Momenta;p_{#pi} (MeV/c);Entries", "bins": 50, "min": 0.0, "max": 1500},
        "initNuc": {"title": "Initial Nucleon Momenta;p (MeV/c);Entries", "bins": 50, "min": 0.0, "max": 1500},
        "pion_angles": {"title": "Pion Scattering Angle;#theta_{#pi} (rad);Entries", "bins": 50, "min": 0.0, "max": 3.2},
        "nucleon_outgoing_momentum": {"title": "Outgoing Nucleon Momentum;p_{N} (MeV/c);Events", "bins": 50, "min": 0.0, "max": 1500},
        "lepton_momentum": {"title": "Outgoing Lepton Momentum;p_{l} (MeV/c);Events", "bins": 10, "min": 0.0, "max": 1000},
        "p_miss": {"title": "Missing Momentum;p_{miss} (MeV/c);Events", "bins": 50, "min": 0.0, "max": 1500},
        "e_miss": {"title": "Missing Energy;E_{miss} (MeV);Area Normalized", "bins": 550, "min": -100, "max": 100},
        "cos_theta_l": {"title": "Lepton Cosine Angle;cos(#theta_{l});Events", "bins": 50, "min": -1.0, "max": 1.0},
        "theta_l_pi": {"title": "Lepton-Pion Opening Angle;#theta_{l#pi} (rad);Events", "bins": 50, "min": 0.0, "max": 3.2},
        "theta_N_pi": {"title": "Nucleon-Pion Opening Angle;#theta_{N#pi} (rad);Events", "bins": 50, "min": 0.0, "max": 3.2},
        "Q2": {"title": "Four-Momentum Transfer;Q^{2} (MeV^{2});Events", "bins": 50, "min": 0.0, "max": 1.0e6},
        "W": {"title": "Hadronic Invariant Mass;W (MeV);Events", "bins": 50, "min": 900.0, "max": 2000.0},
        "nu": {"title": "Energy Transfer;#nu (MeV);Events", "bins": 50, "min": 0.0, "max": 1500.0},
        "delta_pT": {"title": "Transverse Momentum Imbalance;#delta p_{T} (MeV/c);Events", "bins": 50, "min": 0.0, "max": 1000.0},
        "delta_alphaT": {"title": "Transverse Imbalance Angle;#delta#alpha_{T} (rad);Events", "bins": 50, "min": 0.0, "max": 3.2},
        "delta_phiT": {"title": "Transverse Azimuthal Imbalance;#delta#phi_{T} (rad);Events", "bins": 50, "min": 0.0, "max": 3.2},
        "cos_theta_adler": {"title": "Adler Polar Angle;#theta_{A} (rad);Events", "bins": 50, "min": -1.0, "max": 1.0},
        "phi_adler": {"title": "Adler Azimuthal Angle;#phi_{A} (rad);Events", "bins": 50, "min": -3.2, "max": 3.2}
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
        },
        "phi_vs_cos_theta_adler": {
            "x_var": "cos_theta_adler", "y_var": "phi_adler",
            "title": "Adler Angles: #phi_{A} vs #theta_{A};#theta_{A} (rad);#phi_{A} (rad)",
            "xbins": 50, "xmin": 0.0, "xmax": 3.2,
            "ybins": 50, "ymin": 0.0, "ymax": 6.3
        }
    }

    hists = {var: {} for var in var_configs}
    hists_2d = {var_2d: {} for var_2d in var_2d_configs}
    hists_sigma = {}        
    hists_sigma_W = {}      
    hists_sigma_Q2 = {}     
    hists_sigma_thetaA = {} 
    hists_sigma_phiA = {}   
    canvases = {}

    for idx, fname in enumerate(filename_list):
        sample_label = labels[idx]
        
        for var, config in var_configs.items():
            h_name = f"h_{var}_{idx}"
            h = ROOT.TH1D(h_name, config["title"], config["bins"], config["min"], config["max"])
            h.SetLineColor(colors[idx % len(colors)])
            h.SetLineWidth(3)
            h.SetStats(0)
            hists[var][sample_label] = h

        for var_2d, config in var_2d_configs.items():
            h2_name = f"h2_{var_2d}_{idx}"
            h2 = ROOT.TH2D(
                h2_name, config["title"], 
                config["xbins"], config["xmin"], config["xmax"],
                config["ybins"], config["ymin"], config["ymax"]
            )
            h2.SetStats(0)
            hists_2d[var_2d][sample_label] = h2

        f_in = ROOT.TFile.Open(fname)
        t = f_in.Get("neuttree")

        flux_hist = f_in.Get("flux_numu")
        is_monoenergetic = True
        total_flux = 0.0
        yield_hist = None

        if flux_hist and isinstance(flux_hist, ROOT.TH1):
            total_flux = flux_hist.Integral()
            nbins = flux_hist.GetNbinsX()
            filled_flux_bins = sum(1 for i in range(1, nbins + 1) if flux_hist.GetBinContent(i) > 0)
            
            if filled_flux_bins > 1:
                is_monoenergetic = False
                x_edges = [flux_hist.GetXaxis().GetBinLowEdge(i) for i in range(1, nbins + 2)]
                bin_array = np.array(x_edges, dtype=np.float64)
                yield_hist = ROOT.TH1D(f"yield_{idx}", "Event Yield", nbins, bin_array)

        for event in t:
            nvect = event.vectorbranch
            pi_data = PionProcessing(nvect) 

            if pi_data is not None:
                if yield_hist is not None and hasattr(pi_data, 'E_nu') and pi_data.E_nu is not None:
                    yield_hist.Fill(pi_data.E_nu / 1000.0)
                
                for var in var_configs:
                    val = getattr(pi_data, var, None)
                    if val is not None:
                        hists[var][sample_label].Fill(val)

                for var_2d, config in var_2d_configs.items():
                    x_val = getattr(pi_data, config["x_var"], None)
                    y_val = getattr(pi_data, config["y_var"], None)
                    if x_val is not None and y_val is not None:
                        hists_2d[var_2d][sample_label].Fill(x_val, y_val)

        if not is_monoenergetic:
            n_targets = 12.0  
            
            sigma_hist = yield_hist.Clone(f"sigma_{idx}")
            sigma_hist.SetTitle("Total Cross Section;E_{#nu} (GeV);#sigma (cm^{2}/nucleon)")
            sigma_hist.SetLineColor(colors[idx % len(colors)])
            sigma_hist.GetXaxis().SetRangeUser(0.0, 5.0)
            sigma_hist.SetLineWidth(3)
            sigma_hist.SetStats(0)
            sigma_hist.Divide(flux_hist)
            sigma_hist.Scale(1.0 / n_targets)
            sigma_hist.SetDirectory(0)
            hists_sigma[sample_label] = sigma_hist

            sigma_W = hists["W"][sample_label].Clone(f"sigma_W_{idx}")
            sigma_W.SetTitle("Differential Cross Section vs W;W (MeV);d#sigma/dW (cm^{2}/MeV/nucleon)")
            sigma_W.Scale(1.0 / (total_flux * n_targets))
            sigma_W.Scale(1.0, "width")
            sigma_W.SetDirectory(0)
            hists_sigma_W[sample_label] = sigma_W

            sigma_Q2 = hists["Q2"][sample_label].Clone(f"sigma_Q2_{idx}")
            sigma_Q2.SetTitle("Differential Cross Section vs Q^{2};Q^{2} (MeV^{2});d#sigma/dQ^{2} (cm^{2}/MeV^{2}/nucleon)")
            sigma_Q2.Scale(1.0 / (total_flux * n_targets))
            sigma_Q2.Scale(1.0, "width")
            sigma_Q2.SetDirectory(0)
            hists_sigma_Q2[sample_label] = sigma_Q2

            sigma_thetaA = hists["cos_theta_adler"][sample_label].Clone(f"sigma_thetaA_{idx}")
            sigma_thetaA.SetTitle("Differential Cross Section vs #theta_{A};#theta_{A} (rad);d#sigma/d#theta_{A} (cm^{2}/rad/nucleon)")
            sigma_thetaA.Scale(1.0 / (total_flux * n_targets))
            sigma_thetaA.Scale(1.0, "width")
            sigma_thetaA.SetDirectory(0)
            hists_sigma_thetaA[sample_label] = sigma_thetaA

            sigma_phiA = hists["phi_adler"][sample_label].Clone(f"sigma_phiA_{idx}")
            sigma_phiA.SetTitle("Differential Cross Section vs #phi_{A};#phi_{A} (rad);d#sigma/d#phi_{A} (cm^{2}/rad/nucleon)")
            sigma_phiA.Scale(1.0 / (total_flux * n_targets))
            sigma_phiA.Scale(1.0, "width")
            sigma_phiA.SetDirectory(0)
            hists_sigma_phiA[sample_label] = sigma_phiA

        f_in.Close()

    out_file.cd()

    def draw_sigma_canvas(c_name, c_title, hist_dict, var_key):
        c = ROOT.TCanvas(c_name, c_title, 800, 600)
        c.SetGridy()
        legend = build_legend(var_key)

        max_y = max([h.GetMaximum() for h in hist_dict.values()] or [0])
        is_first = True
        for sample_label, h in hist_dict.items():
            if is_first:
                h.SetMaximum(max_y * 1.2)
                h.Draw("HIST")
                is_first = False
            else:
                h.Draw("HIST SAME")
            legend.AddEntry(h, sample_label, "l")
            
        legend.Draw()
        c.Update()
        c.Write()
        return c

    # --- DRAW & SAVE CROSS SECTION CANVASES ---
    if hists_sigma:
        canvases["sigma_Enu"] = draw_sigma_canvas("c_sigma", "Total Cross Section", hists_sigma, "sigma_Enu")
    if hists_sigma_W:
        canvases["sigma_W"] = draw_sigma_canvas("c_sigma_W", "Differential Cross Section dSigma/dW", hists_sigma_W, "sigma_W")
    if hists_sigma_Q2:
        canvases["sigma_Q2"] = draw_sigma_canvas("c_sigma_Q2", "Differential Cross Section dSigma/dQ2", hists_sigma_Q2, "sigma_Q2")
    if hists_sigma_thetaA:
        canvases["sigma_thetaA"] = draw_sigma_canvas("c_sigma_thetaA", "Differential Cross Section dSigma/dTheta_A", hists_sigma_thetaA, "sigma_thetaA")
    if hists_sigma_phiA:
        canvases["sigma_phiA"] = draw_sigma_canvas("c_sigma_phiA", "Differential Cross Section dSigma/dPhi_A", hists_sigma_phiA, "sigma_phiA")

    # --- DRAW & SAVE 1D CANVASES ---
    for var in var_configs:
        c = ROOT.TCanvas(f"c_{var}", f"Comparison of {var}", 800, 600)
        c.SetGridy()
        legend = build_legend(var)

        max_y = max([h.GetMaximum() for h in hists[var].values()] or [0])
        is_first = True
        
        for sample_label, h in hists[var].items():
            if is_first:
                h.SetMaximum(max_y * 1.2)
                h.Draw("HIST")
                is_first = False
            else:
                h.Draw("HIST SAME")
            
            legend.AddEntry(h, sample_label, "l")
            
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

        for pad_idx, (sample_label, h2) in enumerate(hists_2d[var_2d].items()):
            if num_files > 1:
                c2.cd(pad_idx + 1)
            else:
                c2.cd()

            ROOT.gPad.SetRightMargin(0.15)  
            ROOT.gPad.SetGridx()
            ROOT.gPad.SetGridy()

            raw_title = config["title"].split(";")
            h2.SetTitle(f"{raw_title[0]} ({sample_label});{raw_title[1]};{raw_title[2]}")
            h2.Draw("COLZ")

        c2.Update()
        c2.Write()
        canvases[f"2d_{var_2d}"] = c2

    out_file.Close()
    return canvases