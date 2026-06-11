import ROOT

# 1. Configuration
file_path = "kdar_file_comparisons_all_channels_kdar.root"
output_path = "chi2_results_0_to_70.root"
x_min, x_max = 0.0, 70.0  # Set your range here

f_in = ROOT.TFile.Open(file_path)
f_out = ROOT.TFile(output_path, "RECREATE")

for key in f_in.GetListOfKeys():
    if key.GetClassName() == "TCanvas":
        canvas = key.ReadObj()
        canvas_name = canvas.GetName()
        
        if canvas_name == "c3":
            data_hist = None
            mc_hists = []
            
            for prim in canvas.GetListOfPrimitives():
                name = prim.GetName()
                if name.endswith("_data"):
                    data_hist = prim
                elif name == "hs3":
                    hist_list = prim.GetHists()
                    for h in hist_list:
                        mc_hists.append(h)
            
            if not data_hist or not mc_hists:
                continue
            
            mg = ROOT.TMultiGraph()
            mg.SetName(f"mg_chi2_range_{canvas_name}")
            mg.SetTitle(f"Chi2 N-1 ({x_min}-{x_max} MeV);Missing Energy [MeV];#chi^{{2}}_{{N-1}}")

            for h_mc in mc_hists:
                gr = ROOT.TGraph()
                gr.SetName(f"gr_chi2_{h_mc.GetName()}_range")
                
                # --- Step 1: Calculate Total Chi2 only within [0, 70] ---
                full_chi2_vals = []
                valid_bin_indices = [] # Keep track of which bins fall in range
                
                for i in range(1, data_hist.GetNbinsX() + 1):
                    x_val = data_hist.GetBinCenter(i)
                    
                    # Range check
                    if x_val < x_min or x_val > x_max:
                        continue
                        
                    d, m = data_hist.GetBinContent(i), h_mc.GetBinContent(i)
                    err_sq = (data_hist.GetBinError(i)**2) + (h_mc.GetBinError(i)**2)
                    
                    if err_sq > 0:
                        chi2_val = ((d - m)**2) / err_sq
                        full_chi2_vals.append(chi2_val)
                        valid_bin_indices.append(i)
                
                total_chi2 = sum(full_chi2_vals)
                N = len(full_chi2_vals)
                
                if N <= 1: 
                    continue

                # --- Step 2: Calculate Chi2_{N-1} for bins within [0, 70] ---
                point_idx = 0
                for i in valid_bin_indices:
                    x_val = data_hist.GetBinCenter(i)
                    d, m = data_hist.GetBinContent(i), h_mc.GetBinContent(i)
                    err_sq = (data_hist.GetBinError(i)**2) + (h_mc.GetBinError(i)**2)
                    
                    bin_i_contrib = ((d - m)**2) / err_sq
                    chi2_n_minus_1 = (total_chi2 - bin_i_contrib) / (N - 1)
                    
                    gr.SetPoint(point_idx, x_val, chi2_n_minus_1)
                    point_idx += 1
                
                # --- Styling and Saving ---
                gr.SetLineColor(h_mc.GetLineColor())
                gr.SetMarkerColor(h_mc.GetLineColor())
                gr.SetLineWidth(2)
                
                f_out.cd()
                gr.Write() 
                mg.Add(gr, "L")

            mg.Write()

print(f"Analysis complete. Range [{x_min}, {x_max}] saved to {output_path}")
f_out.Close()
f_in.Close()