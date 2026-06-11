import ROOT

def process_and_save_to_root(file1_path, file2_path, output_name="output_plots.root"):

    f1 = ROOT.TFile.Open(file1_path, "READ")
    f2 = ROOT.TFile.Open(file2_path, "READ") 
    if not f1 or f1.IsZombie() or not f2 or f2.IsZombie():
        print("Error: Could not open one or both ROOT files.")
        return

    hist_groups = {}

    def get_category(name):
        if name.endswith("_MC_SHAPE"): return name[:-9], "MC_SHAPE"
        if name.endswith("_MC"): return name[:-3], "MC"
        if name.endswith("_data_SHAPE"): return name[:-11], "data_SHAPE"
        if name.endswith("_data"): return name[:-5], "data"
        return None, None

    def process_file(f):
        for key in f.GetListOfKeys():
            obj = key.ReadObj()
            
            if not isinstance(obj, ROOT.TH1):
                continue
                
            base_name, cat = get_category(obj.GetName())
            
            if not base_name:
                continue 
                
            if base_name not in hist_groups:
                hist_groups[base_name] = {"MC": None, "MC_SHAPE": None, "data": None, "data_SHAPE": None}
                
            if hist_groups[base_name][cat] is None:
                hist_groups[base_name][cat] = obj.Clone(f"{base_name}_{cat}_final")
                hist_groups[base_name][cat].SetDirectory(0)
            else:
                if cat in ["MC", "MC_SHAPE"]:
                    hist_groups[base_name][cat].Add(obj)


    process_file(f1)
    process_file(f2)
    
    f1.Close()
    f2.Close()

    out_file = ROOT.TFile.Open(output_name, "RECREATE")
    
    for base_name, hists in hist_groups.items():
        
        if hists["MC"] and hists["data"]:
            c_nom = ROOT.TCanvas(f"c_{base_name}_nominal", f"{base_name} Nominal", 800, 600)
            mc = hists["MC"]
            data = hists["data"]
            
            x_label = mc.GetXaxis().GetTitle()
            y_label = mc.GetYaxis().GetTitle()
    
            mc.SetLineColor(ROOT.kRed)
            mc.SetLineWidth(2)
            mc.SetFillColorAlpha(ROOT.kRed, 0.3)
            mc.SetStats(0)
            
            mc.SetTitle(f"{base_name} (Nominal);{x_label};{y_label}")
            
            data.SetMarkerStyle(20) 
            data.SetMarkerColor(ROOT.kBlack)
            data.SetLineColor(ROOT.kBlack)
            data.SetStats(0)
            
            max_y = max(mc.GetMaximum(), data.GetMaximum())
            mc.SetMinimum(0.0)
            mc.SetMaximum(max_y * 1.3)
            
            mc.Draw("HIST")
            data.Draw("PE SAME")
            
            leg = ROOT.TLegend(0.65, 0.75, 0.88, 0.88)
            leg.SetBorderSize(0)
            leg.AddEntry(data, "Data", "pe")
            leg.AddEntry(mc, "MC (Summed)", "f")
            leg.Draw()
            
            c_nom.Write()
            
        if hists["MC_SHAPE"] and hists["data_SHAPE"]:
            c_shape = ROOT.TCanvas(f"c_{base_name}_shape", f"{base_name} Shape", 800, 600)
            mc_s = hists["MC_SHAPE"]
            data_s = hists["data_SHAPE"]
            
            x_label_s = mc_s.GetXaxis().GetTitle()
            y_label_s = mc_s.GetYaxis().GetTitle()
            
            mc_s.SetLineColor(ROOT.kBlue)
            mc_s.SetLineWidth(2)
            mc_s.SetFillColorAlpha(ROOT.kBlue, 0.3)
            mc_s.SetStats(0)
            mc_s.SetTitle(f"{base_name} (Shape);{x_label_s};{y_label_s}")
            
            data_s.SetMarkerStyle(24) 
            data_s.SetMarkerColor(ROOT.kBlack)
            data_s.SetLineColor(ROOT.kBlack)
            data_s.SetStats(0)
            
            max_y_s = max(mc_s.GetMaximum(), data_s.GetMaximum())
            mc_s.SetMinimum(0.0)
            mc_s.SetMaximum(max_y_s * 1.3)
            
            mc_s.Draw("HIST")
            data_s.Draw("PE SAME")
            
            leg2 = ROOT.TLegend(0.65, 0.75, 0.88, 0.88)
            leg2.SetBorderSize(0)
            leg2.AddEntry(data_s, "Data Shape", "pe")
            leg2.AddEntry(mc_s, "MC Shape (Summed)", "f")
            leg2.Draw()
        
            c_shape.Write()

        for cat, hist in hists.items():
            if hist:
                hist.Write()
                
    out_file.Close()
    print(f"Success! All canvases and raw summed histograms saved to {output_name}")

# --- Execution ---
if __name__ == "__main__":
    #file_A = "Nuiscomp_T2K_out_IOPfile_DEFAULT_NoCCQE_100k.root"
    #file_B = "Nuiscomp_T2K_IOPfile_INCL_ABLA_newSF_NEUTSRC_100k.root"
    file_A = "Nuiscomp_MINERVA_out_IOPfile_DEFAULT_NoCCQE_100k.root"
    file_B = "Nuiscomp_MINERVA_out_IOPfile_INCL_ABLA_CCQEonly_100k.root"
    #file_B = "Nuiscomp_T2K_out_IOPfile_DEFAULT_CCQEOnly_100k.root"
    
    process_and_save_to_root(file_A, file_B, "Minerva_output_plots_INCL_ABLA.root")