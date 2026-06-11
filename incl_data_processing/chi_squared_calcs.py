import ROOT


file_path = "T2K_IOP_total_data_comparisons_finalised.root"
f = ROOT.TFile.Open(file_path)

for key in f.GetListOfKeys():
    
    if key.GetClassName() == "TCanvas":
        canvas = key.ReadObj()
        canvas_name = canvas.GetName()
        
        data_hist = None
        mc_hists = []
        
        for prim in canvas.GetListOfPrimitives():
            if prim.InheritsFrom("TH1"):
                name = prim.GetName()
                if name.endswith("_data"):
                    data_hist = prim
                elif name.endswith(("_mc_0", "_mc_1", "_mc_2", "_mc2")):
                    mc_hists.append(prim)
        
        print(f"\n{'='*60}")
        print(f"Canvas: {canvas_name}")
        print(f"{'='*60}")
        
        if not data_hist or not mc_hists:
            continue
            
        print(f"{'MC Sample Name':<30} | {'Chi2 / NDOF'}")
        print("-" * 60)
        
        for h_mc in mc_hists:
            name = h_mc.GetName()
            chi2 = 0.0
            ndf = 0
            
            for i in range(1, data_hist.GetNbinsX() + 1):
                d_val = data_hist.GetBinContent(i)
                m_val = h_mc.GetBinContent(i)
                
                d_err = data_hist.GetBinError(i)
                m_err = h_mc.GetBinError(i)
                err_squared = (d_err**2) + (m_err**2)

                if err_squared > 0:
                    chi2 += ((d_val - m_val)**2) / err_squared
                    ndf += 1
            
            print(f"{name:<30} | {chi2:.2f} / {ndf}")

f.Close()