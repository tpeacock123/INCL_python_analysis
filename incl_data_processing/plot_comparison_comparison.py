import ROOT
import os

# --- Configuration ---
#input_files = [
# {
# "path": "NUISCOMP_study_Nuiscomp_T2K_IOPfile_DEFAULT_newSF_NEUTSRC_100k.root",
# "label": "MC (Default Cascade)"
# },
# {
# "path": "NUISCOMP_study_Nuiscomp_T2K_IOPfile_INCL_ABLA_newSF_NEUTSRC_100k.root",
# "label": "MC (INCL++ with ABLA)"
# },
# {
# "path": "NUISCOMP_study_Nuiscomp_T2K_IOPfile_INCL_NoABLA_newSF_NEUTSRC_100k.root",
# "label": "MC (INCL++ standalone)"
# }
#]
#input_files = [
# {
# "path": "NUISCOMP_study_Nuiscomp_MicroBOONe_IOPfile_DEFAULT_newSF_NEUTSRC_100k.root",
# "label": "MC (Default Cascade)"
# },
# {
# "path": "NUISCOMP_study_Nuiscomp_MicroBOONe_IOPfile_INCL_ABLA_newSF_NEUTSRC_100k.root",
# "label": "MC (INCL++ with ABLA)"
# },
#] 
input_files = [
    {
        "path": "Minerva_output_plots_DEFAULT.root",
        "label": "MC (Default Cascade)"
    },
    {
        "path": "Minerva_output_plots_INCL_ABLA.root",
        "label": "MC (INCL++ with ABLA)"
    },
    {
        "path": "Minerva_output_plots_INCL_noABLA.root",
        "label": "MC (INCL++ without ABLA)"
    }
]

output_filename = "MINERVA_IOP_total_data_comparisons_.root"

color_palette = [
    ROOT.kRed, ROOT.kBlue, ROOT.kGreen + 2, ROOT.kMagenta, 
    ROOT.kOrange + 7, ROOT.kCyan + 2, ROOT.kYellow + 2, ROOT.kViolet + 1
]
line_style_palette = [1, 2, 3, 9, 10] 

ROOT.gROOT.SetBatch(True)

# --- Global Style Adjustments ---
# Set the grid color to a light grey
ROOT.gStyle.SetGridColor(ROOT.kGray)
# Optional: Set grid style to dotted (3) or dashed (2) for a lighter feel. 1 is solid.
ROOT.gStyle.SetGridStyle(3) 

def main():
    if not input_files:
        print("Error: No input files specified.")
        return

    # 1. Open all input files
    opened_files = []
    for file_info in input_files:
        f = ROOT.TFile.Open(file_info["path"], "READ")
        if not f or f.IsZombie():
            print(f"Error: Could not open {file_info['path']}")
            return
        opened_files.append({"file": f, "label": file_info["label"]})

    # 2. Create the output ROOT file
    out_file = ROOT.TFile(output_filename, "RECREATE")

    base_file = opened_files[0]["file"]

    for key in base_file.GetListOfKeys():
        plot_name = key.GetName()
        
        data_hist = None
        mc_hists = [] 
        skip_plot = False

        # Extract Data/MC from canvases
        for i, file_dict in enumerate(opened_files):
            f = file_dict["file"]
            label = file_dict["label"]
            
            canvas = f.Get(plot_name)
            if not canvas or not canvas.InheritsFrom("TCanvas"):
                skip_plot = True
                continue

            mc_found = False
            for prim in canvas.GetListOfPrimitives():
                name = prim.GetName().lower()
                
                if i == 0 and "data" in name and not data_hist:
                    data_hist = prim.Clone(f"{plot_name}_data")
                    data_hist.SetDirectory(0)
                
                elif "mc" in name and not mc_found:
                    mc_hist = prim.Clone(f"{plot_name}_mc_{i}")
                    mc_hist.SetDirectory(0)
                    mc_hists.append((mc_hist, label))
                    mc_found = True
            
            if not mc_found:
                skip_plot = True

        if skip_plot or not data_hist or len(mc_hists) != len(opened_files):
            continue

        # 3. Create a canvas for this specific plot
        c_temp = ROOT.TCanvas(plot_name, f"Comparison: {plot_name}", 800, 600)
        
        # Enable grid on the canvas
        c_temp.SetGrid() 
        
        data_hist.SetMarkerStyle(20)
        data_hist.SetMarkerColor(ROOT.kBlack)
        data_hist.SetLineColor(ROOT.kBlack)

        # Force Y-axis to start at 0
        data_hist.SetMinimum(0)

        # Set Y-axis maximum to 30% above the highest data point
        max_data_val = data_hist.GetMaximum()
        data_hist.SetMaximum(max_data_val * 1.3) 
        
        # --- Halve the size of axis titles and labels ---
        # X-axis adjustments
        x_axis = data_hist.GetXaxis()
        x_axis.SetTitleSize(x_axis.GetTitleSize() * 0.5)
        x_axis.SetLabelSize(x_axis.GetLabelSize() * 0.5)
        
        # Y-axis adjustments
        y_axis = data_hist.GetYaxis()
        y_axis.SetTitleSize(y_axis.GetTitleSize() * 0.5)
        y_axis.SetLabelSize(y_axis.GetLabelSize() * 0.5)
        # Optionally, you might need to increase the title offset if the smaller title gets too close
        # y_axis.SetTitleOffset(y_axis.GetTitleOffset() * 1.5)
        
        data_hist.Draw("E1 P")

        leg_height = min(0.4, 0.05 * (len(mc_hists) + 1))
        leg = ROOT.TLegend(0.50, 0.88 - leg_height, 0.88, 0.88)
        leg.SetBorderSize(0)
        leg.AddEntry(data_hist, "Data", "lp")

        for i, (mc_hist, label) in enumerate(mc_hists):
            mc_hist.SetLineColor(color_palette[i % len(color_palette)])
            mc_hist.SetLineStyle(line_style_palette[i % len(line_style_palette)])
            mc_hist.SetLineWidth(3)
            mc_hist.Draw("HIST SAME")
            leg.AddEntry(mc_hist, label, "l")

        leg.Draw()

        # Write the canvas to the output file
        out_file.cd()
        c_temp.Write()
        
        # Cleanup memory for this loop
        c_temp.Close()

    # 4. Finalize
    out_file.Close()
    for file_dict in opened_files:
        file_dict["file"].Close()
        
    print(f"Done! Combined canvases saved to {output_filename}")

if __name__ == "__main__":
    main()