import ROOT
import sys

def main():
    filename = sys.argv[1]
    filename_chunk = filename.split(".")[0]

    input_filename = filename
    output_filename = "NUISCOMP_study_{}.root".format(filename_chunk)

    # 1. Run in batch mode to prevent canvas pop-ups
    ROOT.gROOT.SetBatch(True)

    # 2. Open the input file
    infile = ROOT.TFile.Open(input_filename, "READ")
    if not infile or infile.IsZombie():
        print(f"Error: Could not open {input_filename}")
        sys.exit(1)

    # 3. Create the output file
    outfile = ROOT.TFile(output_filename, "RECREATE")

    print(f"Scanning {input_filename} for matching Data/MC pairs...")
    
    pairs_plotted = 0

    # 4. Loop over all keys in the file
    for key in infile.GetListOfKeys():
        # Check if it's a histogram
        if not key.ReadObj().InheritsFrom("TH1"):
            continue

        data_name = key.GetName()
        mc_name = None
        plot_identifier = ""

        # --- A. Check for standard _data suffix ---
        if data_name.endswith("_data"):
            base_name = data_name[:-5] 
            mc_name = base_name + "_MC"
            plot_identifier = base_name
        
        # --- B. Check for _data_SHAPE suffix ---
        elif data_name.endswith("_data_SHAPE"):
            base_name = data_name[:-11]
            mc_name = base_name + "_mc_SHAPE"
            
            # Fallback in case the MC file uses uppercase "_MC_SHAPE" instead
            if not infile.GetListOfKeys().Contains(mc_name) and infile.GetListOfKeys().Contains(base_name + "_MC_SHAPE"):
                mc_name = base_name + "_MC_SHAPE"
                
            plot_identifier = base_name + "_SHAPE"
        
        # If it doesn't match either pattern, skip it
        else:
            continue

        # 5. Check if the corresponding MC histogram exists
        if infile.GetListOfKeys().Contains(mc_name):
            h_data = infile.Get(data_name)
            h_mc = infile.Get(mc_name)

            # Ensure the MC object is also a histogram
            if not h_mc.InheritsFrom("TH1"):
                continue

            # 6. Create a Canvas (named uniquely so standard and SHAPE don't overwrite)
            canvas_name = f"canvas_{plot_identifier}"
            c1 = ROOT.TCanvas(canvas_name, f"Data vs MC for {plot_identifier}", 800, 600)

            # 7. Style the histograms
            # Data: Black dots with error bars
            h_data.SetMarkerStyle(20)
            h_data.SetMarkerColor(ROOT.kBlack)
            h_data.SetLineColor(ROOT.kBlack)
            h_data.SetStats(0) # Turn off the stats box
            
            # MC: Red solid line
            h_mc.SetLineColor(ROOT.kRed)
            h_mc.SetLineWidth(2)
            h_mc.SetStats(0)

            # 8. Adjust the Y-axis range to fit both histograms
            max_y = max(h_data.GetMaximum(), h_mc.GetMaximum())
            h_data.SetMaximum(max_y * 1.2) # Add 20% headroom for the legend

            # 9. Draw the histograms
            h_data.Draw("PE")
            h_mc.Draw("HIST SAME")

            # 10. Add a Legend
            leg = ROOT.TLegend(0.7, 0.75, 0.88, 0.88)
            leg.SetBorderSize(0)
            leg.AddEntry(h_data, "Data", "pe")
            leg.AddEntry(h_mc, "MC", "l")
            leg.Draw()

            # 11. Write the canvas to the output file
            outfile.cd()
            c1.Write()
            pairs_plotted += 1

    # Clean up
    infile.Close()
    outfile.Close()

    print(f"Done! Successfully paired, plotted, and saved {pairs_plotted} canvases to {output_filename}.")

if __name__ == "__main__":
    main()