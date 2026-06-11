import ROOT
import sys
import array

def main():
    # =========================================================
    # 1. Get the first histogram
    # =========================================================
    f1 = ROOT.TFile.Open("Nuiscomp_T2K_IOPfile_INCL_ABLA_newSF_NuWroSRC_100k.root", "READ")
    if not f1 or f1.IsZombie():
        print("Error: Could not open file1.root")
        sys.exit(1)
        
    hist1_orig = f1.Get("MCStudy_CCQE_TProt")
    if not hist1_orig:
        print("Error: MCStudy_CCQE_TProt not found!")
        sys.exit(1)

    # =========================================================
    # 2. Helper function to scale the X-axis 
    # =========================================================
    def scale_x_axis(hist, scale_factor, new_name):
        xaxis = hist.GetXaxis()
        nbins = xaxis.GetNbins()
        
        # Check if the histogram uses variable bin widths
        xbins = xaxis.GetXbins()
        if xbins.GetSize() > 0:
            new_edges = [xbins[i] * scale_factor for i in range(xbins.GetSize())]
            new_edges_arr = array.array('d', new_edges)
            new_hist = ROOT.TH1D(new_name, hist.GetTitle(), nbins, new_edges_arr)
        else:
            # Uniform bin widths
            xmin = xaxis.GetXmin() * scale_factor
            xmax = xaxis.GetXmax() * scale_factor
            new_hist = ROOT.TH1D(new_name, hist.GetTitle(), nbins, xmin, xmax)
            
        # Copy bin contents and errors (including underflow/overflow bins 0 and nbins+1)
        for i in range(0, nbins + 2):
            new_hist.SetBinContent(i, hist.GetBinContent(i))
            new_hist.SetBinError(i, hist.GetBinError(i))
            
        new_hist.SetDirectory(0)
        return new_hist

    # Scale the X-axis of the MCStudy histogram by 1000
    hist1 = scale_x_axis(hist1_orig, 1000.0, "MCStudy_CCQE_TProt_scaledX")
    hist1.SetLineColor(ROOT.kBlue)
    hist1.SetLineWidth(2)

    # =========================================================
    # 3. Open the second file and extract the canvases
    # =========================================================
    f2 = ROOT.TFile.Open("incl_analysis_IOPfile_INCL_ABLA_newSF_NuWroSRC_100k.root", "READ")
    if not f2 or f2.IsZombie():
        print("Error: Could not open the second file")
        sys.exit(1)

    c1 = f2.Get("c1")
    c6 = f2.Get("c6")

    if not c1 or not c6:
        print("Error: Could not find canvas c1 or c6!")
        sys.exit(1)

    # =========================================================
    # 4. Robust helper function to extract totals
    # =========================================================
    def get_total_from_canvas(canvas, name_suffix):
        hists_found = []

        def search_primitives(primitives):
            for obj in primitives:
                if obj.InheritsFrom(ROOT.THStack.Class()):
                    stack_hists = obj.GetStack()
                    if stack_hists and stack_hists.GetSize() > 0:
                        return stack_hists.Last()
                elif obj.InheritsFrom(ROOT.TPad.Class()):
                    result = search_primitives(obj.GetListOfPrimitives())
                    if result: 
                        return result
                elif obj.InheritsFrom(ROOT.TH1.Class()):
                    hists_found.append(obj)
            return None

        found_total = search_primitives(canvas.GetListOfPrimitives())

        if found_total:
            total_hist = found_total.Clone(f"summed_{name_suffix}")
            total_hist.SetDirectory(0)
            return total_hist
            
        if hists_found:
            total_hist = hists_found[0].Clone(f"summed_{name_suffix}")
            total_hist.SetDirectory(0)
            for i in range(1, len(hists_found)):
                total_hist.Add(hists_found[i])
            return total_hist

        return None

    # =========================================================
    # 5. Extract and process second file
    # =========================================================
    hist_c1 = get_total_from_canvas(c1, "c1")
    hist_c6 = get_total_from_canvas(c6, "c6")

    if not hist_c1 or not hist_c6:
        print("Error: Failed to extract histograms from the canvases.")
        sys.exit(1)

    hist_total = hist_c1.Clone("hist_total")
    hist_total.SetDirectory(0) 
    hist_total.Add(hist_c6)
    hist_total.SetLineColor(ROOT.kRed)
    hist_total.SetLineWidth(2)

    f1.Close()
    f2.Close()

    # =========================================================
    # 6. Normalize Histograms
    # =========================================================
    # Normalize hist_total to an area of 1
    if hist_total.Integral() > 0:
        hist_total.Scale(1.0 / hist_total.Integral())

    # Normalize hist1 to an area of 1
    if hist1.Integral() > 0:
        hist1.Scale(1.0 / hist1.Integral())

    # =========================================================
    # 7. Plot everything on a new canvas
    # =========================================================
    c_final = ROOT.TCanvas("c_final", "Combined Histograms", 800, 600)

    max1 = hist1.GetMaximum()
    max_total = hist_total.GetMaximum()

    # Multiply the highest peak by 1.1 to give some breathing room at the top
    if max_total > max1:
        hist_total.SetMaximum(max_total * 1.1)
        hist_total.Draw("HIST")
        hist1.Draw("HIST SAME")
    else:
        hist1.SetMaximum(max1 * 1.1)
        hist1.Draw("HIST")
        hist_total.Draw("HIST SAME")

    leg = ROOT.TLegend(0.55, 0.75, 0.85, 0.85)
    leg.SetBorderSize(0)
    leg.AddEntry(hist1, "MCStudy_CCQE_TProt (X-axis x1000)", "l")
    leg.AddEntry(hist_total, "Total Stack (c1 + c6)", "l")
    leg.Draw()

    c_final.Update()

    # Save outputs
    c_final.SaveAs("Combined_Histograms.root")
    c_final.SaveAs("Combined_Histograms.pdf")
    input("Press Enter to exit...")

if __name__ == "__main__":
    main()