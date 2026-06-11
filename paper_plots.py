import numpy as np
from array import array
import sys 
import argparse
import ROOT
from incl_analysis_functions import *

ROOT.gSystem.Load("libNEUTROOTClass.so")
ROOT.gSystem.Load("libNEUTOutput.so")
ROOT.gSystem.Load("libNEUTReWeight.so")
ROOT.TH1.AddDirectory(False)


def spectral_function2d(filename,additional_filename):
    import ROOT

    ROOT.gStyle.SetOptStat(0)

    # --- First File ---
    filename_chunk = filename.split(".")[0].split("out_")[1]
    f = ROOT.TFile(filename)
    t = f.Get("neuttree")   

    h_miss = ROOT.TH2D("h_miss", "; P_{N} [MeV/c]; E_{RMV} [MeV/c]", 150, 0, 450, 80, 0, 80)

    for event in t:
        nvect = event.vectorbranch
        nvect_class = nvect_reader(nvect)
        eMiss = nvect_class.E_miss
        pMiss = nvect_class.P_miss
        h_miss.Fill(pMiss, eMiss)


    # --- Second File ---
    f2 = ROOT.TFile(additional_filename)
    t2 = f2.Get("neuttree")  

    # IMPORTANT: Ensure the internal ROOT name ("h2_miss") is unique so it doesn't overwrite h_miss
    h2_miss = ROOT.TH2D("h2_miss", "; P_{N} [MeV/c]; E_{RMV} [MeV/c]", 150, 0, 450, 80, 0, 80)

    for event in t2:
        nvect = event.vectorbranch
        nvect_class = nvect_reader(nvect)
        eMiss = nvect_class.E_miss
        pMiss = nvect_class.P_miss
        h2_miss.Fill(pMiss, eMiss)


    # --- Syncing the Color Scales ---
    # Find the global maximum to ensure the color scales map exactly the same way
    global_max = max(h_miss.GetMaximum(), h2_miss.GetMaximum())

    h_miss.SetMaximum(global_max)
    h2_miss.SetMaximum(global_max)

    h_miss.SetMinimum(0.00001)
    h2_miss.SetMinimum(0.00001)


    # --- Styling ---
    ROOT.gStyle.SetPalette(ROOT.kBlueRedYellow)
    ROOT.gStyle.SetNumberContours(99)

    # Create a wider canvas to accommodate two plots
    c1 = ROOT.TCanvas("c1", "c1", 1600, 700) 
    c1.Divide(2, 1) # Divide into 1 row, 2 columns

    # Helper function to apply the same axis formatting
    def format_axes(hist):
        for axis in [hist.GetXaxis(), hist.GetYaxis()]:
            axis.SetLabelColor(ROOT.kBlack)
            axis.SetTitleColor(ROOT.kBlack)
            axis.SetAxisColor(ROOT.kWhite) 

    # --- Draw Pad 1 ---
    pad1 = c1.cd(1)
    pad1.SetFillColor(ROOT.kWhite) 
    pad1.SetFrameFillColor(ROOT.kBlack)
    pad1.SetRightMargin(0.15) # Give space for the Z-axis color palette
    format_axes(h_miss)
    h_miss.Draw() # Use COLZ to draw the 2D color map

    # --- Draw Pad 2 ---
    pad2 = c1.cd(2)
    pad2.SetFillColor(ROOT.kWhite) 
    pad2.SetFrameFillColor(ROOT.kBlack)
    pad2.SetRightMargin(0.15) # Give space for the Z-axis color palette
    format_axes(h2_miss)
    h2_miss.Draw() # Use COLZ to draw the 2D color map

    # --- Update and Save ---
    c1.Update()
    c1.SaveAs("spectral_function_{}.png".format(filename_chunk))

def SRC_plot(filename,additional_filename):
    ROOT.gStyle.SetOptStat(0)

    # Helper function to apply the same axis formatting
    def format_axes(hist):
        for axis in [hist.GetXaxis(), hist.GetYaxis()]:
            axis.SetLabelColor(ROOT.kBlack)
            axis.SetTitleColor(ROOT.kBlack)
            axis.SetAxisColor(ROOT.kWhite)

    # --- First File ---
    filename_chunk = filename.split(".")[0].split("out_")[1] if "out_" in filename else "output"
    f = ROOT.TFile(filename)
    t = f.Get("neuttree")   

    p_mf_1, e_mf_1 = array('d'), array('d')
    p_src_1, e_src_1 = array('d'), array('d')

    for event in t:
        nvect = event.vectorbranch
        nvect_class = nvect_reader(nvect)
        if nvect_class.eventType == EventType.MF:
            p_mf_1.append(nvect_class.P_miss)
            e_mf_1.append(nvect_class.E_miss)
        elif nvect_class.eventType == EventType.SRC:
            p_src_1.append(nvect_class.P_miss)
            e_src_1.append(nvect_class.E_miss)

    print( "NEUT", len(p_src_1)/(len(p_mf_1) + len(p_src_1)))

    g1_mf = ROOT.TGraph(len(p_mf_1), p_mf_1, e_mf_1)
    g1_src = ROOT.TGraph(len(p_src_1), p_src_1, e_src_1)

    # --- Second File ---
    f2 = ROOT.TFile(additional_filename)
    t2 = f2.Get("neuttree")  

    p_mf_2, e_mf_2 = array('d'), array('d')
    p_src_2, e_src_2 = array('d'), array('d')

    for event in t2:
        nvect = event.vectorbranch
        nvect_class = nvect_reader(nvect)
        if nvect_class.eventType == EventType.MF:
            p_mf_2.append(nvect_class.P_miss)
            e_mf_2.append(nvect_class.E_miss)
        elif nvect_class.eventType == EventType.SRC:
            p_src_2.append(nvect_class.P_miss)
            e_src_2.append(nvect_class.E_miss)

    print( "NuWro", len(p_src_2)/(len(p_mf_2) + len(p_src_2)))

    g2_mf = ROOT.TGraph(len(p_mf_2), p_mf_2, e_mf_2)
    g2_src = ROOT.TGraph(len(p_src_2), p_src_2, e_src_2)

    # --- Styling the Scatter Points ---
    # High contrast colors for black background: Cyan and Magenta
    color_mf = ROOT.kAzure-1
    color_src = ROOT.kRed-4

    for g in [g1_mf, g2_mf]:
        g.SetMarkerStyle(20) 
        g.SetMarkerSize(0.4)
        g.SetMarkerColorAlpha(color_mf,0.5) 

    for g in [g1_src, g2_src]:
        g.SetMarkerStyle(20)
        g.SetMarkerSize(0.4)
        g.SetMarkerColorAlpha(color_src,0.5)

    # --- Dummy Histograms for Axes Formatting ---
    h_dummy1 = ROOT.TH2D("hd1", "; P_{N} [MeV/c]; E_{RMV} [MeV/c]", 150, 0, 600, 150, 0, 250)
    h_dummy2 = ROOT.TH2D("hd2", "; P_{N} [MeV/c]; E_{RMV} [MeV/c]", 150, 0, 600, 150, 0, 250)

    # --- Legend Setup (Clearer & Bigger Markers) ---
    # Create dummy graphs just for the legend so the markers appear larger in the box
    leg_dummy_mf = ROOT.TGraph()
    leg_dummy_mf.SetMarkerStyle(20)
    leg_dummy_mf.SetMarkerColor(color_mf)
    leg_dummy_mf.SetMarkerSize(1.5) # Much larger for visibility

    leg_dummy_src = ROOT.TGraph()
    leg_dummy_src.SetMarkerStyle(20)
    leg_dummy_src.SetMarkerColor(color_src)
    leg_dummy_src.SetMarkerSize(1.5)

    def create_legend():
        # Better positioning and sizing
        leg = ROOT.TLegend(0.60, 0.75, 0.88, 0.88)
        leg.AddEntry(leg_dummy_mf, "Mean Field", "p")
        leg.AddEntry(leg_dummy_src, "SRC", "p")
        leg.SetTextColor(ROOT.kWhite)
        # Add a semi-transparent dark background with a white border
        leg.SetFillColorAlpha(ROOT.kBlack, 0.7) 
        leg.SetLineColor(ROOT.kWhite)
        leg.SetBorderSize(1)
        return leg

    # --- Setup Canvas ---
    c1 = ROOT.TCanvas("c1", "c1", 1600, 700) 
    c1.Divide(2, 1)

    # --- Draw Pad 1 ---
    pad1 = c1.cd(1)
    pad1.SetFillColor(ROOT.kWhite) 
    pad1.SetFrameFillColor(ROOT.kBlack)
    format_axes(h_dummy1)
    
    h_dummy1.Draw()           
    g1_mf.Draw("P SAME")      
    g1_src.Draw("P SAME")

    # --- Draw Pad 2 ---
    pad2 = c1.cd(2)
    pad2.SetFillColor(ROOT.kWhite) 
    pad2.SetFrameFillColor(ROOT.kBlack)
    format_axes(h_dummy2)
    
    h_dummy2.Draw()
    g2_mf.Draw("P SAME")
    g2_src.Draw("P SAME")

    leg2 = create_legend()
    leg2.Draw()

    # --- Update and Save ---
    c1.Update()
    c1.SaveAs("spectral_scatter_{}.png".format(filename_chunk))

def transparency(filename, filename_list=[]):
    # Fix 1: Avoid using mutable default arguments
        
    hist_list = []
    filename_list.append(filename)

    # Setup colors AND line styles to match the plot
    # Line styles in ROOT: 1=Solid, 2=Dash, 3=Dot, 7=Long Dash, 9=Dash-Dot, etc.
    styles = [
        (ROOT.kBlack, 3),    # Dotted Black
        (ROOT.kRed, 7),      # Dashed Red
        (ROOT.kBlue, 9),     # Dash-Dot Blue
        (ROOT.kGreen+2, 5),  # Long-Dash Green
        (ROOT.kMagenta, 2),  # Short-Dash Magenta
        (ROOT.kCyan, 1)      # Solid Cyan (fallback)
    ]

    for filename in filename_list:
        filename_chunk = filename.split(".")[0].split("out_")[-1]

        f = ROOT.TFile(filename)
        t = f.Get("neuttree")   

        # Changed X-axis range to 0.0 -> 1.0 (GeV)
        h_trans = ROOT.TH1D("h_trans_{}".format(filename_chunk), "Numerator", 50, 0.0, 1.0)
        h_total = ROOT.TH1D("h_total_{}".format(filename_chunk), "Denominator", 50, 0.0, 1.0)
        h_trans.Sumw2()
        h_total.Sumw2()

        for event in t:
            nvect = event.vectorbranch
            nvect_class = nvect_reader(nvect) 
            if nvect_class.eventType == EventType.MF: 
                # Calculate FSI Proton KE in MeV
                fsiProtonKE_MeV = np.sqrt(nvect_class.fsiProton**2 + 938.27**2) - 938.27
                
                if nvect_class.intChannel == intChannel_CCQE.qeDeEX or nvect_class.intChannel ==  intChannel_CCQE.noCascadeFSI:
                    h_trans.Fill(fsiProtonKE_MeV)
                    h_total.Fill(fsiProtonKE_MeV)
                else:
                    h_total.Fill(fsiProtonKE_MeV)

        h_ratio = h_trans.Clone("h_ratio_{}".format(filename_chunk))
        # Set titles: Title; X-axis; Y-axis
        h_ratio.SetTitle(";T_{p} [GeV];Transparency p C")
        h_ratio.Divide(h_trans, h_total, 1.0, 1.0, "B")
        
        # Detach histogram from file
        h_ratio.SetDirectory(0) 
        
        hist_list.append((filename_chunk, h_ratio))
        f.Close()

    # --- PLOTTING AND SAVING SECTION ---

    # Suppress the default plot title box at the top
    ROOT.gStyle.SetOptTitle(0)

    out_file = ROOT.TFile("paper_transparency_combined_output.root", "RECREATE")

    # Create Canvas and add top/right ticks to match the reference image
    canvas = ROOT.TCanvas("c1", "Transparency Plot", 800, 600)
    canvas.SetTickx(1)
    canvas.SetTicky(1)
    # Increase bottom margin slightly for the legend
    canvas.SetBottomMargin(0.15) 
    canvas.SetLeftMargin(0.12)

    # Create Legend (x1, y1, x2, y2 in normalized coordinates)
    # Placed at the bottom, spanning across
    legend = ROOT.TLegend(0.15, 0.18, 0.85, 0.35) 
    legend.SetNColumns(2) # Split into two columns
    legend.SetBorderSize(1) # Solid border
    legend.SetTextSize(0.04)

    for i, (name, h) in enumerate(hist_list):
        color, style = styles[i % len(styles)]
        
        h.SetLineColor(color)
        h.SetLineStyle(style)
        h.SetLineWidth(2) # Make lines slightly thicker for visibility
        
        h.SetStats(0) 
        h.SetMinimum(0.0) 
        h.SetMaximum(1.0) # Force Y-axis to 1.0

        # Adjust axis label and title sizes
        h.GetXaxis().SetTitleSize(0.06)
        h.GetXaxis().SetTitleOffset(1.0)
        h.GetYaxis().SetTitleSize(0.06)
        h.GetYaxis().SetTitleOffset(0.9)

        # Draw lines ("HIST" or "L") instead of error bars ("E")
        if i == 0:
            h.Draw("HIST") 
        else:
            h.Draw("HIST SAME")
            
        # Add to legend. "l" indicates a line entry.
        legend.AddEntry(h, name, "l")

    legend.Draw()
    out_file.cd() 
    canvas.Write("transparency_canvas") 

    out_file.Close()
import ROOT

def ccqe_combined_plots(filename, filename_list=None):
    if filename_list is None:
        filename_list = []
        
    # Combine into a single list to iterate over
    files_to_process = [filename] + filename_list
    
    # Store histograms separately
    hist_list_proton = []
    hist_list_ex = []

    styles = [
        (ROOT.kRed, 7),      # Dashed Red
        (ROOT.kBlue, 9),     # Dash-Dot Blue
        (ROOT.kGreen+2, 5),  # Long-Dash Green
        (ROOT.kMagenta, 2),  # Short-Dash Magenta
        (ROOT.kCyan, 1)      # Solid Cyan (fallback)
    ]

    # --- 1. Loop over files and events ONCE ---
    for fname in files_to_process:
        filename_chunk = fname.split(".")[0].split("out_")[-1]

        f = ROOT.TFile(fname)
        t = f.Get("neuttree")   

        # Initialize both histograms
        h_leadin = ROOT.TH1D(f"h_leading_{filename_chunk}", "Leading Proton Distribution", 50, 0.0, 1500)
        h_leadin.Sumw2()
        h_leadin.SetTitle(";P_{N} [GeV];Counts")

        h_ex = ROOT.TH1D(f"h_excitation_E_{filename_chunk}", "Excitation Energy Distribution", 40, 0.0, 100)
        h_ex.Sumw2()
        h_ex.SetTitle(";E_{ex} [MeV];Counts")

        # Single event loop
        for event in t:
            nvect = event.vectorbranch
            # Assuming nvect_reader and EventType are defined elsewhere
            nvect_class = nvect_reader(nvect) 
            
            if nvect_class.eventType == EventType.MF: 
                # Fill Excitation Energy
                h_ex.Fill(nvect_class.excitation_E_CCQE())
                
                # Fill Leading Proton (using 'if' instead of 'continue' so we don't break the loop)
                if nvect_class.HMPMom:
                    h_leadin.Fill(nvect_class.HMPMom)

        h_leadin.SetDirectory(0) 
        h_ex.SetDirectory(0) 
        
        hist_list_proton.append((filename_chunk, h_leadin))
        hist_list_ex.append((filename_chunk, h_ex))
        f.Close()

    ROOT.gStyle.SetOptTitle(0)

    out_file = ROOT.TFile("paper_combined_comparisons.root", "RECREATE")

    # --- 2. Create Canvas and divide it into 2 columns, 1 row ---
    # Width is doubled (1600) to accommodate two 800x600 plots side-by-side
    canvas = ROOT.TCanvas("c_combined", "CCQE Comparisons", 1600, 600)
    canvas.Divide(2, 1)

    latex_font = 132

    # Helper function to style histograms and legends so we don't repeat code
    def draw_pad(pad_num, hist_list, pad_title):
        canvas.cd(pad_num)
        ROOT.gPad.SetTickx(1)
        ROOT.gPad.SetTicky(1)
        ROOT.gPad.SetBottomMargin(0.15) 
        ROOT.gPad.SetLeftMargin(0.12)

        legend = ROOT.TLegend(0.60, 0.65, 0.85, 0.85) 
        legend.SetNColumns(1) 
        legend.SetTextFont(latex_font) 
        legend.SetBorderSize(0) 
        legend.SetFillStyle(0)  
        legend.SetTextSize(0.03) 

        # Find global max for this specific pad
        global_max = 0.0
        for name, h in hist_list:
            if h.GetMaximum() > global_max:
                global_max = h.GetMaximum()

        for i, (name, h) in enumerate(hist_list):
            color, style = styles[i % len(styles)]
            
            h.SetLineColor(color)
            h.SetLineStyle(style)
            h.SetLineWidth(2) 
            h.SetStats(0) 

            h.GetXaxis().SetTitleFont(latex_font)
            h.GetXaxis().SetLabelFont(latex_font)
            h.GetYaxis().SetTitleFont(latex_font)
            h.GetYaxis().SetLabelFont(latex_font)

            h.GetXaxis().SetTitleSize(0.04) 
            h.GetXaxis().SetLabelSize(0.03) 
            h.GetXaxis().SetTitleOffset(1.0)
            
            h.GetYaxis().SetTitleSize(0.04) 
            h.GetYaxis().SetLabelSize(0.03) 
            h.GetYaxis().SetTitleOffset(0.9)

            if i == 0:
                h.SetMaximum(global_max * 1.2)
                h.Draw("HIST") 
            else:
                h.Draw("HIST SAME")
                
            legend.AddEntry(h, name, "l")

        legend.Draw()
        # Return legend so it doesn't get garbage collected by Python before the canvas saves
        return legend 

    # --- 3. Draw Pad 1 (Leading Proton) ---
    leg1 = draw_pad(1, hist_list_proton, "Leading Proton")

    # --- 4. Draw Pad 2 (Excitation Energy) ---
    leg2 = draw_pad(2, hist_list_ex, "Excitation Energy")

    # --- 5. Save output ---
    out_file.cd() 
    canvas.Write("ccqe_combined") 
    out_file.Close()