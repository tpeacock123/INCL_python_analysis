import array
import sys
import numpy as np
import ROOT
from incl_analysis_functions import *

ROOT.gSystem.Load("libNEUTROOTClass.so")
ROOT.gSystem.Load("libNEUTOutput.so")
ROOT.gSystem.Load("libNEUTReWeight.so")

ROOT.TH1.AddDirectory(False)


def _excitation_plots_2d(
    filename, output_filename="fsi_excitation_2d_plots.root"
):
    
    channel_map = {
        intChannel_CCQE.noCascadeFSI: "noCascadeFSI",
        intChannel_CCQE.noCascadeFSIPhoton: "noCascadeFSIPhoton",
        intChannel_CCQE.qeDeEX: "QE_deex",
        intChannel_CCQE.oneProton: "oneProton",
        intChannel_CCQE.multipleNucleon: "multipleNucleon_noCluster",
        intChannel_CCQE.nuclearCluster: "nuclearCluster",
        intChannel_CCQE.protonPion: "protonPion",
        intChannel_CCQE.other: "other",
    }

    nbins_2d = 65
    xmin_2d = -5.0
    xmax_2d = 60.0

    # Initialize 2D histograms per channel
    hists_2d = {}
    for ch_enum, ch_name in channel_map.items():
        h2 = ROOT.TH2D(
            f"h2_pre_vs_post_{ch_name}",
            f"{ch_name};Pre-FSI Excitation Energy [MeV];Post-FSI Excitation Energy [MeV]",
            nbins_2d,
            xmin_2d,
            xmax_2d,
            nbins_2d,
            xmin_2d,
            xmax_2d,
        )
        h2.SetStats(0)
        hists_2d[ch_enum] = h2

    # Inclusive 2D histogram
    h2_inclusive = ROOT.TH2D(
        "h2_pre_vs_post_inclusive",
        "Inclusive;Pre-FSI Excitation Energy [MeV];Post-FSI Excitation Energy [MeV]",
        nbins_2d,
        xmin_2d,
        xmax_2d,
        nbins_2d,
        xmin_2d,
        xmax_2d,
    )
    h2_inclusive.SetStats(0)

    f = ROOT.TFile.Open(filename)
    if not f or f.IsZombie():
        print(f"Error opening file {filename}")
        return None, None

    t = f.Get("neuttree")

    for i, event in enumerate(t):
        if i % 100000 == 0:
            print(f"Event number {i}")
        if i == 500000:
            break

        nvect = event.vectorbranch
        nvect_class = nvect_reader(nvect)

        eventType = nvect_class.eventType
        intChannel = nvect_class.intChannel

        if eventType == EventType.MF:
            if (intChannel == intChannel_CCQE.muOnly) or (
                intChannel == intChannel_CCQE.neutronPion
            ):
                continue

            pre_FSI_E = nvect_class.neutrino_E()
            post_FSI_E = nvect_class.post_casc_E_star()

            if pre_FSI_E is None or post_FSI_E is None:
                continue

            if intChannel in hists_2d:
                hists_2d[intChannel].Fill(pre_FSI_E, post_FSI_E)
                h2_inclusive.Fill(pre_FSI_E, post_FSI_E)

    out_file = ROOT.TFile(output_filename, "RECREATE")
    h2_inclusive.Write()
    for h2 in hists_2d.values():
        h2.Write()
    out_file.Close()

    return hists_2d, h2_inclusive


def excitation_plots_2d(
    filename, output_filename="fsi_excitation_2d_plots_Pn.root"
):
    
    channel_map = {
        intChannel_CCQE.noCascadeFSI: "no Cascade/deex",
#       intChannel_CCQE.noCascadeFSIPhoton: "noCascadeFSIPhoton",
        intChannel_CCQE.qeDeEX: "TNo cascade, deexcitation",
        intChannel_CCQE.oneProton: "one proton",
        intChannel_CCQE.multipleNucleon: "Cascade, nucleons",
        intChannel_CCQE.nuclearCluster: "Cascsade, Clusters",
#        intChannel_CCQE.protonPion: "protonPion",
#        intChannel_CCQE.other: "other",
    }

    nbins_2d = 65
    xmin_2d = -5.0
    xmax_2d = 60.0

    # Initialize 2D histograms per channel
    hists_2d = {}
    for ch_enum, ch_name in channel_map.items():
        h2 = ROOT.TH2D(
            f"h2_pre_vs_post_{ch_name}",
            f"{ch_name};Input proton momentum [MeV];Post-FSI Excitation Energy [MeV]",
            nbins_2d,
            xmin_2d,
            xmax_2d,
            nbins_2d,
            xmin_2d,
            xmax_2d,
        )
        h2.SetStats(0)
        hists_2d[ch_enum] = h2

    # Inclusive 2D histogram
    h2_inclusive = ROOT.TH2D(
        "h2_pre_vs_post_inclusive",
        "Inclusive; Input proton momentum [MeV];Post-FSI Excitation Energy [MeV]",
        nbins_2d,
        xmin_2d,
        xmax_2d,
        nbins_2d,
        xmin_2d,
        xmax_2d,
    )
    h2_inclusive.SetStats(0)

    f = ROOT.TFile.Open(filename)
    if not f or f.IsZombie():
        print(f"Error opening file {filename}")
        return None, None

    t = f.Get("neuttree")

    for i, event in enumerate(t):
        if i % 100000 == 0:
            print(f"Event number {i}")
        if i == 500000:
            break

        nvect = event.vectorbranch
        nvect_class = nvect_reader(nvect)

        eventType = nvect_class.eventType
        intChannel = nvect_class.intChannel

        if eventType == EventType.MF:
            if (intChannel == intChannel_CCQE.muOnly) or (
                intChannel == intChannel_CCQE.neutronPion
            ):
                continue

            pre_FSI_E = nvect_class.PreFSIProtMom()
            post_FSI_E = nvect_class.post_casc_E_star()

            if pre_FSI_E is None or post_FSI_E is None:
                continue

            if intChannel in hists_2d:
                hists_2d[intChannel].Fill(pre_FSI_E, post_FSI_E)
                h2_inclusive.Fill(pre_FSI_E, post_FSI_E)

    out_file = ROOT.TFile(output_filename, "RECREATE")
    h2_inclusive.Write()
    for h2 in hists_2d.values():
        h2.Write()
    out_file.Close()

    return hists_2d, h2_inclusive