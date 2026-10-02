import ROOT
import sys 
import numpy as np
import math
ROOT.gSystem.Load("libNEUTROOTClass.so")
ROOT.gSystem.Load("libNEUTOutput.so")
ROOT.gSystem.Load("libNEUTReWeight.so")
from enum import Enum

def reader(filename):
    
    f = ROOT.TFile(filename)
    t = f.Get("neuttree")   
    channel = []

    for event in t:
        
        nvect = event.vectorbranch
        nopart = nvect.Npart()

        print("~~~~~~~~~Particles~~~~~~~") 
        print("ID | V. ID  | flag | Parent ID | particle | mom | mass")
        for i in range(nopart):
            pinfo = nvect.PartInfo(i)
            print(i, " ", pinfo.fIsAlive,"       " , pinfo.fStatus,"  ",nvect.ParentIdx(i),"              ", pinfo.fPID, "    ",pinfo.fP.X(), pinfo.fP.Y(), pinfo.fP.Z(), pinfo.fMass)
            pinfo.fP.X()**2 + pinfo.fP.X()**2 +pinfo.fP.Y()**2
            P = (pinfo.fP.X()**2 + pinfo.fP.Y()**2 +pinfo.fP.Z()**2)

reader("out_pi_sf_10k.root")