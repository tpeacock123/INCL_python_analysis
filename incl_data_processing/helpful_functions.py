import ROOT
import sys 
import numpy as np
ROOT.gSystem.Load("libNEUTROOTClass.so")
ROOT.gSystem.Load("libNEUTOutput.so")
ROOT.gSystem.Load("libNEUTReWeight.so")
from enum import Enum

class EventType(Enum):
    MF = 0
    SRC = 1
    twop2h = 2
    other = 3

class nvect_reader:
    def __init__(self, nvect_,flavour = 2212):
        self.nvect = nvect_
        self.nopart = self.nvect.Npart()
        self.novert = self.nvect.NnucFsiVert()
        self.nosteps = self.nvect.NnucFsiStep()
        self.beam_flavour = flavour
        self.isnofsi = True if self.nopart <=5 else False

        self.eventType = self.event_type()
        
    def event_type(self):
        if self.twop2h() == True:
            return EventType.twop2h
        elif self.src() == True:
            return EventType.SRC
        else:
            return EventType.MF
        
    def muon_energy(self):
        muon_energy = 0
        for i in range(self.nopart):
            pinfo = self.nvect.PartInfo(i)
            if ( i==2 and pinfo.fPID == 13):
                muon_energy =np.sqrt((pinfo.fP.X()**2 + pinfo.fP.Y()**2 +pinfo.fP.Z()**2)+pinfo.fMass**2)-pinfo.fMass
            elif(i ==2 and pinfo.fPID != 13):
                muon_energy = -999.9

        return muon_energy
    
    def casc_proton_energy(self):
        p_casc_energy = 0
        alive = 2
        pid = 0
        for i in range(self.nopart):
            pinfo = self.nvect.PartInfo(i)

            if (i == 3 and ((pinfo.fPID == abs(2112)) or pinfo.fPID == abs(2212))):
                if (pinfo.fIsAlive == 0):
                    p_casc_energy =np.sqrt((pinfo.fP.X()**2 + pinfo.fP.Y()**2 +pinfo.fP.Z()**2)+pinfo.fMass**2)-pinfo.fMass
                    pid =pinfo.fPID
                    alive=pinfo.fIsAlive
                elif (pinfo.fIsAlive == 1):
                    p_casc_energy =np.sqrt((pinfo.fP.X()**2 + pinfo.fP.Y()**2 +pinfo.fP.Z()**2)+pinfo.fMass**2)-pinfo.fMass
                    pid =pinfo.fPID
                    alive=pinfo.fIsAlive

        return p_casc_energy,alive, pid
    
    def proton_energy(self, flavour):
        p_casc_energy = []
        alive = 2
        pid = 0
        for i in range(self.nopart):
            pinfo = self.nvect.PartInfo(i)
            if ((pinfo.fPID == abs(flavour))):
                if (pinfo.fIsAlive == 1):
                    p_casc_energy.append(np.sqrt((pinfo.fP.X()**2 + pinfo.fP.Y()**2 +pinfo.fP.Z()**2)+pinfo.fMass**2)-pinfo.fMass)
                    pid =pinfo.fPID
                    alive=pinfo.fIsAlive

        return p_casc_energy

    def proton_energy(self):
        p_casc_energy = []
        alive = 2
        pid = 0

        for i in range(self.nopart):

            pinfo = self.nvect.PartInfo(i)

            if pinfo.fPID == abs(2212):
                if (pinfo.fIsAlive == 1):
                    p_casc_energy.append(np.sqrt((pinfo.fP.X()**2 + pinfo.fP.Y()**2 +pinfo.fP.Z()**2)))
                    pid =pinfo.fPID
                    alive=pinfo.fIsAlive

        return p_casc_energy
    
    def proton_momentum_per_channel(self):
        p_casc_energy = []
        proton = False
        nuclear_remnant = False
        nucleonCounter = 0
        clusterCounter = 0
        transparentProton = False
        transparentEvent = False
        pion = False
        photonCounter = False
        proton_mom_prefsi = 0
        deexcitation_event = False
        init_proton_mom = None

        for i in range(self.nopart):


            pinfo = self.nvect.PartInfo(i)

            if pinfo.fStatus == 10 and (pinfo.fPID < 10000):
                deexcitation_event = True
            
            if pinfo.fPID == abs(2212): #getting proton momentum 
                if (pinfo.fIsAlive == 0 and self.nvect.ParentIdx(i)==2):
                    init_proton_mom = np.array([pinfo.fP.X(),pinfo.fP.Y(), pinfo.fP.Z()])
                    proton_mom_prefsi = np.linalg.norm(init_proton_mom)

                elif (pinfo.fIsAlive == 1 and self.nvect.ParentIdx(i)==2 and self.nopart == 4):
                    p_casc_energy.append(np.sqrt((pinfo.fP.X()**2 + pinfo.fP.Y()**2 +pinfo.fP.Z()**2)))
                    transparentProton = True
                    nucleonCounter += 1
                    proton = True
                                              
                elif (pinfo.fIsAlive == 1):
                    p_casc_energy.append(np.sqrt((pinfo.fP.X()**2 + pinfo.fP.Y()**2 +pinfo.fP.Z()**2)))
                    pre_fsi_proton_mom = np.array([pinfo.fP.X(),pinfo.fP.Y(), pinfo.fP.Z()])
                    nucleonCounter += 1
                    proton = True
                     
                    #if (self.nvect.ParentIdx(i)==4):
                    #    pre_fsi_proton_mom = np.array([pinfo.fP.X(),pinfo.fP.Y(), pinfo.fP.Z()])
                    #    #print(init_proton_mom - pre_fsi_proton_mom )
                    #    #print(np.linalg.norm(init_proton_mom - pre_fsi_proton_mom ))
                    #    if (np.linalg.norm(pre_fsi_proton_mom - init_proton_mom) < 0.007):
                    #        transparentProton = True 
                    if (self.nvect.ParentIdx(i) != 2 and init_proton_mom is not None):
                        if (np.linalg.norm(pre_fsi_proton_mom - init_proton_mom) < 0.007):
                            transparentProton = True


            elif pinfo.fPID == abs(2112):
                if (pinfo.fIsAlive == 1):
                    nucleonCounter += 1

            elif pinfo.fPID > abs(10000):
                if (pinfo.fIsAlive == 1):
                    nuclear_remnant = True

            elif 200 < abs(pinfo.fPID) < 250:
                if (pinfo.fIsAlive == 1):
                    pion = True
            
            elif ((10000 > pinfo.fPID > 1000) and (pinfo.fIsAlive == 1)):
                if((pinfo.fPID != 2212) and (pinfo.fPID != 2112) and (pinfo.fPID != 6011)and (pinfo.fPID != 13)):
                   clusterCounter += 1

            elif (pinfo.fPID == 0 and pinfo.fIsAlive == 1):
                photonCounter +=1
                 
        return p_casc_energy,nuclear_remnant, nucleonCounter, clusterCounter, transparentProton, pion, photonCounter, proton, proton_mom_prefsi,deexcitation_event

    def src(self):

        if self.isnofsi == True:
            if self.nopart == 5:
                return True 
            else:
                return False
        else:
            src_counter = 0
            proton_mom2 = 0
            init_proton_mom = 0
            for i in range(self.nopart):
                pinfo = self.nvect.PartInfo(i)
                if pinfo.fPID == abs(2112):
                    if (pinfo.fIsAlive == 0) and (self.nvect.ParentIdx(i)== 0) and (pinfo.fStatus == -1):
                        init_proton_mom = np.array([pinfo.fP.X(),pinfo.fP.Y(), pinfo.fP.Z()])
                        continue

                if pinfo.fPID == abs(2212) or pinfo.fPID == abs(2112):
                    if pinfo.fStatus == 5:
                        proton_mom2 = np.array([-pinfo.fP.X(),-pinfo.fP.Y(), -pinfo.fP.Z()])

                        if np.linalg.norm(init_proton_mom - proton_mom2) < 0.1:
                            return True   
            return False
    
    def twop2h(self):   
        twop2hcounter = 0
        for i in range(self.nopart):
            pinfo = self.nvect.PartInfo(i)
            if pinfo.fPID == abs(2112) or pinfo.fPID == abs(2212):
                 if(self.nvect.ParentIdx(i)== 0 and pinfo.fStatus == -1):
                     twop2hcounter +=1
        if twop2hcounter > 1:
            return True
        return 0

    def missing_E_calc(self):

        E_nu = 0.0
        E_lep = 0.0 
        T_had = 0.0
        for i in range(self.nopart):
            pinfo = self.nvect.PartInfo(i)
            if pinfo.fPID == abs(14):
                E_nu = np.linalg.norm(np.array([pinfo.fP.X(),pinfo.fP.Y(), pinfo.fP.Z()]))
            if pinfo.fPID == abs(13):
                E_lep = np.sqrt(np.linalg.norm(np.array([pinfo.fP.X(),pinfo.fP.Y(), pinfo.fP.Z()]))**2 + 105.0**2)
            if pinfo.fPID == abs(2212) and (self.nvect.ParentIdx(i)==2):
                T_had = np.sqrt(np.linalg.norm(np.array([pinfo.fP.X(),pinfo.fP.Y(), pinfo.fP.Z()]))**2 + 938.00**2) - 938.00

        return E_nu - E_lep - T_had
             
    def neutron_mom(self):
        for i in range(self.nopart):
            pinfo = self.nvect.PartInfo(i)
            if ((pinfo.fPID == abs(2112)) and  (pinfo.fStatus == -1)):
                return np.linalg.norm(np.array([pinfo.fP.X(),pinfo.fP.Y(), pinfo.fP.Z()]))

    def remnant(self):
        A = 0
        Z = 0
        for i in range(self.nopart):
            pinfo = self.nvect.PartInfo(i)
            if ((pinfo.fPID > 1000) and pinfo.fIsAlive == 1):
                if (pinfo.fPID ==2112):
                    A +=1
                elif (pinfo.fPID ==2212):
                    A +=1
                    Z+=1
                else:
                    remnant = str(pinfo.fPID)
                    Z+= int(remnant[0])

                    if remnant[2] != 0:
                        A+= 10*int(remnant[2]) + int(remnant[3])
                    else:
                        A+= int(remnant[3])
        if A != 12:
            self.Print()
        print(A,Z)
        return 0 
    
    def get_deltaPT(self):
        proton_momentum = np.asarray([0,0,0])

        for i in range(self.nopart):

            pinfo = self.nvect.PartInfo(i)
            #get neutrino momentum direction
            if pinfo.fPID == abs(14):
                neutrino_mometum = np.asarray([pinfo.fP.X(), pinfo.fP.Y() ,pinfo.fP.Z()])

            if pinfo.fPID == abs(13):
                lepton_momentum = np.asarray([pinfo.fP.X(), pinfo.fP.Y() ,pinfo.fP.Z()])

            if pinfo.fPID == abs(2212): #getting proton momentum  
                if (pinfo.fIsAlive == 1):
                    proton_mometum_placeholder = np.asarray([pinfo.fP.X(), pinfo.fP.Y() ,pinfo.fP.Z()])
                    if(np.linalg.norm(proton_mometum_placeholder) > np.linalg.norm(proton_momentum)):
                        proton_momentum = proton_mometum_placeholder
        
        if np.linalg.norm(proton_momentum) == 0:
            return -999.99
        else:
            normal_v = neutrino_mometum/np.linalg.norm(neutrino_mometum) 
            proton_momentum_long = np.dot(proton_momentum, normal_v) * normal_v
            proton_momentum_trans = proton_momentum - proton_momentum_long

            lept_momentum_long = np.dot(lepton_momentum, normal_v) * normal_v
            lept_momentum_trans = lepton_momentum - lept_momentum_long
            #print(np.linalg.norm(proton_momentum_trans + lept_momentum_trans))
            DPT = proton_momentum_trans + lept_momentum_trans
            #print(DPT)
            DPT_norm = np.linalg.norm(DPT)
            DaT = np.arccos( ( np.dot(-lept_momentum_trans,DPT) /(  np.linalg.norm(lept_momentum_trans) *  np.linalg.norm(DPT) )   ) )
            #print(DaT)
            #print(np.rad2deg(DaT))
            return DPT_norm, np.rad2deg(DaT)

    def interaction_counter(self):
        noflags = 0 
        interaction = 0
        for i in range(self.nosteps):
            nucfsisinfo = self.nvect.NucFsiStepInfo(i)
            noflags += nucfsisinfo.fPBFlag 
            interaction_flag = nucfsisinfo.fVertFlagStep

            if(-4 < interaction_flag < 4 ):
                interaction += 1 

        return noflags, interaction, 0
    
    def is_pauliblocked(self):
        PB = 0
        for i in range(self.novert):
            nucfsivinfo = self.nvect.NucFsiVertInfo(i)
            flag = nucfsivinfo.fVertFlag
            if flag  < 0:
                flaglist = [-1]
            else:
                flaglist = [int(i) for i in str(flag)]
            if len(flaglist) == 4:
                PB = 1
        return PB

    def neutrino_multiplicity(self):
        nucleon_no = 0
        for i in range(self.nopart):
            pinfo = self.nvect.PartInfo(i)
            if (i >= 2 and (pinfo.fPID == 2212 or pinfo.fPID == 2112) and pinfo.fIsAlive == 1):
                nucleon_no +=1

        if nucleon_no == 0:
            for j in range(self.nopart):
                pinfo = self.nvect.PartInfo(j)
                #print(j, " ", pinfo.fIsAlive,"       ",self.nvect.ParentIdx(j),"              ", pinfo.fPID, "    ",pinfo.fP.X(), pinfo.fP.Y(), pinfo.fP.Z(), pinfo.fMass)
        
        return nucleon_no

    def neutrino_daughter_nuclei_energies(self):
        init_ke =[]
        flavour = []
        for i in range(self.nopart): 
            pinfo = self.nvect.PartInfo(i)
            if i >= 1 and (pinfo.fPID == 2212 or pinfo.fPID == 2112) and pinfo.fIsAlive == 1:
                #print(i, " ", pinfo.fIsAlive,"       ",self.nvect.ParentIdx(i),"              ", pinfo.fPID, "    ",pinfo.fP.X(), pinfo.fP.Y(), pinfo.fP.Z(), pinfo.fMass)
                init_ke.append(np.sqrt((pinfo.fP.X()**2 + pinfo.fP.Y()**2 +pinfo.fP.Z()**2)+pinfo.fMass**2)-pinfo.fMass)
                flavour.append(pinfo.fPID)
        return init_ke, flavour

    def Print(self):
        print("~~~~~~~~~Particles~~~~~~~") 
        print("ID | V. ID  | flag | Parent ID | particle | mom | mass")
        for i in range(self.nopart):
            pinfo = self.nvect.PartInfo(i)
            print(i, " ", pinfo.fIsAlive,"       " , pinfo.fStatus,"  ",self.nvect.ParentIdx(i),"              ", pinfo.fPID, "    ",pinfo.fP.X(), pinfo.fP.Y(), pinfo.fP.Z(), pinfo.fMass)
            pinfo.fP.X()**2 + pinfo.fP.X()**2 +pinfo.fP.Y()**2
            
            P = (pinfo.fP.X()**2 + pinfo.fP.Y()**2 +pinfo.fP.Z()**2)
            print("energy", np.sqrt(P+pinfo.fMass**2)-pinfo.fMass)
    
    def kinetic_energy(self,pinfo):
        return np.sqrt((pinfo.fP.X()**2 + pinfo.fP.Y()**2 +pinfo.fP.Z()**2)+pinfo.fMass**2)-pinfo.fMass


