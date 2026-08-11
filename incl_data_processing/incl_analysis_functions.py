import ROOT
import sys 
import numpy as np
import math
ROOT.gSystem.Load("libNEUTROOTClass.so")
ROOT.gSystem.Load("libNEUTOutput.so")
ROOT.gSystem.Load("libNEUTReWeight.so")
from enum import Enum

class EventType(Enum):
    MF = 0
    SRC = 1
    twop2h = 2
    piProd = 3
    NC = 4
    NC_mf = 5
    NC_SRC = 6
    NC_piProd = 7

class intChannel_CCQE(Enum):
    noCascadeFSI = 0 
    qeDeEX = 1
    oneProton = 2
    multipleNucleon = 3
    nuclearCluster = 4
    protonPion = 5
    muOnly = 6
    neutronPion = 7
    other = 8
    noCascadeFSIPhoton = 9

class intChannel_CC0pi(Enum):
    noCascadeFSI = 0 # two protons leave without changing 
    qeDeEX = 1 # two protons leave without changing, plus deexcitation
    elasticProton = 2 # both protons change energy. 
    multipleNucleons = 3 
    nuclearClusters = 4
    protonPion = 5
    muOnly = 6
    other = 7
    neutronPion = 8

def create_histo(name, title, color, fill_style, data,n_bins,x_min,x_max, line = False):

    h = ROOT.TH1F(name, title, n_bins, x_min, x_max)
    for val in data:
        h.Fill(val)
    
    h.SetLineColor(color-2)  
    h.SetLineWidth(2)
    h.SetFillColor(color)
    h.SetFillStyle(fill_style) 
    
    if line:
        h.SetLineStyle(fill_style) 

    
    # h.SetStats(0)
    return h

import numpy as np
import numpy as np

class PionProcessing:
    def __init__(self, nvect):
        self.nvect = nvect
        self.nopart = self.nvect.Npart()  # Stored locally to loop over particles
        
        # 1. Pion variables
        self.pion_count = 0
        self.pion_momenta = 0.0
        self.pion_angles = 0.0
        self.pion_3mom = np.zeros(3)

        # 2. Nucleon variables
        self.initProton = None
        self.nucleon_outgoing_momentum = None
        self.outgoing_nucleon_3mom = np.zeros(3)

        # 3. Lepton variables
        self.lepton_momentum = None
        self.lepton_3mom = np.zeros(3)
        self.cos_theta_l = None

        # 4. Missing kinematics variables
        self.p_miss = None
        self.e_miss = None

        # 5. Global & Invariant Kinematics
        self.Q2 = None
        self.W = None
        self.nu = None

        # 6. Angular correlations
        self.theta_l_pi = None
        self.theta_N_pi = None

        # 7. Single-Transverse Kinematic Variables (STVs)
        self.delta_pT = None
        self.delta_alphaT = None
        self.delta_phiT = None

        # Execute extraction and calculation methods
        self.extract_pion_data()
        self.extract_nucleon_data()
        self.extract_lepton_data()
        self.calculate_missing_kinematics()
        self.calculate_global_kinematics()
        self.calculate_angular_correlations()
        self.calculate_transverse_variables()

    def extract_pion_data(self):
        """Extracts pion multiplicity, 3-momentum, and scattering angle."""
        pion_pids = {211, -211, 111}  # pi+, pi-, pi0
        
        for i in range(self.nopart):
            pinfo = self.nvect.PartInfo(i)
            if pinfo.fPID in pion_pids and pinfo.fIsAlive == 1:
                self.pion_3mom = np.array([pinfo.fP.X(), pinfo.fP.Y(), pinfo.fP.Z()])
                p_mag = np.linalg.norm(self.pion_3mom)
                
                # Polar angle wrt beam (z-axis)
                angle = np.arccos(np.clip(self.pion_3mom[2] / p_mag, -1.0, 1.0)) if p_mag > 0 else 0.0
                
                self.pion_momenta = p_mag
                self.pion_angles = angle
                self.pion_count += 1

    def extract_nucleon_data(self):
        """Extracts pre-FSI proton momentum and final-state nucleon momentum."""
        for i in range(self.nopart):
            pinfo = self.nvect.PartInfo(i)
            
            if pinfo.fPID in (2212, 2112):  # Protons and neutrons
                vec = np.array([pinfo.fP.X(), pinfo.fP.Y(), pinfo.fP.Z()])
                p_mag = np.linalg.norm(vec)
                
                # Initial pre-FSI proton
                if pinfo.fIsAlive == 0:
                    self.initProton = p_mag
                
                # Primary outgoing final-state nucleon
                if pinfo.fIsAlive == 1:
                    self.nucleon_outgoing_momentum = p_mag
                    self.outgoing_nucleon_3mom = vec

    def extract_lepton_data(self):
        """Extracts outgoing charged lepton momentum and cosine angle."""
        charged_leptons = {11, -11, 13, -13, 15, -15}  # e, mu, tau
        
        for i in range(self.nopart):
            pinfo = self.nvect.PartInfo(i)
            if pinfo.fPID in charged_leptons and pinfo.fIsAlive == 1:
                vec = np.array([pinfo.fP.X(), pinfo.fP.Y(), pinfo.fP.Z()])
                self.lepton_3mom = vec
                self.lepton_momentum = np.linalg.norm(vec)
                
                if self.lepton_momentum > 0:
                    self.cos_theta_l = vec[2] / self.lepton_momentum
                break

    def calculate_missing_kinematics(self):
        """Calculates missing momentum (p_miss) and missing energy (e_miss)."""
        p_nu = np.array([
            self.nvect.PartInfo(0).fP.X(),
            self.nvect.PartInfo(0).fP.Y(),
            self.nvect.PartInfo(0).fP.Z()
        ]) 

        p_miss_vec = p_nu - self.lepton_3mom - self.outgoing_nucleon_3mom  - self.pion_3mom
        self.p_miss = np.linalg.norm(p_miss_vec)
        
        M_N = 938.272
        E_nu = np.linalg.norm(p_nu)
        E_l = np.sqrt(self.lepton_momentum**2 + 105.658**2) if self.lepton_momentum is not None else 0.0
        E_N = np.sqrt(self.nucleon_outgoing_momentum**2 + M_N**2) if self.nucleon_outgoing_momentum is not None else 0.0
        E_pi = np.sqrt(self.pion_momenta**2 + 139.4**2) if self.pion_momenta is not None else 0.0
        
        self.e_miss = E_nu - E_l - (E_N - M_N) - E_pi

        # 3. Missing energy for 1-pion production (Eq. 4 adapted)
        # Note: E_N is total relativistic energy (E_N = T_N + M_N)


    def calculate_global_kinematics(self):
        """Calculates Q2, Energy Transfer (nu), and Hadronic Invariant Mass (W)."""
        if self.lepton_momentum is None:
            return

        p_nu = np.array([
            self.nvect.PartInfo(0).fP.X(),
            self.nvect.PartInfo(0).fP.Y(),
            self.nvect.PartInfo(0).fP.Z()
        ])
        
        E_nu = np.linalg.norm(p_nu)
        E_l = np.sqrt(self.lepton_momentum**2 + 105.658**2)
        
        # 1. Energy Transfer nu (q0)
        self.nu = E_nu - E_l

        # 2. Q^2 = -q^2 = |p_nu - p_l|^2 - (E_nu - E_l)^2
        q_3vec = p_nu - self.lepton_3mom
        self.Q2 = np.dot(q_3vec, q_3vec) - (self.nu)**2

        # 3. Hadronic Invariant Mass W = sqrt((E_N + E_pi)^2 - |p_N + p_pi|^2)
        if self.nucleon_outgoing_momentum is not None and self.pion_momenta is not None:
            M_N = 938.272
            E_N = np.sqrt(self.nucleon_outgoing_momentum**2 + M_N**2)
            E_pi = np.sqrt(self.pion_momenta**2 + 139.4**2)
            
            E_had = E_N + E_pi
            p_had_3vec = self.outgoing_nucleon_3mom + self.pion_3mom
            p_had_sq = np.dot(p_had_3vec, p_had_3vec)
            
            W_sq = E_had**2 - p_had_sq
            self.W = np.sqrt(max(0.0, W_sq))

    def calculate_angular_correlations(self):
        """Calculates opening angles between final-state particles."""
        p_l_mag = self.lepton_momentum
        p_N_mag = self.nucleon_outgoing_momentum
        p_pi_mag = self.pion_momenta

        # Opening angle theta_(lepton, pion)
        if p_l_mag and p_pi_mag:
            cos_l_pi = np.dot(self.lepton_3mom, self.pion_3mom) / (p_l_mag * p_pi_mag)
            self.theta_l_pi = np.arccos(np.clip(cos_l_pi, -1.0, 1.0))

        # Opening angle theta_(nucleon, pion)
        if p_N_mag and p_pi_mag:
            cos_N_pi = np.dot(self.outgoing_nucleon_3mom, self.pion_3mom) / (p_N_mag * p_pi_mag)
            self.theta_N_pi = np.arccos(np.clip(cos_N_pi, -1.0, 1.0))

    def calculate_transverse_variables(self):
        """Calculates Single-Transverse Kinematic Variables (delta_pT, delta_alphaT, delta_phiT)."""
        if self.lepton_momentum is None or self.nucleon_outgoing_momentum is None or self.pion_momenta is None:
            return

        # 2D transverse momentum vectors (x, y components)
        pT_l = self.lepton_3mom[:2]
        pT_N = self.outgoing_nucleon_3mom[:2]
        pT_pi = self.pion_3mom[:2]

        # 1. Transverse momentum imbalance vector
        delta_pT_vec = pT_l + pT_N + pT_pi
        self.delta_pT = np.linalg.norm(delta_pT_vec)

        pT_l_mag = np.linalg.norm(pT_l)

        # 2. Transverse imbalance angle (delta_alphaT)
        if pT_l_mag > 0 and self.delta_pT > 0:
            cos_alphaT = np.dot(-pT_l, delta_pT_vec) / (pT_l_mag * self.delta_pT)
            self.delta_alphaT = np.arccos(np.clip(cos_alphaT, -1.0, 1.0))

        # 3. Transverse azimuthal imbalance (delta_phiT)
        pT_had_vec = pT_N + pT_pi
        pT_had_mag = np.linalg.norm(pT_had_vec)

        if pT_l_mag > 0 and pT_had_mag > 0:
            cos_phiT = np.dot(-pT_l, pT_had_vec) / (pT_l_mag * pT_had_mag)
            self.delta_phiT = np.arccos(np.clip(cos_phiT, -1.0, 1.0))
                
import numpy as np


class CCQEProcessing:

    def __init__(self, nvect):
        self.nvect = nvect
        self.nopart = (
            self.nvect.Npart()
        )  # Stored locally to loop over particles

        # 1. Outgoing & Initial Nucleon variables
        self.initProton = None
        self.nucleon_outgoing_momentum = None
        self.outgoing_nucleon_3mom = np.zeros(3)
        self.nucleon_angle = None
        self.nucleon_pid = None

        # 2. Lepton variables
        self.lepton_momentum = None
        self.lepton_3mom = np.zeros(3)
        self.cos_theta_l = None
        self.lepton_pid = None

        # 3. Missing kinematics variables
        self.p_miss = None
        self.e_miss = None

        # 4. Global & Invariant Kinematics
        self.Q2 = None
        self.W = None
        self.nu = None

        # 5. Angular correlations
        self.theta_l_N = None  # Opening angle between outgoing lepton and nucleon

        # 6. Single-Transverse Kinematic Variables (STVs)
        self.delta_pT = None
        self.delta_alphaT = None
        self.delta_phiT = None

        # Execute extraction and calculation methods
        self.extract_nucleon_data()
        self.extract_lepton_data()
        self.calculate_missing_kinematics()
        self.calculate_global_kinematics()
        self.calculate_angular_correlations()
        self.calculate_transverse_variables()

    def _get_lepton_mass(self, pid):
        """Returns mass in MeV/c^2 based on PID."""
        abs_pid = abs(pid)
        if abs_pid == 11:
            return 0.511  # Electron
        elif abs_pid == 13:
            return 105.658  # Muon
        elif abs_pid == 15:
            return 1776.86  # Tau
        return 105.658  # Default to muon

    def _get_nucleon_mass(self, pid):
        """Returns mass in MeV/c^2 based on PID."""
        if pid == 2212:
            return 938.272  # Proton
        elif pid == 2112:
            return 939.565  # Neutron
        return 938.272  # Default to proton

    def extract_nucleon_data(self):
        """Extracts pre-interaction initial nucleon momentum and primary final-state nucleon momentum."""
        for i in range(self.nopart):
            pinfo = self.nvect.PartInfo(i)

            if pinfo.fPID in (2212, 2112):  # Proton or Neutron
                vec = np.array([pinfo.fP.X(), pinfo.fP.Y(), pinfo.fP.Z()])
                p_mag = np.linalg.norm(vec)

                # Initial pre-FSI target nucleon
                if pinfo.fIsAlive == 0:
                    self.initProton = p_mag

                # Primary outgoing final-state nucleon
                if pinfo.fIsAlive == 1:
                    self.nucleon_outgoing_momentum = p_mag
                    self.outgoing_nucleon_3mom = vec
                    self.nucleon_pid = pinfo.fPID
                    if p_mag > 0:
                        angle = np.arccos(
                            np.clip(vec[2] / p_mag, -1.0, 1.0)
                        )
                        self.nucleon_angle = angle

    def extract_lepton_data(self):
        """Extracts primary outgoing charged lepton momentum and cosine angle."""
        charged_leptons = {11, -11, 13, -13, 15, -15}  # e, mu, tau

        for i in range(self.nopart):
            pinfo = self.nvect.PartInfo(i)
            if pinfo.fPID in charged_leptons and pinfo.fIsAlive == 1:
                vec = np.array([pinfo.fP.X(), pinfo.fP.Y(), pinfo.fP.Z()])
                self.lepton_3mom = vec
                self.lepton_momentum = np.linalg.norm(vec)
                self.lepton_pid = pinfo.fPID

                if self.lepton_momentum > 0:
                    self.cos_theta_l = vec[2] / self.lepton_momentum
                break

    def calculate_missing_kinematics(self):
        """Calculates missing momentum (p_miss) and missing energy (e_miss) for a 2-body final state."""
        if (
            self.lepton_momentum is None
            or self.nucleon_outgoing_momentum is None
        ):
            return

        p_nu = np.array([
            self.nvect.PartInfo(0).fP.X(),
            self.nvect.PartInfo(0).fP.Y(),
            self.nvect.PartInfo(0).fP.Z(),
        ])

        # Missing momentum 3-vector: p_miss = p_nu - p_l - p_N
        p_miss_vec = p_nu - self.lepton_3mom - self.outgoing_nucleon_3mom
        self.p_miss = np.linalg.norm(p_miss_vec)

        # Particle energies
        m_l = self._get_lepton_mass(self.lepton_pid)
        m_N = self._get_nucleon_mass(self.nucleon_pid)

        E_nu = np.linalg.norm(p_nu)
        E_l = np.sqrt(self.lepton_momentum**2 + m_l**2)
        E_N = np.sqrt(self.nucleon_outgoing_momentum**2 + m_N**2)

        # Missing energy: E_miss = E_nu - E_l - T_N (where T_N = E_N - M_N)
        self.e_miss = E_nu - E_l - (E_N - m_N)

    def calculate_global_kinematics(self):
        """Calculates Q2, Energy Transfer (nu), and Hadronic Invariant Mass (W)."""
        if self.lepton_momentum is None:
            return

        p_nu = np.array([
            self.nvect.PartInfo(0).fP.X(),
            self.nvect.PartInfo(0).fP.Y(),
            self.nvect.PartInfo(0).fP.Z(),
        ])

        m_l = self._get_lepton_mass(self.lepton_pid)
        m_N = self._get_nucleon_mass(self.nucleon_pid)

        E_nu = np.linalg.norm(p_nu)
        E_l = np.sqrt(self.lepton_momentum**2 + m_l**2)

        # 1. Energy Transfer nu (q0)
        self.nu = E_nu - E_l

        # 2. Q^2 = -q^2 = |p_nu - p_l|^2 - (E_nu - E_l)^2
        q_3vec = p_nu - self.lepton_3mom
        self.Q2 = np.dot(q_3vec, q_3vec) - (self.nu) ** 2

        # 3. Hadronic Invariant Mass W = sqrt(M_N^2 + 2*M_N*nu - Q^2)
        W_sq = m_N**2 + 2.0 * m_N * self.nu - self.Q2
        self.W = np.sqrt(max(0.0, W_sq))

    def calculate_angular_correlations(self):
        """Calculates opening angle between final-state lepton and nucleon (theta_l_N)."""
        p_l_mag = self.lepton_momentum
        p_N_mag = self.nucleon_outgoing_momentum

        if p_l_mag and p_N_mag:
            cos_l_N = np.dot(self.lepton_3mom, self.outgoing_nucleon_3mom) / (
                p_l_mag * p_N_mag
            )
            self.theta_l_N = np.arccos(np.clip(cos_l_N, -1.0, 1.0))

    def calculate_transverse_variables(self):
        """Calculates Single-Transverse Kinematic Variables (delta_pT, delta_alphaT, delta_phiT) for CCQE."""
        if (
            self.lepton_momentum is None
            or self.nucleon_outgoing_momentum is None
        ):
            return

        # 2D transverse momentum vectors (x, y components)
        pT_l = self.lepton_3mom[:2]
        pT_N = self.outgoing_nucleon_3mom[:2]

        # 1. Transverse momentum imbalance vector: delta_pT = pT_l + pT_N
        delta_pT_vec = pT_l + pT_N
        self.delta_pT = np.linalg.norm(delta_pT_vec)

        pT_l_mag = np.linalg.norm(pT_l)

        # 2. Transverse imbalance angle (delta_alphaT)
        if pT_l_mag > 0 and self.delta_pT > 0:
            cos_alphaT = np.dot(-pT_l, delta_pT_vec) / (
                pT_l_mag * self.delta_pT
            )
            self.delta_alphaT = np.arccos(np.clip(cos_alphaT, -1.0, 1.0))

        # 3. Transverse azimuthal imbalance (delta_phiT)
        pT_N_mag = np.linalg.norm(pT_N)
        if pT_l_mag > 0 and pT_N_mag > 0:
            cos_phiT = np.dot(-pT_l, pT_N) / (pT_l_mag * pT_N_mag)
            self.delta_phiT = np.arccos(np.clip(cos_phiT, -1.0, 1.0))

class nvect_reader:
    def __init__(self, nvect_,flavour = 2212):
        
        self.intChannel = None
        self.HMPMom = None
        self.PreFSIProtMom = 0
        self.nvect = nvect_
        self.nopart = self.nvect.Npart()
        self.novert = self.nvect.NnucFsiVert()
        self.nosteps = self.nvect.NnucFsiStep()
        self.beam_flavour = flavour
        #if self.nopart <=5:
        #    self.isnofsi = True 
        #else:
        self.isnofsi = False
        self.istransparent = False
        self.fsiProton = 0.0
        self.neutrino = self.nu()
        self.NC_flag = self.NC()

        self.nocasc_pi = self.nocasc_pion()
        #self.nocasc_pi = True
        #self.nocasc_piprocessing = None
        if self.nocasc_pi == True:
            self.nocasc_piprocessing = PionProcessing(self.nvect)
            #self.nocasc_piprocessing = CCQEProcessing(self.nvect)

        if self.NC_flag == False:
            self.lepton_mass = 105.00
            if self.nubar == True:
                self.incoming_nucleon = 2112
                self.outgoing_nucleon = 2212
                self.incoming_nu = self.neutrino
                self.outgoing_lep = -(abs(self.neutrino)-1)
                self.outgoing_mass = 939.565
            else:
                self.incoming_nucleon = 2112
                self.outgoing_nucleon = 2212
                self.incoming_nu = self.neutrino
                self.outgoing_lep = self.neutrino - 1
                self.outgoing_mass = 938.272

        else:
            self.lepton_mass = 0.0
            if self.nubar == True:
                self.incoming_nucleon = 2212
                self.outgoing_nucleon = 2212
                self.incoming_nu = self.neutrino
                self.outgoing_lep = self.neutrino
                self.outgoing_mass = 938.272
            else:
                self.incoming_nucleon = 2212
                self.outgoing_nucleon = 2212
                self.incoming_nu = self.neutrino
                self.outgoing_lep = self.neutrino
                self.outgoing_mass = 938.272

        self.eventType = self.event_type()
        if self.eventType == EventType.MF or  self.eventType == EventType.NC_mf :
            self.intChannel = self.interaction_channel_CCQE()
        elif (self.eventType == EventType.SRC or self.eventType == EventType.NC_SRC or  self.eventType == EventType.twop2h):
            self.intChannel = self.interaction_channel_CC0pi()

        if (self.eventType == EventType.MF) or (self.eventType == EventType.NC_mf)  or (self.eventType == EventType.SRC) or (self.eventType == EventType.twop2h):
            self.E_miss = self.missing_E_calc()
            self.P_miss = self.missing_mom()
        if (self.eventType == EventType.NC_mf or self.eventType == EventType.NC_SRC ):
            self.E_miss = self.missing_E_calc_NC()
            print(self.E_miss)
            self.P_miss = self.missing_mom()

    def nocasc_pion(self):
        for i in range(self.nopart):
            pinfo = self.nvect.PartInfo(i)
            if (pinfo.fIsAlive == 1 and pinfo.fPID == 211 and self.nopart == 5):
                return True
    
    def nu(self):
         for i in range(self.nopart):
            pinfo = self.nvect.PartInfo(i)
            if pinfo.fPID == -14:
                self.nubar = True
                return -14
            elif pinfo.fPID == 14:
                self.nubar = False
                return 14
            if pinfo.fPID == -16:
                self.nubar = True
                return -16
            elif pinfo.fPID == 16:
                self.nubar = False
                return 16
            if pinfo.fPID == -12:
                self.nubar = True
                return -12
            elif pinfo.fPID == 12:
                self.nubar = False
                return 12
                
    def event_type(self):

        if self.NC_Pi_prod() == True:
            return EventType.NC_piProd
        elif self.pion_prod() == True:
            return EventType.piProd
        elif self.src() == True & self.NC() == True: 
            return EventType.NC_SRC
        elif self.NC() == True:
            return EventType.NC_mf

        elif self.twop2h() == True:
            return EventType.twop2h
        elif self.src() == True:
            return EventType.SRC
        else:
            return EventType.MF
        
    def multiplicity(self):  
        multiplicity = 0
        particles = []
        for i in range(self.nopart):
            pinfo = self.nvect.PartInfo(i)  
            if pinfo.fIsAlive == 1:
                if pinfo.fPID != self.outgoing_lep:
                    particles.append(pinfo.fPID)
                    

        return particles
           
    def proton_momentum_per_channel_CCQE(self):
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
        deex_counter = 0

        for i in range(self.nopart):

            pinfo = self.nvect.PartInfo(i)
            
            if pinfo.fPID == self.outgoing_nucleon: #getting proton momentum 
                if (pinfo.fIsAlive == 0 and self.nvect.ParentIdx(i)==2):
                    init_proton_mom = np.array([pinfo.fP.X(),pinfo.fP.Y(), pinfo.fP.Z()])
                    proton_mom_prefsi = np.linalg.norm(init_proton_mom)
                    self.fsiProton = proton_mom_prefsi

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
                     
                    if (self.nvect.ParentIdx(i)==4):
                        pre_fsi_proton_mom = np.array([pinfo.fP.X(),pinfo.fP.Y(), pinfo.fP.Z()])
                        

                    if (self.nvect.ParentIdx(i) != 2 and init_proton_mom is not None):
                        if (np.linalg.norm(pre_fsi_proton_mom - init_proton_mom)/np.linalg.norm(init_proton_mom) < 0.005):
                            transparentProton = True


            elif pinfo.fPID == (2112):
                if (pinfo.fIsAlive == 1):
                    nucleonCounter += 1

            elif 200 < abs(pinfo.fPID) < 250:
                if (pinfo.fIsAlive == 1):
                    pion = True

            elif((pinfo.fPID > 1000000) and (pinfo.fIsAlive == 1)):
                   clusterCounter += 1

            elif (pinfo.fPID == 22 and pinfo.fIsAlive == 1):
                photonCounter +=1

            if (pinfo.fStatus == 10 and pinfo.fIsAlive == 1 and pinfo.fPID != 22 ):
                deex_counter +=1

            if deex_counter > 1:
                deexcitation_event = True
                 
        return p_casc_energy,nuclear_remnant, nucleonCounter, clusterCounter, transparentProton, pion, photonCounter, proton, proton_mom_prefsi,deexcitation_event

    def proton_momentum_per_channel_SRC(self):
        p_casc_energy = []

        init_src_nucleon_mom = 0.0
        init_main_proton_mom = 0.0
        fsi_neutron_energy = []
        fsi_proton_energy = []
        src_partner_flavour = 0
        nucleonCounter = 0
        clusterCounter = 0
        nuclear_remnant = False
        deex_counter = 0
        pion = False
        proton = False
        photonCounter = False
        proton_mom_prefsi = 0
        deexcitation_event = False
        init_proton_mom = None

        for i in range(self.nopart):
            pinfo = self.nvect.PartInfo(i)

            if pinfo.fPID ==(self.incoming_nucleon):
                if(pinfo.fIsAlive == 0 and self.nvect.ParentIdx(i)==0 and pinfo.fStatus == 7):
                    init_src_nucleon_mom = np.linalg.norm(np.array([pinfo.fP.X(),pinfo.fP.Y(), pinfo.fP.Z()]))

                elif (pinfo.fIsAlive == 1):
                    #proton = True
                    fsi_neutron_energy.append(np.sqrt((pinfo.fP.X()**2 + pinfo.fP.Y()**2 +pinfo.fP.Z()**2)))
    
                    nucleonCounter += 1
            
            if pinfo.fPID == self.outgoing_nucleon: #getting proton momentum 
                if (pinfo.fIsAlive == 0 and self.nvect.ParentIdx(i)==2 and pinfo.fStatus ==7):
                    init_main_proton_mom = np.linalg.norm(np.array([pinfo.fP.X(),pinfo.fP.Y(), pinfo.fP.Z()]))

                elif(pinfo.fIsAlive == 0 and self.nvect.ParentIdx(i)==0 and pinfo.fStatus ==7):
                    init_src_nucleon_mom = np.linalg.norm(np.array([pinfo.fP.X(),pinfo.fP.Y(), pinfo.fP.Z()]))


                elif (pinfo.fIsAlive == 1):
                    proton = True
                    fsi_proton_energy.append(np.sqrt((pinfo.fP.X()**2 + pinfo.fP.Y()**2 +pinfo.fP.Z()**2)))
    
                    nucleonCounter += 1

            elif 200 < (pinfo.fPID) < 250:
                if (pinfo.fIsAlive == 1):
                    pion = True
            
            elif((pinfo.fPID > 1000000) and (pinfo.fIsAlive == 1)):
                   clusterCounter += 1

            elif (pinfo.fPID == 22 and pinfo.fIsAlive == 1):
                photonCounter +=1

            if (pinfo.fStatus == 10 and pinfo.fIsAlive == 1):
                deex_counter +=1

            if deex_counter > 1:
                deexcitation_event = True
        
        all_fsi_energies = fsi_neutron_energy + fsi_proton_energy

        main_present = any(math.isclose(init_main_proton_mom, energy, abs_tol=0.01) for energy in all_fsi_energies)
        src_present = any(math.isclose(init_src_nucleon_mom, energy, abs_tol=0.01) for energy in all_fsi_energies)

        if main_present and src_present:
            transparentNucleons = 2
        elif main_present or src_present:
            transparentNucleons = 1
        else:
            transparentNucleons = 0
                 
        return fsi_proton_energy, nuclear_remnant, nucleonCounter, clusterCounter, transparentNucleons, pion, photonCounter, deexcitation_event,proton

    def proton_momentum_per_channel_2p2h(self):
        p_casc_energy = []

        init_second_nucleon_mom = 0.0
        init_main_proton_mom = 0.0
        fsi_neutron_energy = []
        fsi_proton_energy = []
        second_partner_flavour = 0
        nucleonCounter = 0
        clusterCounter = 0
        nuclear_remnant = False
        pion = False
        proton = False
        photonCounter = False
        proton_mom_prefsi = 0
        deexcitation_event = False
        init_proton_mom = None
        deex_counter = 0

        for i in range(self.nopart):
            pinfo = self.nvect.PartInfo(i)

            if pinfo.fPID ==2112:
                if(pinfo.fIsAlive == 0 and self.nvect.ParentIdx(i)==2 and pinfo.fStatus == 7):
                    init_main_nucleon_mom = np.linalg.norm(np.array([pinfo.fP.X(),pinfo.fP.Y(), pinfo.fP.Z()]))

                if(pinfo.fIsAlive == 0 and self.nvect.ParentIdx(i)==3 and pinfo.fStatus == 7):
                    init_second_nucleon_mom = np.linalg.norm(np.array([pinfo.fP.X(),pinfo.fP.Y(), pinfo.fP.Z()]))
                    second_partner_flavour = 2212

                elif (pinfo.fIsAlive == 1):
                    #proton = True
                    fsi_neutron_energy.append(np.sqrt((pinfo.fP.X()**2 + pinfo.fP.Y()**2 +pinfo.fP.Z()**2)))
    
                    nucleonCounter += 1
            
            if pinfo.fPID == 2212: #getting proton momentum 
                if (pinfo.fIsAlive == 0 and (self.nvect.ParentIdx(i)==2) and pinfo.fStatus ==7):
                    init_main_proton_mom = np.linalg.norm(np.array([pinfo.fP.X(),pinfo.fP.Y(), pinfo.fP.Z()]))

                elif(pinfo.fIsAlive == 0 and self.nvect.ParentIdx(i)==3 and pinfo.fStatus ==7):
                    init_second_nucleon_mom = np.linalg.norm(np.array([pinfo.fP.X(),pinfo.fP.Y(), pinfo.fP.Z()]))
                    second_partner_flavour = 2212

                elif (pinfo.fIsAlive == 1):
                    proton = True
                    fsi_proton_energy.append(np.sqrt((pinfo.fP.X()**2 + pinfo.fP.Y()**2 +pinfo.fP.Z()**2)))
    
                    nucleonCounter += 1


            elif 200 < (pinfo.fPID) < 250:
                if (pinfo.fIsAlive == 1):
                    pion = True
            
            elif((pinfo.fPID > 1000000) and (pinfo.fIsAlive == 1)):
                clusterCounter += 1

            elif (pinfo.fPID == 22 and pinfo.fIsAlive == 1):
                photonCounter +=1

            if (pinfo.fStatus == 10 and pinfo.fIsAlive == 1):
                deex_counter +=1

            if deex_counter > 1:
                deexcitation_event = True

            elif (pinfo.fPID == 0 and pinfo.fIsAlive == 1):
                photonCounter +=1


        all_fsi_energies = fsi_neutron_energy + fsi_proton_energy
        main_present = any(math.isclose(init_main_proton_mom, energy, abs_tol=0.01) for energy in all_fsi_energies)
        second_present = any(math.isclose(init_second_nucleon_mom, energy, abs_tol=0.01) for energy in all_fsi_energies)

        if main_present and second_present:
            transparentNucleons = 2
        elif main_present or second_present:
            transparentNucleons = 1
        else:
            transparentNucleons = 0
                 
        return fsi_proton_energy, nuclear_remnant, nucleonCounter, clusterCounter, transparentNucleons, pion, photonCounter, deexcitation_event,proton

    def NC(self):
        neutrino_counter = 0
        for i in range(self.nopart):
            pinfo = self.nvect.PartInfo(i)
            if abs(pinfo.fPID) == 14:
                neutrino_counter += 1
                if neutrino_counter == 2:
                    return True
        return False

    def pion_prod(self):
        for i in range(self.nopart):
            pinfo = self.nvect.PartInfo(i)
            if pinfo.fPID == abs(211) or pinfo.fPID == abs(111):
                if(self.nvect.ParentIdx(i)== 2 and pinfo.fStatus == 0):
                    return True
                elif(self.nvect.ParentIdx(i)== 2 and pinfo.fStatus == 7):
                    return True
        return False

    def NC_Pi_prod(self):
        neutrino_counter = 0
        neutrino = False
        pion = False

        for i in range(self.nopart):
            pinfo = self.nvect.PartInfo(i)
            if abs(pinfo.fPID) == 14:
                neutrino_counter += 1
                if neutrino_counter == 2:
                    neutrino = True
            if pinfo.fPID == abs(211) or pinfo.fPID == abs(111):
                 if(self.nvect.ParentIdx(i)== 2 and pinfo.fStatus == 0):
                    pion = True

            if(neutrino == True  and pion) == True:
                return True
        return False
   
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
                if pinfo.fPID == self.incoming_nucleon:
                    if (pinfo.fIsAlive == 0) and (self.nvect.ParentIdx(i)== 0) and (pinfo.fStatus == -1):
                        init_proton_mom = np.array([pinfo.fP.X(),pinfo.fP.Y(), pinfo.fP.Z()])
                        continue

                if pinfo.fPID == (2212) or pinfo.fPID == (2112):
                    if pinfo.fStatus == 7:
                        proton_mom2 = np.array([-pinfo.fP.X(),-pinfo.fP.Y(), -pinfo.fP.Z()])

                        if np.linalg.norm(init_proton_mom - proton_mom2) < 0.1:
                            return True   
            return False
    
    def twop2h(self):   
        twop2hcounter = 0
        for i in range(self.nopart):
            pinfo = self.nvect.PartInfo(i)
            if pinfo.fPID == (2112) or pinfo.fPID == (2212):
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
            if pinfo.fPID == self.incoming_nu:
                E_nu = np.linalg.norm(np.array([pinfo.fP.X(),pinfo.fP.Y(), pinfo.fP.Z()]))
            if pinfo.fPID == self.outgoing_lep:
                E_lep = np.sqrt(np.linalg.norm(np.array([pinfo.fP.X(),pinfo.fP.Y(), pinfo.fP.Z()]))**2 + 105.0**2)
            if pinfo.fPID == self.outgoing_nucleon and (self.nvect.ParentIdx(i)==2):
                T_had = np.sqrt(np.linalg.norm(np.array([pinfo.fP.X(),pinfo.fP.Y(), pinfo.fP.Z()]))**2 + 938.00**2) - 938.00

        return E_nu - E_lep - T_had

    def missing_E_calc_NC(self):

        E_nu = 0.0
        E_lep = 0.0 
        T_had = 0.0
        NC_count = False
        for i in range(self.nopart):
            pinfo = self.nvect.PartInfo(i)
            if pinfo.fPID == self.incoming_nu and NC_count == False:
                E_nu = np.linalg.norm(np.array([pinfo.fP.X(),pinfo.fP.Y(), pinfo.fP.Z()]))
                NC_count = True
            elif pinfo.fPID == self.outgoing_lep :
                E_lep = np.sqrt(np.linalg.norm(np.array([pinfo.fP.X(),pinfo.fP.Y(), pinfo.fP.Z()]))**2)
            if pinfo.fPID == self.outgoing_nucleon and (self.nvect.ParentIdx(i)==2):
                T_had = np.sqrt(np.linalg.norm(np.array([pinfo.fP.X(),pinfo.fP.Y(), pinfo.fP.Z()]))**2 + 938.00**2) - 938.00
        print(E_nu, E_lep, T_had)
        return E_nu - E_lep - T_had

    def excitation_E_CCQE(self):
        E_nu = 0.0
        E_lep = 0.0 
        E_had = 0.0
        p_had = 0.0
        NC_catch = False
        for i in range(self.nopart):
            pinfo = self.nvect.PartInfo(i)
            if pinfo.fPID == self.incoming_nu and NC_catch == False:
                p_nu = np.array([pinfo.fP.X(),pinfo.fP.Y(), pinfo.fP.Z()])
                E_nu = np.linalg.norm(np.array([pinfo.fP.X(),pinfo.fP.Y(), pinfo.fP.Z()]))
                NC_catch = True 
            elif (pinfo.fPID) == self.outgoing_lep:
                p_lep = np.array([pinfo.fP.X(),pinfo.fP.Y(), pinfo.fP.Z()])
                E_lep = np.sqrt(np.linalg.norm(np.array([pinfo.fP.X(),pinfo.fP.Y(), pinfo.fP.Z()]))**2 + self.lepton_mass**2)
            if pinfo.fPID == self.outgoing_nucleon and (self.nvect.ParentIdx(i)==2):
                p_had += np.array([pinfo.fP.X(),pinfo.fP.Y(), pinfo.fP.Z()])
                E_had += np.sqrt(np.linalg.norm(np.array([pinfo.fP.X(),pinfo.fP.Y(), pinfo.fP.Z()]))**2 + self.outgoing_mass**2) 

        print(E_nu + 11174.86 - E_lep - E_had)
        E_star = E_nu + 11174.86 - E_lep - E_had
        p_star = 0 # p_nu - p_lep - p_had

        if self.NC_flag == True:
            E = np.sqrt(E_star**2 - np.linalg.norm(p_star**2)) -  10252.61
            print("HERE")
            return E   

        if self.nubar == True:
            E = np.sqrt(E_star**2 - np.linalg.norm(p_star**2)) - 10252.61
        elif self.nubar == False:
            E = np.sqrt(E_star**2 - np.linalg.norm(p_star**2)) - 10254.27

        return E   
    
    def excitation_E_SRC(self):
        E_nu = 0.0
        E_lep = 0.0 
        E_had = 0.0
        p_had = 0.0
        NC_catch = False
        proton_counter = 0
        for i in range(self.nopart):
            pinfo = self.nvect.PartInfo(i)
            if (pinfo.fPID) == self.incoming_nu and NC_catch == False:
                p_nu = np.array([pinfo.fP.X(),pinfo.fP.Y(), pinfo.fP.Z()])
                E_nu = np.linalg.norm(np.array([pinfo.fP.X(),pinfo.fP.Y(), pinfo.fP.Z()]))
                NC_catch = True
            elif (pinfo.fPID) == self.outgoing_lep:
                p_lep = np.array([pinfo.fP.X(),pinfo.fP.Y(), pinfo.fP.Z()])
                E_lep = np.sqrt(np.linalg.norm(np.array([pinfo.fP.X(),pinfo.fP.Y(), pinfo.fP.Z()]))**2 + self.lepton_mass**2)

            if pinfo.fPID == self.outgoing_nucleon and (self.nvect.ParentIdx(i)==2):
                p_had += np.array([pinfo.fP.X(),pinfo.fP.Y(), pinfo.fP.Z()])
                E_had += np.sqrt(np.linalg.norm(np.array([pinfo.fP.X(),pinfo.fP.Y(), pinfo.fP.Z()]))**2 + self.outgoing_mass**2) 
                if pinfo.fPID == (2212):
                    proton_counter += 1
            if (pinfo.fPID == (2212) or pinfo.fPID == (2112)) and (self.nvect.ParentIdx(i)==0) and (pinfo.fStatus ==7 or pinfo.fStatus ==0):
                p_had += np.array([pinfo.fP.X(),pinfo.fP.Y(), pinfo.fP.Z()])
                if pinfo.fPID == (2212):
                    proton_counter += 1
                    E_had += np.sqrt(np.linalg.norm(np.array([pinfo.fP.X(),pinfo.fP.Y(), pinfo.fP.Z()]))**2 + 938.27**2) 
                elif pinfo.fPID == 2112:
                    E_had += np.sqrt(np.linalg.norm(np.array([pinfo.fP.X(),pinfo.fP.Y(), pinfo.fP.Z()]))**2 + 939.56**2) 

        print(E_nu, E_lep, E_had)

        E_star = E_nu + 11174.86 - E_lep - E_had
        print(E_star)
        p_star = 0.0 #p_nu - p_lep - p_had
        print(proton_counter)
        if (proton_counter == 1):
            Nuclear_mass = (5)*938.27 + (5)*939.57  - 64.75
        elif (proton_counter == 2):
            Nuclear_mass = (4)*938.27 + (6)*939.57  - 64.98
        else:
            Nuclear_mass = (6)*938.27 + (4)*939.57  - 60.32
        print(Nuclear_mass)

        E = np.sqrt(E_star**2 - np.linalg.norm(p_star**2)) - Nuclear_mass
        print(E)
        print(E)
        if E < 0.0:
            print(E)

        return E   
    
    def excitation_E_2p2h(self):

        E_nu = 0.0
        E_lep = 0.0 
        E_had = 0.0
        p_had = 0.0
        proton_counter =0
        neutron_counter = 0
        for i in range(self.nopart):
            pinfo = self.nvect.PartInfo(i)
            if pinfo.fPID == self.incoming_nu:
                p_nu = np.array([pinfo.fP.X(),pinfo.fP.Y(), pinfo.fP.Z()])
                E_nu = np.linalg.norm(np.array([pinfo.fP.X(),pinfo.fP.Y(), pinfo.fP.Z()]))
            if pinfo.fPID == self.outgoing_lep:
                p_lep = np.array([pinfo.fP.X(),pinfo.fP.Y(), pinfo.fP.Z()])
                E_lep = np.sqrt(np.linalg.norm(np.array([pinfo.fP.X(),pinfo.fP.Y(), pinfo.fP.Z()]))**2 + 105.0**2)
            if pinfo.fPID == (2212) and (pinfo.fStatus ==7 or pinfo.fStatus == 0) and ((self.nvect.ParentIdx(i)==2) or (self.nvect.ParentIdx(i)==3)):
                p_had += np.array([pinfo.fP.X(),pinfo.fP.Y(), pinfo.fP.Z()])
                E_had += np.sqrt(np.linalg.norm(np.array([pinfo.fP.X(),pinfo.fP.Y(), pinfo.fP.Z()]))**2 + 938.00**2)
                proton_counter += 1
            if pinfo.fPID == (2112) and (pinfo.fStatus ==7 or pinfo.fStatus == 0) and ((self.nvect.ParentIdx(i)==2) or (self.nvect.ParentIdx(i)==3)):
                p_had += np.array([pinfo.fP.X(),pinfo.fP.Y(), pinfo.fP.Z()])
                E_had += np.sqrt(np.linalg.norm(np.array([pinfo.fP.X(),pinfo.fP.Y(), pinfo.fP.Z()]))**2 + 938.00**2) 
                neutron_counter += 1

        E_star = E_nu + 11174.86 - E_lep - E_had
        p_star = p_nu - p_lep - p_had
        if (neutron_counter == 1 and proton_counter == 1):
            Nuclear_mass = (6-proton_counter)*938.37 + (6-neutron_counter)*939.57  - 64.75
        elif (proton_counter == 2):
            Nuclear_mass = (6-proton_counter)*938.37 + (6-neutron_counter)*939.57  - 64.98
        elif (neutron_counter == 2):
            Nuclear_mass = (6-proton_counter)*938.37 + (6-neutron_counter)*939.57  - 60.32
        

        E = np.sqrt(E_star**2 - np.linalg.norm(p_star**2)) - Nuclear_mass
        return E

    def missing_energy_2p2h(self):
        self.Print()

        E_nu = 0.0
        E_lep = 0.0 
        T_had = 0.0
        for i in range(self.nopart):
            pinfo = self.nvect.PartInfo(i)
            if pinfo.fPID == self.incoming_nu:
                E_nu = np.linalg.norm(np.array([pinfo.fP.X(),pinfo.fP.Y(), pinfo.fP.Z()]))
            if pinfo.fPID == self.outgoing_lep:
                E_lep = np.sqrt(np.linalg.norm(np.array([pinfo.fP.X(),pinfo.fP.Y(), pinfo.fP.Z()]))**2 + 105.0**2)
            
            if pinfo.fPID == (2212)and (pinfo.fStatus ==7 or pinfo.fStatus == 0) and ((self.nvect.ParentIdx(i)==2) or (self.nvect.ParentIdx(i)==3)):
                T_had += np.sqrt(np.linalg.norm(np.array([pinfo.fP.X(),pinfo.fP.Y(), pinfo.fP.Z()]))**2 + 938.00**2) - 938.00
            if pinfo.fPID == (2112) and (pinfo.fStatus ==7 or pinfo.fStatus == 0) and ((self.nvect.ParentIdx(i)==2) or (self.nvect.ParentIdx(i)==3)):
                T_had += np.sqrt(np.linalg.norm(np.array([pinfo.fP.X(),pinfo.fP.Y(), pinfo.fP.Z()]))**2 + 938.00**2) - 938.00

        return E_nu - E_lep - T_had

    def missing_mom(self):
        for i in range(self.nopart):
            pinfo = self.nvect.PartInfo(i)
            if ((pinfo.fPID == self.incoming_nucleon) and  (pinfo.fStatus == -1)):
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
            if pinfo.fPID == self.incoming_nu:
                neutrino_mometum = np.asarray([pinfo.fP.X(), pinfo.fP.Y() ,pinfo.fP.Z()])

            if pinfo.fPID == self.outgoing_lep:
                lepton_momentum = np.asarray([pinfo.fP.X(), pinfo.fP.Y() ,pinfo.fP.Z()])

            if pinfo.fPID == self.outgoing_nucleon: #getting proton momentum  
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

    def Print(self):
        print("~~~~~~~~~Particles~~~~~~~") 
        print("ID | V. ID  | flag | Parent ID | particle | mom | mass")
        for i in range(self.nopart):
            pinfo = self.nvect.PartInfo(i)
            print(i, " ", pinfo.fIsAlive,"       " , pinfo.fStatus,"  ",self.nvect.ParentIdx(i),"              ", pinfo.fPID, "    ",pinfo.fP.X(), pinfo.fP.Y(), pinfo.fP.Z(), pinfo.fMass)
            pinfo.fP.X()**2 + pinfo.fP.X()**2 +pinfo.fP.Y()**2
            P = (pinfo.fP.X()**2 + pinfo.fP.Y()**2 +pinfo.fP.Z()**2)
            print("energy", np.sqrt(P+pinfo.fMass**2)-pinfo.fMass)
        
        print(self.eventType.name)
        if (not (self.eventType == EventType.piProd or self.eventType == EventType.NC_piProd) ):
            print(self.intChannel.name)
        return self.eventType.name
    
    def kinetic_energy(self,pinfo):
        return np.sqrt((pinfo.fP.X()**2 + pinfo.fP.Y()**2 +pinfo.fP.Z()**2)+pinfo.fMass**2)-pinfo.fMass

    def interaction_channel_CCQE(self):

        p_casc_energy,nuclear_remnant, nucleonCounter, clusterCounter, transparentProton, pion, photonCounter,proton,prefsi_proton_mom, deex_event = self.proton_momentum_per_channel_CCQE()

        if len(p_casc_energy) > 0:
            self.HMPMom = np.asarray(p_casc_energy).max()
            if (self.eventType != EventType.NC_mf):
                self.DpT,self.DaT = self.get_deltaPT()
            else:
                self.DpT = 0.0
                self.DaT = 0.0

        self.PreFSIProtMom = prefsi_proton_mom
        if(proton != True):
            
            if  (nucleonCounter > 0) or (clusterCounter > 0) or (pion == True):
                return intChannel_CCQE.neutronPion
            else:
                return intChannel_CCQE.muOnly
        
            
        elif (transparentProton == True) and (clusterCounter <= 1) and (proton == True)  and (pion == False) and (deex_event == False):
            if photonCounter == 0:
                return intChannel_CCQE.noCascadeFSI
            else:
                print("here")
                return intChannel_CCQE.noCascadeFSIPhoton

        elif (transparentProton == True) and (deex_event == True) and (proton == True):  
            return intChannel_CCQE.qeDeEX
        
        elif (transparentProton == False) and (proton == True) and (clusterCounter <= 1) and (nucleonCounter == 1)  and (pion == False):
            return intChannel_CCQE.oneProton
            
        elif (transparentProton == False) and (proton == True) and (clusterCounter <= 1) and (nucleonCounter > 1)  and (pion == False):
            return intChannel_CCQE.multipleNucleon
        
        elif (transparentProton == False) and (proton == True) and (clusterCounter >= 1) and (nucleonCounter >= 1):
            return intChannel_CCQE.nuclearCluster

        elif (transparentProton == False) and (proton == True) and (clusterCounter <= 1) and (pion == True):   
            return intChannel_CCQE.protonPion
        else:
            print( p_casc_energy,nuclear_remnant, nucleonCounter, clusterCounter, transparentProton, pion, photonCounter,proton, deex_event)
            self.Print()
            return intChannel_CCQE.other

    def interaction_channel_CC0pi(self):
        if self.eventType == EventType.SRC or self.eventType == EventType.NC_SRC:
            p_casc_energy, nuclear_remnant, nucleonCounter, clusterCounter, transparentNucleons, pion, photonCounter, deex_event,proton = self.proton_momentum_per_channel_SRC()
        if self.eventType == EventType.twop2h:
            p_casc_energy, nuclear_remnant, nucleonCounter, clusterCounter, transparentNucleons, pion, photonCounter, deex_event,proton = self.proton_momentum_per_channel_2p2h()
        
        if len(p_casc_energy) > 0:
            self.HMPMom = np.asarray(p_casc_energy).max()
            self.DpT,self.DaT = self.get_deltaPT()

        if (transparentNucleons > 0):
            self.istransparent = True

        if(proton != True):
            self.intChannel = intChannel_CC0pi.neutronPion
            #self.Print()

            #need to add event with leading neutrons too
            if  (nucleonCounter > 0) or (clusterCounter > 0) or (pion == True):
                return intChannel_CCQE.neutronPion
            else:
                return intChannel_CCQE.muOnly
        
            
        elif (transparentNucleons == 2) and (clusterCounter <= 1) and (proton == True)  and (pion == False) and (photonCounter == False) and (deex_event == False):
            return intChannel_CC0pi.noCascadeFSI
        
        elif (transparentNucleons == 2) and (deex_event == True) and (proton == True):  
            return intChannel_CC0pi.qeDeEX
        
        elif (transparentNucleons < 2) and (proton == True) and (clusterCounter < 1)  and (nucleonCounter == 2)  and (pion == False):
            return intChannel_CC0pi.elasticProton
        
        elif (transparentNucleons <= 1) and (proton == True) and (clusterCounter < 1)  and (nucleonCounter >= 3)  and (pion == False):
            return intChannel_CC0pi.multipleNucleons
        
        elif (transparentNucleons <= 1) and (proton == True) and (clusterCounter >= 1) and (nucleonCounter >= 1):
            return intChannel_CC0pi.nuclearClusters

        elif (transparentNucleons < 2) and (proton == True) and (clusterCounter <= 0)  and (pion == True):   
            return intChannel_CC0pi.protonPion
        else:
            return intChannel_CC0pi.other

    def particles(self):
        particle_list = []
        energy_list = []
        for i in range(self.nopart):
            pinfo = self.nvect.PartInfo(i)
            if pinfo.fIsAlive == 1 and pinfo.fStatus == 10:
                  particle_list.append(pinfo.fPID)
            P = (pinfo.fP.X()**2 + pinfo.fP.Y()**2 +pinfo.fP.Z()**2)
            energy_list.append(np.sqrt(P))
        
        return particle_list,energy_list




def create_ratio(h_in, h_tot):
    h_ratio = h_in.Clone(h_in.GetName() + "_ratio")
    h_ratio.Divide(h_tot) 


    h_ratio.SetFillStyle(0)
    h_ratio.SetLineWidth(2)
    h_ratio.SetLineColor(h_in.GetLineColor()) 
    return h_ratio