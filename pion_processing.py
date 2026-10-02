from enum import Enum

from incl_analysis_functions import *
from default_processing import *
from kdar_processing import *
from paper_plots import *
import numpy as np
from array import array
import sys 
import argparse

import ROOT
import numpy as np

ROOT.gSystem.Load("libNEUTROOTClass.so")
ROOT.gSystem.Load("libNEUTOutput.so")
ROOT.gSystem.Load("libNEUTReWeight.so")
ROOT.TH1.AddDirectory(False)

class PionProcessing:
    def __init__(self, nvect_):
        self.nvect = nvect_
        self.nopart = self.nvect.Npart()  # Stored locally to loop over particles


    

        self.HMPMom = None
        self.PreFSIProtMom = 0
        self.nvect = nvect_
        self.nopart = self.nvect.Npart()
        self.novert = self.nvect.NnucFsiVert()
        self.nosteps = self.nvect.NnucFsiStep()
        self.E_nu = None
            
        
        # 1. Pion variables
        self.pion_count = 0
        self.pion_flavour = None
        self.pion_momenta = None
        self.pion_angles = None
        self.pion_3mom = np.zeros(3)

        self.initNuc = None
        self.outNuc = None
        self.nucleon_outgoing_momentum = None
        self.outgoing_nucleon_3mom = np.zeros(3)

        self.lepton_momentum = None
        self.lepton_3mom = np.zeros(3)
        self.cos_theta_l = None

        self.p_miss = None
        self.e_miss = None

        self.Q2 = None
        self.W = None
        self.q0 = None

        self.theta_l_pi = None
        self.theta_N_pi = None

        self.delta_pT = None
        self.delta_alphaT = None
        self.delta_phiT = None

        self.cos_theta_adler = None
        self.phi_adler = None

        self.pbed = False
        self.NCCheck = False

        for i in range(self.nopart):
            pinfo = self.nvect.PartInfo(i)
            if pinfo.fStatus == 5 and  pinfo.fIsAlive == 0:
                self.pbed = True
        

        # Execute extraction and calculation methods
        if self.pbed == False:

            self.extract_pion_data()
            self.extract_nucleon_data()
            self.extract_lepton_data()
            self.calculate_missing_kinematics()
            self.calculate_global_kinematics()
            self.calculate_angular_correlations()
            self.calculate_transverse_variables()
            self.calculate_adler_angles()
            self.neutrino_xsec()

    def extract_pion_data(self):
        """Extracts pion multiplicity, 3-momentum, and scattering angle."""
        pion_pids = {211, -211, 111}  # pi+, pi-, pi0
        
        for i in range(self.nopart):
            pinfo = self.nvect.PartInfo(i)
            if pinfo.fPID in pion_pids and pinfo.fIsAlive == 1:
                self.pion_flavour = pinfo.fPID
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
                
                # Initial pre-FSI nucleon
                if pinfo.fIsAlive == 0:
                    self.initNuc = p_mag
                                   
                # Primary outgoing final-state nucleon
                if pinfo.fIsAlive == 1 and self.nvect.ParentIdx(i) == 2:
                    if (p_mag != self.initNuc):

                        self.outNuc = pinfo.fPID
                        
                        self.nucleon_outgoing_momentum = p_mag
                        self.outgoing_nucleon_3mom = vec



    def extract_lepton_data(self):
        """Extracts outgoing charged lepton momentum and cosine angle."""
        charged_leptons = {11, -11, 13, -13, 15, -15}  # e, mu, tau
        neutrinos = {12, -12, 14, -14, 14, -14}  # e, mu, tau
        
        for i in range(self.nopart):
            pinfo = self.nvect.PartInfo(i)
            if pinfo.fPID in charged_leptons or pinfo.fPID in neutrinos and pinfo.fIsAlive == 1:
                vec = np.array([pinfo.fP.X(), pinfo.fP.Y(), pinfo.fP.Z()])
                self.lepton_3mom = vec
                self.lepton_momentum = np.linalg.norm(vec)
                
                if self.lepton_momentum > 0:
                    self.cos_theta_l = vec[2] / self.lepton_momentum

            if  pinfo.fPID in neutrinos and pinfo.fIsAlive == 1:
                print("here")
                self.NCCheck = True 


    def calculate_missing_kinematics(self):

        # 1. Neutrino initial momentum
        p_nu = np.array([
            self.nvect.PartInfo(0).fP.X(),
            self.nvect.PartInfo(0).fP.Y(),
            self.nvect.PartInfo(0).fP.Z()
        ]) 
        E_nu = np.linalg.norm(p_nu)

        # 2. Vector missing momentum
        p_final_vec = self.lepton_3mom + self.outgoing_nucleon_3mom + self.pion_3mom
        p_miss_vec = p_nu - p_final_vec
        self.p_miss = float(np.linalg.norm(p_miss_vec))

        # 3. Nucleon mass (PDG values in MeV)
        # 2212 = proton (938.272 MeV), neutron = 939.565 MeV
        M_N = 938.272 if self.outNuc == 2212 else 939.565
        corr = 0
        #if self.pion_count == 0:
        #    corr = 1.29
            
        print(self.NCCheck)
        # 4. Final state energies
        if self.NCCheck == False:
            m_muon = 105.658

        if self.NCCheck == True:
            m_muon = 0.0
        print(m_muon)
        print(self.pion_flavour)
        m_pion = 139.570 if self.pion_flavour == abs(211) else 134.97

        E_l = np.sqrt(self.lepton_momentum**2 + m_muon**2) if self.lepton_momentum is not None else 0.0
        
        # Kinetic energy of nucleon (T_N = E_N - M_N)
        T_N = 0.0
        if self.nucleon_outgoing_momentum is not None:
            E_N = np.sqrt(self.nucleon_outgoing_momentum**2 + M_N**2)
            T_N = E_N - M_N

        # Total energy of pion(s)
        E_pi = 0.0
        if self.pion_momenta is not None:
            # Sums correctly whether pion_momenta is a scalar float or a 1D array
            E_pi = np.sum(np.sqrt(np.asanyarray(self.pion_momenta)**2 + m_pion**2))
        print(E_nu,  E_l , T_N , E_pi , corr)
        # 5. Missing energy (E_miss = E_nu - E_l - T_N - E_pi)
        self.e_miss = E_nu - E_l - T_N - E_pi + corr


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
        self.q0 = E_nu - E_l

        # 2. Q^2 = -q^2 = |p_nu - p_l|^2 - (E_nu - E_l)^2
        q_3vec = p_nu - self.lepton_3mom
        self.Q2 = np.dot(q_3vec, q_3vec) - (self.q0)**2


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

        if p_l_mag and p_pi_mag:
            cos_l_pi = np.dot(self.lepton_3mom, self.pion_3mom) / (p_l_mag * p_pi_mag)
            self.theta_l_pi = np.arccos(np.clip(cos_l_pi, -1.0, 1.0))

        if p_N_mag and p_pi_mag:
            cos_N_pi = np.dot(self.outgoing_nucleon_3mom, self.pion_3mom) / (p_N_mag * p_pi_mag)
            self.theta_N_pi = np.arccos(np.clip(cos_N_pi, -1.0, 1.0))

    def calculate_transverse_variables(self):
        """Calculates Single-Transverse Kinematic Variables (delta_pT, delta_alphaT, delta_phiT)."""
        if self.lepton_momentum is None or self.nucleon_outgoing_momentum is None or self.pion_momenta is None:
            return

        pT_l = self.lepton_3mom[:2]
        pT_N = self.outgoing_nucleon_3mom[:2]
        pT_pi = self.pion_3mom[:2]

        delta_pT_vec = pT_l + pT_N + pT_pi
        self.delta_pT = np.linalg.norm(delta_pT_vec)

        pT_l_mag = np.linalg.norm(pT_l)

        if pT_l_mag > 0 and self.delta_pT > 0:
            cos_alphaT = np.dot(-pT_l, delta_pT_vec) / (pT_l_mag * self.delta_pT)
            self.delta_alphaT = np.arccos(np.clip(cos_alphaT, -1.0, 1.0))

        pT_had_vec = pT_N + pT_pi
        pT_had_mag = np.linalg.norm(pT_had_vec)

        if pT_l_mag > 0 and pT_had_mag > 0:
            cos_phiT = np.dot(-pT_l, pT_had_vec) / (pT_l_mag * pT_had_mag)
            self.delta_phiT = np.arccos(np.clip(cos_phiT, -1.0, 1.0))


    def neutrino_xsec(self):
        self.E_nu = self.nvect.PartInfo(0).fP.E()

    def calculate_adler_angles(self):
        """
        Calculates the Adler frame angles (cos_theta_A, phi_A) for single pion production.
        - phi_A domain: [-pi, pi]
        - z-axis: Direction of 3-momentum transfer q_hcm
        - y-axis: Normal to lepton scattering plane (p_nu x p_l)
        - x-axis: Completes right-handed system (y_axis x z_axis)
        """
        if self.lepton_momentum is None or self.nucleon_outgoing_momentum is None or self.pion_momenta is None:
            return

        # 1. Kinematics setup
        p_nu_3 = np.array([
            self.nvect.PartInfo(0).fP.X(),
            self.nvect.PartInfo(0).fP.Y(),
            self.nvect.PartInfo(0).fP.Z()
        ])
        p_l_3 = np.array(self.lepton_3mom)
        p_N_3 = np.array(self.outgoing_nucleon_3mom)
        p_pi_3 = np.array(self.pion_3mom)

        # Masses (MeV)
        m_mu = 105.658
        m_N = 938.272 if self.outNuc == 2212 else 939.565
        m_pi = 139.570

        # Energies
        E_nu = np.linalg.norm(p_nu_3)
        E_l = np.sqrt(np.dot(p_l_3, p_l_3) + m_mu**2)
        E_N = np.sqrt(np.dot(p_N_3, p_N_3) + m_N**2)
        E_pi = np.sqrt(np.dot(p_pi_3, p_pi_3) + m_pi**2)

        # 2. Hadronic rest frame boost
        E_had = E_N + E_pi
        p_had_3 = p_N_3 + p_pi_3
        beta_vec = p_had_3 / E_had
        beta2 = np.dot(beta_vec, beta_vec)

        if beta2 >= 1.0 or beta2 == 0.0:
            return

        gamma = 1.0 / np.sqrt(1.0 - beta2)

        def boost_to_hcm(E, p_3):
            bp = np.dot(beta_vec, p_3)
            return p_3 + ((gamma - 1.0) / beta2 * bp - gamma * E) * beta_vec

        p_nu_hcm = boost_to_hcm(E_nu, p_nu_3)
        p_l_hcm = boost_to_hcm(E_l, p_l_3)
        p_pi_hcm = boost_to_hcm(E_pi, p_pi_3)

        # 3. Coordinate axes
        q_hcm = p_nu_hcm - p_l_hcm
        q_mag = np.linalg.norm(q_hcm)
        if q_mag == 0:
            return

        z_axis = q_hcm / q_mag

        cross_y = np.cross(p_nu_hcm, p_l_hcm)
        y_mag = np.linalg.norm(cross_y)
        if y_mag == 0:
            return
        y_axis = cross_y / y_mag

        x_axis = np.cross(y_axis, z_axis)

        # 4. Angle extraction
        p_pi_mag = np.linalg.norm(p_pi_hcm)
        if p_pi_mag == 0:
            return

        # Polar angle cos(theta_A)
        self.cos_theta_adler = np.dot(p_pi_hcm, z_axis) / p_pi_mag

        # Azimuthal angle phi_A in [-pi, pi]
        p_pi_x = np.dot(p_pi_hcm, x_axis)
        p_pi_y = np.dot(p_pi_hcm, y_axis)
        self.phi_adler = np.arctan2(p_pi_y, p_pi_x)