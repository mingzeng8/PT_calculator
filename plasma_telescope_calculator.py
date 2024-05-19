from scipy import pi
from scipy.constants import c, physical_constants

class plasma_telescope:
    def __init__(self, laser_wavelength_um = 0.8, laser_power_PW = 1., w_0_um = 10., np_cm__3 = None, d_um = None, a_2 = None, w_2_um = None):
        self.laser_wavelength_um = laser_wavelength_um
        self.laser_power_PW = laser_power_PW
        self.w_0_um = w_0_um
        if np_cm__3 is not None:
            self.set_np_cm__3(np_cm__3)
        if d_um is not None:
            self.d_um = d_um
        if a_2 is not None:
            self.set_a_2(a_2)
        if w_2_um is not None:
            self.set_w_2_um(w_2_um)

    def get_one_o_k_um(self): return self.laser_wavelength_um/2./pi

    one_o_k_um = property(get_one_o_k_um)

    def get_a0(self): return (self.laser_power_PW/2.1491120853293987e-5)**0.5*self.laser_wavelength_um/self.w_0_um

    a0 = property(get_a0)

    def get_zR_um(self): return pi/self.laser_wavelength_um*self.w_0_um**2

    zR_um = property(get_zR_um)

    def get_np_cm__3(self): return self._np_cm__3

    def set_np_cm__3(self, value):
        self._np_cm__3 = value
        self._one_o_kp_um = (4*pi*physical_constants['classical electron radius'][0]*value)**-0.5*1e3
        self._k_o_kp = self._one_o_kp_um/self.one_o_k_um

    np_cm__3 = property(get_np_cm__3, set_np_cm__3)

    def get_one_o_kp_um(self): return self._one_o_kp_um

    def set_one_o_kp_um(self, value):
        self._one_o_kp_um = value
        self._np_cm__3 = 1e6/(4*pi*physical_constants['classical electron radius'][0]*value*value)
        self._k_o_kp = self._one_o_kp_um/self.one_o_k_um

    one_o_kp_um = property(get_one_o_kp_um, set_one_o_kp_um)

    def get_k_o_kp(self): return self._k_o_kp

    def set_k_o_kp(self, value):
        self._k_o_kp = value
        self._one_o_kp_um = value*self.one_o_k_um
        self._np_cm__3 = 1e6/(4*pi*physical_constants['classical electron radius'][0]*self._one_o_kp_um*self._one_o_kp_um)

    k_o_kp = property(get_k_o_kp, set_k_o_kp)

    def get_kpw_0(self): return self.w_0_um/self.one_o_kp_um

    kpw_0 = property(get_kpw_0)

    def get_a_2(self): return self._a_2

    def set_a_2(self, value):
        self._a_2 = value
        self._w_2_um = self.a0*self.w_0_um/value

    a_2 = property(get_a_2, set_a_2)

    def get_w_2_um(self): return self._w_2_um

    def set_w_2_um(self, value):
        self._w_2_um = value
        self._a_2 = self.a0*self.w_0_um/value

    w_2_um = property(get_w_2_um, set_w_2_um)

    def get_kpw_2(self): return self._w_2_um/self.one_o_kp_um

    def set_kpw_2(self, value):
        self.set_w_2_um(value*self.one_o_kp_um)

    kpw_2 = property(get_kpw_2, set_kpw_2)

    def get_kpd(self): return self.d_um/self.one_o_kp_um

    def set_kpd(self, value):
        self.d_um = value*self.one_o_kp_um

    kpd = property(get_kpd, set_kpd)

    def match_know_a2(self):
        self.k_o_kp = (self.laser_power_PW/2.1775058026562748e-6/self.a_2**3)**0.5
        self.set_kpd_know_kpw_2()

    def match_know_np(self):
        #self.set_np_cm__3() should have been already called before this function
        self.a_2 = (self.laser_power_PW/2.1775058026562748e-6/self.k_o_kp**2)**(1/3)
        self.set_kpd_know_kpw_2()

    def set_kpd_know_kpw_2(self):
        self.set_kpd( ((self.kpw_2/self.kpw_0)**2-1.)**0.5*self.kpzeta() )

    def set_kpw_2_know_kpd(self):
        self.set_w_2_um((1.+(self.kpd/self.kpzeta())**2)**0.5*self.w_0_um)

    def kpzR(self):
        return self.zR_um/self.one_o_kp_um

    def kpdM(self):
        return self.kpzR()*(self.P_o_Pc()-1.)**0.5

    def w1_um(self):
        return self.w_0_um*(1+(self.d_um/self.zR_um)**2)**0.5

    def kpw1(self):
        return self.w1_um()/self.one_o_kp_um

    def a1(self):
        return self.a0*self.w_0_um/self.w1_um()

    def kpzeta(self):
        #return 0.475*(self.kpw_0)**2*self.k_o_kp-1.2*self.k_o_kp-13.
        return 0.95*self.kpzR()-1.2*self.k_o_kp-13.

    def kpd_eff(self):
        tmp = self.kpzR()*((self.a0/self.kpw_0**2)**(2/3)-1.)**0.5
        if tmp>self.kpd: return tmp
        else: return self.kpd

    def d_eff_um(self):
        return self.kpd_eff()*self.one_o_kp_um

    def a0kpw_0(self):
        return self.a0*self.kpw_0

    def P_o_Pc(self):
    #calculate P over P_c
        return self.a0kpw_0()**2/32.

    def dM_um(self):
        return self.kpdM()*self.one_o_kp_um

    def zeta_um(self):
        return self.kpzeta()*self.one_o_kp_um

    def kpl(self):
        return 21.*self.kpd/(self.kpw_0)**2.08

    def l_um(self):
        return self.kpl()*self.one_o_kp_um

    def kp_L_dephasing(self):
    # This only works if matched
        return self.a_2**0.5 * self.k_o_kp**2 *4./3.
        #return self.a_2**0.5 * self.k_o_kp**2 *4

    def L_dephasing_um(self):
    # This only works if matched
        return self.kp_L_dephasing() * self.one_o_kp_um

    def tau_opt_fs(self):
    # This only works if matched
        return self.w_2_um*2.e9/3./c

    def omega_p_tau_opt(self):
    # This only works if matched
        return self.kpw_2*2./3.
    
    def energy_gain_GeV(self):
    # This only works if matched
        return self.a_2*2./3.*self.k_o_kp**2*0.0005109989499961642

    def print_parameters(self):
        column_size = 19
        precision = 5
        print('-'*(column_size*4))
        print("{2:{0}}{3:<{0}.{1}}| {4:{0}}{5:<{0}.{1}}".format(column_size, precision, "n_p [cm^-3]", self.np_cm__3, "P [PW]", self.laser_power_PW))
        print("{2:{0}}{3:<{0}.{1}}| {4:{0}}{5:<{0}.{1}}".format(column_size, precision, "1/k_p [um]", self.one_o_kp_um, "a_0", self.a0))
        print("{2:{0}}{3:<{0}.{1}}| {4:{0}}{5:<{0}.{1}}".format(column_size, precision, "lambda [um]", self.laser_wavelength_um, "k/k_p", self.k_o_kp))
        print("{2:{0}}{3:<{0}.{1}}| {4:{0}}{5:<{0}.{1}}".format(column_size, precision, "w_0 [um]", self.w_0_um, "k_p w_0", self.kpw_0))
        print('-'*(column_size*4))
        print("{2:{0}}{3:<{0}.{1}}| {4:{0}}{5:<{0}.{1}}".format(column_size, precision, "z_R [um]", self.zR_um, "k_p z_R", self.kpzR()))
        print("{2:{0}}{3:<{0}.{1}}| {4:{0}}{5:<{0}.{1}}".format(column_size, precision, "zeta [um]", self.zeta_um(), "k_p zeta", self.kpzeta()))
        print("{2:{0}}{3:<{0}.{1}}| {4:{0}}{5:<{0}.{1}}".format(column_size, precision, "d [um]", self.d_um, "k_p d", self.kpd))
        print("{2:{0}}{3:<{0}.{1}}| {4:{0}}{5:<{0}.{1}}".format(column_size, precision, "d_eff [um]", self.d_eff_um(), "k_p d_eff", self.kpd_eff()))
        print("{2:{0}}{3:<{0}.{1}}| {4:{0}}{5:<{0}.{1}}".format(column_size, precision, "d_M [um]", self.dM_um(), "k_p d_M", self.kpdM()))
        print("{2:{0}}{3:<{0}.{1}}| {4:{0}}{5:<{0}.{1}}".format(column_size, precision, "l [um]", self.l_um(), "k_p l", self.kpl()))
        print('-'*(column_size*4))
        print("{2:{0}}{3:<{0}.{1}}| {4:{0}}{5:<{0}.{1}}".format(column_size, precision, "a_1", self.a1(), "a_2", self.a_2))
        print("{2:{0}}{3:<{0}.{1}}| {4:{0}}{5:<{0}.{1}}".format(column_size, precision, "w_1 [um]", self.w1_um(), "k_p w_1", self.kpw1()))
        print("{2:{0}}{3:<{0}.{1}}| {4:{0}}{5:<{0}.{1}}".format(column_size, precision, "w_2 [um]", self.w_2_um, "k_p w_2", self.kpw_2))
        print("{2:{0}}{3:<{0}.{1}}| {4:{0}}{5:<{0}.{1}}".format(column_size, precision, "P/P_c", self.P_o_Pc(), "a_1/(k_p w_1)^2", self.a1()/self.kpw1()**2))
        print('-'*(column_size*4))
        print("{2:{0}}{3:<{0}.{1}}| {4:{0}}{5:<{0}.{1}}".format(column_size, precision, "L_dephasing [um]", self.L_dephasing_um(), "k_p L_dephasing", self.kp_L_dephasing()))
        print("{2:{0}}{3:<{0}.{1}}| {4:{0}}{5:<{0}.{1}}".format(column_size, precision, "tau_opt [fs]", self.tau_opt_fs(), "omega_p * tau_opt", self.omega_p_tau_opt()))
        print("{2:{0}}{3:<{0}.{1}}| {4:{0}}{5:<{0}.{1}}".format(column_size, precision, "FWHM duration [fs]", self.tau_opt_fs()*1.1774100225154747, "Energy gain [GeV]", self.energy_gain_GeV()))
        print('-'*(column_size*4))

if __name__ == '__main__':
    # Set a_2, match w_2 and n_p
    PT1 = plasma_telescope(laser_wavelength_um = 0.8, laser_power_PW = 5/27.*0.9394372786996513, w_0_um = 5.0, a_2 = 6.2)
    PT1.match_know_a2()
    # Set n_p, d, and calculate other parameters
    #PT1 = plasma_telescope(laser_wavelength_um = 0.8, laser_power_PW = 5/27.*0.9394372786996513, w_0_um = 5.0, np_cm__3 = 5.2e18, d_um = 25., a_2 = None, w_2_um = None)
    #PT1.set_kpw_2_know_kpd()
    PT1.print_parameters()
