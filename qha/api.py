import os
import pathlib
import time

from qha.calculator import Calculator, SamePhDOSCalculator, DifferentPhDOSCalculator
from qha.isotopes import IsotopesCalculator
from qha.basic_io.out import save_x_tp, save_x_tv, save_to_output, make_starting_string, make_tp_info, make_ending_string
from qha.settings import DEFAULT_SETTINGS
import numpy as np

class qha():
    def __init__(self,**kwargs):
        self.settings = kwargs
        self.calc = None
    def init(self,**kwargs):
        """
        Initialize the calculation
        """
        #==== Update settings from kwargs ====#       
        user_settings = self.settings.copy()
        user_settings.update(kwargs)

        for key in ('input', 'calculation',
                'thermodynamic_properties', 'static_only', 'energy_unit',
                'T_MIN', 'NT', 'DT', 'DT_SAMPLE',
                'P_MIN', 'NTV', 'DELTA_P', 'DELTA_P_SAMPLE', 'order', 'p_min_modifier',
                'output_directory', 'high_verbosity'):
            if key not in user_settings:
                user_settings[key] = DEFAULT_SETTINGS[key]
        
        #=== Set critical parameters default values ===#
        if 'input_from_file' not in user_settings:
            user_settings['input_from_file'] = False
        
        calculation_type = user_settings['calculation'].lower()
        if calculation_type == 'single':
            calc = Calculator(user_settings)
            print("You have single-configuration calculation assumed.")
        elif calculation_type == 'same phonon dos':
            calc = SamePhDOSCalculator(user_settings)
            print("You have same-phonon-dos calculation assumed.")
        elif calculation_type == 'different phonon dos':
            calc = DifferentPhDOSCalculator(user_settings)
            print("You have different-phonon-dos calculation assumed.")
        elif calculation_type == 'iso':
            calc = IsotopesCalculator(user_settings)
            print("You have isotopes calculation assumed.")
        else:
            raise ValueError("Calculation type not recognized.")
        
        if user_settings['input_from_file']:
            calc.read_input()
        else:
            calc.set_input(user_settings)        
        self.calc=calc
      
    def check(self,**kwargs):
        """
        From Desired T range to calculate the P range that can be calculate
        """
        #=== Set critical parameters default values ===#
        if self.calc is None:
            self.init(**kwargs)
        self.calc.refine_grid()
        return self.calc.p_tv_gpa[:, 0].max(), self.calc.p_tv_gpa[:, -1].min()

    def run(self,**kwargs):
        """
        Run the calculation
        """
        if self.calc is None:
            self.init(**kwargs)
        user_settings = self.settings.copy()
        P_min, P_max = self.check()
        if  user_settings['P_MIN']<P_min:
            raise ValueError("P_MIN is too small, the minimum P is %f"%P_min)
        if  user_settings['P_MIN']+user_settings['DELTA_P']*user_settings['NTV']>P_max:
            raise ValueError("NTV is too large, the maximum P is %f"%P_max)
        self.output={}
        temperature_array = self.calc.temperature_array
        desired_pressures_gpa = self.calc.desired_pressures_gpa
        temperature_sample = self.calc.temperature_sample_array
        p_sample_gpa = self.calc.pressure_sample_array
        self.output['temperature_array']=temperature_array
        self.output['desired_pressures_gpa']=desired_pressures_gpa
        self.output['temperature_sample']=temperature_sample
        self.output['p_sample_gpa']=p_sample_gpa

        calculation_option = {'F': 'f_tp',
                              'G': 'g_tp',
                              'H': 'h_tp',
                              'U': 'u_tp',
                              'V': 'v_tp',
                              'Cv': 'cv_tp_jmolk',
                              'Cp': 'cp_tp_jmolk',
                              'Bt': 'bt_tp_gpa',
                              'Btp': 'btp_tp',
                              'Bs': 'bs_tp_gpa',
                              'alpha': 'alpha_tp',
                              'gamma': 'gamma_tp',
                              'iso':'isotopes'
                              }
        IsoModes=['iso'] #Only iso first, add other modes later
        if 'iso' in self.calc.settings['thermodynamic_properties']:
            #Check if isotopes calculation is requested. Only in 'Iso' mode such property is supported
            calculation_type = user_settings['calculation'].lower()
            if calculation_type not in IsoModes:
                raise ValueError("Isotopes calculation is only supported in 'Iso' type mode")

        for idx in self.calc.settings['thermodynamic_properties']:
            if idx in ['F', 'G', 'H', 'U']:
                attr_name = calculation_option[idx] + \
                    '_' + self.calc.settings['energy_unit']
                self.output[attr_name]=getattr(self.calc, attr_name)
                self.output[idx]=getattr(self.calc, attr_name)

            if idx == 'V':
                v_bohr3 = calculation_option[idx] + '_' + 'bohr3'
                v_ang3 = calculation_option[idx] + '_' + 'ang3'
                self.output[v_bohr3]=getattr(self.calc, v_bohr3)
                self.output[v_ang3]=getattr(self.calc, v_ang3)
                self.output[idx]=getattr(self.calc, v_ang3)

            if idx in ['Cv', 'Cp', 'Bt', 'Btp', 'Bs', 'alpha', 'gamma']:
                attr_name = calculation_option[idx]
                self.output[attr_name]=getattr(self.calc, attr_name)
                self.output[idx]=getattr(self.calc, attr_name)

            if idx in ['iso']:
            #Add isotopes calculation
            #Todo
                attr_name = calculation_option[idx]
                self.output[attr_name]=getattr(self.calc, attr_name)
                self.output[idx]=getattr(self.calc, attr_name)
        return self.output.copy()
    
    def thermodict(self):
        """
        Return the output dictionary
        """
        return self.output.copy()
    def propertiesTP(self,idx):
        """
        Return the properties
        """
        Table=self.output.copy()[idx]
        def get_value_from_table(T,P):
            #T is in K, P is in GPa
            T_array=self.output['temperature_array']
            P_array=self.output['desired_pressures_gpa']
            if T<T_array[0] or T>T_array[-1]:
                raise ValueError("T is out of range")
            if P<P_array[0] or P>P_array[-1]:
                raise ValueError("P is out of range")
            #Estimate the value from the sorted table using linear interpolation
            T_idx=np.searchsorted(T_array,T)
            P_idx=np.searchsorted(P_array,P)
            if T_idx==0:
                T_idx=1
            if P_idx==0:
                P_idx=1
            T1=T_array[T_idx-1]
            T2=T_array[T_idx]
            P1=P_array[P_idx-1]
            P2=P_array[P_idx]
            V1=Table[T_idx-1,P_idx-1]
            V2=Table[T_idx,P_idx-1]
            V3=Table[T_idx-1,P_idx]
            V4=Table[T_idx,P_idx]
            V12=V1+(V2-V1)*(T-T1)/(T2-T1)
            V34=V3+(V4-V3)*(T-T1)/(T2-T1)
            V=V12+(V34-V12)*(P-P1)/(P2-P1)
            return V
        return get_value_from_table
        
        