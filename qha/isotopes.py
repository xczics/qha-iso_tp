import textwrap
from typing import Dict, Any, Optional

import numpy as np
from lazy_property import LazyProperty

import qha.multi_configurations.different_phonon_dos as different_phonon_dos
import qha.multi_configurations.same_phonon_dos as same_phonon_dos
import qha.tools
from qha.grid_interpolation import FinerGrid
from qha.basic_io.out import save_to_output
from qha.basic_io import read_input
from qha.single_configuration import free_energy
from qha.thermodynamics import *
from qha.type_aliases import Vector
from qha.unit_conversion import gpa_to_ry_b3, ry_b3_to_gpa, b3_to_a3, ry_to_j_mol, ry_to_ev, ry_to_j
from qha.v2p import v2p
from qha.calculator import Calculator
from typing import Iterator, Union, Tuple
from qha.type_aliases import Vector, Array3D,Matrix
import pathlib
from qha.Calclnbeta import GetLnbeta

__all__ = ['IsotopesCalculator']

class IsotopesCalculator(Calculator):
    def __init__(self, user_settings):
        super().__init__(user_settings)
    
    def read_input(self):
        try:
            formula_unit_number, volumes, static_energies, frequencies,q_weights = read_input(
                self.settings['input'])
        except KeyError:
            raise KeyError(
                "The 'input' option must be given in your settings!")
        if 'input_h' not in self.settings and 'input_l' not in self.settings:
            raise KeyError(
                "The 'input_h' and 'input_l' options must be given in your settings!")
        if 'input_h' not in self.settings:
            freqh = frequencies
        else:
            _, _, _, freqh, _ = read_input(self.settings['input_h'])
        if 'input_l' not in self.settings:
            freql = frequencies
        else:
            _, _, _, freql, _ = read_input(self.settings['input_l'])

        if not qha.tools.is_monotonic_decreasing(volumes):
            raise RuntimeError(
                "Check the input file to make sure the volume decreases!")

        self._formula_unit_number: int = formula_unit_number
        self._volumes = volumes
        self._static_energies = static_energies
        self._frequencies = frequencies
        self._q_weights = q_weights
        self._freqh = freqh
        self._freql = freql
        self._an=self.settings['an']
        self._factor_freq_2_THz=self.settings['freq2THz']

    @property
    def freqh(self):
        return self._freqh
    @property
    def freql(self):
        return self._freql
    @property
    def an(self):
        return self._an
    @property
    def factor_freq_2_THz(self):
        return self._factor_freq_2_THz

    @LazyProperty
    def lnbeta_vt(self) -> Matrix:
        """Return the 1000lnbeta at each volumes and temperatures."""
        lnbeta_vt = np.zeros((self.volumes.size, self.temperature_array.size), dtype=float)
        for i in range(self.volumes.size):
            lnbeta_vt[i,:]=GetLnbeta(self._freqh[i,:,:]*self._factor_freq_2_THz,self._freql[i,:,:]*self._factor_freq_2_THz,self.q_weights,self.an,self.temperature_array)
        return lnbeta_vt
    @LazyProperty
    def lnbeta_dense_vt(self) -> Matrix:
        """Return the 1000lnbeta at each densed-grid volumes and temperatures."""
        lnbeta_tv=np.nan_to_num(self.lnbeta_vt.T)
        d = self.settings
        try:
            p_min, p_min_modifier, ntv, order = d['P_MIN'], d['p_min_modifier'], d['NTV'], d['order']
        except KeyError:
            raise KeyError(
                "All the 'P_MIN', 'p_min_modifier', 'NTV', 'order' options must be given in your settings!")
        r = FinerGrid(p_min - p_min_modifier, ntv, order=order)
        _, lnbeta_densed_vt, _ = r.refine_grid(self.volumes,lnbeta_tv ,self._v_ratio)
        return lnbeta_densed_vt
    @LazyProperty
    def isotopes(self) -> Vector:
        return v2p(self.lnbeta_dense_vt, self.p_tv_au, self.desired_pressures)


        
