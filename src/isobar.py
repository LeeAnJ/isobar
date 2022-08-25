import sys
from dataclasses import dataclass
from typing import Literal
import numpy as np
from scipy.interpolate import CubicSpline

import units_literal as units
import converters as conv
from fluid_class import (RP10Fluid)


@dataclass(kw_only=True)
class Parameter:
    value: float
    units: str


@dataclass(kw_only=True)
class Pressure(Parameter):
    units: Literal[units.p]


@dataclass(kw_only=True)
class Temperature(Parameter):
    units: Literal[units.t]


class Isobar:
    def __init__(self,
                 fluid: RP10Fluid,
                 p: Pressure,
                 t_min: Temperature,
                 t_max: Temperature,
                 dt: float) -> None:

        self.fluid = fluid
        self.p_kpa = conv.convert_arg_to_internal_units(p.value, p.units)  # p, [kPa]
        self.t_min_k = conv.convert_arg_to_internal_units(t_min.value, t_min.units)  # t, [k]
        self.t_max_k = conv.convert_arg_to_internal_units(t_max.value, t_max.units)  # t, [k]
        self.dt = dt

        self.array_n = self._calc_array_n() # calc. number of elements in array - n; (index_max = n-1)
        self.t_k = self._set_t_array()
        self.h_jmol = self._set_h_array()

        self.t_bubble_k, self.h_bubble_jmol = self._bubble_point()
        self.t_dew_k, self.h_dew_jmol = self._dew_point()

        self.set_optional_t_for_bubble_point(n_1ph=2, n_2ph=3, dt=3.0)
        self.set_optional_t_for_dew_point(n_1ph=2, n_2ph=3, dt=3.0)
        sys.exit()

        self._insert_bubble_point()
        self._insert_dew_point()

        # cubic spline over (t,h) arrays
        self.cs = CubicSpline(self.t_k, self.h_jmol)

        # for i in range(self.array_n):
        #     print(i, self.t_k[i], self.h_jmol[i])

    def _calc_array_n(self) -> int:
        if self.t_max_k <= self.t_min_k:
            sys.exit('t_max_k <= t_min_k while calc. "_calc_array_n" in Iosbar class__init__ ')
        n = (self.t_max_k - self.t_min_k)/self.dt
        if n <= 2:
            n = 2
        return int(n)+1

    def _set_t_array(self) -> np.array:
        _t_k = np.empty(self.array_n)
        for i in range(self.array_n):
            _t_k[i] = self.t_min_k + (self.t_max_k - self.t_min_k)/(self.array_n-1)*i
        return _t_k

    def _set_h_array(self) -> np.array:
        _h_jmol = np.empty(self.array_n)
        for i in range(self.array_n):
            self.fluid.calc_spec_state(t=(self.t_k[i], 'k'), p=(self.p_kpa, 'kpa'))
            if self.fluid.error.index > 0:
                self.fluid.error.print_and_terminate()
            else:
                _h_jmol[i] = self.fluid.state.get_data(flag='blk', x_symbol='h', x_units='jmol')
        return _h_jmol

    def _bubble_point(self) -> (float, float):   # t_bubble_k, h_bubble_jmol
        self.fluid.calc_sat_state(sat_curve_flag='l', p=(self.p_kpa, 'kpa'))

        if self.fluid.error.index > 0:
            self.fluid.error.print_and_terminate()
        _t_l_k = self.fluid.state.get_data(flag='bubble', x_symbol='t', x_units='k')
        _h_l_jmol = self.fluid.state.get_data(flag='bubble', x_symbol='h', x_units='jmol')
        return _t_l_k, _h_l_jmol

    def _dew_point(self) -> (float, float):  # t_dew_k, h_dew_jmol
        self.fluid.calc_sat_state(sat_curve_flag='v', p=(self.p_kpa, 'kpa'))

        if self.fluid.error.index > 0:
            self.fluid.error.print_and_terminate()
        _t_v_k = self.fluid.state.get_data(flag='dew', x_symbol='t', x_units='k')
        _h_v_jmol = self.fluid.state.get_data(flag='dew', x_symbol='h', x_units='jmol')
        return _t_v_k, _h_v_jmol

    def _insert_bubble_point(self):
        if self.t_min_k < self.t_bubble_k < self.t_max_k:
            # self.t_k[i-1] <= self.t_bubble_k < self.t_k[i]
            index = np.searchsorted(self.t_k, self.t_bubble_k, side='right', sorter=None)
            self.array_n += 1
            self.t_k = np.insert(self.t_k, index, self.t_bubble_k)
            self.h_jmol = np.insert(self.h_jmol, index, self.h_bubble_jmol)
            # print('i = ', index)

    def _insert_dew_point(self):
        if self.t_min_k < self.t_dew_k < self.t_max_k:
            # self.t_k[i-1] <= self.t_bubble_k < self.t_k[i]
            index = np.searchsorted(self.t_k, self.t_dew_k, side='right', sorter=None)
            self.array_n += 1
            self.t_k = np.insert(self.t_k, index, self.t_dew_k)
            self.h_jmol = np.insert(self.h_jmol, index, self.h_dew_jmol)

    def set_optional_t_for_bubble_point(self, n_1ph: int, n_2ph: int, dt: float) -> None:
        t_k = [(lambda i: self.t_bubble_k + dt*i)(i) for i in range(-n_1ph, 0)] + \
              [(lambda i: self.t_bubble_k + dt * i)(i) for i in range(1, n_2ph+1)]
        print(t_k)

    def set_optional_t_for_dew_point(self, n_1ph: int, n_2ph: int, dt: float) -> None:
        t_k = [(lambda i: self.t_dew_k + dt*i)(i) for i in range(-n_2ph, 0)] + \
              [(lambda i: self.t_dew_k + dt * i)(i) for i in range(1, n_1ph+1)]
        print(t_k)

    def get_h_jmol_with_linear_interpolation(self, t_k: float) -> float:
        if t_k < self.t_min_k or t_k > self.t_max_k:
            sys.exit('argument t_k in "get_h_jmol_with_linear_interpolation" is out of acceptable range')
        return np.interp(t_k, self.t_k, self.h_jmol)

    def get_h_jmol_with_refprop10(self, t_k: float) -> float:
        if t_k < self.t_min_k or t_k > self.t_max_k:
            sys.exit('argument t_k in "get_h_jmol_with_refprop10" is out of acceptable range')

        self.fluid.calc_spec_state(t=(t_k, 'k'), p=(self.p_kpa, 'kpa'))
        if self.fluid.error.index > 0:
            self.fluid.error.print_and_terminate()
        else:
            return self.fluid.state.get_data(flag='blk', x_symbol='h', x_units='jmol')

    def get_h_jmol_with_cubic_spline(self, t_k: float) -> float:
        if t_k < self.t_min_k or t_k > self.t_max_k:
            sys.exit('argument t_k in "get_h_jmol_with_cubic_spline" is out of acceptable range')
        return self.cs(t_k)
