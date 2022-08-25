"""
An example of calling the legacy API of REFPROP
By Ian Bell, NIST, 2018, ian.bell@nist.gov
https://github.com/usnistgov/REFPROP-wrappers/blob/master/wrappers/python/ctypes/examples/test_SETREF.py

temperature                     K
pressure, fugacity              kPa
density                         mol/L
composition                     mole fraction
quality                         mole basis (moles vapor/total moles)
enthalpy, internal energy       J/mol
Gibbs, Helmholtz free energy    J/mol
entropy, heat capacity          J/(mol.K)
speed of sound                  m/s
Joule-Thomson coefficient       K/kPa
d(p)/d(rho)                     kPa.L/mol
d2(p)/d(rho)2                   kPa.(L/mol)^2
viscosity                       microPa.s (10^-6 Pa.s)
thermal conductivity            W/(m.K)
dipole moment                   debye
surface tension                 N/m

TRANSPORT PROPERTIES:
   viscosity [microPa.s (10^-6 Pa.s)] and thermal conductivity [W/(m.K)] are calculated for single phase only:
               vis, therm_cond, ierr, herr = r.TRNPRPdll(T,D,z)
   for 2-phase: 0 < q < 1 one should calculate for each phase independently
               vis_liq, therm_cond_liq, ierr, herr = r.TRNPRPdll(T,Dl,x), where x - liq.phase composition
               vis_vap, therm_cond_vap, ierr, herr = r.TRNPRPdll(T,Dv,y), where y - vap.phase composition

   surface tension [N/m] for liquid phase only:
               sigma, ierr, herr = r.SURFTdll(T, Dl, x), where x - liq.phase composition

   Joule-Thomson coefficient [K/kPa]:
               P,e,h,s,Cv,Cp,w,hjt = r.THERMdll(T,D,z)
"""

from __future__ import print_function
import random
import sys
# import timeit

# import cProfile, pstats, io
# from pstats import SortKey

from fluid_class import (RP10Fluid)
from isobar import Isobar, Pressure, Temperature

# input fluid's data:
mixture = RP10Fluid(names=("isobutane", "ethane", "methane"), composition=((0.60, 0.10, 0.30), 'kg/kg'))

# input isobar parameters:
t_min_k = 153.15
t_max_k = 333.15
p_high = Isobar(fluid=mixture,
                p=Pressure(value=20.0, units='bar'),
                t_min=Temperature(value=t_min_k, units='k'),
                t_max=Temperature(value=t_max_k, units='k'),
                dt=10.0)
# sys.exit('ok')

t_random = []
for i in range(1000):
    rand = random.uniform(t_min_k, t_max_k)
    t_random.append(rand)
    # print(i, t_random[i])
# sys.exit()

tau_range = 100
t_range = 100

err_min = float('+inf')
err_max = float('-inf')
err_avr = 0
err_sum = 0

# pr=cProfile.Profile()
# pr.enable()

# for j in range(tau_range):
for i in range(t_range):
    # start = timeit.default_timer()
    _h = p_high.get_h_jmol_with_linear_interpolation(t_k=t_random[i])
    __h = p_high.get_h_jmol_with_refprop10(t_k=t_random[i])
    ___h = p_high.get_h_jmol_with_cubic_spline(t_k=t_random[i])

    err = abs(___h - __h)/__h*100.0

    if err <= err_min:
        err_min = err
    if err >= err_max:
        err_max = err
    err_sum += err

    print(t_random[i],_h)
err_avr = err_sum/t_range
print(err_avr)

#         stop = timeit.default_timer()
#         dtau_interpol = stop-start
#
#         dtau_sum += dtau_interpol
# dtau_avr = dtau_sum/tau_range/t_range

# pr.disable()
# s=io.StringIO()
# sortby=SortKey.CUMULATIVE
# ps=pstats.Stats(pr,stream=s).sort_stats(sortby)
# ps.print_stats()
# print(s.getvalue())

# print('\nВремя расчета avr, сек: ', dtau_avr)  # время расчета см. выше

sys.exit()

# start = timeit.default_timer()
# for i in range(10):
#     p_high.get_h_jmol_with_refprop10(t_k=t_random[i])
# stop = timeit.default_timer()
# dtau_rp10 = stop - start
# print('\nВремя расчета, сек: ', dtau_rp10)  # время расчета см. выше
# print('dtau_rp10/dtau_interpol = ', dtau_rp10/dtau_interpol)
#
# print('ok')
# print(t_random[0])
# print(t_random[0], p_high.get_h_jmol(t_k=t_random[0]))
