"""

"""
import sys
import copy
import converters as conv
from fluid_error_class import FluidRP10Error


def combine_units_props_tag_in_dict(units_internal, units_user, props_symbol, data_tag):
    """
    """
    # combine property's symbol and units in dic: {'t': {'k': 0, 'c': 0}, 'p': {'kpa': 0, 'bar': 0}, ...}
    dic_props = {}      # {'t': {'k': 0, 'c': 0}, 'p': {'kpa': 0, 'bar': 0}, ...}
    for i in range(len(units_internal)):
        # dic_units = {units_internal[i]: 0, units_user[i]: 0}    # {'k': 0, 'c': 0},
        dic_units = {units_internal[i]: None, units_user[i]: None}    # {'k': 0, 'c': 0},
        dic_props[props_symbol[i]] = copy.deepcopy(dic_units)                  # 't': {'k': 0, 'c': 0}

    # combine tags (blk, 'blk', 'blk_liq', 'blk_vap', 'dew', 'bubble') with dict. {'t': {'k': 0, 'c': 0}, ...}
    # data = {'blk': {'t': {'k': 0, 'c': 0}, ...},  'blk_liq': {'t': {'k': 0, 'c': 0}, ...}, ...}
    data = {}
    for item in data_tag:
        data[item] = copy.deepcopy(dic_props)
    return data


class RP10FluidData:
    # properties and units
    props_symbol = ['t', 'p', 'd', 'h', 'e', 's', 'cp', 'cv', 'q', 'w', 'eta', 'tcx', 'mm', 'x']
    units_internal = ['k', 'kpa', 'moll', 'jmol', 'jmol', 'jmolk', 'jmolk', 'jmolk', 'molmol', 'ms', 'upas', 'wmk',
                      'gmol', 'molmol']
    units_user = ['c', 'bar', 'kgm3', 'jkg', 'jkg', 'jkgk', 'jkgk', 'jkgk', 'kgkg', 'ms', 'upas', 'wmk',
                  'gmol', 'kgkg']
    # the following units are used in properties converting from units_internal to units_user
    units_internal_for_converter = ['k', 'kpa', 'moll', 'jmol', 'jmolk', 'molmol', 'ms', 'upas', 'wmk', 'gmol']
    units_user_for_converter = ['c', 'bar', 'kgm3', 'jkg', 'jkgk', 'kgkg']

    data_tag = ['blk', 'blk_liq', 'blk_vap', 'dew', 'bubble']
    property_tag = {'thermal': 0, 'transport': 0}

    data = combine_units_props_tag_in_dict(units_internal, units_user, props_symbol, data_tag)

    def __init__(self, mm_component_gmol, func_dict):
        """
        mm_component_gmol - tuple of mixture's component molar masses
        """
        self.error = FluidRP10Error()   # -> self.error.index; self.error.message

        self.data = copy.deepcopy(RP10FluidData.data)

        # phase= 'l', 'v', 'l_v', depending on q, mol/mol
        self.data_phase_symbol = None

        # data_rp10_routine_status['blk'] = 1/0 1-was called and thermal data are updated; 0 - not called yet
        self.data_rp10_routine_status = {
            'blk': {'thermal': 0, 'transport': 0},
            'blk_liq': {'thermal': 0, 'transport': 0},
            'blk_vap': {'thermal': 0, 'transport': 0},
            'dew': {'thermal': 0, 'transport': 0},
            'bubble': {'thermal': 0, 'transport': 0}
        }
        self.rp10_func = func_dict
        # print(self.rp10_func['therm'])
        # print(self.rp10_func['trnprp'])

        self.mm_component_gmol = mm_component_gmol

    # ------------------------------------------------------------------------------------------------------------------
    # SET FUNCTIONS
    # ------------------------------------------------------------------------------------------------------------------
    def set_data_therm_after_flsh(self, 
                                  flag='blk', 
                                  t_k=None, p_kpa=None, d_moll=None, h_jmol=None, e_jmol=None, s_jmolk=None,
                                  cp_jmolk=None, cv_jmolk=None, q_molmol=None, w_ms=None, mm_gmol=None,
                                  x_molmol=None):

        if flag not in RP10FluidData.data_tag:    # not in ['blk', 'blk_liq', 'blk_vap', 'dew', 'bubble']
            sys.exit('input arg "flag" is out of rang in set_data_therm_after_flsh function of RP10FluidData class')

        self.data[flag]['t']['k'] = t_k
        self.data[flag]['p']['kpa'] = p_kpa
        self.data[flag]['d']['moll'] = d_moll
        self.data[flag]['h']['jmol'] = h_jmol
        self.data[flag]['e']['jmol'] = e_jmol
        self.data[flag]['s']['jmolk'] = s_jmolk
        self.data[flag]['cp']['jmolk'] = cp_jmolk
        self.data[flag]['cv']['jmolk'] = cv_jmolk
        self.data[flag]['q']['molmol'] = q_molmol
        self.data[flag]['w']['ms'] = w_ms
        self.data[flag]['mm']['gmol'] = mm_gmol
        if x_molmol is None:
            self.data[flag]['x']['molmol'] = None
        else:
            self.data[flag]['x']['molmol'] = x_molmol[:]    # array!!! (not list) of 20 elements

        # set phase symbol: 'l', 'v', 'l_v' or None, depending on q, mol/mol
        if flag == 'blk':
            self.data_phase_symbol = set_phase_flag(self.data['blk']['q']['molmol'])

        return

    def set_data_therm_after_therm(self,
                                   flag=None,
                                   h_jmol=None, e_jmol=None, s_jmolk=None,
                                   cp_jmolk=None, cv_jmolk=None,
                                   w_ms=None):

        if flag not in RP10FluidData.data_tag:    # not in ['blk', 'blk_liq', 'blk_vap', 'dew', 'bubble']
            sys.exit('input arg "flag" is out of rang in set_data_therm_after_therm function of RP10FluidData class')

        self.data[flag]['h']['jmol'] = h_jmol
        self.data[flag]['e']['jmol'] = e_jmol
        self.data[flag]['s']['jmolk'] = s_jmolk
        self.data[flag]['cp']['jmolk'] = cp_jmolk
        self.data[flag]['cv']['jmolk'] = cv_jmolk
        self.data[flag]['w']['ms'] = w_ms

        return

    def set_data_wrapper(self,
                         t_k=None, p_kpa=None,
                         d_moll=None, dl_moll=None, dv_moll=None,
                         h_jmol=None, e_jmol=None, s_jmolk=None,
                         cp_jmolk=None, cv_jmolk=None,
                         q_molmol=None, w_ms=None, mm_gmol=None,
                         x_molmol=None, xl_molmol=None, xv_molmol=None,
                         mm_liq_gmol=None, mm_vap_gmol=None):

        # put all separated data in one 'self.state.data' set:
        # blk (rp10 flsh-routines calculate bulk (blk)-set of THERMAL properties (NO TRANSPORT properties)
        self.set_data_therm_after_flsh(flag='blk', t_k=t_k, p_kpa=p_kpa, d_moll=d_moll, h_jmol=h_jmol, e_jmol=e_jmol,
                                       s_jmolk=s_jmolk, cp_jmolk=cp_jmolk, cv_jmolk=cv_jmolk, q_molmol=q_molmol,
                                       w_ms=w_ms,  mm_gmol=mm_gmol, x_molmol=x_molmol)

        self.set_rp10_routine_status(flag='blk', flag_func='thermal', status=1)

        # phase= 'l', 'v', 'l_v', depending on q, mol/mol
        phase_symbol = self.get_phase_symbol()

        # knowing blk-properties and phase state ('l', 'l_v' or 'v') set blk_liq and blk_vap
        if phase_symbol == 'l':
            self.set_data_therm_after_flsh(flag='blk_liq', t_k=t_k, p_kpa=p_kpa, d_moll=d_moll, h_jmol=h_jmol,
                                           e_jmol=e_jmol, s_jmolk=s_jmolk, cp_jmolk=cp_jmolk, cv_jmolk=cv_jmolk,
                                           q_molmol=q_molmol, w_ms=w_ms, mm_gmol=mm_gmol, x_molmol=x_molmol)
            self.set_rp10_routine_status(flag='blk_liq', flag_func='thermal', status=1)

            self.set_data_therm_after_flsh(flag='blk_vap')     # set to None self.state.data['blk_vap'][props][units]
            self.set_rp10_routine_status(flag='blk_vap', flag_func='thermal', status=1)
            self.set_rp10_routine_status(flag='blk_vap', flag_func='transport', status=-1)  # not to be calculated

        elif phase_symbol == 'v':
            self.set_data_therm_after_flsh(flag='blk_liq')     # set to None self.state.data['blk_liq'][props][units]
            self.set_rp10_routine_status(flag='blk_liq', flag_func='thermal', status=1)
            self.set_rp10_routine_status(flag='blk_liq', flag_func='transport', status=-1)  # not to be calculated

            self.set_data_therm_after_flsh(flag='blk_vap', t_k=t_k, p_kpa=p_kpa, d_moll=d_moll, h_jmol=h_jmol,
                                           e_jmol=e_jmol, s_jmolk=s_jmolk, cp_jmolk=cp_jmolk, cv_jmolk=cv_jmolk,
                                           q_molmol=q_molmol, w_ms=w_ms, mm_gmol=mm_gmol, x_molmol=x_molmol)
            self.set_rp10_routine_status(flag='blk_vap', flag_func='thermal', status=1)

        else:   # phase == 'l_v'
            self.set_rp10_routine_status(flag='blk', flag_func='transport', status=-1)  # not to be calculated

            self.set_data_therm_after_flsh(flag='blk_liq', t_k=t_k, p_kpa=p_kpa,
                                           d_moll=dl_moll, q_molmol=q_molmol,
                                           mm_gmol=mm_liq_gmol, x_molmol=xl_molmol)
            self.set_rp10_routine_status(flag='blk_liq', flag_func='thermal', status=0)

            self.set_data_therm_after_flsh(flag='blk_vap', t_k=t_k, p_kpa=p_kpa,
                                           d_moll=dv_moll, q_molmol=q_molmol,
                                           mm_gmol=mm_vap_gmol, x_molmol=xv_molmol)
            self.set_rp10_routine_status(flag='blk_vap', flag_func='thermal', status=0)

    def set_data_to_none(self, flag=None):  # flag in ['dew', 'bubble']
        # props_symbol = ['t', 'p', 'd', 'h', 'e', 's', 'cp', 'cv', 'q', 'w', 'eta', 'tcx', 'mm', 'x']
        # units_internal = ['k', 'kpa', 'moll', 'jmol', 'jmol', 'jmolk', 'jmolk', 'jmolk', 'molmol', 'ms', 'upas', 'wmk',
        #                   'gmol', 'molmol']
        # units_user = ['c', 'bar', 'kgm3', 'jkg', 'jkg', 'jkgk', 'jkgk', 'jkgk', 'kgkg', 'ms', 'upas', 'wmk',
        #               'gmol', 'kgkg']

        if flag not in ['dew', 'bubble', 'blk', 'blk_liq', 'blk_vap']:
            sys.exit('input arg "flag" is out of rang in "set_data_to_none" function of RP10FluidData class')

        n = len(RP10FluidData.props_symbol)
        for i in range(n):
            symbol = RP10FluidData.props_symbol[i]
            units_internal = RP10FluidData.units_internal[i]
            units_user = RP10FluidData.units_user[i]
            self.data[flag][symbol][units_internal] = None
            self.data[flag][symbol][units_user] = None

    def set_data_after_sat_state(self, flag=None, t_k=None, p_kpa=None, d_moll=None,
                                 h_jmol=None, e_jmol=None, s_jmolk=None, cp_jmolk=None,
                                 cv_jmolk=None, w_ms=None, eta_upas=None, tcx_wmk=None,
                                 mm_gmol=None, x_molmol=None):

        if flag not in ['dew', 'bubble']:
            sys.exit('input arg "flag" is out of rang in set_data_after_sat_state function of RP10FluidData class')

        self.data[flag]['t']['k'] = t_k
        self.data[flag]['p']['kpa'] = p_kpa
        self.data[flag]['d']['moll'] = d_moll
        self.data[flag]['h']['jmol'] = h_jmol
        self.data[flag]['e']['jmol'] = e_jmol
        self.data[flag]['s']['jmolk'] = s_jmolk
        self.data[flag]['cp']['jmolk'] = cp_jmolk
        self.data[flag]['cv']['jmolk'] = cv_jmolk

        if flag == "bubble":
            self.data[flag]['q']['molmol'] = 0
        else:
            self.data[flag]['q']['molmol'] = 1

        self.data[flag]['w']['ms'] = w_ms
        self.data[flag]['eta']['upas'] = eta_upas
        self.data[flag]['tcx']['wmk'] = tcx_wmk
        self.data[flag]['mm']['gmol'] = mm_gmol

        if x_molmol is None:
            self.data[flag]['x']['molmol'] = None
        else:
            self.data[flag]['x']['molmol'] = x_molmol[:]    # array!!! (not list) of 20 elements

        return

    def set_rp10_routine_status(self,
                                flag='blk_liq',
                                flag_func='thermal',
                                status=0):
        """

        :param flag: either 'blk_liq' or 'blk_vap'
        :param flag_func: either 'thermal' or 'transport'
        :param status: -1 - not to be called; 0 - not called; 1 - has been called at least once;
        :return:
        """
        self.data_rp10_routine_status[flag][flag_func] = status

    # def set_data_to_zero(self):
    #     # zeroing: self.state.data['blk'/'blk_liq'/'blk_vap'][props][units] = 0
    #     self.set_data_therm_after_flsh(flag='blk')
    #     self.set_data_therm_after_flsh(flag='blk_liq')
    #     self.set_data_therm_after_flsh(flag='blk_vap')

    # ------------------------------------------------------------------------------------------------------------------
    # CALC. PROPERTIES FUNCTIONS
    # ------------------------------------------------------------------------------------------------------------------
    def calc_therm_in_2ph(self, flag):
        p_kpa, e_jmol, h_jmol, s_jmolk, cv_jmolk, cp_jmolk, w_ms, hjt = \
            self.rp10_func['therm'](self.data[flag]['t']['k'],
                                    self.data[flag]['d']['moll'],
                                    self.data[flag]['x']['molmol'])
        # print(flag)
        # print(p_kpa, e_jmol, h_jmol, s_jmolk, cv_jmolk, cp_jmolk, w_ms, hjt)
        self.set_data_therm_after_therm(flag=flag, h_jmol=h_jmol, e_jmol=e_jmol, s_jmolk=s_jmolk, cp_jmolk=cp_jmolk,
                                        cv_jmolk=cv_jmolk, w_ms=w_ms)
        self.set_rp10_routine_status(flag=flag, flag_func='thermal', status=1)
        return

    def calc_trnprp(self, flag):
        # eta, tcx, ierr, herr = TRNPRPdll(T, D, z)
        #  = eta_upas
        # self.data[flag]['tcx']['wmk'] = tcx_wmk
        self.data[flag]['eta']['upas'], self.data[flag]['tcx']['wmk'], ierr, herr = self.rp10_func['trnprp'](
            self.data[flag]['t']['k'],
            self.data[flag]['d']['moll'],
            self.data[flag]['x']['molmol']
        )
        if ierr > 0:
            self.error.index = ierr
            self.error.message = herr
            self.error.set_location('location=self.calc_trnprp(flag) <- self.get_data(flag, x_symbol, x_units)')
            print(self.error.index)
            print(self.error.message)
            sys.exit(self.error.location)
            # return self.error.index, self.error.message, self.error.location

        self.set_rp10_routine_status(flag=flag, flag_func='transport', status=1)
        return

    # ------------------------------------------------------------------------------------------------------------------
    # GET FUNCTIONS
    # ------------------------------------------------------------------------------------------------------------------
    def get_rp10_routine_status(self, flag='blk', flag_props='thermal'):
        # flag E ('blk', 'blk_liq', 'blk_vap')
        # flag_props E ('thermal', 'transport')
        return self.data_rp10_routine_status[flag][flag_props]

    def get_phase_symbol(self):
        return self.data_phase_symbol   # phase symbol: 'l', 'v', 'l_v'

    def get_sat_curve_flag(self):
        if self.data['bubble']['t']['k'] is not None and self.data['dew']['t']['k'] is None:
            return 'bubble'
        elif self.data['bubble']['t']['k'] is None and self.data['dew']['t']['k'] is not None:
            return "dew"
        elif self.data['bubble']['t']['k'] is not None and self.data['dew']['t']['k'] is not None:
            sys.exit("Crit.errror: both 'dew' and 'bubble' curves are calculated. sat_curve_flag is undefined")
        else:
            sys.exit("Crit.errror: you are trying to get sat_curve_flag before sat.state calculations")

    # def get_quality(self, q_units='molmol'):
    #     if q_units == 'molmol':
    #         return self.data['blk']['q']['molmol']
    #     elif q_units == 'kgkg':
    #         return conv.convert_quality_molmol_to_kgkg(q_molmol=self.data['blk']['q']['molmol'],
    #                                                    mm_liq_gmol=self.data['blk_liq']['mm']['gmol'],
    #                                                    mm_vap_gmol=self.data['blk_vap']['mm']['gmol'])

    def get_data(self, flag, x_symbol, x_units):
        # x_symbol = ['t', 'p', 'd', 'h', 'e', 's', 'cp', 'cv', 'q', 'w', 'eta', 'tcx', 'mm', 'x']
        # FluidRP10Properties.units_internal_for_converter = ['k', 'kpa', 'moll', 'jmol', 'jmolk', 'molmol', 'ms',
        #                                                     'upas', 'wmk', 'gmol']
        # FluidRP10Properties.units_user_for_converter = ['c', 'bar', 'kgm3', 'jkg', 'jkgk', 'kgkg']

        x_units = conv.convert_units_string(x_units)

        if flag not in ['blk', 'blk_liq', 'blk_vap', 'dew', 'bubble']:
            print(flag)
            sys.exit('input arg "flag" is out of rang in get_data function of FluidRP10Properties class')

        # for 2ph-state calculations of 'blk_liq' or 'blk_vap' are required for the first time request
        if flag == 'blk_liq' and self.get_rp10_routine_status(flag=flag, flag_props='thermal') == 0:
            self.calc_therm_in_2ph(flag='blk_liq')
        if flag == 'blk_vap' and self.get_rp10_routine_status(flag=flag, flag_props='thermal') == 0:
            self.calc_therm_in_2ph(flag='blk_vap')

        # for transport properties calculations are required for the first time request
        if x_units == 'upas' or x_units == 'wmk':
            if self.get_rp10_routine_status(flag=flag, flag_props='transport') == 0:
                self.calc_trnprp(flag=flag)

        # if x_units are INTERNAL, then just return x_value with no units conversion:
        if x_units in RP10FluidData.units_internal_for_converter:
            return self.data[flag][x_symbol][x_units]

        # if x_units are USER, then conversion from internal to user units is required:
        elif x_units in RP10FluidData.units_user_for_converter:
            # property is already converted
            if self.data[flag][x_symbol][x_units] is not None:
                return self.data[flag][x_symbol][x_units]

            index = RP10FluidData.units_user_for_converter.index(x_units)
            x_units_internal = RP10FluidData.units_internal_for_converter[index]
            x_value = self.data[flag][x_symbol][x_units_internal]

            # possible situation when property['blk_liq'] or property['blk_vap'] is None (single phase case)
            if x_value is None:
                return None

            mm_gmol = self.data[flag]['mm']['gmol']

            if x_units == "c":
                x_value = conv.convert_K_to_oC(t_K=x_value)
            elif x_units == "bar":
                x_value = conv.convert_kpa_to_bar(p_kpa=x_value)
            elif x_units == "kgm3":
                x_value = conv.convert_molL_to_kgm3(d_molL=x_value, mm_gmol=mm_gmol)
            elif x_units == "jkg":
                x_value = conv.convert_Jmol_to_Jkg(x_Jmol=x_value, mm_gmol=mm_gmol)
            elif x_units == "jkgk":
                x_value = conv.convert_JmolK_to_JkgK(x_JmolK=x_value, mm_gmol=mm_gmol)
            elif x_units == "kgkg":
                if x_symbol == "x":    # composition
                    x_value, mm = conv.convert_composition_molmol_to_kgkg(
                        composition_mol=x_value,
                        component_mm_gmol=self.mm_component_gmol
                    )
                else:   # x_symbol == "q"
                    x_value = conv.convert_quality_molmol_to_kgkg(
                        q_molmol=self.data['blk']['q']['molmol'],
                        mm_liq_gmol=self.data['blk_liq']['mm']['gmol'],
                        mm_vap_gmol=self.data['blk_vap']['mm']['gmol']
                    )
            # let's put converted value x_value into corresponding field of the properties dict
            self.data[flag][x_symbol][x_units] = x_value
            return x_value
        else:
            print(x_units)
            sys.exit('input arg "x_units" is out of rang in get_data function of FluidRP10Properties class')


# ---------------------------------------------------------------------------------------------------------------------
# STATIC FUNCTIONS
# ---------------------------------------------------------------------------------------------------------------------
def set_phase_flag(q):  # phase = l, v, l_v
    """
    set phase flag = 'l', 'v', 'l_v', depending on molar quality value
    :param q: molar quality q, [mol/mol'] E [0...1]
    :return: string = 'l' or 'v' or 'l_v'
    """
    if q is None:
        sys.exit("q is None in 'set_phase_flag'; after calc. q should not be None")
    else:
        if q < 0:
            return 'l'
        elif q > 1:
            return 'v'
        else:
            return 'l_v'
