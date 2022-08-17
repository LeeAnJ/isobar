import os
import sys
import array as arr
import copy
import ctREFPROP.ctREFPROP as ct
from prettytable import PrettyTable

from fluid_error_class import FluidRP10Error
from fluid_properties_class import RP10FluidData
# from src.rp10_wrapper import converters as conv
import converters as conv

# max number of components in a mixture of fluids acceptable for RefProp 10
component_number_max = 20

# set the path to the folder with RefProp 10 libs and 2 folders
# os.getcwd() - returns string with the pass to the project folder: "D:\Cabinet\Ongoing\python_proj\refprop_wrapper"
path = os.getcwd() + r"\rp10\refprop_lib"
os.environ['RPPREFIX'] = path
print(os.getcwd())
print(path)
# instantiate the RefProp library REFPROP64.DLL (or REFPROP.DLL)
r = ct.REFPROPFunctionLibrary(os.environ['RPPREFIX'], 'dll')

# tell REFPROP what the root directory is that it should use to get to the folders: FLUIDS and MIXTURES
r.SETPATHdll(os.environ['RPPREFIX'])


class RP10Fluid:
    # class variables:
    fluid_is_active = []  # index=self.personal_index; 1 - activated/ 0 - deactivated
    component_datum_tag = ["TC", "PC", "DC", "TTRP", "PTRP", "MM", "ACF",
                           "DIPOLE", "TNBP", "REOS", "ODP", "GWP", "TMIN",
                           "TMAX", "PMAX", "HEATCOMB", "HFRM"]
    component_datum_key = ["t_cr__k", "p_cr__kpa", "d_cr__moll", "t_trp__k", "p_trp__kpa", "mm__gmol", "acf",
                           "dipole__debye", "t_nbl__k", "reos__jmolk", "odp", "gwp", "t_min__k",
                           "t_max__k", "p_max__kpa", "heatcomb__jmol", "hfrm__jmol"]

    flsh_input_argument_symbol = ["t", "p", "d", "h", "s", "e", "q", "x"]

    def __init__(self, names, composition=None):  # composition_mol=None, composition_kg=None):
        """
        names - tuple of N strings: ("butane",) or ("butane", "ethane", "methane")
        composition = ((0.4, 0.4, 0.2,), 'molmol')) or ((0.4, 0.4, 0.2,), 'kgkg'))
        """
        self.error = FluidRP10Error()  # -> self.error.index; self.error.message

        self.personal_index = None
        self.component_name = None
        self.components_number = None
        self.component_datum = None
        self.component_mm_gmol = None
        self.composition_molmol = None
        self.composition_kgkg = None
        self.mm_g_mol = None

        # set self.component_name and self.components_number
        self.set_components_names_and_number(names)

        # get dict of basic properties for every component:  self.component_datum[i]["prop_name"],
        # i E [0..self.components_number]; "prop_name"=["mm__gmol",... see component_datum_key above
        self.retrieve_component_datum()

        # from self.component_datum retrieve components molar masses as a tuple: self.component_mm_gmol
        self.retrieve_component_mm_gmol()

        # mixture component compositions both mol/mol and kg/kg and mixture molar mass, g/mol
        self.set_component_composition_and_mixture_molar_mass(composition)

        # while calling methods of class specific formats of components names and compositions are required
        self.component_names_as_single_string = self.component_names__to_string()

        # for fluid activation/setup demanded: self.components_number, self.component_names_as_single_string
        # and self.component_composition_mol
        self.activate_fluid_within_rp10_internal_subs()

        # "state" attribute contains dict with a number of thermal and transport parameters
        # to be init (zeroing) before flash subs and fill with data after calculations in flash subs
        # tuple of the components' molar masses [g/mol] - is argument in FluidRP10Data __init__ function
        self.state = RP10FluidData(mm_component_gmol=self.component_mm_gmol,
                                   func_dict={'therm': r.THERMdll, 'trnprp': r.TRNPRPdll})

    # ------------------------------------------------------------------------------------------------------------------
    # FUNCTIONS USED WHILE __INIT__ CALLING
    # ------------------------------------------------------------------------------------------------------------------
    def set_components_names_and_number(self, names):
        if type(names) is not tuple:
            print('mixture input "names" is not a tuple. consider closing COMMA in tuple "names"=("butane",)!!!')
            sys.exit('program terminated in "set_components_names_and_number" <- __init__ RP10Fluid instance')
        self.component_name = names
        self.components_number = len(self.component_name)

    def retrieve_component_datum(self):
        """

        :return:
        """
        z = [1.0] + [0.0] * 19  # pure fluid as a single component mixture
        default = r.GETENUMdll(0, "default").iEnum  # GETENUMdll provide access to the RP10 standard unit systems

        data_list = []  # list of n dict. with each component basic datum
        for name in self.component_name:
            data_dic = {}  # dict. with each component basic datum
            for i, k in enumerate(self.component_datum_tag):
                key = self.component_datum_key[i]
                datum = r.REFPROPdll(name, "", k, default, 0, 0, 0, 0, z)
                data_dic[key] = datum.Output[0]  # add property to dict.
            data_list.append(data_dic)  # add dict. to the list of dict.
        self.component_datum = copy.deepcopy(tuple(data_list))
        return

    def retrieve_component_mm_gmol(self):
        self.component_mm_gmol = tuple([self.component_datum[i]["mm__gmol"] for i in range(self.components_number)])

    def set_pure_fluid_composition_and_molar_mass(self, x):  # x=None or x=((1.0,),'mol/mol' or 'kg/kg')
        """
        check pure fluid input compositions equality to 1.0: sum(x) = 1.0 or None (admissible for pure fluid)
        :param x: tuple of self.components_number figures (in sum = 1) or None
        """
        if x is not None:
            if type(x[0]) is not tuple and type(x[0]) is not arr.array:
                print('mixture input composition is not a tuple. consider closing COMMA in tuple x=(1.0,)!!!')
                sys.exit('program terminated in "set_pure_fluid_composition_and_molar_mass" <- '
                         'set_component_composition_and_mixture_molar_mass <- '
                         '__init__ RP10Fluid instance')
        # pure fluid
        if x is None or x[0][0] == 1.0:
            self.composition_molmol = arr.array('d', [1] + [0] * 19)
            self.composition_kgkg = arr.array('d', [1] + [0] * 19)
            self.mm_g_mol = self.component_datum[0]["mm__gmol"]
        # multi component mixture
        else:
            print('mixture input composition: neither None nor sum(x) != 1 for pure fluid')
            sys.exit('program terminated in "set_pure_fluid_composition_and_molar_mass" <- '
                     'set_component_composition_and_mixture_molar_mass <- '
                     '__init__ RP10Fluid instance')

    def set_mixture_composition_and_molar_mass(self, x):  # x=None or x=((1.0,),'mol/mol' or 'kg/kg')
        check_composition_correctness(x)  # sum(x) != 1 => sys.exit
        x_values_array = arr.array('d', list(x[0]) + [0] * (20 - self.components_number))
        x_units = conv.convert_units_string(x[1])  # input units -> 'molmol' or 'kgkg'

        if x_units == 'molmol':
            self.composition_molmol = x_values_array[:]
            self.composition_kgkg, self.mm_g_mol = conv.convert_composition_molmol_to_kgkg(
                self.composition_molmol,
                component_mm_gmol=tuple([self.component_datum[i]["mm__gmol"] for i in range(self.components_number)])
            )
        elif x_units == 'kgkg':
            self.composition_kgkg = x_values_array[:]
            self.composition_molmol, self.mm_g_mol = conv.convert_arg_to_internal_units(
                x=self.composition_kgkg,
                x_units='kgkg',
                component_mm_gmol=tuple(
                    [self.component_datum[i]["mm__gmol"] for i in range(self.components_number)])
            )
        else:
            print('mixture composition  input units: neither "molmol" nor "kgkg"')
            sys.exit('program terminated in "set_mixture__component_composition_and_molar_mass" <- '
                     '__init__ RP10Fluid instance')

    def set_component_composition_and_mixture_molar_mass(self, x):  # x=((0.6, 0.3, 0.1),'mol/mol' or 'kg/kg') or x=None
        """
        check multi component mixture input compositions equality to 1.0: sum(x) = 1.0
        :param x: tuple of self.components_number figures (in sum = 1)
        """
        if self.components_number == 1:  # case of pure fluid
            self.set_pure_fluid_composition_and_molar_mass(x)
        else:  # case of mixture with n > 1
            self.set_mixture_composition_and_molar_mass(x)
        return

    def component_names__to_string(self):
        suffix = '.FLD'
        prefix = '|'
        names_string = self.component_name[0] + suffix  # name.FLD
        for name in self.component_name[1:]:
            names_string += prefix + name + suffix  # name.FLD|name.FLD|name.FLD
        return names_string.upper()  # NAME.FLD|NAME.FLD|NAME.FLD

    def activate_fluid_within_rp10_internal_subs(self):
        # setup RP 10 for new instance of FluidRP10 class: self.names and self.component_composition_mol
        self.setup_rp10_dlls()  # init RP10 internal subs
        # deactivate previous active fluid (element of FluidRP10.fluid_is_active = 1)
        active_fluid_deactivation()  # find list element = 1 and set its value to 0
        # add current fluid to the list as an activated one (i.e. 1)
        RP10Fluid.fluid_is_active.append(1)  # put current fluid into the list and set it as active: 1
        self.personal_index = len(RP10Fluid.fluid_is_active) - 1

    def setup_rp10_dlls(self):
        """
        calling SETUP (r.SETUPdll(...)) to initialize the RefProp program and set the pure fluid component
         name. it is possible to call SETREF as well to set the reference state in RefProp subs. if ref.state
         is other than 'DEF': 'NBP', 'ASH', 'IIR', ... use  SETREF
        :return: error codes ierr, herr
        """
        ierr, herr = r.SETUPdll(self.components_number, self.component_names_as_single_string, 'HMX.BNC', 'DEF')
        if ierr != 0:
            self.error.index = ierr
            self.error.message = herr
            self.error.location = 'setup_rp10_dlls <- activate/reactivate _fluid_within_rp10_internal_subs'
            self.error.print_and_terminate()
        # if ierr <= 0:
        #     # ierr, herr = r.SETREFdll("DEF", 1, self.component_composition_aggregated_mol, 0, 0, 0, 0)
        # return ierr, herr

    # ------------------------------------------------------------------------------------------------------------------
    # CALC. STATE (FLSH/SAT) - FUNCTIONS
    # ------------------------------------------------------------------------------------------------------------------
    def reactivate_fluid_within_rp10_internal_subs(self):
        # setup RP 10 for current instance of FluidRP10 class: self.names and self.component_composition_mol
        self.setup_rp10_dlls()  # init RP10 internal subs
        active_fluid_deactivation()  # find list element = 1 and set its value to 0
        # set current fluid as an activated one
        RP10Fluid.fluid_is_active[self.personal_index] = 1  # current fluid is set as an active one (= 1 in list)

    def control_input_composition(self, x_value_tuple_or_array, x_units):
        # input composition control: if tuple -> turn to array, if kgkg -> molmol
        x_molmol = None
        mm_gmol = None
        if x_value_tuple_or_array is None:  # in this case units (molmol or kgkg) matters not
            x_molmol = self.composition_molmol[:]  # array.array
            mm_gmol = self.mm_g_mol
        else:
            if type(x_value_tuple_or_array) is tuple:
                x_ = arr.array('d', list(x_value_tuple_or_array) + [0] * (20 - self.components_number))
            else:
                x_ = x_value_tuple_or_array

            if x_units == 'molmol':
                x_molmol = x_[:]  # array.array
                mm_gmol = calc_molar_mass(composition_molmol=x_molmol)
            if x_units == 'kgkg':
                # x_molmol, mm_gmol = conv.convert_composition_kgkg_to_molmol(x_, self.component_mm_gmol)
                x_molmol, mm_gmol = conv.convert_arg_to_internal_units(x=x_,
                                                                       x_units=x_units,
                                                                       component_mm_gmol=self.component_mm_gmol)
        return x_molmol, mm_gmol

    def input_kwargs_unpacking(self, **kwargs):
        argument = []
        input_argument_symbol = []
        for arg_symbol, arg_attributes in kwargs.items():

            arg_symbol = arg_symbol.lower()
            if arg_symbol not in RP10Fluid.flsh_input_argument_symbol:
                sys.exit('argument ' + arg_symbol.lower() + 'in fun.: calc_spec_state is unacceptable. Crit.Error!!!')

            input_argument_symbol.append(arg_symbol)
            _units = conv.convert_units_string(arg_attributes[1])  # input units -> 'molmol' or 'kgkg'
            argument.append({'symbol': arg_symbol.lower(), 'value': arg_attributes[0], 'units': _units})

        if 'x' not in input_argument_symbol:
            x_molmol, mm_gmol = self.control_input_composition(x_value_tuple_or_array=None,
                                                               x_units='molmol')
        else:
            x_index = input_argument_symbol.index('x')
            x_molmol, mm_gmol = self.control_input_composition(x_value_tuple_or_array=argument[x_index]['value'],
                                                               x_units=argument[x_index]['units'])
        return input_argument_symbol, argument, x_molmol, mm_gmol

    # ------------------------------------------------------------------------------------------------------------------
    # FLSH - FUNCTIONS WRAPPER: calc_spec_state(self, **kwargs)
    # ------------------------------------------------------------------------------------------------------------------
    def calc_spec_state(self, **kwargs):

        # check current fluid (RP10Fluid class instance) status of activity
        if RP10Fluid.fluid_is_active[self.personal_index] == 0:
            self.reactivate_fluid_within_rp10_internal_subs()

        # turn input dict-kwargs into lists
        input_argument_symbol, argument, x_molmol, mm_gmol = self.input_kwargs_unpacking(**kwargs)

        # init/zeroing errors and flsh-func. outputdata
        self.error.initiate()  # index = 0, message = 'ok', location = ' '
        ierr = 0
        herr = "ok"
        p = t = d = dl = dv = x = y = q = e = h = s = cv = cp = w = None

        for tag in ['blk', 'blk_liq', 'blk_vap']:
            # set to None: self.state.data['blk'/'blk_liq'/'blk_vap'][props][units] = None
            self.state.set_data_to_none(flag=tag)
            # print(tag)
            # print(self.state.data[tag])
            # flag_func='thermal' or 'transport'
            self.state.set_rp10_routine_status(flag=tag, flag_func='thermal', status=0)  # not calculated yet
            self.state.set_rp10_routine_status(flag=tag, flag_func='transport', status=0)  # not calculated yet
        # sys.exit('fuck')

        # state calculation with RP10 internal FLSH routines
        # f(t,p), f(t,d), f(t,h), f(t,s), f(t,e), f(t,q) ---------------------------------------------------------------
        if "t" in input_argument_symbol:
            index = input_argument_symbol.index('t')
            t = conv.convert_arg_to_internal_units(argument[index]['value'], argument[index]['units'])  # t, [K]
            if "p" in input_argument_symbol:  # tp - arguments + bulk composition, [mol/mol]
                index = input_argument_symbol.index('p')
                p = conv.convert_arg_to_internal_units(argument[index]['value'], argument[index]['units'])  # p, [kPa]
                d, dl, dv, x, y, q, e, h, s, cv, cp, w, ierr, herr = r.TPFLSHdll(t, p, x_molmol)
            elif "d" in input_argument_symbol:  # td - arguments + bulk composition, [mol/mol]
                index = input_argument_symbol.index('d')
                d = conv.convert_arg_to_internal_units(
                    argument[index]['value'],
                    argument[index]['units'],
                    mm_gmol=mm_gmol)  # d, [mol/L]
                p, dl, dv, x, y, q, e, h, s, cv, cp, w, ierr, herr = r.TDFLSHdll(t, d, x_molmol)
            elif "h" in input_argument_symbol:  # th - arguments + bulk composition, [mol/mol]
                index = input_argument_symbol.index('h')
                h = conv.convert_arg_to_internal_units(
                    argument[index]['value'],
                    argument[index]['units'],
                    mm_gmol=mm_gmol)  # h, [J/mol]
                kr = 1
                p, d, dl, dv, x, y, q, e, s, cv, cp, w, ierr, herr = r.THFLSHdll(t, h, x_molmol, kr)
                if ierr > 0:
                    kr = 2
                    p, d, dl, dv, x, y, q, e, s, cv, cp, w, ierr, herr = r.THFLSHdll(t, h, x_molmol, kr)
            elif "s" in input_argument_symbol:  # ts - arguments + bulk composition, [mol/mol]
                index = input_argument_symbol.index('s')
                s = conv.convert_arg_to_internal_units(
                    argument[index]['value'],
                    argument[index]['units'],
                    mm_gmol=mm_gmol)  # s,[J/molK]
                kr = 2
                p, d, dl, dv, x, y, q, e, h, cv, cp, w, ierr, herr = r.TSFLSHdll(t, s, x_molmol, kr)
                if ierr > 0:
                    kr = 1
                    p, d, dl, dv, x, y, q, e, h, cv, cp, w, ierr, herr = r.TSFLSHdll(t, s, x_molmol, kr)
            elif "e" in input_argument_symbol:  # te - arguments + bulk composition, [mol/mol]
                """
                flash calculation given temperature, bulk enthalpy, and bulk composition.
                Often in the liquid, two solutions exist, one of them in the two phase.
                If this is the case, call THFLSH with kr=2 to get the single-phase state.
                    kr - -Flag specifying desired root for multi-valued inputs:
                    kr = 1 - return lower density root
                    kr = 2 - return higher density root
                """
                index = input_argument_symbol.index('e')
                e = conv.convert_arg_to_internal_units(
                    argument[index]['value'],
                    argument[index]['units'],
                    mm_gmol=mm_gmol)  # e, [J/mol]
                kr = 2
                p, d, dl, dv, x, y, q, h, s, cv, cp, w, ierr, herr = r.TEFLSHdll(t, e, x_molmol, kr)
                if ierr > 0:
                    kr = 1
                    p, d, dl, dv, x, y, q, h, s, cv, cp, w, ierr, herr = r.TEFLSHdll(t, e, x_molmol, kr)
            elif "q" in input_argument_symbol:  # tq - arguments + bulk composition, [mol/mol]
                # q, [mol/mol] -> kq = 1, q, [kg/kg] -> kq = 2
                index = input_argument_symbol.index('q')
                q_value = argument[index]['value']
                q_units = conv.convert_units_string(argument[index]['units'])
                if q_units == 'molmol':
                    kq = 1
                elif q_units == 'kgkg':
                    kq = 2
                else:
                    sys.exit('q_units is our of range / neither molmol nor kgkg in input for TQFLSHdll')
                p, d, dl, dv, x, y, e, h, s, cv, cp, w, ierr, herr = r.TQFLSHdll(t, q_value, x_molmol, kq)
                if q_units == 'kgkg':
                    q = conv.convert_arg_to_internal_units(x=q_value,
                                                           x_units=q_units,
                                                           mm_liq_gmol=r.WMOLdll(x),
                                                           mm_vap_gmol=r.WMOLdll(y))
                else:
                    q = q_value

        # f(p,d), f(p,h), f(p,s), f(p,e), f(p,q) -----------------------------------------------------------------------
        elif "p" in input_argument_symbol:
            index = input_argument_symbol.index('p')
            p = conv.convert_arg_to_internal_units(argument[index]['value'], argument[index]['units'])  # p, [kPa]
            if "d" in input_argument_symbol:  # pd - arguments + bulk composition, [mol/mol]
                index = input_argument_symbol.index('d')
                d = conv.convert_arg_to_internal_units(
                    argument[index]['value'],
                    argument[index]['units'],
                    mm_gmol=mm_gmol)  # d, [mol/L]
                t, dl, dv, x, y, q, e, h, s, cv, cp, w, ierr, herr = r.PDFLSHdll(p, d, x_molmol)
            elif "h" in input_argument_symbol:  # ph - arguments + bulk composition, [mol/mol]
                index = input_argument_symbol.index('h')
                h = conv.convert_arg_to_internal_units(
                    argument[index]['value'],
                    argument[index]['units'],
                    mm_gmol=mm_gmol)  # h, [J/mol]
                t, d, dl, dv, x, y, q, e, s, cv, cp, w, ierr, herr = r.PHFLSHdll(p, h, x_molmol)
            elif "s" in input_argument_symbol:  # ps - arguments + bulk composition, [mol/mol]
                index = input_argument_symbol.index('s')
                s = conv.convert_arg_to_internal_units(
                    argument[index]['value'],
                    argument[index]['units'],
                    mm_gmol=mm_gmol)  # s,[J/molK]
                t, d, dl, dv, x, y, q, e, h, cv, cp, w, ierr, herr = r.PSFLSHdll(p, s, x_molmol)
            elif "e" in input_argument_symbol:  # pe - arguments + bulk composition, [mol/mol]
                index = input_argument_symbol.index('e')
                e = conv.convert_arg_to_internal_units(
                    argument[index]['value'],
                    argument[index]['units'],
                    mm_gmol=mm_gmol)  # e - energy, [J/mol]
                t, d, dl, dv, x, y, q, h, s, cv, cp, w, ierr, herr = r.PEFLSHdll(p, e, x_molmol)
            elif "q" in input_argument_symbol:  # pq - arguments + bulk composition, [mol/mol]
                # q, [mol/mol] -> kq = 1, q, [kg/kg] -> kq = 2
                index = input_argument_symbol.index('q')
                q_value = argument[index]['value']
                q_units = conv.convert_units_string(argument[index]['units'])
                if q_units == 'molmol':
                    kq = 1
                elif q_units == 'kgkg':
                    kq = 2
                else:
                    sys.exit('q_units is our of range / neither molmol nor kgkg in input for TQFLSHdll')
                t, d, dl, dv, x, y, e, h, s, cv, cp, w, ierr, herr = r.PQFLSHdll(p, q_value, x_molmol, kq)
                if q_units == 'kgkg':
                    q = conv.convert_arg_to_internal_units(x=q_value,
                                                           x_units=q_units,
                                                           mm_liq_gmol=r.WMOLdll(x),
                                                           mm_vap_gmol=r.WMOLdll(y))
                else:
                    q = q_value

        # f(d,h), f(d,s), f(d,e) ---------------------------------------------------------------------------------------
        elif "d" in input_argument_symbol:
            index = input_argument_symbol.index('d')
            d = conv.convert_arg_to_internal_units(
                argument[index]['value'],
                argument[index]['units'],
                mm_gmol=mm_gmol)  # d, [mol/L]
            if "h" in input_argument_symbol:  # dh - arguments + bulk composition, [mol/mol]
                index = input_argument_symbol.index('h')
                h = conv.convert_arg_to_internal_units(
                    argument[index]['value'],
                    argument[index]['units'],
                    mm_gmol=mm_gmol)  # h, [J/mol]
                t, p, dl, dv, x, y, q, e, s, cv, cp, w, ierr, herr = r.DHFLSHdll(d, h, x_molmol)
            elif "s" in input_argument_symbol:  # ds - arguments + bulk composition, [mol/mol]
                index = input_argument_symbol.index('s')
                s = conv.convert_arg_to_internal_units(
                    argument[index]['value'],
                    argument[index]['units'],
                    mm_gmol=mm_gmol)  # s,[J/molK]
                t, p, dl, dv, x, y, q, e, h, cv, cp, w, ierr, herr = r.DSFLSHdll(d, s, x_molmol)
            elif "e" in input_argument_symbol:  # pe - arguments + bulk composition, [mol/mol]
                index = input_argument_symbol.index('e')
                e = conv.convert_arg_to_internal_units(
                    argument[index]['value'],
                    argument[index]['units'],
                    mm_gmol=mm_gmol)  # e, [J/mol]
                t, p, dl, dv, x, y, q, h, s, cv, cp, w, ierr, herr = r.DEFLSHdll(d, e, x_molmol)

        # f(h,s) -------------------------------------------------------------------------------------------------------
        # [HSFLSH error 260] h-s inputs are two-phase or out of bounds, iterative routine is not available to
        # find a solution. Test calc. failed in 2ph envelope?!!
        # очень странно. когда вход. данные были в user_units: j/kg и j/kg.k тест проходил нормально?!!
        elif "h" in input_argument_symbol:
            index = input_argument_symbol.index('h')
            h = conv.convert_arg_to_internal_units(
                argument[index]['value'],
                argument[index]['units'],
                mm_gmol=mm_gmol)  # h, [J/mol]
            index = input_argument_symbol.index('s')
            s = conv.convert_arg_to_internal_units(
                argument[index]['value'],
                argument[index]['units'],
                mm_gmol=mm_gmol)  # s, [J/mol.K]
            t, p, d, dl, dv, x, y, q, e, cv, cp, w, ierr, herr = r.HSFLSHdll(h, s, x_molmol)
        else:
            sys.exit("critical error in state-func: input_argument_symbol is out of acceptable values: 't','p','d',...")

        if ierr > 0:
            self.error.index = ierr
            self.error.message = herr
            self.error.set_location(location='r.??FLSHdll(...) <- calc_spec_state')
            return  # self.error.index, self.error.message, self.error.location

        # post calculations data/results processing
        xl_molmol = x[:]
        xv_molmol = y[:]
        mm_liq_gmol = r.WMOLdll(xl_molmol)
        mm_vap_gmol = r.WMOLdll(xv_molmol)

        if 0 < q < 1:
            cv, cp, w = None, None, None

        # after spec.state calc. some properties have SPECIFIC (non-physical) VALUES !!!:
        #       subcooled liq.phase:    q = -998 (mol/mol)
        #       2ph.:                   cv, cp, w, eta, tcx = None
        #       super heated vap.phase: q = 998 (mol/mol)

        # in fun. 'self.state.set_data_wrapper' the calculated properties (fun. input arguments)
        # will be put into 'self.state.data['blk']'.
        # while setting up 'self.state.data['blk']' the phase symbol to be retrieved:
        #          self.data_phase_symbol = set_phase_flag(self.data['blk']['q']['molmol'])
        # depending on the value of self.data_phase_symbol: 'l', 'v', 'l_v', the data-sets
        # self.state.data['blk_liq'] or self.state.data['blk_vap'] will be assigned for single phase ('l' or 'v').
        # for 2ph case ('l_v') only t_k, p_kpa, dl_moll, dv_moll, q_molmol, mm_gmol, xl_molmol, xv_molmol
        # will be put into self.state.data['blk_liq'] and self.state.data['blk_vap']. all the rest props. are None
        # and to be calc. while calling 'get_data()' func.
        # transport data: eta, tcx are None in all data-sets and to be calc. while calling 'get_data()' func as well.

        self.state.set_data_wrapper(t_k=t, p_kpa=p, d_moll=d, dl_moll=dl, dv_moll=dv, h_jmol=h, e_jmol=e,
                                    s_jmolk=s, cp_jmolk=cp, cv_jmolk=cv, q_molmol=q, w_ms=w,
                                    mm_gmol=mm_gmol, x_molmol=x_molmol, xl_molmol=xl_molmol,
                                    xv_molmol=xv_molmol, mm_liq_gmol=mm_liq_gmol, mm_vap_gmol=mm_vap_gmol)

        return

    # ------------------------------------------------------------------------------------------------------------------
    # SAT - FUNCTIONS WRAPPER: calc_sat_state(self, sat_curve_flag=None, **kwargs)
    # ------------------------------------------------------------------------------------------------------------------
    def calc_sat_state(self, sat_curve_flag=None, **kwargs):
        # sat_curve_flag - mandatory, t, p, x - optional
        # sat_curve_flag = bubble, l, liq, dew, v, vap
        # kwargs = (t=(300.0, 'K'), p=(100.0, 'kPa'), x=((0.6,0.3,0.1), 'mol/mol'))

        if RP10Fluid.fluid_is_active[self.personal_index] == 0:
            self.reactivate_fluid_within_rp10_internal_subs()

        flag = set_sat_curve_symbol(sat_curve_flag=sat_curve_flag)  # 'l', 'liq', 'bubble' -> 'bubble'...

        # turn input dict-kwargs into lists
        input_argument_symbol, argument, x_molmol, mm_gmol = self.input_kwargs_unpacking(**kwargs)

        # initial zeroing
        self.error.initiate()  # index = 0, message = 'ok', location = ' '
        ierr = 0
        herr = "ok"
        # set to None: self.state.data['blk'/'blk_liq'/'blk_vap'][props][units] = None
        self.state.set_data_to_none(flag=flag)

        # print(self.state.data[flag])
        # sys.exit()
        # self.state.set_data_to_none(flag=flag)  # flag in ['dew', 'bubble']

        # self.state.set_data_after_sat_state(flag=flag)
        t_k, p_kpa, dl, dv, xl, xv = None, None, None, None, None, None

        # bubble (sat.liquid): index=1; dew (sat.vapor): index=2
        sat_line_index = get_sat_line_index(flag=flag)

        if "t" in input_argument_symbol:
            index = input_argument_symbol.index('t')
            t_k = conv.convert_arg_to_internal_units(argument[index]['value'], argument[index]['units'])  # t, [K]
            p_kpa, dl, dv, xl, xv, ierr, herr = r.SATTdll(t_k, x_molmol, sat_line_index)
        elif "p" in input_argument_symbol:
            index = input_argument_symbol.index('p')
            p_kpa = conv.convert_arg_to_internal_units(argument[index]['value'], argument[index]['units'])  # p, [kPa]
            t_k, dl, dv, xl, xv, ierr, herr = r.SATPdll(p_kpa, x_molmol, sat_line_index)

        if ierr > 0:
            self.error.index = ierr
            self.error.message = herr
            self.error.set_location(location='r.SATPdll(...) or r.SATTdll(...)  <- calc_sat_state')
            return

            # input args. for r.THERMdll
        t_k_ = t_k  # keep value of t_k for set_data_after_sat_state

        if sat_line_index == 1:  # sat.liquid
            d_ = dl
            x_ = xl[:]
            mm_gmol = r.WMOLdll(x_)
        else:  # sat_line_index = 2: sat.vapor
            d_ = dv
            x_ = xv[:]
            mm_gmol = r.WMOLdll(x_)

        # call  r.THERMdll(t_k_, d_, x_), r.TRNPRPdll(t_k_, d_, x_)
        p_kpa_, e_jmol, h_jmol, s_jmolk, cv_jmolk, cp_jmolk, w_ms, hjt = r.THERMdll(t_k_, d_, x_)
        eta_upas, tcx_wmk, ierr, herr = r.TRNPRPdll(t_k_, d_, x_)
        if ierr > 0:
            self.error.index = ierr
            self.error.message = herr
            self.error.set_location(location='r.TRNPRPdll(...)  <- calc_sat_state')
            return

        # mark functions usage for given flag = 'dew' or 'bubble'
        self.state.set_rp10_routine_status(flag=flag, flag_func='thermal', status=1)
        self.state.set_rp10_routine_status(flag=flag, flag_func='transport', status=1)

        # update therm and transport data in data-set
        self.state.set_data_after_sat_state(flag=flag, t_k=t_k, p_kpa=p_kpa, d_moll=d_, h_jmol=h_jmol, e_jmol=e_jmol,
                                            s_jmolk=s_jmolk, cp_jmolk=cp_jmolk, cv_jmolk=cv_jmolk, w_ms=w_ms,
                                            eta_upas=eta_upas, tcx_wmk=tcx_wmk, mm_gmol=mm_gmol, x_molmol=x_)
        # after sat.state calc. some properties have SPECIFIC VALUES:
        #   sat.curve bubble/liq: q = 0 (mol/mol,kg/kg)
        #   sat.curve dew/vap: q = 1 (mol/mol,kg/kg)

        return

    # ------------------------------------------------------------------------------------------------------------------
    # PRINT FUNCTIONS
    # ------------------------------------------------------------------------------------------------------------------
    def print_spec_state(self, units_tag='internal'):
        if units_tag not in ['internal', 'user']:
            print('units_tag input argument is neither "internal" nor "user"')
            sys.exit('program was terminated in "print_spec_state" ')

        print(" ")
        if units_tag == 'internal':
            print("                     SPECIFIED STATE POINT /internal units/")
            composition = self.state.get_data(flag='blk', x_symbol='x', x_units='molmol')
            composition_name__units = "composition [mol/mol]"
            _quality = self.state.get_data(flag='blk', x_symbol='q', x_units='molmol')
            if 0 <= _quality <= 1:
                quality = '{:.5f}'.format(_quality)
            elif _quality <= 0:
                quality = 'sub cooled liquid'
            else:
                quality = 'super heated vapour'
            quality_name__units = "quality, [mol/mol]"
        else:
            print("                     SPECIFIED STATE POINT /user units/")
            composition = self.state.get_data(flag='blk', x_symbol='x', x_units='kgkg')
            composition_name__units = "composition [kg/kg]"
            _quality = self.state.get_data(flag='blk', x_symbol='q', x_units='kgkg')
            if 0 <= _quality <= 1:
                quality = '{:.5f}'.format(_quality)
            elif _quality <= 0:
                quality = 'sub cooled liquid'
            else:
                quality = 'super heated vapour'

            quality_name__units = "quality, [kg/kg]"

        table1 = PrettyTable(header=False, padding_width=3)
        list_ = []
        for i in range(self.components_number):
            list_.append('{:.4f}'.format(composition[i]))
        table1.add_rows([
            ["Fluid:", self.component_name],
            [composition_name__units, list_],
            ["molar mass [g/mol]", '{:.4f}'.format(self.state.get_data(flag='blk', x_symbol='mm', x_units='gmol'))]
        ])
        table1.align = "l"
        print(table1)

        table2 = PrettyTable(border=False, header=False)
        table2.header = False
        table2.add_rows([
            ["phase state:", self.state.get_phase_symbol()],
            [quality_name__units, quality]
        ])

        table2.align = "l"
        print(table2)

        table3 = PrettyTable(padding_width=2)
        table3.field_names = ["Property", "Units", "liq", "bulk", "vap"]

        parameter = [
            {'symbol': 't', 'units': None, 'name': 'Temperature'},
            {'symbol': 'p', 'units': None, 'name': 'Pressure'},
            {'symbol': 'd', 'units': None, 'name': 'Density'},
            {'symbol': 'h', 'units': None, 'name': 'Enthalpy'},
            {'symbol': 's', 'units': None, 'name': 'Entropy'},
            {'symbol': 'e', 'units': None, 'name': 'Energy'},
            {'symbol': 'cp', 'units': None, 'name': 'Cp'},
            {'symbol': 'cv', 'units': None, 'name': 'Cv'},
            {'symbol': 'w', 'units': None, 'name': 'w'},
            {'symbol': 'mm', 'units': None, 'name': 'Molar mass'},
            {'symbol': 'eta', 'units': None, 'name': 'Viscosity'},
            {'symbol': 'tcx', 'units': None, 'name': 'Therm.cond.'}
        ]
        if units_tag == 'internal':
            parameter[0]['units'] = 'K'
            parameter[1]['units'] = 'kPa'
            parameter[2]['units'] = 'mol/L'
            parameter[3]['units'] = 'J/mol'
            parameter[4]['units'] = 'J/mol.K'
            parameter[5]['units'] = 'J/mol'
            parameter[6]['units'] = 'J/mol.K'
            parameter[7]['units'] = 'J/mol.K'
            parameter[8]['units'] = 'm/s'
            parameter[9]['units'] = 'g/mol'
            parameter[10]['units'] = 'uPa.s'
            parameter[11]['units'] = 'W/m.K'
        else:
            parameter[0]['units'] = 'C'
            parameter[1]['units'] = 'bar'
            parameter[2]['units'] = 'kg/m3'
            parameter[3]['units'] = 'J/kg'
            parameter[4]['units'] = 'J/kg.K'
            parameter[5]['units'] = 'J/kg'
            parameter[6]['units'] = 'J/kg.K'
            parameter[7]['units'] = 'J/kg.K'
            parameter[8]['units'] = 'm/s'
            parameter[9]['units'] = 'g/mol'
            parameter[10]['units'] = 'uPa.s'
            parameter[11]['units'] = 'W/m.K'

        print_table_list = []
        for _dict in parameter:
            symbol = _dict['symbol']
            units = _dict['units']
            name = _dict['name']
            _list = list()
            _list.append(name)
            _list.append(units)
            for _flag in ['blk_liq', 'blk', 'blk_vap']:
                prop = self.state.get_data(flag=_flag, x_symbol=symbol, x_units=units)
                # prop = '{:.4f}'.format(prop) if prop != 0 else None
                prop = '{:.4f}'.format(prop) if prop is not None else None
                _list.append(prop)
            print_table_list.append(_list)
        table3.add_rows(print_table_list)

        table3.align["Property"] = "l"
        table3.align["Units"] = "r"
        table3.align["liq"] = "r"
        table3.align["bulk"] = "r"
        table3.align["vap"] = "r"

        print(table3)

        if units_tag == 'internal':
            composition_l = self.state.get_data(flag='blk_liq', x_symbol='x', x_units='molmol')
            composition_v = self.state.get_data(flag='blk_vap', x_symbol='x', x_units='molmol')
            list_data_1 = ["liq.phase comp., [mol/mol]"]
            list_data_2 = ["vap.phase comp., [mol/mol]"]
        else:
            composition_l = self.state.get_data(flag='blk_liq', x_symbol='x', x_units='kgkg')
            composition_v = self.state.get_data(flag='blk_vap', x_symbol='x', x_units='kgkg')
            list_data_1 = ["liq.phase comp., [kg/kg]"]
            list_data_2 = ["vap.phase comp., [kg/kg]"]

        table4 = PrettyTable(border=False, header=False)

        _phase_symbol = self.state.get_phase_symbol()
        if _phase_symbol == 'l_v':
            for i in range(self.components_number):
                list_data_1.append('{:.4f}'.format(composition_l[i]))
                list_data_2.append('{:.4f}'.format(composition_v[i]))
        elif _phase_symbol == 'l':
            for i in range(self.components_number):
                list_data_1.append('{:.4f}'.format(composition_l[i]))
                list_data_2.append('none')
        else:
            for i in range(self.components_number):
                list_data_1.append('none')
                list_data_2.append('{:.4f}'.format(composition_v[i]))

        table4.add_rows([list_data_1, list_data_2])
        table1.align = "l"
        print(table4)
        print(" ")

    def print_sat_state(self, sat_curve_symbol, units_tag='internal'):
        if units_tag not in ['internal', 'user']:
            print('units_tag input argument is neither "internal" nor "user"')
            sys.exit('program was terminated in "print_sat_state" ')

        sat_curve_symbol = set_sat_curve_symbol(sat_curve_flag=sat_curve_symbol)  # 'l', 'liq', 'bubble' -> 'bubble'...

        print(" ")
        if units_tag == 'internal':
            print("                     SATURATED STATE POINT /internal units/")
            if sat_curve_symbol == 'dew':
                print("                          dew /saturated vapor/ line")
            else:
                print("                        bubble /saturated liquid/ line")
            composition = self.state.get_data(flag=sat_curve_symbol, x_symbol='x', x_units='molmol')
            composition_name__units = "composition [mol/mol]"
        else:
            print("                     SPECIFIED STATE POINT /user units/")
            if sat_curve_symbol == 'dew':
                print("                          dew /saturated vapor/ line")
            else:
                print("                        bubble /saturated liquid/ line")
            composition = self.state.get_data(flag=sat_curve_symbol, x_symbol='x', x_units='kgkg')
            composition_name__units = "composition [kg/kg]"

        table1 = PrettyTable(header=False, padding_width=3)
        list_ = []
        for i in range(self.components_number):
            list_.append('{:.4f}'.format(composition[i]))
        table1.add_rows([
            ["Fluid:", self.component_name],
            [composition_name__units, list_],
            ["molar mass [g/mol]", '{:.4f}'.format(self.state.get_data(flag=sat_curve_symbol,
                                                                       x_symbol='mm',
                                                                       x_units='gmol'))]
        ])
        table1.align = "l"
        print(table1)

        table3 = PrettyTable(padding_width=2)
        table3.field_names = ["Property", "Units", "State on " + sat_curve_symbol + "-line"]

        parameter = [
            {'symbol': 't', 'units': None, 'name': 'Temperature'},
            {'symbol': 'p', 'units': None, 'name': 'Pressure'},
            {'symbol': 'd', 'units': None, 'name': 'Density'},
            {'symbol': 'h', 'units': None, 'name': 'Enthalpy'},
            {'symbol': 's', 'units': None, 'name': 'Entropy'},
            {'symbol': 'e', 'units': None, 'name': 'Energy'},
            {'symbol': 'cp', 'units': None, 'name': 'Cp'},
            {'symbol': 'cv', 'units': None, 'name': 'Cv'},
            {'symbol': 'w', 'units': None, 'name': 'w'},
            {'symbol': 'mm', 'units': None, 'name': 'Molar mass'},
            {'symbol': 'eta', 'units': None, 'name': 'Viscosity'},
            {'symbol': 'tcx', 'units': None, 'name': 'Therm.cond.'}
        ]
        if units_tag == 'internal':
            parameter[0]['units'] = 'K'
            parameter[1]['units'] = 'kPa'
            parameter[2]['units'] = 'mol/L'
            parameter[3]['units'] = 'J/mol'
            parameter[4]['units'] = 'J/mol.K'
            parameter[5]['units'] = 'J/mol'
            parameter[6]['units'] = 'J/mol.K'
            parameter[7]['units'] = 'J/mol.K'
            parameter[8]['units'] = 'm/s'
            parameter[9]['units'] = 'g/mol'
            parameter[10]['units'] = 'uPa.s'
            parameter[11]['units'] = 'W/m.K'
        else:
            parameter[0]['units'] = 'C'
            parameter[1]['units'] = 'bar'
            parameter[2]['units'] = 'kg/m3'
            parameter[3]['units'] = 'J/kg'
            parameter[4]['units'] = 'J/kg.K'
            parameter[5]['units'] = 'J/kg'
            parameter[6]['units'] = 'J/kg.K'
            parameter[7]['units'] = 'J/kg.K'
            parameter[8]['units'] = 'm/s'
            parameter[9]['units'] = 'g/mol'
            parameter[10]['units'] = 'uPa.s'
            parameter[11]['units'] = 'W/m.K'

        print_table_list = list()
        for _dict in parameter:
            symbol = _dict['symbol']
            units = _dict['units']
            name = _dict['name']
            _list = list()
            _list.append(name)
            _list.append(units)
            prop = self.state.get_data(flag=sat_curve_symbol, x_symbol=symbol, x_units=units)
            prop = '{:.4f}'.format(prop) if prop is not None else None
            _list.append(prop)
            print_table_list.append(_list)

        table3.add_rows(print_table_list)

        table3.align["Property"] = "l"
        table3.align["Units"] = "r"
        table3.align["liq"] = "r"
        table3.align["bulk"] = "r"
        table3.align["vap"] = "r"

        print(table3)

    def convert_dataset_to_user_units(self, flag_data: str):

        if flag_data.lower() in ['l', 'liq', 'bubble']:
            flag = 'bubble'
        elif flag_data.lower() in ['v', 'vap', 'dew']:
            flag = 'dew'
        elif flag_data.lower() in ['blk', 'blk_vap', 'blk_liq']:
            flag = flag_data.lower()
        else:
            sys.exit('flag_data is out of range of admissible values in "convert_dataset_to_user_units"')

        parameter = [
            {'symbol': 't', 'units': 'C'},
            {'symbol': 'p', 'units': 'bar'},
            {'symbol': 'd', 'units': 'kg/m3'},
            {'symbol': 'h', 'units': 'J/kg'},
            {'symbol': 's', 'units': 'J/kg.K'},
            {'symbol': 'e', 'units': 'J/kg'},
            {'symbol': 'cp', 'units': 'J/kg.K'},
            {'symbol': 'cv', 'units': 'J/kg.K'},
            {'symbol': 'w', 'units': 'm/s'},
            {'symbol': 'mm', 'units': 'g/mol'},
            {'symbol': 'eta', 'units': 'uPa.s'},
            {'symbol': 'tcx', 'units': 'W/m.K'}
        ]

        for item in parameter:
            prop = self.state.get_data(flag=flag, x_symbol=item['symbol'], x_units=item['units'])
            prop = prop if prop is not None else None

        if flag in ["bubble", "dew"]:
            # since at sat.curve liq. q = 0 and vap. q = 1 => q_kgkg = q_molmol
            self.state.data[flag]['q']['kgkg'] = self.state.data[flag]['q']['molmol']
        else:   # flag in ['blk', 'blk_vap', 'blk_liq']
            q_kgkg = self.state.get_data(flag=flag, x_symbol='q', x_units='kgkg')

        # converted data will be saved while calling self.state.get_data(...) in func. self.get_data(...)
        x_kgkg = self.state.get_data(flag=flag, x_symbol='x', x_units='kgkg')

    # ------------------------------------------------------------------------------------------------------------------
    # GET FUNCTIONS
    # ------------------------------------------------------------------------------------------------------------------
    def get_name(self):
        prefix = ' / '
        names_string = self.component_name[0]
        for name in self.component_name[1:]:
            names_string += prefix + name   # 'name / name / ... name'
        return names_string.upper()  # 'NAME / NAME / ... NAME'

    def get_composition_kgkg(self):
        prefix = ' / '
        composition_string = str(self.composition_kgkg[0])
        for composition in self.composition_kgkg[1:self.components_number]:
            composition_string += prefix + str(composition)   # 'composition / composition / ... composition'
        return composition_string

    def get_composition_molmol(self):
        prefix = ' / '
        composition_string = str(self.composition_molmol[0])
        for composition in self.composition_molmol[1:self.components_number]:
            composition_string += prefix + str(composition)   # 'composition / composition / ... composition'
        return composition_string

    def get_molar_mass_gmol(self):
        return self.mm_g_mol


# ----------------------------------------------------------------------------------------------------------------------
# STATIC FUNCTIONS
# ----------------------------------------------------------------------------------------------------------------------
def set_sat_curve_symbol(sat_curve_flag: str = '') -> str:
    if sat_curve_flag.lower() in ['bubble', 'l', 'liq']:
        return "bubble"
    elif sat_curve_flag.lower() in ['dew', 'v', 'vap']:
        return "dew"
    else:
        sys.exit("sat_curve_flag is out of range !=['bubble', 'l', 'liq', 'dew', 'v', 'vap']")


def check_composition_correctness(x):
    """
    check whether sum(x) == 1
    :param x: mixture composition  x=((0.6, 0.3, 0.1),'mol/mol' or 'kg/kg') or x=None
    :return:
    """
    if x is None:
        print('mixture input composition is None while n_component > 1')
        sys.exit('program terminated in "check_composition_correctness" <- __init__ RP10Fluid instance')
    else:  # x is not None
        if abs(sum(list(x[0])) - 1.0) > 1.0e-3:
            print(sum(list(x[0])))
            print(abs(sum(list(x[0])) - 1.0))
            print('mixture input composition: sum(x) != 1')
            sys.exit('program terminated in "check_composition_correctness" <- __init__ RP10Fluid instance')


def is_composition_ok(n, comp_mol, comp_kg):
    flag_mol, flag_kg = False, False  # supposedly  there are errors in molar and in mass compositions
    error_message_mol, error_message_kg = "ok", "ok"

    if comp_mol is None:
        error_message_mol = "there is no data in (composition_mol)"
    elif len(comp_mol) != n:
        error_message_mol = "len(composition_mol) != comp.number: crit.error"
    elif abs(1.0 - sum(comp_mol)) > 1.0e-6:
        error_message_mol = "sum(composition_mol) != 1: crit.error"
    else:
        flag_mol = True

    if comp_kg is None:
        error_message_kg = "there is no data in (composition_kg)"
    elif len(comp_kg) != n:
        error_message_kg = "len(composition_kg) != comp.number: crit.error"
    elif abs(1.0 - sum(comp_kg)) > 1.0e-6:
        error_message_kg = "sum(composition_kg) != 1: crit.error"
    else:
        flag_kg = True

    return flag_mol, flag_kg, error_message_mol, error_message_kg


def active_fluid_deactivation():
    # deactivate previous active fluid (element of FluidRP10.fluid_is_active = 1)
    try:
        # find index of the list element whose value is 1 (list - FluidRP10.fluid_is_active)
        index_of_active_fluid = RP10Fluid.fluid_is_active.index(1)  # index of the element = 1
        RP10Fluid.fluid_is_active[index_of_active_fluid] = 0  # set element with index to 0
    except ValueError:
        pass


def calc_molar_mass(composition_molmol=None):
    return r.WMOLdll(composition_molmol)  # 'kg/kmol'


def get_sat_line_index(flag=None):
    if flag == 'bubble':
        return 1  # sat.vapor line
    elif flag == 'dew':
        return 2  # sat.liquid line
    else:
        sys.exit('flag != dew or bubble in get_sat_line_index(flag) <- self.calc_sat_state(flag, p, p_units, ...')
