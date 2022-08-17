"""

"""
import sys
import array as arr


# ----------------------------------------------------------------------------------------------------------------------
# temperature: [oC] <-> [K]
# ----------------------------------------------------------------------------------------------------------------------
def convert_oC_to_K(t_oC):
    """
    convert t, [oC] to t, [K]
    :param t_oC: t, [oC]
    :return: t, [K]
    """
    return t_oC + 273.15


def convert_K_to_oC(t_K):
    """
    convert t, [K] to t, [oC]
    :param t_K: t, [K]
    :return: t, [oC]
    """
    if t_K != 0:
        return t_K - 273.15
    else:
        return t_K


# ----------------------------------------------------------------------------------------------------------------------
# pressure: [bar] <-> [kPa]
# ----------------------------------------------------------------------------------------------------------------------
def convert_bar_to_kpa(p_bar):
    """
    convert p, [bar] -> p, [kPa]
    :param p_bar: p, [bar]
    :return: p, [kPa]
    """
    return p_bar*100.0

def convert_kpa_to_bar(p_kpa):
    """
    convert p, [kPa] -> p, [bar]
    :param p_kpa: p, [kPa]
    :return: p, [bar]
    """
    return p_kpa/100.0


# ----------------------------------------------------------------------------------------------------------------------
# specific energy: [J/mol] <-> [J/kg]
# ----------------------------------------------------------------------------------------------------------------------
def convert_Jkg_to_Jmol(x_Jkg, mm_gmol):
    """
    x [J/kg] * m [g/mol] = x [J/kg] * m [kg/(1000*mol)] = x*m/1000 [J/mol]
    :param x_Jkg: argument in [J/kg]
    :param mm_gmol: molar mass, [g/mol]
    :return: x, [J/mol]
    """
    if x_Jkg == 0:
        return 0
    else:
        return x_Jkg * mm_gmol / 1000.0


def convert_Jmol_to_Jkg(x_Jmol, mm_gmol):
    """
    x [J/mol] / m [g/mol] = x [J/mol] / m [kg/(1000*mol)] = x*m/1000 [J/mol]
    :param x_Jmol: argument in [J/mol]
    :param mm_gmol: molar mass, [g/mol]
    :return: x, [J/kg]
    """
    if x_Jmol == 0:
        return 0
    else:
        return x_Jmol / mm_gmol * 1000.0

# ----------------------------------------------------------------------------------------------------------------------
# density: [mol/L] <-> [kg/m3]
# ----------------------------------------------------------------------------------------------------------------------
def convert_molL_to_kgm3(d_molL, mm_gmol):
    """
    x [mol/L] * m [g/mol] = x * m [g/L] = x * m [0.001 kg / 0.001 m3]
    :param d_molL: argument in [mol/L]
    :param mm_gmol: molar mass, [g/mol]
    :return: d, [kg/m3]
    """
    return d_molL*mm_gmol


def convert_kgm3_to_molL(d_kgm3, mm_gmol):
    """
    d [kg/m3] / m [g/mol] = x * m [(1000 g * mol) / (1000 L) * g] = x / m [mol/L]
    :param d_kgm3: argument in [kg/m3]
    :param mm_gmol: molar mass, [g/mol]
    :return: d, [mol/L]
    """
    if d_kgm3 == 0:
        return 0
    else:
        return d_kgm3 / mm_gmol


# ----------------------------------------------------------------------------------------------------------------------
# mass rate: [mol/s] <-> [kg/s]
# ----------------------------------------------------------------------------------------------------------------------
def convert_mols_to_kgs(m_mols, mm_gmol):
    """
    x [mol/s] * mm [g/mol] = x * mm [g/s] = x * m [0.001 kg / s]
    :param m_mols: mass rate in [mol/s]
    :param mm_gmol: molar mass, [g/mol]
    :return: m, [kg/s]
    """
    return m_mols*mm_gmol*0.001


def convert_kgs_to_mols(m_kgs, mm_gmol):
    """
    x [kg/s] / mm [g/mol] = x [1000 g/s] / mm [g/mol] = x / mm [1000 mol/s]
    :param m_kgs: mass rate in [kg/s]
    :param mm_gmol: molar mass, [g/mol]
    :return: m, [mol/s]
    """
    return m_kgs / mm_gmol*1000.0


def convert_mass_rate_to_internal_units(m, m_units, mm_gmol):
    if m_units == 'mols':
        return m
    elif m_units == 'kgs':
        return convert_kgs_to_mols(m, mm_gmol)
    else:
        sys.exit('unsupported units for mass rate in "convert_mass_rate_to_internal_units". program terminated')


def convert_mass_rate_to_user_units(m, m_units, mm_gmol):
    if m_units == 'kgs':
        return m
    elif m_units == 'mols':
        return convert_mols_to_kgs(m, mm_gmol)
    else:
        sys.exit('unsupported units for mass rate in "convert_mass_rate_to_user_units". program terminated')


# ----------------------------------------------------------------------------------------------------------------------
# specific entropy/heat: [J/mol.K] <-> [J/kg.K]
# ----------------------------------------------------------------------------------------------------------------------
def convert_JkgK_to_JmolK(x_JkgK, mm_gmol):
    """
    x [J/kg.K] * m [g/mol] = x [J/kg.K] * m [kg/(1000*mol)] = x*m/1000 [J/mol.K]
    :param x_JkgK: argument in [J/kg.K]
    :param mm_gmol: molar mass, [g/mol]
    :return: x, [J/mol.K]
    """
    return x_JkgK * mm_gmol / 1000.0


def convert_JmolK_to_JkgK(x_JmolK, mm_gmol):
    """
    x [J/mol.K] / m [g/mol] = x [J/mol.K] / m [kg/(1000*mol)] = x*m/1000 [J/mol.K]
    :param x_JmolK: argument in [J/mol.K]
    :param mm_gmol: molar mass, [g/mol]
    :return: x, [J/kg.K]
    """
    if x_JmolK == 0:
        return 0
    else:
        return x_JmolK / mm_gmol * 1000.0


# ----------------------------------------------------------------------------------------------------------------------
# mixture composition: [mol/mol] <-> [kg/kg]
# ----------------------------------------------------------------------------------------------------------------------
def convert_composition_molmol_to_kgkg(composition_mol, component_mm_gmol):
    """ mixture components compositions conversion: composition_mol -> composition_kg
        ref.: Kirillin...Sheindlin Tech.termod., 4th ed., 1983, p.17 eq(1.42), 19 eq(1.55)
        c_kg[i] = c_mol[i]*m_g_mol[i] / sum(c_mol[i]*m_g_mol[i]), i E (1..Ncomp)
        mixture molar mass: mm_g_mol = sum(c_mol[i]*m_g_mol[i])
        units: composition_mol: [mol/mol]; molar_mass: [g/mol] -> composition_kg: [g/g] = [kg/kg]
        :return: composition_kg: array of 20 elements, mm_g_mol
    """
    n = len(component_mm_gmol)
    a = [composition_mol[i] * component_mm_gmol[i] for i in range(n)]
    mm_g_mol = sum(a)
    if mm_g_mol == 0: # in case composition_mol = 0 /blk_liq or blk_vap after data zeroing/
        return composition_mol, mm_g_mol
    else:
        x_kgkg = arr.array('d', [a_i / mm_g_mol for a_i in a] + [0]*(20-n))
        return x_kgkg, mm_g_mol


def convert_composition_kgkg_to_molmol(composition_kg, component_mm_gmol):
    """ mixture components compositions conversion: composition_kg -> composition_mol
        ref.: Kirillin...Sheindlin Tech.termod., 4th ed., 1983, p.17 eq(1.43), p.19 eq(1.54)
        c_mol[i] = c_kg[i]/m_g_mol[i] / sum(c_kg[i]/m_g_mol[i]), i E (1..Ncomp)
        mixture molar mass: mm_g_mol = 1/sum(c_kg[i]*m_g_mol[i])
        units: composition_kg: [kg/kg]; molar_mass: [g/mol] -> composition_mol: [mol/mol]
        :return: composition_mol: tuple of Ncomp figures (in sum = 1), mm_g_mol
    """
    n = len(component_mm_gmol)
    a = [composition_kg[i] / component_mm_gmol[i] for i in range(n)]
    sum_a = sum(a)
    if sum_a == 0:  # in case composition_kg = 0 /blk_liq or blk_vap after data zeroing/
        mm_g_mol = 0
        return composition_kg, mm_g_mol
    else:
        mm_g_mol = 1.0 / sum_a
        x_molmol = arr.array('d', [a_i * mm_g_mol for a_i in a] + [0]*(20-n))
        return x_molmol,mm_g_mol


# ----------------------------------------------------------------------------------------------------------------------
# mixture quality: [mol/mol] <-> [kg/kg]
# ----------------------------------------------------------------------------------------------------------------------
def convert_quality_molmol_to_kgkg(q_molmol, mm_liq_gmol, mm_vap_gmol):
    # """
    # conversion of molar based quality, [mol/mol] into mass based quality, [kg/kg]
    # see RP10\FORTRAN\UTILITY -> subroutine QMASS (line 266)
    # :param q_molmol: [mol/mol];
    # :param mm_liq_gmol:
    # :param mm_vap_gmol:
    # :return: q_kgkg, [kg/kg]
    # """
    if 0 < q_molmol < 1:
        return q_molmol*mm_vap_gmol / ((1-q_molmol)*mm_liq_gmol + q_molmol*mm_vap_gmol)   # q_kgkg, [kg/kg]
    else:
        return q_molmol


def convert_quality_kgkg_to_molmol(q_kgkg, mm_liq_gmol, mm_vap_gmol):
    # """ convertion of mass based quality, [kg/kg] into molar based quality, [mol/mol]
    #     see RP10\FORTRAN\UTILITY -> subroutine QMOLE (line 310)
    #     units:  q_kgkg, [kg/kg]; molar_mass: mm_liq_gmol, mm_vap_gmol, [g/mol];
    #     :return: q_molmol, [mol/mol]
    # """
    if 0 < q_kgkg < 1:
        return q_kgkg/mm_vap_gmol / ((1-q_kgkg)/mm_liq_gmol + q_kgkg/mm_vap_gmol)   # q_molmol, [mol/mol]
    else:
        return q_kgkg


#                  t     p      d       h,e      s        q,x
internal_units = ["k", "kpa", "moll", "jmol", "jmolk", "molmol"]
internal_units_1 = ["k", "kpa"]                 # t, p               - conv(arg)
internal_units_2 = ["moll", "jmol", "jmolk"]        # d, h, e, s, cp, cv - conv(arg, mm_gmol)

user_units = ["c", "bar", "kgm3", "jkg", "jkgk", "kgkg"]
user_units_1 = ["c", "bar"]                   # t, p               - conv(arg)
user_units_2 = ["kgm3", "jkg", "jkgk"]        # d, h, e, s, cp, cv - conv(arg, mm_gmol)

convertors_user_internal_1 = [convert_oC_to_K,          # t          [k]
                              convert_bar_to_kpa]       # p          [bar]
convertors_user_internal_2 = [convert_kgm3_to_molL,     # d          [kgm3]
                              convert_Jkg_to_Jmol,      # h, e       [jkg]
                              convert_JkgK_to_JmolK]    # s, cp, cv  [jkgk]

convertors_internal_user_1 = [convert_K_to_oC,          # t          [k]
                              convert_kpa_to_bar]       # p          [bar]
convertors_internal_user_2 = [convert_molL_to_kgm3,     # d          [kgm3]
                              convert_Jmol_to_Jkg,      # h, e       [jkg]
                              convert_JmolK_to_JkgK]    # s, cp, cv  [jkgk]


# ----------------------------------------------------------------------------------------------------------------------
# wrappers for for RP10 FLSH-routines arguments' units convertors listed above
# ----------------------------------------------------------------------------------------------------------------------
def convert_units_string(units):
    """
    the purpose is to put input arg.: units-string into lower case and remove symbols: /, . from theat string
    :param units: string like  "J/mol.K", "J/molK", "K", "kPa", "J/kg.K", "k"
    :return: string units in lower case and without "/", ".": "jmolk", "k", "kpa", etc.
    """
    _units = units.lower()
    for character in "./":
        _units = _units.replace(character, "")
    return _units


def convert_arg_to_internal_units(x, x_units, mm_gmol=None, mm_liq_gmol=None, mm_vap_gmol=None, component_mm_gmol = None):
    """
    this is a wrapper for a number of functions converting arguments in user units (see list: user_units above)
    into RP10 internal units (see list: internal_units above)
    if argument is already in internal units the function just return x value without any convertions.
    :param x: argument's (t,p,d,h,s,...) value,
    :param x_units: argument's units, typically user's units: [C, bar, Jkg, JkgK, kg/m3, kg/kg]
    :optional param mm_gmol: mixture molar mass required for specific values conversion, like J/kg -> J/mol, etc.
    :optional param  mm_liq_gmol, mm_vap_gmol: liq., vap.-phase molar masses required for quality conversion kg/kg -> mol/mol
    :optional param  component_mm_gmol: mixture components molar masses required for composition conversion kg/kg -> mol/mol
    """

    # convert J/mol.K, J/molK, K, kPa --> jmolk, jmolk, k, kpa, etc.
    units = convert_units_string(x_units)

    if units in internal_units:
        return x
    elif units in user_units_1: # ["c", "bar", ] t, p - conv(x, x_units)
        index = user_units_1.index(units)
        return convertors_user_internal_1[index](x)
    elif units in user_units_2: # d, h, e, s, cp, cv - conv(x, x_units, mm_gmol)
        if mm_gmol is None:
            print("Critical error : no molar mass among input parameters")
            sys.exit("program terminated at  : convert_arg_to_internal_units(J/kg -> J/mol)")
        index = user_units_2.index(units)
        return convertors_user_internal_2[index](x, mm_gmol)
    elif units == 'kgkg':   # q, x
        if type(x) == tuple or type(x) == arr.array:    # x_composition - conv(x, x_units, component_mm_gmol)
            if component_mm_gmol is None:
                print("Critical error : no all components molar masses 'component_mm_gmol' among input parameters")
                sys.exit("program terminated at  : convert_arg_to_internal_units(x, kg/kg -> x, mol/mol)")
            return convert_composition_kgkg_to_molmol(x, component_mm_gmol)
        elif type(x) == float or type(x) == int:    # q - conv(x, x_units, mm_liq_gmol, mm_vap_gmol)
            if mm_liq_gmol is None or mm_vap_gmol is None:
                print("Critical error : no liq or vap molar mass among input parameters")
                sys.exit("program terminated at  : convert_arg_to_internal_units(q, kg/kg -> q, mol/mol)")
            return convert_quality_kgkg_to_molmol(x, mm_liq_gmol, mm_vap_gmol)
        else:
            sys.exit('arg with units of "kg/kg" neither tuple/array nor float/int in convert_arg_to_internal_units ?!!')
    else:
        sys.exit('units are out of range in convert_arg_to_internal_units')


def convert_arg_to_user_units(x, x_units, mm_gmol=None, mm_liq_gmol=None, mm_vap_gmol=None, component_mm_gmol = None):
    # convert J/kg.K, J/kgK, C, bar --> jkgk, jkgk, c, bar, etc.
    units = convert_units_string(x_units)

    if units in internal_units_1: # ["k", "kpa", ] t, p - conv(x, x_units)
        index = internal_units_1.index(units)
        return convertors_internal_user_1[index](x)
    elif units in internal_units_2: # d, h, e, s, cp, cv - conv(x, x_units, mm_gmol)
        if mm_gmol is None:
            print("Critical error : no molar mass among input parameters")
            sys.exit("program terminated at  : convert_arg_to_user_units(J/mol -> J/kg)")
        index = internal_units_2.index(units)
        return convertors_internal_user_2[index](x, mm_gmol)
    elif units == 'molmol':   # q, x
        if type(x) == tuple or type(x) == arr.array:    # x_composition - conv(x, x_units, component_mm_gmol)
            if component_mm_gmol is None:
                print("Critical error : no all components molar masses 'component_mm_gmol' among input parameters")
                sys.exit("program terminated at  : convert_arg_to_user_units(x, mol/mol -> x, kg/kg)")
            return convert_composition_molmol_to_kgkg(x, component_mm_gmol)
        elif type(x) == float or type(x) == int:    # q - conv(x, x_units, mm_liq_gmol, mm_vap_gmol)
            if mm_liq_gmol is None or mm_vap_gmol is None:
                print("Critical error : no liq or vap molar mass among input parameters")
                sys.exit("program terminated at  : convert_arg_to_user_units(q, mol/mol -> q, kg/kg)")
            return convert_quality_molmol_to_kgkg(x, mm_liq_gmol, mm_vap_gmol)
        else:
            sys.exit('arg with units of "mol/mol" neither tuple/array nor float/int in convert_arg_to_user_units ?!!')
    else:
        sys.exit('units are out of range in convert_arg_to_user_units')
