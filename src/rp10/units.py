from enum import Enum
"""
1. UnitsInternal is an iterable object:
        
        for units in UnitsInternal:
            print units.name, units.value
            
        t, 'k'
        p, 'kpa'
        ...    
2. you can pass an argument into a function as func_name(units: UnitsInternal)
                                                   print(units.name, units.value)
                                                               
    func_name(units=UnitsInternal.t)
    will be printed t, c 

3. Pattern Matching!!!
    func_name(units: UnitsInternal) -> None:
        match units:
            case units.t
                print(units.t.value)
            case units.p
                print(units.p.value)
            case _:
                print('no matches')
                            
4. if inherit from both Enum and str: class UnitsInternal(str, Enum):, then one can use 'value' 
as a string without referring directly to .value attrib.:
    UnitsInternal.t.upper() # C
    UnitsInternal.t == 'c'  # True
            
"""


class UnitsInternal(str, Enum):
    t = 'k'
    p = 'kpa'
    d = 'moll'
    h = 'jmol'
    e = 'jmol'
    s = 'jmolk'
    cp = 'jmolk'
    cv = 'jmolk'
    q = 'molmol'
    w = 'ms'
    eta = 'upas'    # viscosity
    tcx = 'wmk'     # heat conductivity, [W/m.K]
    mm = 'gmol'     # molar mass, [g/mol]
    m = 'mols'      # mass rate, [mol/sec]
    x = 'molmol'    # composition


class UnitsUser(str, Enum):
    t = 'c'
    p = 'bar'
    d = 'kgm3'
    h = 'jkg'
    e = 'jkg'
    s = 'jkgk'
    cp = 'jkgk'
    cv = 'jkgk'
    q = 'kgkg'
    w = 'ms'
    eta = 'upas'    # viscosity
    tcx = 'wmk'     # heat conductivity, [W/m.K]
    mm = 'gmol'     # molar mass, [g/mol]
    m = 'kgs'       # mass rate, [kg/sec]
    x = 'kgkg'      # composition

