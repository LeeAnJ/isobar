указываем явно путь к папке, в кот. хранятся библ. REFPROP.DLL, REFPROP64.DLL и 2 папки:
FLUIDS и MIXTURES (больше и не надо см [1]
os.environ['RPPREFIX'] = r'D:\Cabinet\Ongoing\python_proj\refprop_wrapper\refprop_lib'

команды иниц. dll-библ. refprop и указ.пути к папкам Fluids и Mixtures:
r = ct.REFPROPFunctionLibrary(os.environ['RPPREFIX'], 'dll') 
r.SETPATHdll(os.environ['RPPREFIX'])

теперь надо указать параметры смеси (состав смеси ТОЛЬКО на молярной основе, т.к. во всех
внутр. процедурах расчета св-в используется только лишь мол. состав;
состав смеси - список из 20 эл-тов. поэтому составы всех "лишних" (20-N) позиций в списке 
заполняем нулями):
component_number_max = 20 - это константа для всех смесей рефпроп
component_number = 3
component_name = r'NITROGEN.FLD|OXYGEN.FLD|WATER.FLD'
component_composition_mol = [0.788840469789739,
                             0.209691770450437,
                             0.001467759759824] + [0]*(component_number_max-component_number)






















[1] на этой стр. выложен файл Tutorial.ipynb - подробное описание того, как использовать Python wrapper of NIST REFPROP
https://github.com/usnistgov/REFPROP-wrappers/blob/master/wrappers/python/notebooks/Tutorial.ipynb