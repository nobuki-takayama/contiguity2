from sage.all import *
from chamber_lattice_utils import *
from func_contiguity import *

from pathlib import Path
target_file = Path("get_isom_class.py")
if not target_file.exists():
    raise FileNotFoundError(f"Error: '{target_file}' が存在しない、またはリンク先が切れています。 2026-09-09-rep.py (またはより新しいバージョン) を symbolic link して下さい.")
from get_isom_class import *

load_func_contiguity('2026-09-25-horn-by-gkz.rr')
contiguity_verbose(0)
R2 = PolynomialRing(QQ,'x4,dx4,a1,a2,a3')
dic2 = R2.gens_dict()
a1,a2,a3 = dic2['a1'], dic2['a2'], dic2['a3']
set_constant_part_of_param_vector(0)
set_A([[1,0,0,-1],[0,1,0,1],[0,0,1,1]])  # A matrix
print(show_globals())
data=representatives_of_arrangement([],hg='horn_contiguity2.horn',ring=R2,dic=dic2,param=[a1,a2,a3])
isom_class_simplified=list(map(simplify_face,data[0]))

print('Merge isom class')
data_merged=merge_isom_class(data,func='horn_contiguity2.horn',ring=R2,dic=dic2,param=[a1,a2,a3])

#data, data_berged は Tash/2026-09-25-*.sobj に保存

