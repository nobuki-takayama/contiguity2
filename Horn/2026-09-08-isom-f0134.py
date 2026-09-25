from sage.all import *
from chamber_lattice_utils import *
from func_contiguity import *

from pathlib import Path
target_file = Path("get_isom_class.py")
if not target_file.exists():
    raise FileNotFoundError(f"Error: '{target_file}' が存在しない、またはリンク先が切れています。 2026-09-09-rep.py (またはより新しいバージョン) を symbolic link して下さい.")
from get_isom_class import *

load_func_contiguity('2026-08-28-f0134.rr')
contiguity_verbose(0)
R2 = PolynomialRing(QQ,'x3,x4,dx3,dx4,b1,b2')
dic2 = R2.gens_dict()
b1,b2 = dic2['b1'], dic2['b2']
set_constant_part_of_param_vector(0)
print('\nNote: A を [[1,0,-2,-3],[0,1,3,4]] としているので, hole は (-1,2)')
print('   つまり [[1,1,1,1],[0,1,3,4]] の param を beta とすれば, b1=beta1-beta2, b2=beta2\n\n')
data=representatives_of_arrangement([],hg='f0134',ring=R2,dic=dic2,param=[b1,b2])
save(data,'Trash/2026-09-08-f0134-data.sobj')
#data=load('Trash/2026-09-08-f0134-data.sobj')
isom_class_simplified=list(map(simplify_face,data[0]))

print('Merge isom class')
data_merged=merge_isom_class(data,func='f0134',ring=R2,dic=dic2,param=[b1,b2])
save(data_merged,'Trash/2026-09-08-f0134-merged-data.sobj')

print('2026-09-09  2026-09-09-rep.py 修正版で再度実行')
