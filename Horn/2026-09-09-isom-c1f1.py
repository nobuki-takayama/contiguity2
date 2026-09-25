from sage.all import *
from chamber_lattice_utils import *
from func_contiguity import *

from pathlib import Path
target_file = Path("get_isom_class.py")
if not target_file.exists():
    raise FileNotFoundError(f"Error: '{target_file}' が存在しない、またはリンク先が切れています。 2026-09-09-rep.py (またはより新しいバージョン) を symbolic link して下さい.")
from get_isom_class import *

load_func_contiguity('2026-08-12-1f1-by-gkz.rr')
contiguity_verbose(0)
R = PolynomialRing(QQ,'x3,dx3,a,c')
dic = R.gens_dict()
a,c = dic['a'], dic['c']
set_constant_part_of_param_vector(0)
data=representatives_of_arrangement([],hg='c1f1',ring=R,dic=dic,param=[a,c])
save(data,'Trash/2026-09-09-c1f1-data.sobj')
isom_class_simplified=list(map(simplify_face,data[0]))

print('Merge isom class')
data_merged=merge_isom_class(data,func='c1f1',ring=R,dic=dic,param=[a,c])
save(data_merged,'Trash/2026-09-09-c1f1-merged-data.sobj')
print('昨日までの 2026-09-08-rep.py のbugを修正したので再度実行. 結果はまったく同じだった.')



