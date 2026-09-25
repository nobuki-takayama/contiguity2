from sage.all import *
from chamber_lattice_utils import *
from func_contiguity import *
import sys

from pathlib import Path
target_file = Path("get_isom_class.py")
if not target_file.exists():
    raise FileNotFoundError(f"Error: '{target_file}' が存在しない、またはリンク先が切れています。 2026-09-09-rep.py (またはより新しいバージョン) を symbolic link して下さい.")
from get_isom_class import *   # execute ./link.sh

load_func_contiguity('2026-08-07-banana-one-loop-by-gkz.rr')
contiguity_verbose(0)
R1 = PolynomialRing(QQ,'x4,x5,x6,dx4,dx5,dx6,d,nu1,nu2')
dic1 = R1.gens_dict()
d,nu1,nu2 = dic1['d'], dic1['nu1'], dic1['nu2']
set_constant_part_of_param_vector([4,None,None])
#data=load('Trash/2026-09-08-data-banana1-d-4.sobj')
data=representatives_of_arrangement([],hg='banana1',ring=R1,dic=dic1,param=[nu1,nu2])
#sys.exit()
save(data,'Trash/2026-09-08-data-banana1-d-4.sobj') #保存.
F2=restrict_faces_to_2dim(data[0],[[4,-1,0,0]],ambient_dim=3)

set_constant_part_of_param_vector(0)
print('Merge isom class')
data_merged=merge_isom_class(data,func='banana1',ring=R1,dic=dic1,param=[nu1,nu2])
save(data_merged,'Trash/2026-09-08-banana1-d-4-merged-data.sobj')

data_merged=load('Trash/2026-09-08-banana1-d-4-merged-data.sobj')
F3=restrict_faces_to_2dim(data_merged[0],[[4,-1,0,0]],ambient_dim=3)
show_face(F3)
plot_arrangement_faces(F3,-1,5,-1,5)
#isom_class_simplified=list(map(simplify_face,data[0]))
#plot_arrangement_faces_3d(isom_class_simplified,3,5,-10,10,-10,10)
#set_constant_part_of_param_vector(0)


