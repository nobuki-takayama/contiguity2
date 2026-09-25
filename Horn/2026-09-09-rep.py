from sage.all import *
from chamber_lattice_utils import (
    list_all_arrangement_faces,
    find_chamber_bases_and_shifts,
    intersect_face_with_linear_ideal,
    is_bounded,
    strict_interior_face_for_lattice,
    lattice_points_in_face,
    plot_arrangement_faces,
    plot_arrangement_faces_3d,
    get_equality_constraints,
    get_affine_hull_equations,
    get_faces_not_in_affine_hull,
    list_faces_of_arrangement_on_subspace,
    linear_parametrization_of_affine_space,
    reduce_affine_space_expression,
    get_affine_spaces_from_linear_ideal,
    intersection_of_face_and_affine_space,
    show_face,
    is_same_affine_space,
    simplify_face,
    restrict_faces_to_2dim,
    is_face_contained_with_codim1,
    find_nearby_point_in_face,
    segment_to_face,
    is_affine_space_included_in_arrangement,
    get_new_hypersurfaces
)

from func_contiguity import (
    load_func_contiguity,
    contiguity_verbose,
    func_isom,
    parse_contiguity,
    func_isom_sp,
    func_contiguity,
    func_contiguity_b_ideal,
    func_get_chi,
    bf_factor_to_affine_space,
    b_ideal_to_affine_spaces,
    set_constant_part_of_param_vector,
    get_constant_part_of_param_vector,
    extend_param_vector_by_constant_part,
    extention_position,
    linear_space_by_extension_param
)
import sys
DEPTH_MAX=10

## param は [a,c] のように要素が変数名だけであること.
def starting_wall(hg,ring,dic,param):
    Old=param
    h= vector([1 if j == 0 else 0 for j in range(len(Old))])
    New=list(vector(Old)+vector(h))
    print(f'Old={Old}, New={New}')
    Bf=func_isom_sp(Old,New,func=hg,ring=ring,dic=dic)
    print(f'starting wall Bf={Bf}')
    Bf=b_ideal_to_affine_spaces(Old,Bf[0])
    Eq=[]
    for i in range(len(Bf)):
        for j in range(len(Bf[i])):
            Eq.append(Bf[i][j])
        if (len(Bf[i])>1):
            print(f'Non-principal b-ideal: each of {Bf[i]} might be a fake wall')
    print(f'Starting wall: Eq={Eq}')
    return Eq


def representatives_of_arrangement(eq_list,hg,ring,dic,param):
    dim=len(param)
    if eq_list==[]:
        eq_list=starting_wall(hg,ring,dic,param)
    L=linear_space_by_extension_param() # 0 の時は L=[]
    result = representatives0(L,eq_list,hg,ring,dic,param,0)
    return result

# arrangement の dim最大の部分のみ調べる.
# L   affine space
# Eq  Equations for walls, L での arrangement を決める
# hg  hypergeometric system name, e.g., 'c1f1', 'banana1', 'f0134'
#      対応する asir ファイルは load_func_contiguity で load しておくこと
# param は [a,c] などのパラメータ a,c は ring に属する
def representatives0(L,Eq,hg,ring,dic,param,depth):
    if depth > DEPTH_MAX:
        print('Error: to deep depth. sys.exit()')
        sys.exit()
    print(f'representatives0(L={L}, Eq={Eq}, hg={hg},...,depth={depth})')
    if L==[]:
        Tring,S=linear_parametrization_of_affine_space(L,ring=ring,v=param,ambient_dim=len(param))
    else:
        Tring,S=linear_parametrization_of_affine_space(L,ring=ring,v=param)
    print(f'S={S}')
    Isom_class=[]
    Contiguity=[]
    Isom_class_closure=[]
    Arr=list_faces_of_arrangement_on_subspace(Eq,L)
    dim=Arr[0][0]
    print('dim = ',dim)
    if dim < 1: 
        Isom_class=Arr
        Isom_class_closure=Arr
        Contiguity=[[],[],'0-dim']
        return [Isom_class,Contiguity,Isom_class_closure]
    for C in Arr:
        if C[0] < dim: continue  # 最大次元のもののみ調べる.
        print('In arr, C=',C)
        print('Checking the face ', show_face(C))
        if is_bounded(C):
            C_interior=strict_interior_face_for_lattice(C)
            if C_interior == None: continue
            Isom_class.append(C_interior)  # todo, C の内点を列挙.
            Contiguity.append([[],[],'bounded']) # 
            Isom_class_closure.append(C)
            continue
        # 内点のみの face C_interior を求める.
        C_interior=strict_interior_face_for_lattice(C)
        print('C_interior=',C_interior)
        # 内点がなければ何もしない.
        if C_interior == None: continue
        # shifts が Hilbert basis の集合
        bases, shifts=find_chamber_bases_and_shifts(C_interior[1])
        C_contiguity=[]
        for h in shifts:
            print(f'\nHilbert basis h={h} での同型を調べる.')
            Old=S
            New=list(vector(S)+vector(h))
            print(f'Old={Old}, New={New}')
            # Hilbert basis での shift での同型性を調べる. 
            # S は (sub)affine space のパラメータ付
            print(f'func_isom_sp({Old},{New},..), S={S}')
            Bf=func_isom_sp(Old,New,func=hg,ring=ring,dic=dic)
            print(f'Bf={Bf}, S={S}')
            Bf=b_ideal_to_affine_spaces(S,Bf[0])
            C_contiguity.append([func_contiguity(Old,New,func=hg,ring=ring,dic=dic),func_contiguity(New,Old,func=hg,ring=ring,dic=dic),h])
            if len(Bf)==0:  #todo この時の contiguity は? [] が起きるのは?
                print(f'b_ideal is []. It means no isom between {Old} and {New}')
            for W in Bf:
                print(f'W={W} New wall の係数表現')
                if len(W) > 1:
                    print('codim of W > 1 (non-principal b-ideal is found)!')
                Intersection=intersection_of_face_and_affine_space(W,C_interior)
                print(f'Intersection={Intersection}')
                if Intersection != []:
                    # todo fake wall の除外
                    # A が normal なら初期配置を supp hyperplane にすればWWは生じない
                    if not is_affine_space_included_in_arrangement(W,L,Eq):
                        WW=get_new_hypersurfaces(W,L,Eq)
                        print(f'{WW} by get_new_hypersufaces')
                        Eq2=Eq + WW
                        print(f'Eq={Eq} is updated to {Eq2}. restart.\n\n')
                        return representatives0(L,Eq2,hg,ring,dic,param,depth+1)
                    else:
                        continue
        Isom_class.append(C_interior)
        Contiguity.append(C_contiguity)
        Isom_class_closure.append(C)
    print(f'\n\nChecking lower dim={dim-1} faces -----------')
    for C in Arr:
        if C[0] < dim-1: continue  # 一つ低い次元のもの
        if C[0] == dim: continue   # 一番次元の高いもの
        L = get_affine_hull_equations(C)
        isom2, contiguity2,isom2_closure=representatives0(L,Eq,hg,ring,dic,param,0)
        Isom_class=Isom_class+isom2
        Contiguity=Contiguity+contiguity2
        Isom_class_closure=Isom_class_closure+isom2_closure
        
    return [Isom_class,Contiguity,Isom_class_closure]

# data は representatives_of_arrangement の戻り値
# merge の時には set_constant_part_of_param_vector(0) へ.
from sage.all import *

def merge_isom_class(data, func, ring, dic, param):
    """
    同型類 (isom_class) をマージし、隣接するチャンバー同士が同型である場合、
    それらを結ぶ新しい Face (線分) と Contiguity を追加する。
    ※ func_isom_sp の計算結果をキャッシュし、重複実行を防ぐ。
    """
    isom_class, contiguity, isom_class_closure = data
    
    merged_isom_class = [simplify_face(f) for f in isom_class]
    merged_contiguity = list(contiguity)
    
    n = len(isom_class)
    
    print("\n=== Start Merging Isomorphism Classes ===")
    
    # 実行結果を保存するためのキャッシュ辞書
    isom_cache = {}
    
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
                
            if is_face_contained_with_codim1(isom_class_closure[i], isom_class[j]):
                dim_j, ieqs_j = isom_class[j]
                try:
                    bases, shifts = find_chamber_bases_and_shifts(ieqs_j)
                    if not bases:
                        continue
                    p = bases[0]
                except Exception as e:
                    print(f"Warning: Failed to find base point for face {j}: {e}")
                    continue
                    
                q = find_nearby_point_in_face(p, isom_class[i])
                if q is None:
                    continue
                    
                # 4. キャッシュを利用して p と q の同型性を判定する
                # p と q をタプルに変換して辞書のキーにする
                cache_key = (tuple(p), tuple(q))
                
                if cache_key in isom_cache:
                    Bf = isom_cache[cache_key]
                else:
                    Bf = func_isom_sp(p, q, func=func, ring=ring, dic=dic)
                    isom_cache[cache_key] = Bf
                    
                b_ideal = Bf[0]
                
                is_isom = False
                for gen in b_ideal:
                    if gen != 0:
                        is_isom = True
                        break
                        
                if is_isom:
                    seg_face = segment_to_face(p, q)
                    
                    if seg_face not in merged_isom_class:
                        merged_isom_class.append(seg_face)
                        
                        # (注意: 実行時間がかかる場合は func_contiguity も同様にキャッシュ可能です)
                        c_pq = func_contiguity(p, q, func=func, ring=ring, dic=dic)
                        c_qp = func_contiguity(q, p, func=func, ring=ring, dic=dic)
                        h = [q[k] - p[k] for k in range(len(p))]
                        
                        merged_contiguity.append([[c_pq, c_qp, h]])
                        
                        print(f"Merged: Face {j} and Face {i} via p={p} -> q={q}")
                        
    print("=== Merging Completed ===")
    return [merged_isom_class, merged_contiguity]

# data は representatives_of_arrangement の戻り値
# merge の時には set_constant_part_of_param_vector(0) へ.
# new version が安定して動きだしたら不要.
def merge_isom_class_orig1(data, func, ring, dic, param):
    """
    同型類 (isom_class) をマージし、隣接するチャンバー同士が同型である場合、
    それらを結ぶ新しい Face (線分) と Contiguity を追加する。
    
    入力:
        data: representatives_of_arrangement の戻り値 [isom_class, contiguity, isom_class_closure]
        func: 超幾何系の名前 (例: 'c1f1')
        ring: 多項式環
        dic:  変数の辞書
        param: パラメータ変数のリスト
    出力:
        [merged_isom_class, merged_contiguity]
    """
    if get_constant_part_of_param_vector() != 0:
        raise ValueError('Execute set_constant_part_of_param_vector(0)')
    isom_class, contiguity, isom_class_closure = data
    
    # 既存の Face は簡略化して引き継ぐ
    merged_isom_class = [simplify_face(f) for f in isom_class]
    merged_contiguity = list(contiguity)
    
    n = len(isom_class)
    
    print("\n=== Start Merging Isomorphism Classes ===")
    
    # すべての isom_class のペア (i, j) を探索
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            
            # 1. isom_class_closure[i] に codim 1 で isom_class[j] が含まれるか判定
            if is_face_contained_with_codim1(isom_class_closure[i], isom_class[j]):
                
                # 2. isom_class[j] の hilbert basis の起点を p とする
                dim_j, ieqs_j = isom_class[j]
                try:
                    bases, shifts = find_chamber_bases_and_shifts(ieqs_j)
                    if not bases:
                        continue
                    p = bases[0]
                except Exception as e:
                    print(f"Warning: Failed to find base point for face {j}: {e}")
                    continue
                
                # 3. p の近傍で isom_class[i] の点 q を探す
                q = find_nearby_point_in_face(p, isom_class[i])
                if q is None:
                    continue
                
                # 4. p と q が同じ isom class か判定する
                # func_isom_sp を用いて b-ideal を計算
                Bf = func_isom_sp(p, q, func=func, ring=ring, dic=dic)
                b_ideal = Bf[0]
                #print('p,q=',p,q,' b-ideal=',b_ideal) #for debug
                
                # p, q は整数座標のため、b-ideal は定数のリストとなる。
                # 0 以外の定数が含まれていれば (イデアルが 0 でなければ) 同型と判定。
                is_isom = False
                for gen in b_ideal:
                    if gen != 0:
                        is_isom = True
                        break
                    
                if is_isom:
                    # 5. 線分 [p, q] を標準の Face 表現に変換して追加
                    seg_face = segment_to_face(p, q)
                    
                    # 既に同じ Face が追加されていないか確認 (重複防止)
                    if seg_face not in merged_isom_class:
                        merged_isom_class.append(seg_face)
                        
                        # 6. merged_contiguity に [p, q] の contiguity を追加
                        c_pq = func_contiguity(p, q, func=func, ring=ring, dic=dic)
                        c_qp = func_contiguity(q, p, func=func, ring=ring, dic=dic)
                        h = [q[k] - p[k] for k in range(len(p))]
                        
                        # 既存のフォーマットに合わせてリストのリストとして追加
                        merged_contiguity.append([[c_pq, c_qp, h]])
                        
                        print(f"Merged: Face {j} and Face {i} via p={p} -> q={q}")
                        
    print("=== Merging Completed ===")
    return [merged_isom_class, merged_contiguity]
        
"""
load_func_contiguity('2026-08-12-1f1-by-gkz.rr')
R = PolynomialRing(QQ,'x3,dx3,a,c')
dic = R.gens_dict()
a,c = dic['a'], dic['c']
set_constant_part_of_param_vector(0)
#starting_wall(hg='c1f1',ring=R,dic=dic,param=[a,c])
#sys.exit()
# 0+a+0*b=0, 1+a-c=0
isom_class, contiguity, isom_class_closure=representatives_of_arrangement([[0,1,0],[1,1,-1]],hg='c1f1',ring=R,dic=dic,param=[a,c])
plot_arrangement_faces(isom_class,-3,3,-3,3)
#cl=representatives_of_arrangement([],hg='c1f1',ring=R,dic=dic,param=[a,c])
"""

"""
# 2026.09.04--2026.09.05
# Trash/2026-09-05-f0134-log-a.txt に結果を保存.  bounded の bug fix 前. bug あり
load_func_contiguity('2026-08-28-f0134.rr')
contiguity_verbose(0)
R2 = PolynomialRing(QQ,'x3,x4,dx3,dx4,b1,b2')
dic2 = R2.gens_dict()
b1,b2 = dic2['b1'], dic2['b2']
set_constant_part_of_param_vector(0)
#B=func_isom_sp([b1,b2],[b1+1,b2],func='f0134',ring=R2,dic=dic2)
#Bf=b_ideal_to_affine_spaces([b1,b2],B[0])
#sys.exit()
print('\nNote: A を [[1,0,-2,-3],[0,1,3,4]] としているので, hole は (-1,2)')
print('   つまり [[1,1,1,1],[0,1,3,4]] の param を beta とすれば, b1=beta1-beta2, b2=beta2\n\n')
isom_class, contiguity, isom_class_closure=representatives_of_arrangement([],hg='f0134',ring=R2,dic=dic2,param=[b1,b2])
isom_class_simplified=list(map(simplify_face,isom_class))
save(isom_class_simplified,'Trash/2026-09-05-f0134.sobj') # bounded bug fixed
save(contiguity,'Trash/2026-09-05-f0134-contiguity.sobj')
# plot_arrangement_faces(isom_class_simplified,-4,4,-4,4)
"""

"""
# 2026.09.05 パラメータを一部固定した場合のテスト.
load_func_contiguity('2026-08-12-1f1-by-gkz.rr')
R = PolynomialRing(QQ,'x3,dx3,a,c')
dic = R.gens_dict()
a,c = dic['a'], dic['c']
set_constant_part_of_param_vector([1,None])
isom_class, contiguity, isom_class_closure=representatives_of_arrangement([],hg='c1f1',ring=R,dic=dic,param=[c])
plot_arrangement_faces(isom_class,-3,3,-3,3)
isom_class_simplified=list(map(simplify_face,isom_class))
set_constant_part_of_param_vector(0)
"""

"""
# 2026.09.05 banana1 d を fix したテスト.
load_func_contiguity('2026-08-07-banana-one-loop-by-gkz.rr')
contiguity_verbose(0)
R1 = PolynomialRing(QQ,'x4,x5,x6,dx4,dx5,dx6,d,nu1,nu2')
dic1 = R1.gens_dict()
d,nu1,nu2 = dic1['d'], dic1['nu1'], dic1['nu2']
set_constant_part_of_param_vector([4,None,None])
isom_class, contiguity, isom_class_closure=representatives_of_arrangement([],hg='banana1',ring=R1,dic=dic1,param=[nu1,nu2])
isom_class_simplified=list(map(simplify_face,isom_class))
save(isom_class_simplified,'Trash/2026-09-05-banana1-d-4-b.sobj') #保存.
set_constant_part_of_param_vector(0)
F3=restrict_faces_to_2dim(isom_class_simplified,[[4,-1,0,0]],ambient_dim=3)
#show_face(F3)
#plot_arrangement_faces_3d(isom_class_simplified,3,5,-10,10,-10,10)
#plot_arrangement_faces(F3,-1,5,-1,5)
"""

"""
# 2026.09.07 merge_isom_class 用のデータ作成.
load_func_contiguity('2026-08-12-1f1-by-gkz.rr')
R = PolynomialRing(QQ,'x3,dx3,a,c')
dic = R.gens_dict()
a,c = dic['a'], dic['c']
set_constant_part_of_param_vector(0)
data=representatives_of_arrangement([[0,1,0],[1,1,-1]],hg='c1f1',ring=R,dic=dic,param=[a,c])
save(data,'Trash/2026-09-07-data-rep-c1f1.sobj')
"""

"""
# 2026.09.07 merge_isom_class のテスト.
load_func_contiguity('2026-08-12-1f1-by-gkz.rr')
R = PolynomialRing(QQ,'x3,dx3,a,c')
dic = R.gens_dict()
a,c = dic['a'], dic['c']
set_constant_part_of_param_vector(0)
data=load('Trash/2026-09-07-data-rep-c1f1.sobj')
data_merged=merge_isom_class(data,func='c1f1',ring=R,dic=dic,param=[a,c])
save(data_merged,'Trash/2026-09-07-data-merged-rep-c1f1.sobj')
merged_isom_class,merged_contiguity=data_merged
#plot は Notes/Figs/merged_1F1_isom.png
"""
