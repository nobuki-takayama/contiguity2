"""
2026.08.23
chamber_lattice_utils.py
超平面配置、多面体（チェンバー）の Face 操作、格子点列挙、および図示を行うためのユーティリティ関数群。
SageMath 環境での実行を前提としています。
"""

from sage.geometry.polyhedron.constructor import Polyhedron
from sage.rings.rational_field import QQ
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.geometry.cone import Cone
import sage.plot.colors as colors
from sage.plot.plot3d.shapes2 import point3d
from sage.all import LCM
from sage.matrix.constructor import matrix
from sage.modules.free_module_element import vector
#from sage.symbolic.ring import SR  # <--- シンボリック環を追加. 後ろでやってる.
#
from sage.plot.graphics import Graphics
from sage.plot.point import point
#
import itertools
import math

# =============================================================================
# 1. 超平面配置とFaceの抽出
# =============================================================================

def list_all_arrangement_faces(ieqs_list):
    """
    超平面配置のすべての次元の face を符号ベクトルの探索で抽出する。
    戻り値: [dim, eqs_list] のリスト（次元の降順）
    """
    valid_faces = []
    dim = len(ieqs_list[0]) - 1
    
    def dfs(index, current_ieqs, current_eqns, sigma):
        if current_ieqs or current_eqns:
            try:
                P = Polyhedron(ieqs=current_ieqs, eqns=current_eqns, base_ring=QQ)
                if P.is_empty():
                    return
            except Exception:
                return
        else:
            P = Polyhedron(ambient_dim=dim, base_ring=QQ)
            
        if index == len(ieqs_list):
            pt = P.representative_point()
            for i, sign in enumerate(sigma):
                h = ieqs_list[i]
                val = h[0] + sum(h[j+1] * pt[j] for j in range(dim))
                if sign == '+':
                    if val <= 0: return
                elif sign == '-':
                    if val >= 0: return
                elif sign == '0':
                    if val != 0: return
                    
            eqs_list = []
            for i, sign in enumerate(sigma):
                h = ieqs_list[i]
                if sign == '+':
                    eqs_list.append(h)
                elif sign == '-':
                    eqs_list.append([-x for x in h])
                elif sign == '0':
                    eqs_list.append(h)
                    eqs_list.append([-x for x in h])
                    
            valid_faces.append([P.dim(), eqs_list])
            return
        
        h = ieqs_list[index]
        h_neg = [-x for x in h]
        
        dfs(index + 1, current_ieqs + [h], current_eqns, sigma + ['+'])
        dfs(index + 1, current_ieqs + [h_neg], current_eqns, sigma + ['-'])
        dfs(index + 1, current_ieqs, current_eqns + [h], sigma + ['0'])
        
    dfs(0, [], [], [])
    valid_faces.sort(key=lambda x: x[0], reverse=True)
    return valid_faces

# =============================================================================
# 2. 格子点の網羅（ヒルベルト基底と基点）
# =============================================================================

def find_chamber_bases_and_shifts(ieqs):
    """
    チェンバー内のすべての格子点を網羅するための「基点(Base Points)」と
    「退避錐のヒルベルト基底(Safe Shifts)」を同次化を用いて計算する。
    """
    dim = len(ieqs[0]) - 1
    hom_ieqs = []
    for ieq in ieqs:
        hom_ieqs.append([0, ieq[0]] + list(ieq[1:]))
        
    z_ge_0 = [0, 1] + [0] * dim
    hom_ieqs.append(z_ge_0)
    
    P_hom = Polyhedron(ieqs=hom_ieqs, base_ring=QQ)
    
    rays = []
    for r in P_hom.Vrepresentation():
        if r.is_ray():
            rays.append(r.vector())
        elif r.is_line():
            rays.append(r.vector())
            rays.append(-r.vector())
            
    if rays:
        C_hom = Cone(rays)
        H_hom = C_hom.Hilbert_basis()
    else:
        H_hom = []
        
    base_points = []
    safe_shifts = []
    
    for h in H_hom:
        z = h[0]
        pt = list(h[1:])
        if z == 1:
            base_points.append(pt)
        elif z == 0:
            safe_shifts.append(pt)
            
    base_points.sort()
    safe_shifts.sort()
    return base_points, safe_shifts

# =============================================================================
# 3. 相対的内部と有界性の解析
# =============================================================================

def is_bounded(face):
    """Face が有界 (コンパクト) かどうかを判定する。"""
    dim, ieqs = face
    try:
        P = Polyhedron(ieqs=ieqs, base_ring=QQ)
        return P.is_compact()
    except Exception:
        return False

def get_equality_constraints(face):
    """Face からアフィン包を構成する等式制約を抽出する。"""
    if len(face) == 2 and isinstance(face[1], list) and (len(face[1]) > 0) and isinstance(face[1][0], list):
        ieqs = face[1]
    else:
        ieqs = face[1:]
        
    try:
        P = Polyhedron(ieqs=ieqs, base_ring=QQ)
        eqs_list = []
        for eq in P.equations():
            vec = list(eq)
            lcm = lcm_list([x.denominator() for x in vec if x != 0] + [1])
            vec_int = [int(x * lcm) for x in vec]
            eqs_list.append(vec_int)
        return eqs_list
    except Exception:
        return []

def strict_interior_face_for_lattice(face):
    """
    相対的内部にある格子点のみを含むように Face を縮小する。
    領域が潰れた場合は警告を出力し None を返す。
    """
    dim, ieqs = face
    try:
        P = Polyhedron(ieqs=ieqs, base_ring=QQ)
    except Exception:
        return None
        
    new_ieqs = []
    for eq in P.equations():
        vec = list(eq)
        lcm = lcm_list([x.denominator() for x in vec if x != 0] + [1])
        vec_int = [int(x * lcm) for x in vec]
        g = math.gcd(*[abs(x) for x in vec_int if x != 0]) if any(vec_int) else 1
        if g > 1:
            vec_int = [x // g for x in vec_int]
        new_ieqs.append(vec_int)
        new_ieqs.append([-x for x in vec_int])
        
    for ieq in P.inequalities():
        vec = list(ieq)
        lcm = lcm_list([x.denominator() for x in vec if x != 0] + [1])
        vec_int = [int(x * lcm) for x in vec]
        g = math.gcd(*[abs(x) for x in vec_int if x != 0]) if any(vec_int) else 1
        if g > 1:
            vec_int = [x // g for x in vec_int]
        vec_int[0] -= 1  # 内部へのシフト
        new_ieqs.append(vec_int)
        
    try:
        P_interior = Polyhedron(ieqs=new_ieqs, base_ring=QQ)
        if P_interior.is_empty():
            print(f"Warning: The relative interior of face (dim={dim}) is geometrically empty after shrinking.")
            return None
    except Exception:
        return None
        
    return [dim, new_ieqs]

def lattice_points_in_face(face):
    """
    有界な face の相対的内部にある格子点を列挙する。
    (strict_interior_face_for_lattice と組み合わせて使用可能)
    """
    if not is_bounded(face):
        raise ValueError("Face is not bounded.")
        
    dim, ieqs = face
    P = Polyhedron(ieqs=ieqs, base_ring=QQ)
    all_pts = P.integral_points()
    relative_interior_pts = []
    
    for pt in all_pts:
        is_inside = True
        for ieq in P.inequalities():
            val = ieq[0] + sum(ieq[i+1] * pt[i] for i in range(len(pt)))
            if val == 0:
                is_inside = False
                break
        if is_inside:
            relative_interior_pts.append(list(pt))
            
    return relative_interior_pts

# =============================================================================
# 4. イデアルとFaceの交差
# =============================================================================

def lcm_list(numbers):
    res = 1
    for n in numbers:
        if n != 0:
            res = (res * n) // math.gcd(res, n)
    return res

def extract_linear_coeffs(poly, ring_vars):
    d = poly.dict()
    coeffs = [0] * (len(ring_vars) + 1)
    for exp, c in d.items():
        degree = sum(exp)
        if degree == 0:
            coeffs[0] = c
        elif degree == 1:
            for i, e in enumerate(exp):
                if e == 1:
                    coeffs[i+1] = c
    return coeffs

def intersect_face_with_linear_ideal(face, ideal_gens, ring):
    """face C と ideal I の零点集合 V(I) の共通部分(極大なFaceリスト)を求める。"""
    ring_vars = ring.gens()
    factors_per_gen = []
    for f in ideal_gens:
        factors = [fac[0] for fac in f.factor() if fac[0].degree() > 0]
        factors_per_gen.append(factors)
        
    subspaces_candidates = list(itertools.product(*factors_per_gen))
    base_ieqs = face[1]
    polyhedra_list = []
    
    for combo in subspaces_candidates:
        current_ieqs = list(base_ieqs)
        for L in combo:
            vec = extract_linear_coeffs(L, ring_vars)
            current_ieqs.append(vec)
            current_ieqs.append([-x for x in vec])
            
        try:
            P = Polyhedron(ieqs=current_ieqs, base_ring=QQ)
            if not P.is_empty():
                polyhedra_list.append(P)
        except Exception:
            continue
            
    unique_polyhedra = []
    for P in polyhedra_list:
        if P not in unique_polyhedra:
            unique_polyhedra.append(P)
            
    maximal_polyhedra = []
    for P1 in unique_polyhedra:
        is_maximal = True
        for P2 in unique_polyhedra:
            if P1 == P2: continue
            if P1.intersection(P2) == P1:
                is_maximal = False
                break
        if is_maximal:
            maximal_polyhedra.append(P1)
            
    intersected_faces = []
    for P in maximal_polyhedra:
        eqs_list = []
        for h in P.Hrepresentation():
            vec = list(h)
            lcm = lcm_list([x.denominator() for x in vec if x != 0] + [1])
            vec_int = [int(x * lcm) for x in vec]
            if h.is_equation():
                eqs_list.append(vec_int)
                eqs_list.append([-x for x in vec_int])
            else:
                eqs_list.append(vec_int)
        intersected_faces.append([P.dim(), eqs_list])
        
    intersected_faces.sort(key=lambda x: x[0], reverse=True)
    return intersected_faces

# =============================================================================
# 5. 描画 (2D & 3D)
# =============================================================================

def plot_arrangement_faces(faces, xmin=-3, xmax=3, ymin=-3, ymax=3):
    """2次元の超平面配置の face リストと格子点を図示する。"""
    num_faces = len(faces)
    color_list = colors.rainbow(num_faces)
    G = Graphics()
    
    bbox_ieqs = [
        [-xmin, 1, 0], [xmax, -1, 0],
        [-ymin, 0, 1], [ymax, 0, -1]
    ]
    
    for i, face_data in enumerate(faces):
        original_dim, ieqs = face_data
        try:
            P = Polyhedron(ieqs=ieqs + bbox_ieqs, base_ring=QQ)
            if P.is_empty(): continue
        except Exception:
            continue
            
        c = color_list[i]
        if original_dim == 2:
            G += P.plot(polygon={'color': c, 'alpha': 0.5, 'zorder': 1}, line=False, point=False)
        elif original_dim == 1:
            G += P.plot(line={'color': c, 'thickness': 4, 'zorder': 2}, polygon=False, point=False)
        elif original_dim == 0:
            G += P.plot(point={'color': c, 'size': 60, 'zorder': 3}, polygon=False, line=False)

    lattice_pts = [(x, y) for x in range(int(xmin), int(xmax) + 1) for y in range(int(ymin), int(ymax) + 1)]
    G += point(lattice_pts, color='black', size=15, zorder=4)
    return G

def plot_arrangement_faces_3d(faces, xmin=-2, xmax=2, ymin=-2, ymax=2, zmin=-2, zmax=2):
    """3次元の超平面配置の face リストと格子点を図示する。"""
    num_faces = len(faces)
    color_list = colors.rainbow(num_faces)
    
    bbox_ieqs = [
        [-xmin, 1, 0, 0], [xmax, -1, 0, 0],
        [-ymin, 0, 1, 0], [ymax, 0, -1, 0],
        [-zmin, 0, 0, 1], [zmax, 0, 0, -1]
    ]
    
    plots = []
    for i, face_data in enumerate(faces):
        original_dim, ieqs = face_data
        try:
            P = Polyhedron(ieqs=ieqs + bbox_ieqs, base_ring=QQ)
            if P.is_empty(): continue
        except Exception:
            continue
            
        c = color_list[i]
        if original_dim == 3:
            plots.append(P.plot(polygon={'color': c, 'opacity': 0.15}, line=False, point=False))
        elif original_dim == 2:
            plots.append(P.plot(polygon={'color': c, 'opacity': 0.55}, line=False, point=False))
        elif original_dim == 1:
            plots.append(P.plot(line={'color': c, 'thickness': 5}, polygon=False, point=False))
        elif original_dim == 0:
            plots.append(P.plot(point={'color': c, 'size': 18}, polygon=False, line=False))

    lattice_pts = [(x, y, z) for x in range(int(xmin), int(xmax) + 1) 
                               for y in range(int(ymin), int(ymax) + 1) 
                               for z in range(int(zmin), int(zmax) + 1)]
    plots.append(point3d(lattice_pts, color='black', size=8))
    return sum(plots)

def get_equality_constraints(face):
    if len(face)==0:
        return []
    if face[0]<0:
        return []
    """
    face データから等式制約 (方程式) をすべて抽出する。
    """
    # 入力フォーマットの揺れに対応
    # [dim, [eq1, eq2, ...]] と [dim, eq1, eq2, ...] の両方を許容
    if len(face) == 2 and isinstance(face[1], list) and isinstance(face[1][0], list):
        ieqs = face[1]
    else:
        ieqs = face[1:]
        
    try:
        # 有理数体上で Polyhedron を構成
        P = Polyhedron(ieqs=ieqs, base_ring=QQ)
        
        eqs_list = []
        # equations() はアフィン包を定義する等式(Equationオブジェクト)を返す
        for eq in P.equations():
            # リスト [b, a1, a2, ...] 形式に変換
            vec = list(eq)
            
            # (オプション) 見やすくするために分母を払って整数化する処理
            # SageMathが内部で [0, 1/2, -1/2] のように正規化した場合への対策
            lcm = LCM([x.denominator() for x in vec])
            vec_int = [int(x * lcm) for x in vec]
            
            eqs_list.append(vec_int)
            
        return eqs_list
        
    except Exception as e:
        print(f"Error extracting equations: {e}")
        return []

### 2026-08-30-affine_hull_eq.py より
def get_affine_hull_equations(face):
    """
    与えられた face のアフィン包 (Affine hull) を表現する等式リストを求める。
    戻り値は [c0, c1, ..., cd] のリストのリスト (c0 + c1*x1 + ... = 0)。
    """
    dim, ieqs = face
    try:
        P = Polyhedron(ieqs=ieqs, base_ring=QQ)
        eqs_list = []
        
        # P.equations() はアフィン包を構成する等式オブジェクトを返す
        for eq in P.equations():
            vec = list(eq)
            # 有理数係数の分母を払い、整数化する
            lcm = lcm_list([x.denominator() for x in vec if x != 0] + [1])
            vec_int = [int(x * lcm) for x in vec]
            
            # 全体の最大公約数で割って既約にする
            g = 0
            for x in vec_int:
                g = math.gcd(g, abs(x))
            if g > 1:
                vec_int = [x // g for x in vec_int]
                
            eqs_list.append(vec_int)
            
        return eqs_list
    except Exception as e:
        print(f"Error computing affine hull: {e}")
        return []

def get_faces_not_in_affine_hull(faces, affine_hull_eqs):
    """
    faces のリストのうち、指定された affine_hull_eqs に
    「属さない (完全に含まれていない)」 face のみを列挙する。
    """
    if not affine_hull_eqs:
        # アフィン包の制約が無い(全空間)場合、全てのfaceは空間内に含まれるため
        # 「属さない」face は存在しない
        return []
        
    result = []
    
    for face in faces:
        dim, ieqs = face
        try:
            P = Polyhedron(ieqs=ieqs, base_ring=QQ)
            is_contained = True
            ambient_dim = P.ambient_dim()
            
            for eq in affine_hull_eqs:
                # 頂点 (Vertex) および 半直線/直線 (Ray/Line) をチェック
                for v in P.Vrepresentation():
                    if v.is_vertex():
                        # 頂点の場合: c0 + c1*x1 + c2*x2 + ... = 0 かどうか
                        val = eq[0] + sum(eq[i+1] * v[i] for i in range(ambient_dim))
                    else:
                        # Ray / Line の場合: 方向ベクトルなので定数項 c0 は無視 (同次部分を評価)
                        # c1*d1 + c2*d2 + ... = 0 かどうか
                        val = sum(eq[i+1] * v[i] for i in range(ambient_dim))
                        
                    # 厳密な有理数計算 (QQ) なので、!= 0 の判定が安全に行える
                    if val != 0:
                        is_contained = False
                        break
                        
                if not is_contained:
                    break
                    
            if not is_contained:
                result.append(face)
                
        except Exception as e:
            print(f"Error processing face: {e}")
            
    return result

## 2026-08-30-arrangement_on_subspace.py より
def list_faces_of_arrangement_on_subspace(A, L):
    """
    d次元の超平面配置 A と、アフィン部分空間 L の共通部分 A cap L の
    すべての face を列挙する。
    A, L ともに係数リスト [c0, c1, ..., cd] のリストとして与える。
    出力フォーマットは list_all_arrangement_faces と同じ [dim, eqs_list]。
    """
    valid_faces = []
    
    # 次元を取得 (A か L の空でない方から推定)
    if not A and not L:
        return []
    dim_ambient = len(A[0]) - 1 if A else len(L[0]) - 1
    
    def dfs(index, current_ieqs, current_eqns, sigma):
        # L の方程式は、常に等式として課し続ける
        combined_eqns = L + current_eqns
        
        if current_ieqs or combined_eqns:
            try:
                # QQ上で厳密計算
                P = Polyhedron(ieqs=current_ieqs, eqns=combined_eqns, base_ring=QQ)
                if P.is_empty():
                    return
            except Exception:
                return
        else:
            P = Polyhedron(ambient_dim=dim_ambient, base_ring=QQ)
            
        # 葉 (Leaf) に到達: Aの全ての超平面に対する符号が決定
        if index == len(A):
            pt = P.representative_point()
            # 相対的内部の点で厳密な符号チェック
            for i, sign in enumerate(sigma):
                h = A[i]
                val = h[0] + sum(h[j+1] * pt[j] for j in range(dim_ambient))
                if sign == '+':
                    if val <= 0: return
                elif sign == '-':
                    if val >= 0: return
                elif sign == '0':
                    if val != 0: return
                    
            eqs_list = []
            
            # 1. 出力のために、L の等式 E=0 を E>=0 と -E>=0 のペアとして追加
            for eq in L:
                eqs_list.append(eq)
                eqs_list.append([-x for x in eq])
                
            # 2. A の条件を追加
            for i, sign in enumerate(sigma):
                h = A[i]
                if sign == '+':
                    eqs_list.append(h)
                elif sign == '-':
                    eqs_list.append([-x for x in h])
                elif sign == '0':
                    eqs_list.append(h)
                    eqs_list.append([-x for x in h])
                    
            valid_faces.append([P.dim(), eqs_list])
            return
        
        h = A[index]
        h_neg = [-x for x in h]
        
        # A の超平面に対してのみ分岐を行う
        dfs(index + 1, current_ieqs + [h], current_eqns, sigma + ['+'])
        dfs(index + 1, current_ieqs + [h_neg], current_eqns, sigma + ['-'])
        dfs(index + 1, current_ieqs, current_eqns + [h], sigma + ['0'])
        
    # 深さ優先探索の開始
    dfs(0, [], [], [])
    
    # 次元の大きい順にソート
    valid_faces.sort(key=lambda x: x[0], reverse=True)
    return valid_faces

## 2026-09-03-linear_parametrization.py
def linear_parametrization_of_affine_space(eqs, ambient_dim=None, ring=None, v=None):
    """
    アフィン空間を定義する等式リストから、パラメータ表示 P(t) を生成する。
    
    入力:
        eqs: [[c0, c1, ..., cd], ...] の形式の等式リスト
             (c0 + c1*x1 + ... + cd*xd = 0 を表す)
        ambient_dim: eqs が空リスト(全空間)の場合のみ指定が必要
        ring: (オプション) 戻り値の座標が属する多項式環を指定する
        v: (オプション) パラメータとして使用する ring の元のリスト
        
    出力:
        R: 指定された ring、または生成された多項式環 QQ[t1, ..., tm]
        P: パラメータ表示された座標のリスト [P1(t), ..., Pd(t)]
    """
    if not eqs:
        if ambient_dim is None:
            raise ValueError("等式リストが空の場合、ambient_dim を指定してください。")
        d = ambient_dim
    else:
        d = len(eqs[0]) - 1
        
    A_rows = []
    b_vec = []
    
    for eq in eqs:
        b_vec.append(-eq[0])
        A_rows.append(eq[1:])
        
    A = matrix(QQ, A_rows) if A_rows else matrix(QQ, 0, d)
    b = vector(QQ, b_vec) if b_vec else vector(QQ, 0)
    
    try:
        x0 = A.solve_right(b)
    except ValueError:
        print("Warning: The given affine space is empty (no solutions).")
        return None, None
        
    K = A.right_kernel()
    m = K.dimension()
    basis = K.basis()
    
    # 次元 m に応じてパラメータを設定
    if m > 0:
        if ring is not None and v is not None:
            if len(v) < m:
                raise ValueError(f"指定された変数リスト v の長さ ({len(v)}) が、必要なパラメータ数 ({m}) より不足しています。")
            R = ring
            T = v[:m]  # 変数リストが多めにある場合は先頭から必要な数だけ使用
        else:
            var_names = [f't{i}' for i in range(1, m + 1)]
            R = PolynomialRing(QQ, var_names)
            T = R.gens()
    else:
        # 0次元の場合
        if ring is not None:
            R = ring
        else:
            R = QQ
        T = []
        
    P = []
    for i in range(d):
        expr = R(x0[i])
        for j in range(m):
            expr += T[j] * basis[j][i]
        P.append(expr)
        
    return R, P

from sage.rings.rational_field import QQ
from sage.matrix.constructor import matrix

def reduce_affine_space_expression(E):
    """
    アフィン空間を表現する等式リストから、線形従属（冗長）な等式を削除し、
    独立な等式のみからなる reduced な表現を返す。
    
    入力:
        E: [[c0, c1, ..., cd], ...] の形式の等式リスト
    出力:
        reduced_E: 冗長な式が削除された等式リスト
    """
    if not E:
        return []
        
    # 与えられた方程式リストを QQ 上の行列とする
    M = matrix(QQ, E)
    
    # 転置行列を作成し、ピボット（線形独立な列のインデックス）を取得
    # M^T の独立な列は、元の行列 M の独立な行に対応する
    M_T = M.transpose()
    independent_indices = M_T.pivots()
    
    # 独立と判定されたインデックスの行（元のリストの要素）だけを抽出
    reduced_E = [E[i] for i in independent_indices]
    
    return reduced_E


## 2026-09-04-get_affine.py (.sage)
def get_affine_spaces_from_linear_ideal(ideal_gens, ring):
    """
    一次式の積で生成されるイデアルの生成元リストから、
    零点集合 V(I) を極大なアフィン空間のリストとして抽出する。
    
    入力:
        ideal_gens: イデアルの生成元のリスト
        ring: 多項式環
    出力:
        affine_spaces: アフィン空間を定義する等式リストのリスト
                       [ [[c0, c1, ...], ...], ... ]
    """
    ring_vars = ring.gens()
    factors_per_gen = []
    
    for f in ideal_gens:
        # 定数以外の一次因数を抽出
        factors = [fac[0] for fac in f.factor() if fac[0].degree() > 0]
        if factors:
            factors_per_gen.append(factors)
            
    # 各生成元から1つずつ因子を選ぶ組み合わせの直積
    if factors_per_gen:
        subspaces_candidates = list(itertools.product(*factors_per_gen))
    else:
        subspaces_candidates = [()]
        
    polyhedra_list = []
    
    for combo in subspaces_candidates:
        eqs = []
        for L in combo:
            vec = extract_linear_coeffs(L, ring_vars)
            eqs.append(vec)
            
        try:
            # 方程式系が矛盾していないかチェック
            if eqs:
                P = Polyhedron(eqns=eqs, base_ring=QQ)
            else:
                P = Polyhedron(ambient_dim=len(ring_vars), base_ring=QQ)
                
            if not P.is_empty():
                polyhedra_list.append(P)
        except Exception:
            continue
            
    # 重複の削除と包含関係の整理
    unique_polyhedra = []
    for P in polyhedra_list:
        if P not in unique_polyhedra:
            unique_polyhedra.append(P)
            
    maximal_polyhedra = []
    for P1 in unique_polyhedra:
        is_maximal = True
        for P2 in unique_polyhedra:
            if P1 == P2: continue
            # P1 が P2 に真に含まれる(より小さな部分空間である)場合、P1を破棄
            if P1.intersection(P2) == P1:
                is_maximal = False
                break
        if is_maximal:
            maximal_polyhedra.append(P1)
            
    # 出力フォーマット [[c0, c1, ...], ...] のリストに変換
    affine_spaces = []
    for P in maximal_polyhedra:
        eqs_list = []
        for eq in P.equations():
            vec = list(eq)
            lcm = lcm_list([x.denominator() for x in vec if x != 0] + [1])
            vec_int = [int(x * lcm) for x in vec]
            
            g = 0
            for x in vec_int:
                g = math.gcd(g, abs(x))
            if g > 1:
                vec_int = [x // g for x in vec_int]
                
            eqs_list.append(vec_int)
            
        affine_spaces.append(eqs_list)
        
    return affine_spaces

## 2026-09-04-intersection.py
def intersection_of_face_and_affine_space(L, F):
    """
    アフィン空間 L と Face F の共通部分 (L cap F) を計算する。
    
    入力:
        L: アフィン空間を定義する等式のリスト [[c0, c1, ...], ...]
        F: Face の構造 [dim, [[b, a1, ...], ...]]
    出力:
        共通部分がある場合: [new_dim, eqs_list] (標準のFace形式)
        共通部分がない場合: []
    """
    if not F:
        return []
        
    dim_F, F_ieqs = F
    combined_ieqs = list(F_ieqs)
    
    # アフィン空間 L の等式 (E = 0) を、2つの不等式 (E >= 0, -E >= 0) として追加
    for eq in L:
        combined_ieqs.append(eq)
        combined_ieqs.append([-x for x in eq])
        
    try:
        # 共通部分の多面体を生成
        P = Polyhedron(ieqs=combined_ieqs, base_ring=QQ)
        
        # 共通部分が空集合の場合は空リストを返す
        if P.is_empty():
            return []
            
        eqs_list = []
        
        # 1. 新しい Face の等式制約をペアにして抽出
        for eq in P.equations():
            vec = list(eq)
            # 分母を払う
            lcm = 1
            for x in vec:
                if x != 0:
                    lcm = (lcm * x.denominator()) // math.gcd(lcm, x.denominator())
            vec_int = [int(x * lcm) for x in vec]
            
            # 最大公約数で割って既約にする
            g = 0
            for x in vec_int:
                g = math.gcd(g, abs(x))
            if g > 1:
                vec_int = [x // g for x in vec_int]
                
            eqs_list.append(vec_int)
            eqs_list.append([-x for x in vec_int])
            
        # 2. 新しい Face の不等式制約を抽出
        for ieq in P.inequalities():
            vec = list(ieq)
            lcm = 1
            for x in vec:
                if x != 0:
                    lcm = (lcm * x.denominator()) // math.gcd(lcm, x.denominator())
            vec_int = [int(x * lcm) for x in vec]
            
            g = 0
            for x in vec_int:
                g = math.gcd(g, abs(x))
            if g > 1:
                vec_int = [x // g for x in vec_int]
                
            eqs_list.append(vec_int)
            
        return [P.dim(), eqs_list]
        
    except Exception as e:
        print(f"Error computing intersection: {e}")
        return []

# 2026.09.04
from sage.geometry.polyhedron.constructor import Polyhedron
from sage.rings.rational_field import QQ
from sage.symbolic.ring import SR

def simplify_face(face):
    """
    Face の冗長な不等式・重複する等式を削除し、
    最小限の表現 (Minimal H-representation) に簡略化する。
    """
    if not face:
        return []
        
    dim, ieqs = face
    try:
        # Polyhedron を生成した時点で、内部的に冗長な条件はすべて破棄される
        P = Polyhedron(ieqs=ieqs, base_ring=QQ)
        
        if P.is_empty():
            return []
            
        new_ieqs = []
        
        # 1. 最小化された等式制約を抽出
        for eq in P.equations():
            vec = list(eq)
            lcm = 1
            for x in vec:
                if x != 0:
                    lcm = (lcm * x.denominator()) // math.gcd(lcm, x.denominator())
            vec_int = [int(x * lcm) for x in vec]
            
            g = 0
            for x in vec_int:
                g = math.gcd(g, abs(x))
            if g > 1:
                vec_int = [x // g for x in vec_int]
                
            new_ieqs.append(vec_int)
            new_ieqs.append([-x for x in vec_int])
            
        # 2. 最小化された不等式制約 (ファセット) を抽出
        for ieq in P.inequalities():
            vec = list(ieq)
            lcm = 1
            for x in vec:
                if x != 0:
                    lcm = (lcm * x.denominator()) // math.gcd(lcm, x.denominator())
            vec_int = [int(x * lcm) for x in vec]
            
            g = 0
            for x in vec_int:
                g = math.gcd(g, abs(x))
            if g > 1:
                vec_int = [x // g for x in vec_int]
                
            new_ieqs.append(vec_int)
            
        return [P.dim(), new_ieqs]
        
    except Exception as e:
        print(f"Error simplifying face: {e}")
        return face


def show_face(F, ring=None):
    """
    Face を簡略化した上で、数式のリストとして表示する。
    複数のFaceのリストが与えられた場合は、各Faceに対して再帰的に処理を行う。
    ring が指定されていない場合は、内部で自動的に変数 (x0, x1, ...) を持つ多項式環を定義する。
    """
    # 1. リストのリスト (複数のFace) が与えられた場合の自動マップ処理
    if isinstance(F, list) and len(F) > 0 and isinstance(F[0], list):
        return [show_face(f, ring) for f in F]
        
    # 2. 不正なフォーマットの場合はそのまま返す
    if not (isinstance(F, list) and len(F) == 2 and isinstance(F[0], int)):
        return F
        
    simplified_face = simplify_face(F)
    
    if not simplified_face:
        print("Empty Face")
        return []
        
    dim, ieqs = simplified_face
    
    if not ieqs:
        return []
        
    ambient_dim = len(ieqs[0]) - 1
    
    if ring is None:
        if ambient_dim > 0:
            var_names = [f'x{i}' for i in range(ambient_dim)]
            ring = PolynomialRing(QQ, var_names)
        else:
            ring = QQ
            
    if ring == QQ:
        vars = []
    else:
        vars = ring.gens()
        if len(vars) < ambient_dim:
            raise ValueError(f"指定された環の変数 ({len(vars)}個) が、空間の次元 ({ambient_dim}次元) より不足しています。")
            
    P = Polyhedron(ieqs=ieqs, base_ring=QQ)
    expressions = []
    
    # 等式の表示 (E == 0)
    for eq in P.equations():
        vec = list(eq)
        expr = vec[0] + sum(vec[i+1] * vars[i] for i in range(ambient_dim))
        expressions.append(SR(expr) == 0)  # SR にキャストして方程式オブジェクト化
        
    # 不等式の表示 (I >= 0)
    for ieq in P.inequalities():
        vec = list(ieq)
        expr = vec[0] + sum(vec[i+1] * vars[i] for i in range(ambient_dim))
        expressions.append(SR(expr) >= 0)  # SR にキャストして不等式オブジェクト化
        
    return expressions


def is_same_affine_space(eq_list1, eq_list2, ambient_dim=None):
    """
    二つのアフィン空間が幾何学的に同一であるかを判定する。
    
    入力:
        eq_list1, eq_list2: アフィン空間を定義する等式のリスト [[c0, c1, ...], ...]
        ambient_dim: 両方のリストが空の場合のみ必要
    出力:
        同一の空間 (または共に空集合) であれば True、異なれば False
    """
    # 空間の次元 d を特定
    if eq_list1:
        d = len(eq_list1[0]) - 1
    elif eq_list2:
        d = len(eq_list2[0]) - 1
    else:
        if ambient_dim is None:
            # 両方空リストで次元指定もない場合は、自明に同一(全空間)とする
            return True
        d = ambient_dim
        
    def get_space_components(eqs):
        """方程式系から特解と零空間を取得する。解なし(空集合)の場合は None を返す。"""
        if not eqs:
            # 全空間の場合
            A = matrix(QQ, 0, d)
            b = vector(QQ, 0)
        else:
            A = matrix(QQ, [eq[1:] for eq in eqs])
            b = vector(QQ, [-eq[0] for eq in eqs])
            
        try:
            x0 = A.solve_right(b)
        except ValueError:
            # 解が存在しない (空集合)
            return None
            
        K = A.right_kernel()
        return K, x0

    comp1 = get_space_components(eq_list1)
    comp2 = get_space_components(eq_list2)
    
    # 1. 共に空集合(解なし)の場合は同一とみなす
    if comp1 is None and comp2 is None:
        return True
        
    # 2. 一方だけが空集合の場合は異なる
    if comp1 is None or comp2 is None:
        return False
        
    K1, x0_1 = comp1
    K2, x0_2 = comp2
    
    # 3. 零空間(方向)が異なる場合は異なる
    if K1 != K2:
        return False
        
    # 4. 特解の差が零空間に属していれば完全に同一の空間
    return (x0_1 - x0_2) in K1

# 2026-09-05-2d-section.py
def restrict_faces_to_2dim(F, L, ambient_dim=None):
    """
    d 次元空間内の Face リスト F を、2次元平面 L で切断し、
    2次元空間内の Face リスト F2 として返す。
    
    入力:
        F: d 次元空間の Face のリスト
        L: 切断する 2次元平面を定義する等式のリスト (d-2 個の等式が必要)
        ambient_dim: (オプション) 空間の次元 d。省略した場合は F から推測する。
    出力:
        F2: 2次元平面 (t1, t2) 上にマッピング・簡略化された Face のリスト
    """
    if not F:
        return []

    # 空間の次元 d を特定
    if ambient_dim is None:
        for f in F:
            if f and len(f) == 2 and f[1]:
                ambient_dim = len(f[1][0]) - 1
                break
        if ambient_dim is None:
            raise ValueError("ambient_dim を特定できませんでした。")
    d = ambient_dim

    # 1. L をパラメータ化して、2次元平面 (t1, t2) の座標系を作る
    R, P = linear_parametrization_of_affine_space(L, ambient_dim=d)

    if R == QQ or R.ngens() != 2:
        raise ValueError(f"指定されたアフィン空間 L の次元が 2 ではありません。(現在のパラメータ数: {R.ngens() if R != QQ else 0})")

    t1, t2 = R.gens()
    F2 = []

    for f in F:
        if not f:
            continue
            
        dim_f, ieqs = f
        mapped_ieqs = []

        # 2. 各不等式を (t1, t2) 空間へ射影
        for c in ieqs:
            # c0 + c1*P1 + ... + cd*Pd の代入計算
            expr = R(c[0])
            for i in range(d):
                expr += R(c[i+1]) * P[i]

            # 1次式なので、t1, t2 で微分することで安全に係数を抽出できる
            c_const = expr.subs({t1: 0, t2: 0})
            c_t1 = expr.derivative(t1)
            c_t2 = expr.derivative(t2)

            mapped_ieqs.append([c_const, c_t1, c_t2])

        # 3. 2次元空間上の Face として簡略化し、空集合でなければ追加
        # (ここで simplify_face 内の Polyhedron が次元などを自動で計算・整理する)
        f2_clean = simplify_face([2, mapped_ieqs])

        if f2_clean and f2_clean not in F2:
            F2.append(f2_clean)

    return F2

## 2026.09.07
def is_face_contained_with_codim1(F1, F2):
    """
    Face F1 (閉包) に Face F2 が余次元 1 で含まれるか判定する。
    
    入力:
        F1: 基準となる大きな Face (isom_class_closure[i] に相当)
        F2: 境界となる小さな Face (isom_class[j] に相当)
    出力:
        True: F1 の境界として F2 が含まれ、かつ dim(F2) == dim(F1) - 1 の場合
        False: それ以外
    """
    if not F1 or not F2:
        return False
        
    dim1, ieqs1 = F1
    dim2, ieqs2 = F2
    
    # 余次元1のチェック
    if dim1 - 1 != dim2:
        return False
        
    try:
        P1 = Polyhedron(ieqs=ieqs1, base_ring=QQ)
        P2 = Polyhedron(ieqs=ieqs2, base_ring=QQ)
        # P2 が P1 に完全に包含されているか (共通部分が P2 自身になるか)
        return P1.intersection(P2) == P2
    except Exception as e:
        print(f"Error in inclusion check: {e}")
        return False

def find_nearby_point_in_face(p, F, max_dist=10):
    """
    点 p の近傍で、Face F に属する格子点 q を探す。
    
    入力:
        p: 起点となる座標リスト (例: [0, 1])
        F: 探索対象の Face (isom_class[i] に相当する内部表現)
        max_dist: 探索を打ち切る最大マンハッタン距離
    出力:
        q: 条件を満たす最初の格子点の座標リスト。見つからない場合は None。
    """
    if not F:
        return None
        
    dim, ieqs = F
    d = len(p)
    
    # 与えられた点 pt が Face のすべての制約を満たすか
    def is_in_face(pt):
        for eq in ieqs:
            # val = c0 + c1*x1 + ... + cd*xd
            val = eq[0] + sum(eq[i+1] * pt[i] for i in range(d))
            if val < 0:
                return False
        return True

    # 距離 1 から順に BFS で探索
    for dist in range(1, max_dist + 1):
        # マンハッタン距離が dist となるような d 次元のオフセットベクトルをすべて生成
        for offset in itertools.product(range(-dist, dist + 1), repeat=d):
            if sum(abs(x) for x in offset) == dist:
                q = [p[i] + offset[i] for i in range(d)]
                if is_in_face(q):
                    return q
    return None

def segment_to_face(p, q):
    """
    点 p と点 q を結ぶ線分を標準の Face 表現に変換する。
    
    入力:
        p, q: 座標のリスト (例: [0, 1], [1, 2])
    出力:
        [dim, ieqs_list] の形式の Face
    """
    P = Polyhedron(vertices=[p, q], base_ring=QQ)
    ieqs_list = []
    
    # 等式制約の抽出
    for eq in P.equations():
        vec = list(eq)
        lcm = 1
        for x in vec:
            if x != 0:
                lcm = (lcm * x.denominator()) // math.gcd(lcm, x.denominator())
        vec_int = [int(x * lcm) for x in vec]
        
        g = 0
        for x in vec_int:
            g = math.gcd(g, abs(x))
        if g > 1:
            vec_int = [x // g for x in vec_int]
            
        ieqs_list.append(vec_int)
        ieqs_list.append([-x for x in vec_int])
        
    # 不等式制約の抽出
    for ieq in P.inequalities():
        vec = list(ieq)
        lcm = 1
        for x in vec:
            if x != 0:
                lcm = (lcm * x.denominator()) // math.gcd(lcm, x.denominator())
        vec_int = [int(x * lcm) for x in vec]
        
        g = 0
        for x in vec_int:
            g = math.gcd(g, abs(x))
        if g > 1:
            vec_int = [x // g for x in vec_int]
            
        ieqs_list.append(vec_int)
        
    return [P.dim(), ieqs_list]

## helper test は 2026-09-07-for-merge.py に

## 2026.09.09
from sage.all import matrix, QQ, VectorSpace, vector

def is_affine_space_included_in_arrangement(W, L, Eq):
    """
    L に制限されたアフィン空間 W (\cap L) が、L 上の超平面配置 Eq (の和集合) に
    完全に含まれているかを判定する。
    
    入力:
        W: アフィン空間を定義する等式のリスト
        L: 制限するアフィン空間(ベースとなる空間)の等式のリスト
        Eq: 超平面配置を定義する等式のリスト
    出力:
        含まれていれば True, そうでなければ False
    """
    # 次元 d の特定
    d = -1
    for lst in [L, Eq, W]:
        if lst and len(lst) > 0:
            d = len(lst[0]) - 1
            break
            
    if d < 0:
        return False  # すべて空の場合は判定不能として False

    # L の行空間
    L_mat = matrix(QQ, L) if L else matrix(QQ, 0, d+1)
    M_L = L_mat.row_space()
    
    # W \cap L の行空間
    WL = (W if W else []) + (L if L else [])
    WL_mat = matrix(QQ, WL) if WL else matrix(QQ, 0, d+1)
    M_WL = WL_mat.row_space()
    
    target_empty = vector(QQ, [1] + [0] * d)
    
    # W \cap L が空集合の場合は、任意の部分集合となるため True を返す
    if target_empty in M_WL:
        return True
        
    if not Eq:
        return False
        
    # W \cap L が、Eq のいずれかの超平面(ただし L 自身を含まないもの)に含まれるか
    for eq in Eq:
        v_eq = vector(QQ, eq)
        if v_eq in M_L:
            # v_eq が L の行空間に属する場合、L 上では 0=0 となり壁として機能しない
            continue 
        if v_eq in M_WL:
            # W \cap L 上で v_eq = 0 が成り立つ
            return True
            
    return False

def get_new_hypersurfaces(W, L, Eq):
    """
    W \cap L が Eq on L に含まれない場合、Eq に加えるべき W の超平面のリスト WW を返す。
    L 上での線形独立性を考慮し、冗長な式を省いた「小さい集合」を抽出する。
    """
    WW = []
    if not W:
        return WW
        
    # 次元 d の特定
    d = -1
    for lst in [W, L, Eq]:
        if lst and len(lst) > 0:
            d = len(lst[0]) - 1
            break
    if d < 0:
        return WW
        
    L_mat = matrix(QQ, L) if L else matrix(QQ, 0, d+1)
    M_L = L_mat.row_space()
    
    valid_eqs = []
    for eq in W:
        v_eq = vector(QQ, eq)
        # L 上で自明な式 (0=0 になるもの) を除外
        if v_eq in M_L:
            continue
        # 1=0 のような式を除外
        if eq[0] != 0 and all(c == 0 for c in eq[1:]):
            continue
        valid_eqs.append(eq)
        
    if not valid_eqs:
        return WW
        
    # より係数がスパース（シンプル）な式を優先するためソート
    valid_eqs.sort(key=lambda eq: sum(1 for c in eq if c != 0))
    
    # L をベースにして、ランクが増加する（L 上で独立な）W の式のみを WW に追加
    current_matrix_rows = list(L) if L else []
    current_rank = matrix(QQ, current_matrix_rows).rank() if current_matrix_rows else 0
    
    for eq in valid_eqs:
        M_aug = matrix(QQ, current_matrix_rows + [eq])
        new_rank = M_aug.rank()
        if new_rank > current_rank:
            current_matrix_rows.append(eq)
            WW.append(eq)
            current_rank = new_rank
            
    return WW

