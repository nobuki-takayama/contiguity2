
# `chamber_lattice_utils.py` ユーザーマニュアル

本ライブラリは、SageMath 環境において超平面配置（Hyperplane Arrangement）の Face の網羅的抽出、有理多面体とイデアルの交差計算、相対的内部の格子点列挙、および図示を行うためのユーティリティ群です。

## 1. 動作環境とインポート

* **前提環境:** SageMath 9.x 以上
* **インポート:** 対象のスクリプトと同じディレクトリに `chamber_lattice_utils.py` を配置し、必要な関数をインポートして使用します。

```python
from sage.all import *
from chamber_lattice_utils import (
    list_all_arrangement_faces,
    intersect_face_with_linear_ideal,
    strict_interior_face_for_lattice,
    lattice_points_in_face,
    find_chamber_bases_and_shifts,
    plot_arrangement_faces
)

```

## 2. データ構造の定義

本ライブラリでは、多面体（Face）を以下のリスト形式で表現します。
`[次元 (dim), [不等式1, 不等式2, ...]]`

* **不等式の形式:** $b + a_1 x_1 + a_2 x_2 + \dots \ge 0$ を `[b, a1, a2, ...]` のリストで表現します。
* **等式の表現:** $E = 0$ は、2つの不等式 $E \ge 0$ (`[b, a1, ...]`) と $-E \ge 0$ (`[-b, -a1, ...]`) のペアとして格納されます。

## 3. 主要関数のリファレンスと使用例

### 3.1. 超平面配置からすべての Face を列挙する

**`list_all_arrangement_faces(ieqs_list)`**

* **概要:** 与えられた超平面のリストによって分割される空間の、すべての次元の Face を抽出します。符号ベクトル探索により、幾何学的に空となる領域は自動で枝刈りされます。
* **入力:** 超平面（方程式）を定義する係数リストのリスト。

```python
# x=0, y=0, x-y=0 の3つの超平面による分割
ieqs = [
    [0, 1, 0],   # x = 0
    [0, 0, 1],   # y = 0
    [0, 1, -1]   # x - y = 0
]
faces = list_all_arrangement_faces(ieqs)
# 戻り値: 2次元Chamber(6個), 1次元Edge(6個), 0次元Vertex(1個) のリスト

```

### 3.2. Face と線形イデアルの交差を求める

**`intersect_face_with_linear_ideal(face, ideal_gens, ring)`**

* **概要:** 多面体 $C$ と、一次式の積で生成されるイデアル $I$ の零点集合 $V(I)$ との共通部分を求め、極大な Face のリストを返します。
* **入力:**
* `face`: 対象となる Face 構造
* `ideal_gens`: イデアルの生成元（多項式のリスト）
* `ring`: 変数を定義した多項式環



```python
R = PolynomialRing(QQ, 's1, s2')
s1, s2 = R.gens()

face_C = [2, [[0, 1, 0], [0, 0, 1]]] # 第1象限 (s1>=0, s2>=0)
ideal_gens = [
    (s1 - 1) * (s1 + s2 - 1),
    (s2 - 2) * (s1 + s2 - 1)
]

# 共通部分の計算
intersections = intersect_face_with_linear_ideal(face_C, ideal_gens, R)
# 戻り値: 線分 s1+s2=1 (s1>=0, s2>=0) と 点 (1, 2) の2つのFaceリスト

```

### 3.3. 格子点を抽出するために Face を縮小する

**`strict_interior_face_for_lattice(face)`**

* **概要:** 境界上の点を除外し、相対的内部（Relative Interior）の格子点のみを含むように領域を縮小（不等式の定数項を -1）します。縮小後に領域が潰れる場合は `None` を返します。

```python
face_example = [2, [[0, 1, 0], [0, 0, 1], [3, -1, -1]]] # x>=0, y>=0, x+y<=3
interior_face = strict_interior_face_for_lattice(face_example)
# 戻り値: [2, [[-1, 1, 0], [-1, 0, 1], [2, -1, -1]]] (x>=1, y>=1, x+y<=2)

```

### 3.4. 有界な Face 内部の格子点を列挙する

**`lattice_points_in_face(face)`**

* **概要:** 有界（コンパクト）な Face の相対的内部にある格子点をすべて列挙します。

```python
# 先ほど縮小した interior_face を渡すか、元の face を渡すことができます
pts = lattice_points_in_face(face_example)
# 戻り値: [[1, 1]]

```

### 3.5. チェンバーの基点とヒルベルト基底を計算する

**`find_chamber_bases_and_shifts(ieqs)`**

* **概要:** 同次化（Cayley Trick）を用いて、非有界なチェンバー内のすべての格子点を網羅するための「最小の基点（Base Points）」と「退避錐のヒルベルト基底（Safe Shifts）」を計算します。

```python
# -1 -s1 -2s2 >= 0, -1 +s1 -s2 >= 0 のチェンバー
ieqs_chamber = [[-1, -1, -2], [-1, 1, -1]]
bases, shifts = find_chamber_bases_and_shifts(ieqs_chamber)
# bases: [[0, -1], [1, -1]], shifts: [[-1, -1], [0, -1], [1, -1], [2, -1]]

```

### 3.6. 描画機能 (2D / 3D)

**`plot_arrangement_faces(faces, xmin, xmax, ymin, ymax)`**
**`plot_arrangement_faces_3d(faces, xmin, xmax, ymin, ymax, zmin, zmax)`**

* **概要:** Face のリストと格子点を指定された描画範囲のバウンディングボックスで切り取り、次元ごとに異なるスタイル（透明度、太さなど）で図示します。

```python
# faces は list_all_arrangement_faces などで取得したリスト
G = plot_arrangement_faces(faces, xmin=-3, xmax=3, ymin=-3, ymax=3)
G.show(aspect_ratio=1)

```

### 3.7. face が span する affine 空間の定義式.
**`get_equality_constrants(face)`**

* **概要:** face が span する affine 空間の定義式を出力する.


```python
face_example = [1, [0,1,0], [0,1,-1], [0,0,1], [0,-2,2]]
eqs = get_equality_constraints(face_example)
```

### 3.8. face の affine hull を求める. affine space に属さない face の列挙
**`get_affine_hull_equations(face)`**
**`get_faces_not_in_affine_hull(faces,L)`**


* **概要:** get_affine_hull_equations は face が span する affine 空間の定義式を出力する. 3.7 と同様な機能(のはず).
get_faces_not_in_affine_hull(faces,L) は affine space L に属さない faces を列挙する.


```python
face_example = [1, [[0,1,0,0], [0,0,1,0],[0,0,0,1],[-1,1,0,0],[1,-1,0,0],[-1,0,1,0],[1,0,-1,0],[-1,0,0,1],[1,0,0,-1]]]
get_affine_hull_equations(face_example)
[[-1, 0, 0, 1], [-1, 0, 1, 0], [-1, 1, 0, 0]]

# テスト用の 3次元 Face 群
# 空間: x, y, z
faces_example = [
    # 2D Face: z = 0, x >= 0, y >= 0 (xy平面上の第1象限)
    [2, [[0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1], [0, 0, 0, -1]]],
    
    # 3D Chamber: z >= 0, x >= 0, y >= 0 (第1八分空間)
    [3, [[0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]],
    
    # 1D Edge: x = 0, z = 0, y >= 0 (y軸の正の部分)
    [1, [[0, 1, 0, 0], [0, -1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1], [0, 0, 0, -1]]]
]

print("--- 1. Affine Hull の抽出 ---")
target_face = faces_example[0]  # xy平面上の2D Face
print(f"Target Face (dim={target_face[0]}): {target_face[1]}")
hull_eqs = get_affine_hull_equations(target_face)
print(f"Affine Hull Equations: {hull_eqs}")
# 期待出力: z = 0 を表す [[0, 0, 0, 1]] が抽出される
--- 1. Affine Hull の抽出 ---
Target Face (dim=2): [[0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1], [0, 0, 0, -1]]
Affine Hull Equations: [[0, 0, 0, 1]]


print("\n--- 2. Affine Hull に属さない Face の列挙 ---")
filtered_faces = get_faces_not_in_affine_hull(faces_example, hull_eqs)

print(f"Faces not in Affine Hull (Equations: {hull_eqs}):")
for f in filtered_faces:
    print(f"dim={f[0]}, ieqs={f[1]}")
    
# 期待出力:
# 3D Chamber は z方向の厚みがあるため抽出される
# 1D Edge (y軸) は z=0 平面に完全に含まれるため抽出されない (除外される)
--- 2. Affine Hull に属さない Face の列挙 ---
Faces not in Affine Hull (Equations: [[0, 0, 0, 1]]):
dim=3, ieqs=[[0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]

```

### 3.9. 
**`list_faces_of_arrangement_on_subspace(A_arrangement,L_subspace)`**

* **概要:**  affine subspace L_subspace で定義された A_arrangement の すべての face を列挙する. A_arrangement は全体空間の座標で定義定義.
list_faces_of_arrangement_on_subspace(A_arrangement,[])
は list_all_arrangement_faces(A_arrangement) と同様の結果を与える.

```python
f1=list_all_arrangement_faces([[0,1,0],[0,1,1],[0,-1,-1]])
f2=list_faces_of_arrangement_on_subspace([[0,1,0],[0,1,1],[0,-1,-1]],[])
f1==f2
True
```

```python
print("=== Intersection of Arrangement A and Subspace L ===")

# 3次元空間 (x, y, z) の設定
# L: z = 1 の平面 (すなわち 1 - z = 0 -> [1, 0, 0, -1])
L_subspace = [
    [1, 0, 0, -1]
]

# A: x = 0 と y = 0 の2つの超平面
# これにより、z=1 平面上で4つの象限(2次元Face)が形成されるはず
A_arrangement = [
    [0, 1, 0, 0],  # x = 0
    [0, 0, 1, 0]   # y = 0
]

print(f"Subspace L: {L_subspace}")
print(f"Arrangement A: {A_arrangement}\n")

faces = list_faces_of_arrangement_on_subspace(A_arrangement, L_subspace)

dim_counts = {}
for f in faces:
    d = f[0]
    dim_counts[d] = dim_counts.get(d, 0) + 1
    print(f)
    
print("\n=== Summary ===")
for d, count in sorted(dim_counts.items(), reverse=True):
    print(f"{d}-dimensional faces: {count}")
```

結果
```python
=== Intersection of Arrangement A and Subspace L ===
Subspace L: [[1, 0, 0, -1]]
Arrangement A: [[0, 1, 0, 0], [0, 0, 1, 0]]

[2, [[1, 0, 0, -1], [-1, 0, 0, 1], [0, 1, 0, 0], [0, 0, 1, 0]]]
[2, [[1, 0, 0, -1], [-1, 0, 0, 1], [0, 1, 0, 0], [0, 0, -1, 0]]]
[2, [[1, 0, 0, -1], [-1, 0, 0, 1], [0, -1, 0, 0], [0, 0, 1, 0]]]
[2, [[1, 0, 0, -1], [-1, 0, 0, 1], [0, -1, 0, 0], [0, 0, -1, 0]]]
[1, [[1, 0, 0, -1], [-1, 0, 0, 1], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, -1, 0]]]
[1, [[1, 0, 0, -1], [-1, 0, 0, 1], [0, -1, 0, 0], [0, 0, 1, 0], [0, 0, -1, 0]]]
[1, [[1, 0, 0, -1], [-1, 0, 0, 1], [0, 1, 0, 0], [0, -1, 0, 0], [0, 0, 1, 0]]]
[1, [[1, 0, 0, -1], [-1, 0, 0, 1], [0, 1, 0, 0], [0, -1, 0, 0], [0, 0, -1, 0]]]
[0, [[1, 0, 0, -1], [-1, 0, 0, 1], [0, 1, 0, 0], [0, -1, 0, 0], [0, 0, 1, 0], [0, 0, -1, 0]]]

=== Summary ===
2-dimensional faces: 4
1-dimensional faces: 4
0-dimensional faces: 1
```


### 3.10. アフィン部分空間のパラメータ表示を生成する

**`linear_parametrization_of_affine_space(eqs, ambient_dim=None, ring=None, v=None)`**

* **概要:** 連立一次方程式で定義される $m$ 次元のアフィン部分空間を、特解と零空間の基底を用いてパラメータ表示 $(P_1, \dots, P_d)$ に変換します。パラメータとして利用する変数を既存の多項式環から指定することが可能です。
* **入力:**
* `eqs`: アフィン空間を定義する等式リスト。各等式は $c_0 + c_1 x_1 + \dots + c_d x_d = 0$ を表す `[c0, c1, ..., cd]` の形式で与えます（`get_affine_hull_equations` の出力をそのまま渡すことができます）。
* `ambient_dim`: `eqs` が空リスト（全空間）の場合のみ、元の空間の次元を指定します。
* `ring`: （オプション）戻り値の座標が属する多項式環。指定しない場合は $QQ[t_1, \dots, t_m]$ が自動生成されます。
* `v`: （オプション）パラメータとして使用する `ring` の元のリスト。部分空間の次元 $m$ に対して指定された変数が多めにある場合は、先頭から $m$ 個を使用します。不足している場合はエラーになります。


* **出力:**
* `R`: 指定された `ring`、または自動生成された多項式環。部分空間が 0次元（点）で `ring` の指定もない場合は有理数体 $QQ$ を返します。
* `P`: パラメータ表示された座標のリスト `[P1, ..., Pd]`。



```python
# 例: 特定の多項式環の変数をパラメータとして割り当てる
R1 = PolynomialRing(QQ, 'x4,x5,x6,dx4,dx5,dx6,d,nu1,nu2')
dic1 = R1.gens_dict()

# パラメータの候補リストを作成
d_var, nu1, nu2 = dic1['d'], dic1['nu1'], dic1['nu2']

# 対象のアフィン空間 (3次元空間内の直線)
# 1 + 1*x1 = 0, 1 + 1*x2 = 0
eqs = [
    [1, 1, 0, 0],
    [1, 0, 1, 0]
]

# 変数リスト [d_var, nu1, nu2] を渡し、必要な数(この場合は1つ)だけ使わせる
R_out, P_out = linear_parametrization_of_affine_space(eqs, ring=R1, v=[d_var, nu1, nu2])

print(R_out) # Multivariate Polynomial Ring in x4, x5, x6, dx4, dx5, dx6, d, nu1, nu2 over Rational Field
print(P_out) # [-1, -1, d]  (※先頭の変数 'd' が自由変数のパラメータとして使われる)

```



```python
# 3次元空間内の直線 (x=1, y=2) をパラメータ化する
# 方程式: -1 + 1*x + 0*y + 0*z = 0, -2 + 0*x + 1*y + 0*z = 0
eqs_line = [
    [-1, 1, 0, 0],
    [-2, 0, 1, 0]
]

R, P = linear_parametrization_of_affine_space(eqs_line)

print(R) # Univariate Polynomial Ring in t1 over Rational Field
print(P) # [1, 2, t1]

#R.gens_dict() で dic を取り出せる.  R.gens() は変数リスト

```

### 3.11. アフィン部分空間の冗長な等式を削除する
**`reduce_affine_space_expression(E)`**
*   **概要:** アフィン部分空間を定義する等式リストの中から、他の式の線形結合で表せる線形従属（冗長）な等式を削除し、独立な等式のみを残した（縮約された）リストを返します。リスト内の元の式（係数）の形はそのまま保持されます。
*   **入力:**
    *   `E`: アフィン空間を定義する等式リスト `[[c0, c1, ..., cd], ...]`。
*   **出力:**
    *   `reduced_E`: 冗長な式が除外された等式リスト。

```python
# eq3 は eq1 と eq2 の和になっているため冗長
E = [
    [-1, 1, 0, 0],
    [-1, 0, 1, 0],
    [-2, 1, 1, 0]
]

reduced_E = reduce_affine_space_expression(E)
print(reduced_E)
# 戻り値: [[-1, 1, 0, 0], [-1, 0, 1, 0]]
```


### 3.12. 線形イデアルからアフィン空間の和集合（極大コンポーネント）を抽出する

**`get_affine_spaces_from_linear_ideal(ideal_gens, ring)`**

* **概要:** 一次式の積で生成されるイデアル $I$ の生成元リストから、その零点集合 $V(I)$ を構成する極大なアフィン空間（既約成分）のリストを抽出します。生成元の因数分解の組み合わせから方程式系を構築し、解を持たない空集合や、他の空間に完全に包含される小さな部分空間（冗長な成分）を幾何学的に判定して自動で除外します。
* **入力:**
* `ideal_gens`: イデアルの生成元のリスト。各多項式は一次式の積（例: $x(x-1)$）で構成されていることを前提とします。
* `ring`: 対象の変数を定義した多項式環。


* **出力:**
* `affine_spaces`: 極大なアフィン空間を定義する等式リストのリスト。各アフィン空間は、等式 $c_0 + c_1 x_1 + \dots = 0$ を表す `[[c0, c1, ...], ...]` の形式で表現されます。



```python
R = PolynomialRing(QQ, 'x, y, z')
x, y, z = R.gens()

# I = < x*(x-1), y >
# V(I) は (x=0, y=0) の直線と (x=1, y=0) の直線の和集合
ideal_gens = [
    x * (x - 1),
    y
]

spaces = get_affine_spaces_from_linear_ideal(ideal_gens, R)

for space in spaces:
    print(space)

# 出力結果:
# [[0, 1, 0, 0], [0, 0, 1, 0]]    (x=0, y=0 を表す)
# [[-1, 1, 0, 0], [0, 0, 1, 0]]   (-1+x=0, y=0 を表す)

```

### 3.13. アフィン空間と Face の共通部分を求める
**`intersection_of_face_and_affine_space(L, F)`**
*   **概要:** アフィン空間 $L$ と Face $F$ の共通部分（交差）を計算します。交差が存在する場合は新しい次元と制約式を計算して標準の Face 形式で返し、幾何学的に交差しない（空集合となる）場合は空リストを返します。
*   **入力:**
    *   `L`: アフィン空間を定義する等式のリスト `[[c0, c1, ...], ...]`。
    *   `F`: 対象となる Face 構造 `[dim, [[b, a1, ...], ...]]`。
*   **出力:**
    *   共通部分がある場合: 新しい Face `[new_dim, eqs_list]`。
    *   共通部分がない場合: 空リスト `[]`。

```python
# F: x >= 0, y >= 0 (2次元の第1象限)
F = [2, [[0, 1, 0], [0, 0, 1]]]

# L: x = 1 (-1 + x = 0)
L = [[-1, 1, 0]]

intersected = intersection_of_face_and_affine_space(L, F)
print(intersected)
# 出力結果: [1, [[-1, 1, 0], [1, -1, 0], [0, 0, 1]]]
# (x=1 の等式と y>=0 の不等式を持つ1次元の半直線)

# 交差しない場合 (例: x = -1)
L_empty = [[1, 1, 0]]
print(intersection_of_face_and_affine_space(L_empty, F))
# 出力結果: []
```

### 3.14. Face を人間が読みやすい形式で表示する
**`show_face(F)`**
*   **概要:** リストのリストとして表現されている Face を、解釈しやすい座標リストや SageMath の関係式オブジェクト（`==` や `>=`）に変換します。デバッグや対話モードでの結果確認に最適です。
*   **入力:** `[dim, [[b, a1, ...], ...]]` 形式の単一の Face、またはそのリスト。
*   **出力:**
    *   `dim == 0` の場合: 頂点座標のリスト（例: `[1, 2]`）を返します。
    *   `dim >= 1` の場合: `x, y, z` などの変数を用いた SageMath の関係式オブジェクトのリスト（例: `[-x + 2 >= 0, y + 1 == 0]`）を返します。等式制約は自動的に判定されます。
    *   リストが入力された場合: 自動的に `map` が適用され、変換された結果のリストが返ります。

```python
F_list = [
    [0, [[1, -1, 0], [2, 0, -1]]],
    [1, [[-1, 1, 0], [1, -1, 0], [0, 0, 1]]]
]

# リストを渡すとまとめて見やすく変換
print(show_face(F_list))
# 出力結果:
# [
#   [1, 2],
#   [x - 1 == 0, y >= 0]
# ]
```

### 3.15. 二つのアフィン空間が同一かどうか判定する
**`is_same_affine_space(eq_list1, eq_list2, ambient_dim=None)`**
*   **概要:** 二つのアフィン空間を定義する等式リストが、幾何学的に全く同じ空間を表しているかどうかを判定します。方程式の係数倍や、線形従属な冗長な式の有無にかかわらず厳密に判定します。また、両者が「空集合（解を持たない矛盾した系）」である場合も同一とみなします。
*   **入力:**
    *   `eq_list1`, `eq_list2`: 比較するアフィン空間を定義する等式のリスト `[[c0, c1, ...], ...]`。
    *   `ambient_dim`: （オプション）両方のリストが空 `[]` の場合のみ指定が必要です。
*   **出力:**
    *   幾何学的に同一（または共に空集合）であれば `True`、そうでなければ `False`。

```python
# L1: x = 1 (3次元空間)
L1 = [[-1, 1, 0, 0]]

# L2: 2x = 2, 0 = 0 (冗長な式を含むが幾何学的には同じ x=1)
L2 = [[-2, 2, 0, 0], [0, 0, 0, 0]]

# L3: x = 2
L3 = [[-2, 1, 0, 0]]

print(is_same_affine_space(L1, L2)) # True
print(is_same_affine_space(L1, L3)) # False
```


### 3.16. Face の冗長な条件を削除し簡略化する

**`simplify_face(face)`**

* **概要:** イデアルの交差計算などで肥大化した Face の制約式リストから、幾何学的に無意味な（他の厳しい条件に吸収される）冗長な不等式や、重複する等式を完全に削除します。SageMath の `Polyhedron` が持つ最小 H表現 (Minimal H-representation) の計算を利用しており、極限まで簡略化・既約化された標準の Face 形式を返します。
* **入力:**
* `face`: 簡略化したい Face 構造 `[dim, [[b, a1, ...], ...]]`。


* **出力:**
* 冗長な条件が削ぎ落とされた新しい Face `[dim, new_ieqs_list]`。空間が空集合になる矛盾した制約の場合は空リスト `[]` を返します。



```python
# テスト用データ (2次元空間 x, y)
# 条件1: -y + 2 >= 0
# 条件2: y - 2 >= 0   (条件1と合わせて y = 2 を意味する)
# 条件3: -y + 2 >= 0  (条件1の重複)
# 条件4: -4x - 3y - 3 >= 0
# 条件5: -4x - 3y - 2 >= 0 (条件4より緩いため冗長)

F_redundant = [1, [
    [2, 0, -1],
    [-2, 0, 1],
    [2, 0, -1],
    [-3, -4, -3],
    [-2, -4, -3]
]]

F_clean = simplify_face(F_redundant)

print(F_clean)
# 出力結果: 
# [1, [[-2, 0, 1], [2, 0, -1], [-9, -4, 0]]]
# (y - 2 == 0 を表す等式ペアと、y=2 を代入して整理された最小の不等式 -4x - 9 >= 0 のみが残る)

```

### 3.17. パラメータベクトルの定数部分をテンプレートで設定・拡張する
**`set_constant_part_of_param_vector(template)`**
**`extend_param_vector_by_constant_part(values)`**
**`linear_space_by_extension_param()`**

*   **概要:** 変数を含むパラメータベクトルを生成する際に、値が固定されている部分（定数）と変動する部分（変数）を管理します。`set_constant_part_of_param_vector` で設定したテンプレートはモジュール内に記憶され、以降の `extend_param_vector_by_constant_part` の呼び出し時に `None` の部分が与えられた変数リストで自動的に埋められます。

なお, この関数は func_contiguity.py で定義されている.

*   **入力:**
    *   `template`: 定数と変数スロット `None` を含むリスト（例: `[1, 2, None, None, 3]`）。`0` を指定するとテンプレート機能が無効になります（初期値は `0`）。
    *   `values`: テンプレートの `None` の部分に順番に挿入する変数のリスト `[a, b, ...]`。
*   **出力:**
    *   `extend_param_vector_by_constant_part` は、テンプレートの指定位置に変数を挿入した新しいリストを返します。テンプレートが `0` の場合は `values` をそのまま返します。
    *   `linear_space_by_extension_param` は 数の部分で定義された linear space の定義を戻します. 上の例では [[1,-1,0,0,0],[2,0,-1,0,0,0],[3,0,0,0,0,-1]].
テンプレートが0の時は [] を戻す.

```python
# テンプレートを記憶させる
set_constant_part_of_param_vector([1, 2, None, None, 3])

# 2つの変数 [a, b] を与えて拡張する
expanded = extend_param_vector_by_constant_part(['a', 'b'])
print(expanded)
# 出力結果: [1, 2, 'a', 'b', 3]

# テンプレートを初期状態(無効)に戻す
set_constant_part_of_param_vector(0)
print(extend_param_vector_by_constant_part(['a', 'b']))
# 出力結果: ['a', 'b']
```

### 3.18. 多次元の Face リストを 2次元平面で切断・射影する
**`restrict_faces_to_2dim(F, L, ambient_dim=None)`**
*   **概要:** $d$ 次元空間内の Face リスト $F$ を、指定した 2次元平面 $L$ で切断し、共通部分を 2次元パラメータ空間 $(t_1, t_2)$ 上の Face リストとして返します。内部で `linear_parametrization_of_affine_space` と `simplify_face` を組み合わせることで、高次元の制約を 2D プロット用のデータ形式（`plot_arrangement_faces` 関数等の入力）へ自動的に変換します。
*   **入力:**
    *   `F`: 対象となる $d$ 次元空間の Face のリスト。
    *   `L`: 切断する 2次元平面を定義する等式のリスト。$d$ 次元空間内で 2次元平面を定義するためには、通常 $d-2$ 個の独立な等式が必要です。
    *   `ambient_dim`: （オプション）空間の次元 $d$。省略した場合は $F$ の構造から自動推測されます。
*   **出力:**
    *   `F2`: 2次元空間上にマッピングされ、冗長な式が除外された Face のリスト。交差しない Face は自動的に破棄されます。

```python
# 3次元空間の Face リスト F_3d があるとする
# L: x = 1 で切断する (-1 + 1*x + 0*y + 0*z = 0)
L = [[-1, 1, 0, 0]]

F2 = restrict_faces_to_2dim(F_3d, L, ambient_dim=3)

# この F2 はそのまま 2D描画関数へ渡すことができる
# plot_arrangement_faces(F2)
```



### 3.19. Face が別の Face の余次元 1 の境界として含まれるか判定する

**`is_face_contained_with_codim1(F1, F2)`**

* **概要:** 大きな Face `F1`（の閉包）の境界として、小さな Face `F2` が完全に包含されており、かつその次元の差（余次元）がちょうど 1 であるかどうかを幾何学的に判定します。隣接するチャンバー同士の接続判定などに利用します。
* **入力:**
* `F1`: 基準となる大きな Face 構造 `[dim, ieqs_list]`。
* `F2`: 境界候補となる小さな Face 構造 `[dim, ieqs_list]`。


* **出力:**
* 条件を満たす（包含され、かつ `dim(F2) == dim(F1) - 1` である）場合は `True`、それ以外は `False`。



```python
# F1: x >= 0, y >= 0 (2次元の第1象限)
F1 = [2, [[0, 1, 0], [0, 0, 1]]]

# F2: x == 0, y >= 0 (1次元のy軸上の半直線)
F2 = [1, [[0, 1, 0], [0, -1, 0], [0, 0, 1]]]

print(is_face_contained_with_codim1(F1, F2))
# 出力結果: True

```

### 3.20. 指定した点の近傍にある Face 内の格子点を探索する

**`find_nearby_point_in_face(p, F, max_dist=10)`**

* **概要:** 起点となる座標 `p` から出発し、Face `F` の内部（制約式をすべて満たす領域）に含まれる最初の格子点 `q` をマンハッタン距離による幅優先探索（BFS）で見つけ出します。境界上の点からチャンバー内部の代表点を探す際などに活用します。
* **入力:**
* `p`: 起点となる座標のリスト（例: `[0, 1]`）。
* `F`: 探索対象となる Face 構造 `[dim, ieqs_list]`。
* `max_dist`: （オプション）探索を打ち切る最大のマンハッタン距離。デフォルトは 10。


* **出力:**
* Face 内の格子点が見つかった場合はその座標リスト `q`。見つからなかった場合は `None`。



```python
# 境界上の点 p=(0, 2)
p = [0, 2]

# F: x >= 1, y >= 1 (内部にシフトした Face を想定)
F_int = [2, [[-1, 1, 0], [-1, 0, 1]]]

q = find_nearby_point_in_face(p, F_int)
print(q)
# 出力結果: [1, 2] など、条件を満たす最も近い格子点

```

### 3.21. 2点を結ぶ線分を Face 表現に変換する

**`segment_to_face(p, q)`**

* **概要:** 与えられた2つの点 `p` と `q` を両端とする線分（1次元の Polyhedron）を構築し、それを標準の Face 形式（次元と等式・不等式制約のリスト）に変換して返します。隣接関係（Contiguity）を表すエッジを他の Face と同じデータ構造で扱うために利用します。
* **入力:**
* `p`, `q`: 線分の両端となる座標のリスト（例: `[0, 2]`, `[1, 2]`）。


* **出力:**
* 線分を表す標準の Face 構造 `[dim, ieqs_list]`。



```python
p = [0, 2]
q = [1, 2]

seg_face = segment_to_face(p, q)
print(seg_face)
# 出力結果: [1, [[-2, 0, 1], [2, 0, -1], [0, 1, 0], [1, -1, 0]]]
# (y == 2, x >= 0, -x + 1 >= 0 に相当する制約のリスト)

```


### 3.22. アフィン空間が制限空間上の超平面配置に含まれるか判定する

**`is_affine_space_included_in_arrangement(W, L, Eq)`**

* **概要:** アフィン空間 $L$ に制限された空間 $W \cap L$ が、同じく $L$ 上で定義された超平面配置 $Eq$ （の和集合）に完全に包含されているかを判定します。各等式のリストから行空間（Row space）を生成し、ランク計算を用いて高速かつ厳密に幾何学的な包含関係を評価します。再帰処理における不要な「壁」の追加判定（終了条件）として機能します。
* **入力:**
* `W`: 判定対象のアフィン空間を定義する等式のリスト。
* `L`: 制限のベースとなるアフィン空間（現在のチャンバーを含む空間）を定義する等式のリスト。
* `Eq`: 超平面配置を定義する等式のリスト。


* **出力:**
* 完全に包含されていれば `True`、そうでなければ `False`。



```python
# L: ベースとなる空間
L = [[-4, 1, 0, 0], [3, 0, -1, -1]]

# Eq: L 上の超平面配置
Eq = [[2, 0, -1, -1], [-4, 1, 0, 0], [3, 0, -1, -1], [-4, 1, 0, 0], 
      [-1, 0, 1, 0], [-4, 1, 0, 0], [0, 0, 1, 0], [-4, 1, 0, 0]]

# W: 新たに見つかった壁の候補
W = [[-1, 0, 0, 1], [2, 0, -1, 0], [-4, 1, 0, 0]]

print(is_affine_space_included_in_arrangement(W, L, Eq))
# 出力結果: False (W \cap L は Eq on L に完全には含まれない)

```

### 3.23. 超平面配置に加えるべき新しい独立な超平面を抽出する

**`get_new_hypersurfaces(W, L, Eq)`**

* **概要:** $W \cap L$ が $Eq$ に含まれない場合、$W$ の中から $Eq$ に追加すべき最小限の超平面（方程式）のリストを抽出します。$L$ の行空間を法として線形独立な基底のみを選び出し、かつ非ゼロの係数が少ない（よりシンプルな）式を優先的に選択することで、計算量の爆発を防ぎます。
* **入力:**
* `W`, `L`, `Eq`: `is_affine_space_included_in_arrangement` と同一。


* **出力:**
* 追加すべき線形独立でシンプルな超平面（方程式）のリスト `WW`。



```python
# 上記の L, Eq, W を用いて抽出
WW = get_new_hypersurfaces(W, L, Eq)

print(WW)
# 出力結果: [[-1, 0, 0, 1]] 
# (W の中から L 上で独立かつ最もシンプルな式が選ばれる)

# 抽出した WW を Eq に加えることで再計算へ進む
# Eq_new = Eq + WW

```


## 4. asir との interface: func_contiguity.py


### 4.1 Contiguity の計算

**`func_contiguity(Old,New,ring=Ring, dic=DIC,fuc=FUC)`**

* **概要:** load_func_contiguity で読み込まれた asir 関数の FUNC を呼び出して contiguity を計算する. $D/{\tt New} \rightarrow D/{\tt Old}$ なる左D加群の射, b function, その他 を戻す.
関数に作用させる時は 
$Hom_D(D/{\tt Old},{\cal O}) \rightarrow Hom_D(D/{\tt New},{\cal O})$ なる写像である.


```python
# サンプルプログラムの初期化
load_func_contiguity('2026-08-12-1f1-by-gkz.rr')
R = PolynomialRing(QQ,'x3,dx3,a,c')
dic = R.gens_dict()
a,c = dic['a'], dic['c']
```


```python
sage: load("func_contiguity.py")

# [a,c] と [a+1,c-1] の同型が成立しないかもしれない b-ideal を求める.
# c1f1_isom(Chi) なる asir 関数が必要. 
# これが定義された asir のファイルは load_func_contiguity(file_name) で読み込む.
sage: b=func_isom([1,-1],func='c1f1')
#func_isom calls:  Ans_tmp = c1f1_isom([1, -1]);

# contiguity を求める.
# c1f1_contiguity(Old,New) なる asir 関数が必要.
# 多項式環の定義もしておく.
sage: R=PolynomialRing(QQ,'x3,dx3,a,c')
sage: dic = R.gens_dict()
sage: a,c = dic['a'], dic['c']
sage: c1=func_contiguity([a+1,c-1],[a,c],ring=R,dic=dic,func='c1f1')
sage: c1
[x3*dx3^3 + 3*x3*dx3^2 + dx3^2*c + 3*x3*dx3 + 2*dx3*c + x3 + c,
 -a^3 + 2*a^2*c - a*c^2 - 3*a^2 + 3*a*c - 2*a,
 '[[a+1,c-1], -> ,[a,c]]',
 '[(x1^3*dx2+x1^2*x3)*dx3^2+(2*x1^2*x2*dx2+2*x1*x2*x3+2*x1^2)*dx3+x1*x2^2*dx2+x2^2*x3+2*x1*x2,-a^3+(2*c-3)*a^2+(-c^2+3*c-2)*a,dx1*dx2,1]',
 1,
 Multivariate Polynomial Ring in x3, dx3, a, c over Rational Field]


```

その他の contiguity function
```python
# banana_one_loop case.  
load_func_contiguity('2026-08-07-banana-one-loop-by-gkz.rr')
contiguity_verbose(0)
R1 = PolynomialRing(QQ,'x4,x5,x6,dx4,dx5,dx6,d,nu1,nu2')
dic1 = R1.gens_dict()
d,nu1,nu2 = dic1['d'], dic1['nu1'], dic1['nu2']
```

```python
# f0134 のテスト. A=[[1,0,-2,-3],[0,1,3,4]]
load_func_contiguity('2026-08-28-f0134.rr')
contiguity_verbose(0)
R2 = PolynomialRing(QQ,'x3,x4,dx3,dx4,b1,b2')
dic2 = R2.gens_dict()
b1,b2 = dic2['b1'], dic2['b2']
```


### 4.2. 新しいインタフェース _sp 関数.

**`func_isom_sp(Old,New,ring=Ring, dic=DIC,fuc=FUC)`**

* **概要:** $D/{\tt Old}$ と $D/{\tt New}$ が同型とならないパラメータの条件を戻す.

```python
load_func_contiguity('2026-08-12-1f1-by-gkz.rr')
R = PolynomialRing(QQ,'x3,dx3,a,c')
dic = R.gens_dict()
a,c = dic['a'], dic['c']
# 3. Contiguity の計算 (例: a -> a + 1)
Chi = func_get_chi([a,c],[a+1,c],func='c1f1') # 斉藤論文の \chi を求める
L1 = func_contiguity([a,c], [a+1,c], ring=R, dic=dic,func='c1f1')
L2 = func_contiguity([a+1,c], [a,c], ring=R, dic=dic,func='c1f1')
#  sp は specific parameter Old=[a,c] と New=[a+1,c] の間の同型が壊れる条件.
L3 = func_isom_sp([a,c],[a+1,c], ring=R, dic=dic,func='c1f1')

```

### 4.3. asir.py 
```python
f=asir.asir('[[1,2,3],[4,5,1/2]]')
# type(f) は <class 'asir.AsirElement'>
g=sage_eval(f'{f})
# sage が扱える数字のリストに.
```

注意:   python や sagemath の print は長い出力には自動的に改行を入れる.
  asir.py の asir.eval は長い引数入力は一時ファイルに書き出し load で ox_asir に評価させている.
  load する限りは改行がはいっていても asir 側の問題はない.
  (asir.py の _read_in_file_command(self, filename):
  に end$ を追加した.  end-of-file detected の error packet を作らないため. 2026.09.04)
  したがって pop_string は load の戻り値を戻すのみである.
  結果は変数に格納してから, asir.eval('変数名') で取り出さないといけない.

```python
result=asir.eval('print("aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa")')
print(result)
result=asir.eval(str(2^1000)+';')
print(result)  # 0 しか戻らない
result=asir.eval('Asir_tmp='+str(2^1000)+';')
print(asir.eval('Asir_tmp'))

```


### 4.4. 

**`bf_factor_to_affine_space(S, bfunction_factor)`**

* **概要:**  S(t) で parametrize された affine 空間で b-function を求めると t の関数になる.  b-function b(t) の linear factor を bfunction_factor とした時, これが 0 となる affine 空間の表現を求める. 

```python
load_func_contiguity('2026-08-12-1f1-by-gkz.rr')
R1 = PolynomialRing(QQ,'x4,x5,x6,dx4,dx5,dx6,d,nu1,nu2')
dic1 = R1.gens_dict()
d,nu1,nu2 = dic1['d'], dic1['nu1'], dic1['nu2']
print(bf_factor_to_affine_space([2*nu1+2*1,nu1,1],nu1-3))
## 結果 [[-1, 0, 0, 1], [-3, 0, 1, 0], [-8, 1, 0, 0]]
```

### 4.5. 

**`b_ideal_to_affine_spaces(S, b_ideal)`**

* **概要:**  S(t) で parametrize された affine 空間で b_ideal を求めると t 変数のイデアルとなる.  V(b_ideal) は affine 空間達の和集合になるが, その affine 空間を求める.

戻り値は
[aff1, aff2, ...].
ここで,
aff1=[[1,-1,0,0]] (codim 1)
や
aff2=[[1,-1,0,0],[0,1,1,1]] (codim 2)
などの形式.  (テストは 2026-08-15-uv_contiguity.rr, func_contiguity.py を参照)