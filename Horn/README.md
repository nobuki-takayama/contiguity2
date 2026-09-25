

# Horn型超幾何系 同型分類パッケージ マニュアル Ver 2.

本パッケージは、Risa/Asir の $D$-加群計算機能と SageMath の多面体・整数計画計算機能を連携させ、Horn型超幾何系のパラメータ空間（$\mathbb{Z}^m$ の格子空間）全体における同型分類を自動遂行するための関数群を提供します。

## 1. Risa/Asir コア関数 (GKZ系を経由した Contiguity 導出)

GKZ 超幾何系から制限アルゴリズムを用いて、Horn 系の contiguity relation（隣接関係式）を計算します。

### `horn_contiguity2.set_A(A)`

* **概要:** 行列 $A$ に付随する GKZ 系と、その制限として得られる Horn 型超幾何方程式系を初期化し、必要な module 大域変数（`H_horn_str`, `H_a` など）を設定します。
* **引数:**
* `A`: 左側が単位行列となっている整数要素の行列（リストのリスト）。



### `horn_contiguity2.horn_contiguity(A_start, A_end)`

* **概要:** Horn 系のパラメータを `A_start` から `A_end` へシフトさせた際の contiguity relation（微分作用素）と超幾何 $b$-関数を計算します。
* **引数:**
* `A_start`: シフト前のパラメータのリスト（例: `[a, b, c]`）。退化パラメータ（例: `[p-1, p, q]`）も指定可能です。
* `A_end`: シフト後のパラメータのリスト（例: `[a+1, b, c]`）。


* **戻り値:** `[Numerator, Denominator, Info, UV, Common_factor]` の 5 要素リスト。
* `Numerator` / `Denominator`: 正規化された contiguity 作用素の分子と分母。
* `Info`: シフトの方向を示す文字列情報。
* `UV`: GKZ系での計算に用いられた最短シフトベクトル $u, v$ を用いた作用素。
* `Common_factor`: 可換 GCD によって抽出・約分された、余分なパラメータ因子（Torsion成分）。


* **特徴:** 内部で「最短シフト探索（ヒューリスティクス）」を行い、計算時間の爆発と冗長な因子の発生を防ぎます。

### `horn_contiguity2.horn_isom_sp(Od, New)`

* **概要:** Horn 系のパラメータ `Old` と `New` の間の同型が成立しない条件を出力
* **引数:**
* `Old`: シフト前のパラメータのリスト（例: `[a1, a2, a3]`）。退化パラメータ（例: `[a1-1, a1, a2]`）も指定可能です。
* `New`: シフト後のパラメータのリスト（例: `[a1+1, a2, a3]`）。


* **戻り値:** `[Numerator, Denominator, Info, UV, Common_factor]` の 5 要素リスト。

---

### 準備
このHornフォルダのファイルは tk_horn_contiguity2 フォルダの下へ
コピーしておく. tk_horn_contiguity2 フォルダは作業中のフォルダかまたは ASIRLOADPATH に含まれるフォルダにおいておく.

### 例 (Gauss ${}_2F_1$).
```C
//start openxm asir,  A は Gauss 超幾何の A.  
// F[4] が Gauss 2F1 の方程式. ここからparameter は [c-1,a,b] とわかる.
// a についての contiguity を求める.  
[3009] import("2026-09-25-horn-by-gkz.rr");
tk_uv_contiguity() [saito isom] is used. 2026-08-05.
[3107] A=[[1,0,0,-1],[0,1,0,1],[0,0,1,1]];  
[[1,0,0,-1],[0,1,0,1],[0,0,1,1]]
[3108] F=horn_contiguity2.set_A(A);
[[[a1,b1],[a2,-b2],[a3,-b3]],[[dx1,x4*dx4+a1],[dx2,-x4*dx4-a2],[dx3,-x4*dx4-a3]],[x4],[dx1*dx4-dx2*dx3],[(-x4*dx4-a3)(-x4*dx4-a2) - dx4(x4*dx4+a1)],[(x4^2-x4)*dx4^2+((a2+a3+1)*x4-a1-1)*dx4+a3*a2],[1*(-x4*dx4-a3)*(-x4*dx4-a2) - 1*dx4*(x4*dx4+a1)]]
[3109] horn_contiguity2.show_globals();
A matrix (H_A)=
[ 1 0 0 -1 ]
[ 0 1 0 1 ]
[ 0 0 1 1 ]
A is 3 * 4 matrix.
GKZ parameter (H_beta)=[b1,b2,b3]
GKZ x variables (H_gkz_xvar)=[x1,x2,x3,x4]
GKZ system (H_gkz)=[[-x4*dx4+x1*dx1-b1,x4*dx4+x2*dx2-b2,x4*dx4+x3*dx3-b3,-dx1*dx4+dx2*dx3],[x1,x2,x3,x4]]
------------------------
Horn system (H_horn_str)=
 (-x4*dx4-a3)(-x4*dx4-a2) - dx4(x4*dx4+a1)
Horn sytem is stored in H_horn
Horn sytem xvar(H_horn_xvar) = [x4]
Horn system parameter (H_a) = [a1,a2,a3]
a parameter expressed by beta (H_a_rule)=[[a1,b1],[a2,-b2],[a3,-b3]]
beta parameter expressed by a (H_b_rule)=[[b1,a1],[b2,-a2],[b3,-a3]]
Other globals:
   H_horn_str2=[1*(-x4*dx4-a3)*(-x4*dx4-a2) - 1*dx4*(x4*dx4+a1)]
   H_set_A_data = gkz_natural_params return value
   H_horn_gb (GB of Horn)
   Horn_contiguity
------------------------------

[3110] L=horn_contiguity2.horn_contiguity([c-1,a,b],[c-1,a+1,b]);
UV=[[0,0,0,0],[0,1,0,0]]
Executing tk_uv_contiguity_v1(UV=[[0,0,0,0],[0,1,0,0]],A=[[1,0,0,-1],[0,1,0,1],[0,0,1,1]],Beta=[c-1,-a,-b])
Computing M_chi with iteration for U=[0,0,0,0], V=[0,1,0,0]...
b-ideal B_chi is [1]
U=[0,0,0,0], V=[0,1,0,0], B_poly=1
Executing Saito (2001) Algorithm 4.2...
I_A is homogeneous. Using Block Grevlex order.
Resulting Symmetry Operator E:
dx2
Checking the result by check_symalg.
check_symalg_correct: PASSED
---------- Done checking the result of symalg.
Diff_nf=[0,1]
check_symalg_with_beta: PASSED
---------- Done checking the result of symalg_with_beta.
assoc(S_vars,Beta)=[[s1,c-1],[s2,-a],[s3,-b]]
horn_contiguity2: T0=[dx2,dx2,x1,x2,x3,x4], tk_uv_contiguity(UV=[[0,0,0,0],[0,1,0,0]],A=[[1,0,0,-1],[0,1,0,1],[0,0,1,1]],Base_Beta=[c-1,-a,-b]);
B_dxv=[[1,1],[dx2,1]]
Aop=[x4*dx4-x1*dx1+c-1,-x4*dx4-x2*dx2-a,-x4*dx4-x3*dx3-b]
lr_red1(L=dx2)
lr_red1_top(L=dx2)
lr_red1_top(dx2)-->-x4*dx4-a
lr_red1(L=0)
check_witness: Xvar=[x1,x2,x3], Xvar_all=[x1,x2,x3,x4]
 L[dx2](orig) = L1[-x4*dx4-a](reduced)+Witness[[ 0 1 0 ]] [dxi-Aop[i]] mod Xvar*D
 check_witness: OK.

[-x4*dx4-a,1,[[c-1,-a,-b], -> ,[c-1,-a-1,-b]],[dx2,dx2,1,dx2],1]


[3111] T2=horn_contiguity2.horn_isom_sp([c-1,a,b],[c-1,a+1,b]);
Chi=[0,-1,0], UV=[[0,0,0,0],[0,1,0,0]]
Executing tk_uv_contiguity_v2(UV=[[0,0,0,0],[0,1,0,0]],A=[[1,0,0,-1],[0,1,0,1],[0,0,1,1]],Beta=[c-1,-a,-b])
Computing M_chi with iteration for U=[0,0,0,0], V=[0,1,0,0]...
Generators of M_chi:
[1]
b-ideal B_chi is [1]
U=[0,0,0,0], V=[0,1,0,0], B_poly=1
Executing Saito (2001) Algorithm 4.2...
I_A is homogeneous. Using Block Grevlex order.

Chi=[0,1,0], UV=[[0,1,0,0],[0,0,0,0]]
Executing tk_uv_contiguity_v2(UV=[[0,1,0,0],[0,0,0,0]],A=[[1,0,0,-1],[0,1,0,1],[0,0,1,1]],Beta=[c-1,-a,-b])
Computing M_chi with iteration for U=[0,1,0,0], V=[0,0,0,0]...
Generators of M_chi:
[dx2,dx1*dx4]
b-ideal B_chi is [s2*s1+s2^2]
U=[0,1,0,0], V=[0,0,0,0], B_poly=s2*s1+s2^2
Executing Saito (2001) Algorithm 4.2...
I_A is homogeneous. Using Block Grevlex order.
Vd_perm=[x1,x3,x4,x2,dx1,dx3,dx4,dx2]
I_A=[(-1)*<<1,0,1,0>>+(1)*<<0,1,0,1>>]
I=0; assoc(S_vars,Beta)=[[s1,c-1],[s2,-a],[s3,-b]]
T0=[[a^2+(-c+1)*a,x2*x4*dx4+x1*x4*dx3+x2^2*dx2+x1*x2*dx1+x2,x1,x2,x3,x4]], tk_uv_contiguity([[0,1,0,0],[0,0,0,0]],[[1,0,0,-1],[0,1,0,1],[0,0,1,1]],[c-1,-a,-b]);
length(T0)=1 > 1 ?

[[a^2+(-c+1)*a],[[dx2,dx2,x1,x2,x3,x4]],[[a^2+(-c+1)*a,x2*x4*dx4+x1*x4*dx3+x2^2*dx2+x1*x2*dx1+x2,x1,x2,x3,x4]]]
[3112] fctr(T2[0][0]);
[[1,1],[a,1],[a-c+1,1]]

```

---

### 例 (Appell $F_1$).
```C
// parameter は [a,b,bp,c-1], a,b,$b'$ (bp), c は通常の Appell F1 のパラメータ.
[2701] A=[[1,0,0,0,1,1],
          [0,1,0,0,1,0],
          [0,0,1,0,0,1],
          [0,0,0,1,-1,-1]];;
[2703] F=horn_contiguity2.set_A(A);
[[[a1,-b1],[a2,-b2],[a3,-b3],[a4,b4]],[[dx1,-x6*dx6-x5*dx5-a1],[dx2,-x5*dx5-a2],[dx3,-x6*dx6-a3],[dx4,x6*dx6+x5*dx5+a4]],[x5,x6],[-dx2*dx6+dx3*dx5,dx4*dx6-dx1*dx3,dx4*dx5-dx1*dx2],[dx5(-x6*dx6-a3) - dx6(-x5*dx5-a2),(-x6*dx6-a3)(-x6*dx6-x5*dx5-a1) - dx6(x6*dx6+x5*dx5+a4),(-x5*dx5-a2)(-x6*dx6-x5*dx5-a1) - dx5(x6*dx6+x5*dx5+a4)],[((-x6+x5)*dx5+a2)*dx6-a3*dx5,(x6^2-x6)*dx6^2+((x5*x6-x5)*dx5+(a1+a3+1)*x6-a4-1)*dx6+a3*x5*dx5+a3*a1,((x5-1)*x6*dx5+a2*x6)*dx6+(x5^2-x5)*dx5^2+((a1+a2+1)*x5-a4-1)*dx5+a2*a1],[1*dx5*(-x6*dx6-a3) - 1*dx6*(-x5*dx5-a2),1*(-x6*dx6-a3)*(-x6*dx6-x5*dx5-a1) - 1*dx6*(x6*dx6+x5*dx5+a4),1*(-x5*dx5-a2)*(-x6*dx6-x5*dx5-a1) - 1*dx5*(x6*dx6+x5*dx5+a4)]]
[2704] F[4];
[dx5(-x6*dx6-a3) - dx6(-x5*dx5-a2),(-x6*dx6-a3)(-x6*dx6-x5*dx5-a1) - dx6(x6*dx6+x5*dx5+a4),(-x5*dx5-a2)(-x6*dx6-x5*dx5-a1) - dx5(x6*dx6+x5*dx5+a4)]
// 上記は Appell F1 の方程式. x5, x6 が変数.

```


---

## 2. SageMath 連携インターフェース (`func_contiguity.py`, `asir.sage` / `asir.py`)

Risa/Asir の計算エンジンを SageMath から透過的に呼び出し、数式オブジェクトとして安全に操作するためのラッパー関数群です。

### `set_horn_A(A_matrix)`, `set_A(A_matrix)`

* **概要:** 行列 $A$ を Risa/Asir に送信して初期化し、SageMath 側の多項式環（`PolynomialRing`）と変数辞書を自動生成します。
* **引数:**
* `A_matrix`: 初期化する行列 $A$（例: `[[1,0,0,-1], [0,1,0,1], [0,0,1,1]]`）。

<!--
* **戻り値:** `(ring, dic)` のタプル。
* `ring`: 生成された多項式環（例: `x4, dx4, a1, a2, a3` をジェネレータに持つ）。
* `dic`: 変数名（文字列）からジェネレータへの変換辞書。
-->

### `show_globals()`

* **概要:** 現在セットされている Horn 型超幾何方程式の定義方程式の情報を文字列で表示.

<!--
### `get_horn_system(ring, dic)`

* **概要:** 現在セットされている Horn 型超幾何方程式の定義方程式を取得します。
* **引数:** `set_horn_A` で取得した `ring` と `dic`。
* **戻り値:** SageMath の Symbolic Ring (`SR`) オブジェクトのリスト。
* **特徴:** 完全に展開された多項式ではなく、因数分解された構造（例: `(dx4*x4 + a2)*(dx4*x4 + a3) - ...`）を保ったまま保持されるため、パラメータの代入や解析が容易です。

### `horn_contiguity_by_GKZ(old, new, ring, dic)`

* **概要:** `horn_contiguity_from_a_shift` を SageMath から安全に呼び出し、結果を Sage の多項式/有理式オブジェクトとして返します。
* **引数:** `old` (シフト前リスト), `new` (シフト後リスト), `ring`, `dic`。
* **戻り値:** `[Numerator, Denominator, Info, UV, Common_factor, ring]` のリスト。各要素は Sage のオブジェクトとしてパース済みです。
-->

---

## 3. 同型分類の自動遂行関数 (`get_isom_class.py`)

パラメータ空間の次元を再帰的に落としながら、超幾何系の同型分類（セル分解）を実行します。

### `representatives_of_arrangement(L,hg='horn_contiguity2.horn',ring,dic,param)`

* **概要:** 与えられた格子空間上で contiguity relation を網羅的に計算し、$b$-関数の一次因子を「壁（超平面）」とみなして空間を分割、代表元（representative）を抽出します。
* **引数:**
* `L`: 現在探索している affine 空間の定義式
* `hg`: Asirを呼び出して contiguity を計算する関数（例: `gauss_contiguity`, `f1_contiguity`）。
* `ring`: 環
* `dic`: パラメータ変数とその文字列表現の対応辞書。
* `param`: パラメータ変数


* **戻り値:** 再帰的なリスト構造を持つ分類ログ。各次元（セル）における壁の方程式、計算された作用素、到達した次元（`0-dim face` など）の情報が格納されます。

### `merge_isom_class(data,func='horn_contiguity2.horn',ring,dic,param)`

* **概要:** data の中の同型な cell を検出してまとめる
* **引数:**
* `data`: representtives_of_arrangement の出力.


---

## 4. 基本的な実行ワークフロー (SageMath)

<!--Prog: 2026-09-25-isom-horn.py -->

```python
from sage.all import *
from chamber_lattice_utils import *
from func_contiguity import *
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

plot_arrangement_faces_3d(data_merged[0],-4,4,-4,4,-4,4)

```

---
## 5. サンプル入力ファイル
* 2026-09-24-isom-horn.rr,  A を与えて計算.
* 2026-09-09-isom-c1f1.rr,  1F1
* 2026-09-08-f0134.rr,  non-Cohen-Macaulay example
* 2026-09-08-isom-banana1-no-d.py

## 6. ファイル

## SageMath toplevel
* [chamber_lattice_utils.py](chamber_lattice_utils.py) 
* [func_contiguity.py](func_contiguity.py)  
* [get_isom_class.py](get_isom_class.py) (<2026-09-09-rep.py>)  分類関数本体

## Risa/Asir, 計算アルゴリズムの実装本体
* <2026-08-05-contiguity.rr>
* <2026-08-15-uv_contiguity.rr>
* <tk_weyl.rr>
* <symeq_one_loop.rr>
* <linear_poly.rr>

## Horn 型の contiguity. 一般の A 用 interface
* <2026-09-24-chg-by-gkz.rr>
* <2026-09-25-horn-by-gkz.rr>  

## 1F1, f0134, banana 専用の horn-by-gkz
* <2026-08-12-1f1-by-gkz.rr>
* <2026-08-07-banana-one-loop-by-gkz.rr>
* <2026-08-28-f0134.rr>

## サンプル入力
* [2026-09-25-isom-horn.py](2026-09-25-isom-horn.py)  一般の A
* [2026-09-08-isom-banana1-no-d.py](2026-09-08-isom-banana1-no-d.py)
* [2026-09-08-isom-f0134.py](2026-09-08-isom-f0134.py)
* [2026-09-09-isom-c1f1.py](2026-09-09-isom-c1f1.py)
 

## 7. 参考文献
* <https://github.com/nobuki-takayama/contiguity> : ver 1.
* <https://arxiv.org/abs/2510.05104> : Algorithm についての論文
* [man-chamber-lattice-utils.md](man-chamber-lattice-utils.md) : chamber_lattice_utils.py のマニュアル
