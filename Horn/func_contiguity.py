## Example: L = func_contiguity([a+1,c], [a,c], ring=R, dic=dic,func='c1f1')
##              func_isom(chi=[-1,0],func='c1f1,dic=dic)
# by gemini, 2026.03.07.  Modified 2026.08.08 for banana*. For 1F1, 2026.08.24
from sage.misc.sage_eval import sage_eval
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.rational_field import QQ
from sage.symbolic.ring import SR
from sage.all import load

def load_func_contiguity(file_name):
    global Asir_running
    try:
        print(Asir_running)
    except NameError:
        # SageMathの組み込みasirを利用するか、外部のasir.pyを利用する設定
        load('./asir.py') 
        Asir_running = 'asir is already running'
    
    # Risa/Asir側の最新ファイルをロード
    asir.eval('load("'+file_name+'");')

def contiguity_verbose(level):
    asir.eval(f'Ans_tmp = contiguity_verbose({level});');
    result=asir.eval('Ans_tmp')
    id = sage_eval(result)
    return id

## 2026.09.05
# テンプレートを記憶するためのモジュールレベル変数（初期値は 0）
_GLOBAL_PARAM_TEMPLATE = 0

def set_constant_part_of_param_vector(template):
    """
    パラメータベクトルの定数部分を含むテンプレートを設定・記憶する。
    
    入力:
        template: 定数と None (変数を埋め込む場所) が混在したリスト。
                  (例: [1, 2, None, None, 3])
                  0 を指定すると初期状態に戻る。
    """
    global _GLOBAL_PARAM_TEMPLATE
    _GLOBAL_PARAM_TEMPLATE = template

def get_constant_part_of_param_vector():
    global _GLOBAL_PARAM_TEMPLATE
    return _GLOBAL_PARAM_TEMPLATE

def len_global_param_template():
    global _GLOBAL_PARAM_TEMPLATE
    if _GLOBAL_PARAM_TEMPLATE==0:
        return 0
    else:
        return len(_GLOBAL_PARAM_TEMPLATE)

def extention_position():
    global _GLOBAL_PARAM_TEMPLATE
    if _GLOBAL_PARAM_TEMPLATE == 0:
        return 0
    for i,item in enumerate(_GLOBAL_PARAM_TEMPLATE):
        if item is None:
            return i
    return -1

def linear_space_by_extension_param():
    global _GLOBAL_PARAM_TEMPLATE
    if _GLOBAL_PARAM_TEMPLATE == 0:
        return []
    dim=len(_GLOBAL_PARAM_TEMPLATE)
    L=[]
    for i,item in enumerate(_GLOBAL_PARAM_TEMPLATE):
        if item is None:
            continue
        eq=[-1 if j==i+1 else 0 for j in range(dim+1)]
        eq[0]=item
        L.append(eq)
    return L

    

def extend_param_vector_by_constant_part(values):
    """
    記憶されたテンプレートの None の位置に、与えられた values の要素を順に埋め込んだリストを返す。
    テンプレートが 0 の場合は values をそのまま返す。
    
    入力:
        values: テンプレートの None の部分に埋め込む変数のリスト [a, b, ...]
    出力:
        定数と変数が合成された新しいリスト
    """
    global _GLOBAL_PARAM_TEMPLATE
    
    # テンプレートが 0 の場合は、入力された values をそのまま返す
    if _GLOBAL_PARAM_TEMPLATE == 0:
        return list(values)
        
    result = []
    val_idx = 0
    
    for item in _GLOBAL_PARAM_TEMPLATE:
        if item is None:
            # None の部分に values の要素を順番に割り当てる
            if val_idx < len(values):
                result.append(values[val_idx])
                val_idx += 1
            else:
                raise ValueError(f"提供された values の要素数 ({len(values)}) が、テンプレート内の None の数より不足しています。")
        else:
            # None でない定数部分はそのまま引き継ぐ
            result.append(item)
            
    return result

def func_isom(chi,func,dic):
    print('Warning. func_isom is obsoleted. Do not use it.')
    print('func_isom calls: ', f'Ans_tmp = {func}_isom({chi});');
    asir.eval(f'Ans_tmp = {func}_isom({chi});');
    result=asir.eval('Ans_tmp')
    id = sage_eval(result, locals=dic)
    return id

def parse_contiguity(str,ring,dic):
    result=asir.eval(f'Ans_tmp={str};')
    #print('type(result)=', type(result))
    if result=='0':
        print(f'Warning: str={str} is invalid. Return []')
        return []  # contiguity is not found.
    
    # 各要素を個別に文字列として取得
    c0_str = asir.eval(f'Ans_tmp[0]')
    b0_str = asir.eval(f'Ans_tmp[1]')
    info_str = asir.eval(f'Ans_tmp[2]')
    uv_str = asir.eval(f'Ans_tmp[3]')
    
    # 要素数を確認して Common_factor を取得
    len_ans = int(asir.eval(f'length(Ans_tmp)'))
    if len_ans >= 5:
        cf_str = asir.eval(f'Ans_tmp[4]')
    else:
        cf_str = '1'
    
    # Sage の多項式/有理式に変換
    C0 = sage_eval(c0_str, locals=dic)
    B0 = sage_eval(b0_str, locals=dic)
    CF = sage_eval(cf_str, locals=dic) if cf_str != '1' else ring(1)
    
    # 情報リストはそのまま文字列として保持
    return [C0, B0, info_str, uv_str, CF, ring]

def func_isom_sp(old, new, ring, dic,func):
    global _GLOBAL_PARAM_TEMPLATE
    #print('old=',old,' new=',new)
    """
    パラメータのシフトからcontiguityを計算し、Sageのオブジェクトとして返す。
    old, new: [a1, a2, ...] のようなリスト形式
    """
    if len(old) < len_global_param_template():
        old_orig=old
        old= extend_param_vector_by_constant_part(old)
        new= extend_param_vector_by_constant_part(new)
        print('func_isom_sp: from {old_orig} to {old}, ...')
    oo = str(old).replace("'", "")
    nn = str(new).replace("'", "")
    
    # Asir側で計算し、一時変数 Ans_tmp に格納する
    # ※直接結果文字列を受け取るとカンマ分割が難しいため
    print(f'Ans_tmp = {func}_isom_sp({oo}, {nn});') # for debug
    asir.eval(f'Ans_tmp = {func}_isom_sp({oo}, {nn});')
    c_up=asir.eval('sprintf("%a",Ans_tmp[1])')
    c_down=asir.eval('sprintf("%a",Ans_tmp[2])') # string として受け取る. 必要なら parse_contiguity
    result=asir.eval('Ans_tmp[0]')
    id = sage_eval(result, locals=dic)
    return [id,c_up,c_down]
    

def func_contiguity(old, new, ring, dic,func):
    #print('old=',old,' new=',new)
    if len(old) < len_global_param_template():
        old_orig=old
        old= extend_param_vector_by_constant_part(old)
        new= extend_param_vector_by_constant_part(new)
        print('func_contiguity: from {old_orig} to {old}, ...')
    """
    パラメータのシフトからcontiguityを計算し、Sageのオブジェクトとして返す。
    old, new: [a1, a2, ...] のようなリスト形式
    """
    oo = str(old).replace("'", "")
    nn = str(new).replace("'", "")
    
    # Asir側で計算し、一時変数 Ans_tmp に格納する
    # ※直接結果文字列を受け取るとカンマ分割が難しいため
    asir.eval(f'Ans_tmp = {func}_contiguity({oo}, {nn});')
    result=asir.eval('Ans_tmp')
    if result=='0':
        print('Warning: contiguity is not found.')
        return []  # contiguity is not found.
    #print('type(result) is ',type(result))
    #return parse_contiguity(result,ring,dic)  #Todo, info_str で syntax error
    
    # 各要素を個別に文字列として取得
    c0_str = asir.eval('Ans_tmp[0]')
    b0_str = asir.eval('Ans_tmp[1]')
    info_str = asir.eval('Ans_tmp[2]')
    uv_str = asir.eval('Ans_tmp[3]')
    
    # 要素数を確認して Common_factor を取得
    len_ans = int(asir.eval('length(Ans_tmp)'))
    if len_ans >= 5:
        cf_str = asir.eval('Ans_tmp[4]')
    else:
        cf_str = '1'
    
    # Sage の多項式/有理式に変換
    C0 = sage_eval(c0_str, locals=dic)
    B0 = sage_eval(b0_str, locals=dic)
    CF = sage_eval(cf_str, locals=dic) if cf_str != '1' else ring(1)
    
    # 情報リストはそのまま文字列として保持
    return [C0, B0, info_str, uv_str, CF, ring]

# option b_ideal=1 の時用.
def func_contiguity_b_ideal(old, new, ring, dic,func):
    #print('old=',old,' new=',new)
    if len(old) < len_global_param_template():
        old_orig=old
        old= extend_param_vector_by_constant_part(old)
        new= extend_param_vector_by_constant_part(new)
        print('func_contiguity_b_ideal: from {old_orig} to {old}, ...')
    oo = str(old).replace("'", "")
    nn = str(new).replace("'", "")
    
    # Asir側で計算し、一時変数 Ans_tmp に格納する
    # ※直接結果文字列を受け取るとカンマ分割が難しいため
    asir.eval(f'Ans_tmp = {func}_contiguity({oo}, {nn} | b_ideal=1);')
    size=sage_eval(asir.eval('length(Ans_tmp)'))
    if size=='0':
        print('Warning: contiguity is not found.')
        return []  # contiguity is not found.
    asir.eval('Ans_tmp_orig=Ans_tmp;')
    ans=[]
    for i in range(size):
        print(f'Ans_tmp=Ans_tmp_orig[{i}];')
        asir.eval(f'Ans_tmp=Ans_tmp_orig[{i}];')
        # 各要素を個別に文字列として取得
        c0_str = asir.eval('Ans_tmp[0]')
        b0_str = asir.eval('Ans_tmp[1]')
        info_str = asir.eval('Ans_tmp[2]')
        uv_str = asir.eval('Ans_tmp[3]')
    
        # 要素数を確認して Common_factor を取得
        len_ans = int(asir.eval('length(Ans_tmp)'))
        if len_ans >= 5:
            cf_str = asir.eval('Ans_tmp[4]')
        else:
            cf_str = '1'
        # Sage の多項式/有理式に変換
        C0 = sage_eval(c0_str, locals=dic)
        B0 = sage_eval(b0_str, locals=dic)
        CF = sage_eval(cf_str, locals=dic) if cf_str != '1' else ring(1)
        # 情報リストはそのまま文字列として保持
        ans.append([C0, B0, info_str, uv_str, CF, ring])
    return ans

def func_get_chi(old,new,func):
    print('func_get_chi calls: ', f'Ans_tmp = {func}_get_chi({old},{new});');
    asir.eval(f'Ans_tmp = {func}_get_chi({old},{new});');
    result=asir.eval('Ans_tmp')
    print(result)
    id = sage_eval(result, locals=dic)
    return id

def bf_factor_to_affine_space(param_expr,bfac):
    if len(param_expr) < len_global_param_template():
        param_expr_orig=param_expr
        param_expr= extend_param_vector_by_constant_part(param_expr) 
        print('bf_factor_to_affine_spaces: from {param_expr_orig} to {param_expr}, ...')
    result=asir.eval(f'Asir_tmp=bf_factor_to_affine_space({param_expr},{bfac})')
    result=asir.eval('Asir_tmp')
    return sage_eval(f'{result}')

"""
load_func_contiguity('2026-08-12-1f1-by-gkz.rr')
R1 = PolynomialRing(QQ,'x4,x5,x6,dx4,dx5,dx6,d,nu1,nu2')
dic1 = R1.gens_dict()
d,nu1,nu2 = dic1['d'], dic1['nu1'], dic1['nu2']
print(bf_factor_to_affine_space([2*nu1+2*1,nu1,1],nu1-3))
## 結果 [[-1, 0, 0, 1], [-3, 0, 1, 0], [-8, 1, 0, 0]]
"""

def b_ideal_to_affine_spaces(param_expr,b_ideal):
    if len(param_expr) < len_global_param_template():
        param_expr_orig=param_expr
        param_expr= extend_param_vector_by_constant_part(param_expr) 
        print('b_ideal_to_affine_spaces: from {param_expr_orig} to {param_expr}, ...')
    # 改行コード(\n)をすべて消去して1行にする. 勝手な改行を抑制. end-of-file error がある時はこれを疑え.
#    param_expr = str(param_expr).replace('\n', ' ')
#    b_ideal = str(b_ideal).replace('\n', ' ')
    print(f'b_ideal_to_affine_spaces({param_expr},{b_ideal})')
    result=asir.eval(f'Asir_tmp=b_ideal_to_affine_spaces({param_expr},{b_ideal});')
    result=asir.eval('Asir_tmp')
    return sage_eval(f'{result}')

def set_horn_A(Amat):
    return set_A(Amat)

def set_A(Amat):
    if type(Amat)!=type([]):
        Amat=[list(row) for row in Amat]
    result=asir.eval(f'Asir_tmp=horn_contiguity2.set_A({Amat});')
    result=asir.eval('Asir_tmp')
    return result
def show_globals():
    result=asir.eval(f'Asir_tmp=horn_contiguity2.show_globals();')
    result=asir.eval('Asir_tmp')
    return result
    

"""
load("func_contiguity.py")
load_func_contiguity('2026-08-12-1f1-by-gkz.rr')
R3 = PolynomialRing(QQ, 'b1,b2')
b1,b2 = R3.gens()
ideal3=[256*b1^5+1152*b1^4+(-1440*b2^2+4896*b2-2128)*b1^3+(-2160*b2^3+7128*b2^2-5448*b2+1032)*b1^2+(-1215*b2^4+4374*b2^3-4581*b2^2+1638*b2-144)*b1-243*b2^5+972*b2^4-1269*b2^3+648*b2^2-108*b2,(-256*b2+512)*b1^4+(-768*b2^2+1920*b2-768)*b1^3+(-864*b2^3+2592*b2^2-1904*b2+352)*b1^2+(-432*b2^4+1512*b2^3-1560*b2^2+552*b2-48)*b1-81*b2^5+324*b2^4-423*b2^3+216*b2^2-36*b2]
T1=b_ideal_to_affine_spaces([b1,b2],ideal3)
"""
    

"""
# 1. モジュールのロード
load_func_contiguity('2026-08-12-1f1-by-gkz.rr')

R = PolynomialRing(QQ,'x3,dx3,a,c')
dic = R.gens_dict()
a,c = dic['a'], dic['c']


# 3. Contiguity の計算 (例: a -> a + 1)
Chi = func_get_chi([a,c],[a+1,c],func='c1f1')
I = func_isom(Chi,func='c1f1',dic=dic)
L1 = func_contiguity([a,c], [a+1,c], ring=R, dic=dic,func='c1f1')
L2 = func_contiguity([a+1,c], [a,c], ring=R, dic=dic,func='c1f1')
# 2026.08.27 sp は specific parameter
L3 = func_isom_sp([a,c],[a+1,c], ring=R, dic=dic,func='c1f1')


# 結果の表示
print("chi:",Chi)
print("B-ideal:", I)
print("分子 (Numerator):", L2[0])
print("分母 (Denominator):", L2[1])
print("共通因子 (Common Factor):", L2[4])

print('L3=',L3)
R_gkz = PolynomialRing(QQ,'x1,x2,x3,dx1,dx2,dx3,a,c')
dic_gkz = R_gkz.gens_dict()
a,c = dic_gkz['a'], dic_gkz['c']
L4=parse_contiguity(L3[1],ring=R_gkz,dic=dic_gkz) #todo 最初の二つのみでOK.
print('L4=',L4)
"""




"""
# banana_one_loop case.  sage -n jupyterlab
load_func_contiguity('2026-08-07-banana-one-loop-by-gkz.rr')
contiguity_verbose(0)
R1 = PolynomialRing(QQ,'x4,x5,x6,dx4,dx5,dx6,d,nu1,nu2')
dic1 = R1.gens_dict()
d,nu1,nu2 = dic1['d'], dic1['nu1'], dic1['nu2']
Chi = func_get_chi([d,nu1,nu2],[d,nu1+1,nu2],func='banana1')
#II = func_isom(Chi,func='banana1',dic=dic1)
II = func_isom_sp([d,nu1,nu2],[d,nu1+1,nu2],func='banana1',ring=R1,dic=dic1)
LL2 = func_contiguity([d,nu1+1,nu2], [d,nu1,nu2], ring=R1, dic=dic1,func='banana1')
factor(II[0][0])
factor(LL2[1])

"""


"""
# f0134 のテスト.
load_func_contiguity('2026-08-28-f0134.rr')
contiguity_verbose(0)
R2 = PolynomialRing(QQ,'x3,x4,dx3,dx4,b1,b2')
dic2 = R2.gens_dict()
b1,b2 = dic2['b1'], dic2['b2']
II = func_isom_sp([b1,b2],[b1+1,b2],func='f0134',ring=R2,dic=dic2)
LL2 = func_contiguity([b1,b2], [b1+1,b2], ring=R2, dic=dic2,func='f0134')
LL3= func_contiguity_b_ideal([b1,b2], [b1+1,b2], ring=R2, dic=dic2,func='f0134')
"""

# 2026.09.05
"""
load_func_contiguity('2026-08-12-1f1-by-gkz.rr')
R = PolynomialRing(QQ,'x3,dx3,a,c')
dic = R.gens_dict()
a,c = dic['a'], dic['c']
set_constant_part_of_param_vector([1,None])
T1=func_isom_sp([c],[c+1],func='c1f1',ring=R,dic=dic)
T2=func_contiguity([c],[c+1],func='c1f1',ring=R,dic=dic)
T3=func_contiguity([c+1],[c],func='c1f1',ring=R,dic=dic)
"""

