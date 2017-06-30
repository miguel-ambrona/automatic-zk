#!/usr/bin/env sage

import sys
from sage.all import *

var('g h')

def look_for_non_zero_in_col(M,row,col):
    for i in range(row,len(M.rows())):
        if not M[i][col] == 0:  return i

def gauss_elim(A,b):
    row_cnt = 0
    A_b = A.augment(b)
    for j in range(len(A.columns())):
        row = look_for_non_zero_in_col(A_b,row_cnt,j)
        if row == None:  continue
        A_b.swap_rows(row,row_cnt)
        A_b = A_b.with_row_set_to_multiple_of_row(row_cnt,row_cnt,1/A_b[row_cnt][j])
        for i in range(len(A_b.rows())):
            if i == row_cnt:  continue
            if A_b[i][j] != 0:
                A_b.add_multiple_of_row(i,row_cnt,-A_b[i][j])
        row_cnt += 1
    return A_b

def split_if_DLOG_is_hard(variables, eqs):
    new_vars = []
    for v in variables:
        vg = eval("var('" + str(v) + "_g')")
        vh = eval("var('" + str(v) + "_h')")
        eqs2 = []
        for eq in eqs:
            new_eq = eval("eq.subs(" + str(v) + "=" + str(vg*g+vh*h) + ")")
            eqs2.append(new_eq)
        eqs = eqs2
        new_vars += [vg, vh]
    return [eq.coefficient(g) for eq in eqs] + [eq.coefficient(h) for eq in eqs], new_vars

def create_system(variables, eqs):
    A = []
    b = []
    for eq in eqs:
        residue = eq
        row = []
        for v in variables:
            t = eval("residue.coefficient(" + str(v) + ")")
            residue = residue - v*t
            row.append(t)
        A.append(row)
        b.append(-residue)
    return matrix(A), vector(b)

def ideal_fun(xs,fs,ideal_terms, P):
    X = PolynomialRing(QQ,xs)
    R = PolynomialRing(X,fs)
    exec("(" + xs + ") = X.gens()")
    exec("(" + fs + ") = R.gens()")

    loc_dict = "{"
    for v in (xs + "," + fs).split(","):
        loc_dict += "'" + v + "':" + v + ","
    loc_dict = eval(loc_dict[:-1] + "}")

    ideal_eqs = [sage_eval(preparse(str(t)), locals=loc_dict) for t in ideal_terms]
    polys = [sage_eval(preparse(str(p)), locals=loc_dict) for p in P]

    I = ideal(ideal_eqs[1])
    
    print
    print I
    B = I.groebner_basis()
    print
    for p in polys:
        print p.reduce(B) == 0

    print "Are the same ideal?"
    print I == ideal(p)
        
def example1():
    xs = 'x,xp,xpp'
    fs = 'f1,f2,f3,f1p,f2p,f3p,f1pp,f2pp,f3pp'
    var('a1,b1,b2,b3')
    var(xs)
    var(fs)
    
    verif_eqs = [ f1*g+f2*h - a1*x - b1, f1*(x-f1)*g + f3*h - b2*x - b3 ] + \
                [ f1p*g+f2p*h - a1*xp - b1, f1p*(xp-f1p)*g + f3p*h - b2*xp - b3 ] + \
                [ f1pp*g+f2pp*h - a1*xpp - b1, f1pp*(xpp-f1pp)*g + f3pp*h - b2*xpp - b3 ]

    variables = [a1,b1,b2,b3]
    new_equations, new_variables = split_if_DLOG_is_hard(variables, verif_eqs)

    A, b = create_system(new_variables, new_equations)
    print A
    print b
    print
    S = gauss_elim(A,b)
    new_A = transpose(transpose(S)[:-1])
    new_b = vector([a.full_simplify() for a in transpose(S)[-1]])

    print new_A
    for a in new_b:  print a
    ideal_terms = []
    for i in range(len(new_A.rows())):
        all_zeros = True
        for j in range(len(new_A.columns())):
            if not new_A[i][j] == 0:  all_zeros = False; break
        if all_zeros:
            t = new_b[i]
            ideal_terms.append(t)

            
    a = new_b[0]
    P = (a*(1-a))

    ideal_fun(xs,fs,ideal_terms,[P])

def example2():
    xs = 'x,xp,xpp,xp3,xp4,xp5'
    fs = 'f1,f1p,f1pp,f1p3,f1p4,f1p5,f2,f2p,f2pp,f2p3,f2p4,f2p5,f3,f3p,f3pp,f3p3,f3p4,f3p5'
    var('a1 a2 b1 m1 m0 mm1')
    var(xs)
    var(fs)

    verif_eqs = [ f1*g+f2*h - a1/x - b1 - x*a2, f1*(x+1/x-f1)*g + f3*h - mm1/x - m0 - m1*x ] + \
                [ f1p*g+f2p*h - a1/xp - b1 - xp*a2, f1p*(xp+1/xp-f1p)*g + f3p*h - mm1/xp - m0 - m1*xp ] + \
                [ f1pp*g+f2pp*h - a1/xpp - b1 - xpp*a2, f1pp*(xpp+1/xpp-f1pp)*g + f3pp*h - mm1/xpp - m0 - m1*xpp ] + \
                [ f1p3*g+f2p3*h - a1/xp3 - b1 - xp3*a2, f1p3*(xp3+1/xp3-f1p3)*g + f3p3*h - mm1/xp3 - m0 - m1*xp3 ] + \
                [ f1p4*g+f2p4*h - a1/xp4 - b1 - xp4*a2, f1p4*(xp4+1/xp4-f1p4)*g + f3p4*h - mm1/xp4 - m0 - m1*xp4 ] + \
                [ f1p5*g+f2p5*h - a1/xp5 - b1 - xp5*a2, f1p5*(xp5+1/xp5-f1p5)*g + f3p5*h - mm1/xp5 - m0 - m1*xp5 ]

    variables = [a1,a2,b1,m1,m0,mm1]
    new_equations, new_variables = split_if_DLOG_is_hard(variables, verif_eqs)

    A, b = create_system(new_variables, new_equations)
    print A
    print b
    print
    S = gauss_elim(A,b)
    new_A = transpose(transpose(S)[:-1])
    new_b = vector([a.full_simplify() for a in transpose(S)[-1]])

    print new_A
    for a in new_b:  print a
    ideal_terms = []
    for i in range(len(new_A.rows())):
        all_zeros = True
        for j in range(len(new_A.columns())):
            if not new_A[i][j] == 0:  all_zeros = False; break
        if all_zeros:
            t = new_b[i]
            ideal_terms.append(t)

    print(len(ideal_terms))
    a1v = new_b[0]
    P1 = (a1v*(1-a1v))
    a2v = new_b[2]
    P2 = (a2v*(1-a2v))

    ideal_fun(xs,fs,ideal_terms,[P1,P2])

def example3():
    xs = 'x,xp,xpp,xp3,xp4,xp5'
    fs = 'f1,f1p,f1pp,f1p3,f1p4,f1p5,f2,f2p,f2pp,f2p3,f2p4,f2p5,f3,f3p,f3pp,f3p3,f3p4,f3p5'
    var('a a1 a2 b1 m1 m0 mm1')
    var(xs)
    var(fs)
    
    verif_eqs = [ f1*g+f2*h - a1/x - b1 - x*a2, f1*(x+1/x-f1)*g + f3*h - mm1/x - m0 - m1*x, a-2*a2-a1 ] + \
                [ f1p*g+f2p*h - a1/xp - b1 - xp*a2, f1p*(xp+1/xp-f1p)*g + f3p*h - mm1/xp - m0 - m1*xp ] + \
                [ f1pp*g+f2pp*h - a1/xpp - b1 - xpp*a2, f1pp*(xpp+1/xpp-f1pp)*g + f3pp*h - mm1/xpp - m0 - m1*xpp ] + \
                [ f1p3*g+f2p3*h - a1/xp3 - b1 - xp3*a2, f1p3*(xp3+1/xp3-f1p3)*g + f3p3*h - mm1/xp3 - m0 - m1*xp3 ] + \
                [ f1p4*g+f2p4*h - a1/xp4 - b1 - xp4*a2, f1p4*(xp4+1/xp4-f1p4)*g + f3p4*h - mm1/xp4 - m0 - m1*xp4 ] + \
                [ f1p5*g+f2p5*h - a1/xp5 - b1 - xp5*a2, f1p5*(xp5+1/xp5-f1p5)*g + f3p5*h - mm1/xp5 - m0 - m1*xp5 ]

    variables = [a,a1,a2,b1,m1,m0,mm1]
    new_equations, new_variables = split_if_DLOG_is_hard(variables, verif_eqs)

    A, b = create_system(new_variables, new_equations)
    print A
    print b
    print
    S = gauss_elim(A,b)
    new_A = transpose(transpose(S)[:-1])
    new_b = vector([t.full_simplify() for t in transpose(S)[-1]])

    print new_A
    for t in new_b:  print t
    ideal_terms = []
    for i in range(len(new_A.rows())):
        all_zeros = True
        for j in range(len(new_A.columns())):
            if not new_A[i][j] == 0:  all_zeros = False; break
        if all_zeros:
            t = new_b[i]
            ideal_terms.append(t)

    print(len(ideal_terms))
    av = new_b[0]
    P1 = av*(1-av)*(2-av)*(3-av)
    P2 = av
    P3 = av*(1-av)*(2-av)

    ideal_fun(xs,fs,ideal_terms,[P1,P2,P3])
    
def main():
    example3()

if __name__ == "__main__":
    main()
