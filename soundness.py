#!/usr/bin/env sage

import sys
from sage.all import *

def varlist(name,n):
    return [name+str(k+1) for k in range(n)]

def list2string(l):
    output = ""
    for a in l:
        output += str(a)+","
    return output[:-1]

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

def create_system(eqs):
    variables = []
    for eq in eqs:
        for v in eq.variables():
            if str(v)[0] in ["a","b"] and not v in variables:  variables.append(v)

#    print variables
    A = []
    b = []
    for eq in eqs:
        residue = eq
        row = []
        for v in variables:
            t = eq.coefficient(v)
            residue = residue - v*t
            row.append(t)
        A.append(row)
        b.append(-residue)
    return matrix(A), vector(b), variables

def analyze(fname):
    f = open(fname,'r')
    lines = f.readlines()
    f.close()
    for i in range(len(lines)):
        if "equations" in lines[i]:  break
        exec(lines[i])

    gh = ["g","h"]
    avs = varlist("ag",nA) + varlist("ah",nA)
    bvs = varlist("b",nB) + varlist("bg",nB) + varlist("bh",nB)
    fs = varlist("f",nf)
    fvs = []
    for k in range(rep):
        fvs += [str(a)+"p"+str(k+1) for a in fs]
    xvs = varlist("x",rep)

    X = PolynomialRing(QQ,xvs)
    R = PolynomialRing(QQ,gh+avs+bvs+fvs)
    exec(list2string(xvs) + " = X.gens()")
    exec(list2string(gh+avs+bvs+fvs) + " = R.gens()")

    # Randomized x's:
    for x in xvs:
        exec(x + " = randint(1,100000)")

    local_dict = {}
    for a in gh+avs+bvs+fvs+xvs:
        exec("local_dict['" + a + "'] = " + a)

    equations = []
    for j in range(i+1,len(lines)):
        line = preparse(lines[j].strip())
        if len(line) == 0:  continue
        if "Polynomial" in line:  break
        for k in range(nA,0,-1):
            line = line.replace("a"+str(k), "(g*ag"+str(k)+"+h*ah"+str(k)+")")
        for k in range(nB,0,-1):
            line = line.replace("b"+str(k), "(g*bg"+str(k)+"+h*bh"+str(k)+")")
        for k in range(rep):
            this_line = line.replace("x","x"+str(k+1))
            for a in fs:
                this_line = this_line.replace(a,a.upper()+"p"+str(k+1))
            eq = sage_eval(this_line.replace("F","f"), locals=local_dict)
            equations.append(eq.coefficient(g))
            equations.append(eq.coefficient(h))

    polys = []
    for i in range(j+1,len(lines)):
        polys.append(sage_eval(preparse(lines[i].strip()), locals=local_dict))

    A, b, variables = create_system(equations)
    S = gauss_elim(A,b)
    new_A = transpose(transpose(S)[:-1])
    new_b = vector(transpose(S)[-1])

    vdic = {}
    for k in range(len(variables)):
        vdic[str(variables[k])] = new_b[k]

#    for a in new_b:  print a

    ideal_terms = []
    for i in range(len(new_A.rows())):
        all_zeros = True
        for j in range(len(new_A.columns())):
            if not new_A[i][j] == 0:  all_zeros = False; break
        if all_zeros:
            t = new_b[i]
            ideal_terms.append(t)

    I = ideal(ideal_terms[0])
    B = I.groebner_basis()

    for v in variables:
        aux_polys = []
        for P in polys:
            exec("P = P.subs(" + str(v) + " = vdic['" + str(v) + "'])")
            aux_polys.append(P)
        polys = aux_polys

    for P in polys:
        print P.reduce(B) == 0


if __name__ == '__main__':
    if len(sys.argv) >= 2:
        analyze(sys.argv[1])
