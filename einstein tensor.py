# -*- coding: utf-8 -*-
"""
Created on Thu Apr  9 00:55:44 2020

@author: Matthew Jankowski
"""


from sympy import *

init_printing(use_latex='mathjax')

t, r, theta, phi = symbols('t r theta phi')

dims = {t:0, r:1, theta:2, phi:3}
dims_list = [t, r, theta, phi]

f = Function("f")(r)

g = Function("g")(r)

metric = [[-exp(2*f), 0, 0, 0],
          [0, exp(2*g), 0, 0],
          [0, 0, r**2, 0],
          [0, 0, 0, (r**2)*(sin(theta)**2)]]

metric_inv = [[-1*(exp(-2*f)), 0, 0, 0],
              [0, (exp(-2*g)), 0, 0],
              [0, 0, r**-2, 0],
              [0, 0, 0, (r**-2)*(sin(theta)**-2)]]

def get_christoffel(a, b, c):
    christoffel = 0
    for i in range(4):
        christoffel += (.5)*metric_inv[i][a]*(diff(metric[i][b], dims_list[c]) +
                                              diff(metric[i][c], dims_list[b]) -
                                              diff(metric[b][c], dims_list[i]))

    return christoffel

def get_riemann(a, b, c, d, christoffels):
    riemann = diff(christoffels[a][b][d], dims_list[c]) - diff(christoffels[a][b][c], dims_list[d])
    for i in range(4):
        riemann += christoffels[a][i][c] * christoffels[i][b][d]
        riemann -= christoffels[a][i][d] * christoffels[i][b][c]

    return riemann

def get_ricci(a, b, riemanns):
    ricci = 0
    for i in range(4):
        ricci += riemanns[i][a][i][b]

    return ricci

def flip_ricci(riccis):
    for i in range(4):
        for j in range(4):
            riccis[i][j] *= (metric_inv[i][j] * metric_inv[i][j])
    return riccis

def get_ricci_scalar(riccis):
    ricci_scalar = 0
    for i in range(4):
        for j in range(4):
            ricci_scalar += metric_inv[i][j] * riccis[i][j]

    return ricci_scalar

def get_einstein(riccis, ricci_scalar, a, b):
    return riccis[a][b] - (.5 * metric_inv[a][b] * ricci_scalar)

def main():

    christoffels = []
    for i in range(4):
        christoffels.append([])
        for j in range(4):
            christoffels[i].append([])
            for k in range(4):
                christoffels[i][j].append(get_christoffel(i, j, k))
                # if christoffels[i][j][k] != 0:
                #     print("The Christoffel for {}, {}, {} is {}\n".format(dims_list[i],
                #                                                         dims_list[j],
                #                                                         dims_list[k],
                #                                                         christoffels[i][j][k]))


    riemanns = []
    for i in range(4):
        riemanns.append([])
        for j in range(4):
            riemanns[i].append([])
            for k in range(4):
                riemanns[i][j].append([])
                for l in range(4):
                    riemanns[i][j][k].append(get_riemann(i, j, k, l, christoffels))
                    # if riemanns[i][j][k][l] != 0:
                    #     print("The Riemann for {}, {}, {}, {} is {}\n".format(dims_list[i],
                    #                                                         dims_list[j],
                    #                                                         dims_list[k],
                    #                                                         dims_list[l],
                    #                                                         riemanns[i][j][k][l]))
                    
    riccis = []
    for i in range(4):
        riccis.append([])
        for j in range(4):
            riccis[i].append(get_ricci(i, j, riemanns))
            # if riccis[i][j] != 0:
            #     print("The Ricci for {}, {} is {}\n".format(dims_list[i],
            #                                                 dims_list[j],
            #                                                 latex(riccis[i][j])))

    ricci_scalar = get_ricci_scalar(riccis)
    
    # print("R = {}\n".format(latex(simplify(ricci_scalar))))
    
    flip_ricci(riccis)
    
    einstein = []
    for i in range(4):
        einstein.append([0, 0, 0, 0])
        einstein[i][i] = get_einstein(riccis, ricci_scalar, i, i)
            
    for i in range(4):
        for j in range(4):
            if einstein[i][j] != 0:
                print("G^{{ {} {} }} &= {} \\\\".format(dims_list[i],
                                                        dims_list[i],
                                                        latex(simplify(einstein[i][j]))))
                


if __name__ == "__main__":
    main()





