import math
def checkUni(cords): #cords = [(n1, d1), (n2, d2), ...] list of vertices
    res = [False] * (len(cords)-1)
    for i in range(len(cords)-1):
        n1, d1 = cords[i]
        n2, d2 = cords[i+1]
        if abs(n1*d2 - n2*d1) > 1:
            res[i] = False #the simplex with vertices (n1,d1) and (n2,d2) is not unimodular
        else: res[i] = True #the simplex with vertices (n1,d1) and (n2,d2) is unimodular
    return res

def interpolate(xl, yl, xr, yr, mid):
    #res = A/B+C
    A = (yr-yl)*mid[0]*xr[1]*xl[1] - mid[1]*xl[0]*xr[1]*(yr-yl)
    B = xr[0]*xl[1] - xl[0]*xr[1]
    C = yl*mid[1]
    return int(A/B+C)


def simpDenom(cords): # simplify the denominator
    for i in range(len(cords)):
        n, d = cords[i]
        if n == 0:
            cords[i] = (n, 1)
            continue
        if math.gcd(n, d) > 1:
            b = math.gcd(n, d)
            n = n // b
            d = d // b
            cords[i] = (n, d)
    return cords



def searchMinkowski(cord1, cord2): #cord1 = (n1, d1) cord1 = (n2, d2)
    n1, d1 = cord1
    n2, d2 = cord2
    a, b, c, d = 0, 1, 1, 1
    u = a + c
    v = b + d
    while u/v <= n1/d1 or u/v >= n2/d2:
        if u <= v * (n1 / d1):
            a = u
            b = v
            u = a + c
            v = b + d
        elif u >= v * (n2 / d2):
            c = u
            d = v
            u = a + c
            v = b + d
    return (u, v)