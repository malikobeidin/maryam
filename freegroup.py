
from sage.all import FreeGroup, LaurentPolynomialRing, ZZ

def snappy_string_to_relator(s, F):
    a,b = F.gens()
    R = a*a.inverse()
    for g in s:
        if g == 'a':
            R = R*a
        elif g == 'A':
            R = R*a.inverse()
        elif g == 'b':
            R = R*b
        elif g == 'B':
            R = R*b.inverse()
        else:
            raise Exception("Word must be in a and b and their inverses.")
    return R
    
def change_coordinates(phi_a, phi_b):

        
    G = FreeGroup(['t','u'])
    t, u = G.gens()

    F = FreeGroup(['a','b'])
    a, b = F.gens()

    g_a = a
    g_b = b

    
    nielsen_ops = []
    if phi_a < 0:
        g_a = g_a.inverse()
        phi_a = -phi_a
        nielsen_ops.append('t=t.inverse()')

    if phi_a == 0:
        if phi_b == 0:
            raise Exception("phi is 0")
        else:
            return tuple(reversed(change_coordinates(phi_b, phi_a)))
        
    if phi_b < 0:
        g_b = g_b.inverse()
        phi_b = -phi_b
        nielsen_ops.append('u=u.inverse()')
    while phi_b>0:
        if phi_a>phi_b:
            phi_a = phi_a-phi_b
            g_a = g_a*g_b.inverse()
            nielsen_ops.append('t = t*u')
        else:
            phi_b = phi_b-phi_a
            g_b = g_b*g_a.inverse()
            nielsen_ops.append('u = u*t')

    for op in reversed(nielsen_ops):
        exec(op)
    return (t,u)

def abelianized(R, F):
    a, b = F.gens()
    h_a = 0
    h_b = 0
    S = R.syllables()
    for s in S:
        if s[0] == a:
           h_a += s[1] 
        if s[0] == b:
           h_b += s[1]
    return (h_a, h_b)

def alexander_polynomial(R, F):
    a, b = F.gens()
    L = LaurentPolynomialRing(ZZ, ['t'])
    t = L.gens()[0]
    S = R.syllables()
    h_a = 0
    h_b = 0
    for s in S:
        if s[0] == a:
           h_a += s[1] 
        if s[0] == b:
           h_b += s[1]
    new_a, new_b = change_coordinates(-h_b,h_a)
 #   print(new_a, new_b)
    Rsub = R.subs(a=new_a, b=new_b)
    height = 0
    poly = t-t
#    print(Rsub.Tietze())
    for i in Rsub.Tietze():
  #      print(poly)
        if i in (-1,1):
            height += i
        elif i == 2:
            poly += t**height
        else:
            poly -= t**height
    return poly




def ut_walk_just_a(R, F):
    a, b = F.gens()
    S = R.syllables()
    h_a = 0
    h_b = 0
    for s in S:
        if s[0] == a:
           h_a += s[1] 
        if s[0] == b:
           h_b += s[1]
    new_a, new_b = change_coordinates(-h_b,h_a)

    Rsub = new_a
#    print(Rsub)
    height = 0
    heights = [0]

    for i in Rsub.Tietze():
        if i in (-1,1):
            height += i
        heights.append(height)

    return heights



def u_heights(R, F):
    a, b = F.gens()
    S = R.syllables()
    h_a = 0
    h_b = 0
    for s in S:
        if s[0] == a:
           h_a += s[1] 
        if s[0] == b:
           h_b += s[1]
    new_a, new_b = change_coordinates(-h_b,h_a)

    Rsub = R.subs(a=new_a, b=new_b)
    height = 0
    heights = []

    for i in Rsub.Tietze():
        if i in (-1,1):
            height += i
        else:
            if i>0:
                heights.append((height, 1))
            elif i<0:
                heights.append((height,-1))
            else:
                raise Exception()

    return heights


def ut_walk(R, F):
    a, b = F.gens()
    S = R.syllables()
    h_a = 0
    h_b = 0
    for s in S:
        if s[0] == a:
           h_a += s[1] 
        if s[0] == b:
           h_b += s[1]
    new_a, new_b = change_coordinates(-h_b,h_a)

    Rsub = R.subs(a=new_a, b=new_b)
#    print(Rsub)
    height = 0
    heights = [0]

    for i in Rsub.Tietze():
        if i in (-1,1):
            height += i
        heights.append(height)

    return heights


def ut_walk_with_orientations(R, F):
    a, b = F.gens()
    S = R.syllables()
    h_a = 0
    h_b = 0
    for s in S:
        if s[0] == a:
           h_a += s[1] 
        if s[0] == b:
           h_b += s[1]
    new_a, new_b = change_coordinates(-h_b,h_a)

    Rsub = R.subs(a=new_a, b=new_b)
    return Rsub.Tietze()


def plot_walk(walk, filename):
    import matplotlib.pyplot as plt
    plt.clf()
    plt.plot(range(len(walk)), walk)
    plt.savefig(filename+'.png')

def plot_walk_with_orientations(tietze, filename):
    import matplotlib.pyplot as plt
    plt.clf()
    height = 0
    for i, step_type in enumerate(tietze):
        if step_type == 1:
            plt.plot((i, i+1), (height,height+1), color='black')
            height += 1
        elif step_type == -1:
            plt.plot((i, i+1), (height,height-1), color='black')
            height -= 1
        elif step_type == 2:
            plt.plot((i, i+1), (height,height),color = 'blue')
        elif step_type == -2:
            plt.plot((i, i+1), (height,height), color= 'red')
        else:
            raise Exception("Step type messed up.")            
    plt.savefig(filename+'.png')


def unique_max_and_min(walk):
    max_height = max(walk)
    min_height = min(walk)
    max_count = 0
    min_count = 0
    for i in range(len(walk)-1):
        height, next_height = walk[i], walk[i+1]
        if height == next_height == max_height:
            max_count += 1
        elif height == next_height == min_height:
            min_count += 1
    return (max_count < 2) and (min_count < 2)


def relator_fibers(R, F):
    return unique_max_and_min(ut_walk(R,F))


def monic_alexander_polynomial(R, F):
    A = alexander_polynomial(R,F)
    return A.coefficients()[-1] in (-1, 1)

