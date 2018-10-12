from orbits import *
from random import randint

def test_transmit():
    P1 = Flip([1,5],[11,15], a*b)
    P2 = Flip([7,9],[12,14], c*d)
    P1tP2 = P1.transmit(P2)
    P2tP1 = P2.transmit(P1)
    print(P1tP2)
    print(P2tP1)


def test_transmit2():
    P1 = Flip([1,5],[11,15])
    P2 = Flip([7,9],[12,14])
    P1i = MonoidElement(P1.isometry, IsometryMonoid)
    P2i = MonoidElement(P2.isometry, IsometryMonoid)
    P1 = Flip([1,5],[11,15], P1i)
    P2 = Flip([7,9],[12,14], P2i)
    
    P1tP2 = P1.transmit(P2)
    P2tP1 = P2.transmit(P1)
    return P1tP2, P2tP1
    print(P1tP2)
    print(P2tP1)


                
def random_curve(n,e,a,b):
    while True:
        alpha = randint(1,n-1)
        #    beta = randint(1,n)
        beta = randint(1,n-1)
        gamma = randint(1, min(2*alpha,2*beta)/2)*2
        alpha_twist = randint(0,alpha-1)
        beta_twist = randint(0,beta-1)
        gamma_twist = randint(0, gamma-1)
        G = Pseudogroup(genus_2_curve_pairings(alpha,gamma,beta,alpha_twist,gamma_twist,beta_twist,e,a,b),WordMonoid)
        Gc = G.copy()
        if Gc.reduce() == 1:
            return G


def relator_walk(relator):
    x,y = (0,0)
    walk = [(x,y)]
    for letter in relator:
        if letter == 'a':
            x,y = x+1,y
        if letter == 'A':
            x,y = x-1,y
        if letter == 'b':
            x,y = x,y+1
        if letter == 'B':
            x,y = x,y-1
        walk.append((x,y))
    return walk

def draw_relator_walk(filename, relator):
    import matplotlib.pyplot as plt
    plt.clf()
    fig, ax = plt.subplots()
#        fig.patch.set_visible(False)
#    ax.axis('off')
    walk = relator_walk(relator)
    for i in range(len(walk)-1):
        x1, y1 = walk[i]
        x2, y2 = walk[i+1]
        ax.plot([x1,x2],[y1,y2], c='black')
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_axis_bgcolor('white')
    with open(filename,'w') as outfile:
        fig.canvas.print_png(outfile)

def random_homology_stats(sizes, num_samples_per_size):
    return [[random_curve(size, eh, ah, bh).reduce_to_single_pairing().data for i in range(num_samples_per_size)] for size in sizes]
