from orbits import *
from monoids import *
from freegroup import *
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


                
def random_curve(n):
    zero  = TrivialMonoid.identity
    while True:
        alpha = randint(1,n-1)
        #    beta = randint(1,n)
        beta = randint(1,n-1)
        gamma = randint(1, min(2*alpha,2*beta)/2)*2
        alpha_twist = randint(0,alpha-1)
        beta_twist = randint(0,beta-1)
        gamma_twist = randint(0, gamma-1)
        G = Pseudogroup(genus_2_curve_pairings(alpha,gamma,beta,alpha_twist,gamma_twist,beta_twist,zero,zero,zero),TrivialMonoid)
        if G.reduce() == 1:
            return (alpha,gamma,beta,alpha_twist,gamma_twist,beta_twist)


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
    return [[random_curve(size, eh, ah, bh)[0].reduce_to_single_pairing().data for i in range(num_samples_per_size)] for size in sizes]


def monic_alexander_polynomial_vs_fibers(size, num_samples, e, a, b, F):
    monic_alex_vs_fibers = []
    for i in range(num_samples):
        G, dt_coords = random_curve(size,e,a,b)
        R = snappy_string_to_relator(G.reduce_to_single_pairing().data, F)
        try:
            ma = monic_alexander_polynomial(R,F)
            fi = fibers(R,F)
            monic_alex_vs_fibers.append( (ma, fi, dt_coords, R ) )
            if ma and not fi:
                print(dt_coords, R)
            
        except:
            continue
    return monic_alex_vs_fibers

def monic_alexander_polynomial_but_does_not_fiber(size, num_samples, e, a, b, F):

    exceptional_count = 0
    total_count = 0
    for i in range(num_samples):
        G, dt_coords = random_curve(size,e,a,b)
        R = snappy_string_to_relator(G.reduce_to_single_pairing().data, F)
        try:
            ma = monic_alexander_polynomial(R,F)
            fi = fibers(R,F)

            if ma and not fi:
                exceptional_count += 1
            total_count += 1
            
        except:
            continue
    return exceptional_count, total_count

def fibering_probability(size, num_samples):
    return float(sum(fibers(*random_curve(size))for i in range(num_samples)))/num_samples


def fibering_and_alex_top_coeff_counts(size, num_samples):
    num_fiber = 0
    num_alex_top_coeff_cancels = 0
    num_phi_is_zero = 0
    for i in range(num_samples):
#        print('{}/{}'.format(i,num_samples))
        c = random_curve(size)
        try:
            if fibers(*c):
                num_fiber += 1
            if alexander_poly_top_coeff_cancels(*c):
                num_alex_top_coeff_cancels += 1
        except:
            num_phi_is_zero += 1
    return  num_fiber, num_alex_top_coeff_cancels , num_phi_is_zero

def fibering_and_alex(sizes, num_samples, filename):
    with open(filename+'.txt', 'w') as f:
        f.write('samples per size: {}'.format(num_samples))
        for i, size in enumerate(sizes):
            print(i, len(sizes))
            fi, c, z = fibering_and_alex_top_coeff_counts(size, num_samples)
            f.write('{},{},{},{}\n'.format(size,fi,c,z))

            
        

def fibering_and_signed_fibering_boxes(size, num_samples):
    fibering_boxes = []
    signed_fibering_boxes = []
    for i in range(num_samples):
        print('{}/{}'.format(i,num_samples))
        c = random_curve(size)
        try:
            fibering_boxes.append(fibering_box(*c).data)
            signed_fibering_boxes.append(signed_fibering_box(*c).data)
        except:
            continue
    return  fibering_boxes, signed_fibering_boxes
