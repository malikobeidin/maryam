from types import IntType, LongType
from random import randint
Illegal = "Illegal Operation"

def gcd(x, y):
    if x == 0:
        if y == 0:
            raise ValueError, "gcd(0,0) is undefined."
        else:
            return abs(y)
    x = abs(x)
    y = abs(y)
    while y != 0:
        r = x%y
        x = y
        y = r
    return x

class Interval:
    """
    A finite subinterval of the integers.
    """
    def __init__(self, a, b):
        self.start = min(a,b)
        self.end = max(a,b)
        self.width = self.end - self.start + 1
        
    def __repr__(self):
        return '[%d, %d]'%(self.start, self.end)

    def __cmp__(self, other):
        """
        Intervals are ordered by (start, end) in lex order.
        """
        if self.start != other.start:
            return cmp(self.start, other.start)
        return cmp(self.end, other.end)
    
    def __contains__(self, x):
        """
        True if the Interval contains the integer or Interval argument.
        """
        if type(x) == IntType or type(x) == LongType :
            return self.start <= x <= self.end
        else:
            return self.start <= x.start and x.end <= self.end

    def __xor__(self, other):
        """
        Intersection of two intervals
        """
        start = max(self.start, other.start)
        end = min(self.end, other.end)
        if end < start:
            return None
        else:
            return Interval(start, end)

    def set_end(self, end):
        self.end = end
        self.width = self.end - self.start + 1

    def set_start(self, start):
        self.start = start
        self.width = self.end - self.start + 1
        
def ToInterval(x):
    """
    Converts an integer or a 2-tuple to an Interval.
    """
    if x.__class__ == Interval:
        return x
    if type(x) == IntType or type(x) == LongType :
        return Interval(x,x)
    else:
        return Interval(x[0],x[1])

    
class Isometry:
    """
    An element of the infinite dihedral group acting on the integers.
    """
    def __init__(self, shift, flip=0):
        self.shift, self.flip = shift, flip

    def __repr__(self):
        if self.flip:
            return 'x -> -x + %d'%self.shift
        else:
            return 'x -> x + %d'%self.shift
        

    def __mul__(self, other):
        """
        Composition operator for Isometries.
        """
        flip = self.flip ^ other.flip
        if self.flip:
            shift = self.shift - other.shift
        else:
            shift = other.shift + self.shift
        return Isometry(shift, flip)

    def __pow__(self, n):
        """
        Power operator for Isometries.
        """
        if self.flip:
            if n%2 != 0:
                return Isometry(self.shift, self.flip)
            else:
                return Isometry(0,0)
        else:
            return Isometry(n*self.shift, self.flip)

    def __invert__(self):
        """
        Inversion operator for Isometries.
        """
        if self.flip:
            return Isometry(self.shift, self.flip)
        else:
            return Isometry(-self.shift, self.flip)
    
    def __call__(self, x):
        """
        An Isometry as a mapping (of an integer or an interval).
        """
        if type(x) == IntType or type(x) == LongType:
            if self.flip:
                return -x + self.shift
            else:
                return x + self.shift
        if x.__class__ == Interval:
            return Interval(self(x.start), self(x.end))


class Monoid(object):
    def __init__(self, identity, product_map, involution):
        self.product_map = product_map
        self.involution = involution
        self.identity = MonoidElement(identity, self)

class MonoidElement(object):
    def __init__(self, data, monoid):
        self.data = data
        self.monoid = monoid

    def __mul__(self, other):
        if self.monoid != other.monoid:
            raise TypeError("Elements are not in the same monoid")
        else:
            new_data = self.monoid.product_map(self.data, other.data)
            return MonoidElement(new_data, self.monoid)

    def __repr__(self):
        return str(self.data)

    def __eq__(self, other):
        return self.data == other.data and self.monoid == other.monoid


    def __pow__(self, n):
        p = self.monoid.identity
        if n == 0:
            return p
        elif n>0:            
            for i in range(n):
                p = p*self
        else:
            op = self.op()
            for i in range(abs(n)):
                p = p*op
        return p
             
    
    def op(self):
        return MonoidElement(self.monoid.involution(self.data), self.monoid)



StringMonoid = Monoid( '' , lambda x,y: x+y, lambda t: ''.join(reversed(t)))

def cat(s1, s2):
    return s1+s2

def inv(s):
    rev =  ''.join(reversed(s))
    return rev.swapcase()
                        
StringMonoid = Monoid( '' , lambda x,y: x+y, lambda t: ''.join(reversed(t)))
WordMonoid = Monoid( '' , lambda x,y: x+y, lambda t: ''.join(reversed(t)).swapcase())
IntMonoid = Monoid( 0, lambda x,y: x+y, lambda x: -x)
PositiveIntMonoid = Monoid( 0, lambda x,y: x+y, lambda x: x)
IsometryMonoid = Monoid(Isometry(0), lambda x,y: x*y, lambda x: x**(-1))


class Pairing:
    """
    The restriction of an isometry to a finite interval.
    """
    def __init__(self, domain, isometry, label = None):
        self.domain, self.isometry = domain, isometry
        self.range = self(self.domain)
        if label:
            self.label = label
        else:
            self.label = MonoidElement(1, IntMonoid)
        
        
    def __repr__(self):
        if self.isometry.flip:
            op = ' ~> '
        else:
            op = ' -> '
        return str(self.domain) + op + str(self.range) +', ' + str(self.label)

    def __call__(self, x):
        """
        A Pairing as a mapping.
        """
        if not x in self.domain:
            raise Illegal, "Operand is not contained in domain."
        else:
            return self.isometry(x)

    def __cmp__(self, other):
        """
        Linear ordering of Pairings.
        """
        if self.range.end != other.range.end:
            return cmp(other.range.end, self.range.end)
        if self.domain.width != other.domain.width:
            return cmp(other.domain.width, self.domain.width)
        if self.domain.start != other.domain.start:
            return cmp(self.domain.start, other.domain.start)
        return cmp(self.isometry.flip, other.isometry.flip)

    def __contains__(self,x):
        """
        True if the argument is contained in either the domain or range.
        """
        return x in self.domain or x in self.range

    def is_preserving(self):
        """
        True if the Pairing is orientation preserving.
        """
        return self.isometry.flip == 0 or self.domain.width == 1

    def is_periodic(self):
        """
        True if the Pairing is orientation preserving, and
        its domain and range meet.
        """
        return self.is_preserving() and self.domain ^ self.range

    def is_trivial(self):
        """
        True if the Pairing is  restriction of the identity map.
        """
        return (self.is_preserving and self.isometry.shift == 0 or
	       self.domain.width == 1 and self.domain == self.range)
    
    def contract(self,I):
        """
        Adjust the Pairing to account for removal of a static interval.
        """
        I = ToInterval(I)
        if I ^ self.domain or I ^ self.range:
            raise Illegal, "Contraction interval is not static."
        shift = Isometry( -I.width )
        if I.end < self.domain.start:
            return Pairing(shift(self.domain), shift * self.isometry * ~shift, label = self.label)
        elif I.end < self.range.start:
            return Pairing(self.domain, shift * self.isometry, label = self.label)
        else:
            return self

    def trim(self):
        """
        Trim an orientation reversing pairing so that its domain and
        range become disjoint.
        """
        if self.is_preserving():
            return self
        else:
            intersection = self.domain ^ self.range
            if intersection:
                middle = (self.domain.start + self.range.end - 1)/2
                domain = Interval(self.domain.start, middle)
                return Pairing(domain, self.isometry, label = self.label)
            else:
                return self

    def merge(self, other):
        """
        Merge a periodic Pairing with an overlapping orientation 
        preserving Pairing.
        """
        if self.is_periodic() and other.is_preserving():
            R = Interval(self.domain.start, self.range.end)
            I = R ^ other.domain
            shift = gcd(self.isometry.shift, other.isometry.shift)
            if (other(I) ^ R).width >= self.isometry.shift :
                domain = Interval(R.start, R.end - shift)
                isometry = Isometry(shift)
                return Pairing(domain, isometry)
            else:
                return None
        else:
            raise Illegal, "Pairing cannot be merged."

    def transmit(self, other):
        """
        Left shift the domain and range of another Pairing as far as possible.
        """
        trim = self.trim()
        post = None
        pre = None
        if other.range not in trim.range:
            return other
        domain = other.domain
        if not trim.is_preserving():
#            print('reversing')
            isometry = (trim.isometry**(-1)) * other.isometry
            new_label =  other.label * (self.label**(-1))
            if domain in trim.range:
                isometry = isometry * trim.isometry
                new_label = self.label * new_label
                domain = trim.isometry(domain)
        else:
            shift = trim.isometry.shift
            post = -(1 + (other.range.start - trim.range.start)/shift)
            isometry = (trim.isometry**post) * other.isometry
            new_label = other.label * (self.label**post)
            if domain in trim.range:
                pre = 1 + (other.domain.start - trim.range.start)/shift
                isometry = isometry * (trim.isometry**pre)
                new_label =   (self.label**pre) * new_label
                domain = (trim.isometry**(-pre))(domain)
        range = isometry(domain)
        if range.start < domain.start:
#            print('flipping')
            isometry = isometry**(-1)
            new_label = new_label**(-1)
            domain = range
        return Pairing(domain, isometry, label = new_label)
            
def Shift(domain, range, label = None):
    """
    Constructor for an orientation preserving pairing, given the domain
    and range.
    """
    if domain.__class__ != Interval:
        domain = Interval(domain[0], domain[1])
    if range.__class__ != Interval:
        range = Interval(range[0], range[1])
    if domain.width != range.width:
        raise Illegal, "The domain and range must have the same width."
    if range.start < domain.start:
        domain, range = range, domain
    isometry = Isometry(range.start - domain.start)
    return Pairing(domain, isometry, label = label)

def Flip(domain, range, label=None):
    """
    Constructor for an orientation reversing pairing, given the domain
    and range.
    """
    if domain.__class__ != Interval:
        domain = Interval(domain[0], domain[1])
    if range.__class__ != Interval:
        range = Interval(range[0], range[1])
    if domain.width != range.width:
        raise Illegal, "The domain and range must have the same width."
    if range.start < domain.start:
        domain, range = range, domain
    isometry = Isometry(range.end + domain.start, 1)
    return Pairing(domain, isometry, label = label)

class Pseudogroup:
    """
    Pseudogroup(P,U) is the pseudogroup of maps of the interval U which
    is generated by the Pairings in the list P.  
    """
    def __init__(self, pairings, label_monoid, universe=None):
        self.pairings = pairings
        start = min([p.domain.start for p in self.pairings])
        end = max([p.range.end for p in self.pairings])
	if universe: 
            universe = ToInterval(universe)
	    if start < universe.start or end > universe.end:
                raise ValueError, 'Universe must contain all domains and ranges.'
            self.universe = universe
	else:
            self.universe = Interval(start, end)

        self.label_monoid = label_monoid
        
        
    def __repr__(self):
        result = 'Pseudogroup on %s:\n'%str(self.universe)
	if self.pairings:
            self.pairings.sort()
            for pairing in self.pairings:
                result += str(pairing) + '\n'
	return result 

    def clean(self):
        """
        Get rid of trivial Pairings.
        """
        self.pairings = [p for p in self.pairings if not p.is_trivial()]

    def trim(self):
        """
        Trim all orientation reversing pairings.
        """
        self.pairings = [p.trim() for p in self.pairings]

    def static(self):
        """
        Find a static interval.
        """
        if len(self.pairings) == 0:
            return self.universe
        intervals = [p.domain for p in self.pairings]
        intervals += [p.range for p in self.pairings]
        intervals.sort()
	I = intervals.pop(0)
        start, end = I.start, I.end
        if start > 1:
            return Interval(1, start - 1)
        for interval in intervals:
            if end < interval.start - 1:
                return Interval(end+1, interval.start-1)
            end = max(end,interval.end)
        if end < self.universe.end:
            return Interval(end + 1, self.universe.end)
        return None

    def contract(self):
        """
        Remove all static intervals.  Return the total size.
        """
        result = 0
        I = self.static()
        while I:
            result += I.width
            if I.end != self.universe.end:
                self.pairings =  [p.contract(I) for p in self.pairings ]
            self.universe.set_end(self.universe.end - I.width)
            I = self.static()
        return result
    
    def merge(self):
        """
        Merge periodic pairing whenever possible.
        """
        if len(self.pairings) < 2:
             return
        done=0
        while not done:
            periodics = [p for p in self.pairings if p.is_periodic()]
            done=1
            for p in periodics[:-1]:
                for q in periodics[1+periodics.index(p):]:
                    g = None 
                    try:
                        g = p.merge(q)
	            except: pass
                    if g:
                        self.pairings.remove(p) 
                        self.pairings.remove(q)
                        self.pairings.append(g)
                        done=0
                        break

    def transmit(self):
        """
        Use the largest Pairing to transmit others.
        """
        self.pairings.sort()
        g = self.pairings[0]
        self.pairings = [g] + [ g.transmit(p) for p in self.pairings[1:] ]

    def truncate(self):
        """
        Truncate the largest pairing.
        """
        self.pairings.sort()
        g = self.pairings.pop(0)
        if len(self.pairings) > 0:
            support_end = self.pairings[0].range.end
        else:
            support_end = g.range.start - 1
        if support_end < g.range.start:
            self.universe.set_end(support_end)
            return
        if not g.is_preserving():
            g.trim()
        self.universe.set_end(support_end)
        range = Interval(g.range.start, support_end)
        domain = (~g.isometry)(range)
        self.pairings = [Pairing(domain, g.isometry, label = g.label)] + self.pairings

        
    def simplify(self):
        """
        Do one cycle of the orbit counting reduction algorithm due
        to Agol, Hass and Thurston.
        """
#        print(len(self.pairings))
#        print('cleaning')
        self.clean()
#        print(len(self.pairings))
        if len(self.pairings) == 0:
            self.pairings = None
            return self.universe.width
#        print "cleaned\n", self
        count = self.contract()
#        print "contracted\n", self
        self.trim()
#        print "trimmed\n", self
        self.merge()
#        print "merged\n", self
        self.transmit()
#        print "transmitted\n", self
        self.truncate()
#        print "truncated\n", self
#	print 'count = ', count
        return count


    def copy(self):
        return Pseudogroup([Pairing(P.domain, P.isometry, P.label) for P in self.pairings], self.label_monoid)
    
    def reduce(self):
        """
        Reduce the pseudogroup to nothing.  Return the number of orbits.
        """
        count = 0
        while self.pairings != None:
            count += self.simplify()
        return count

    def reduce_to_single_pairing(self):
        while len(self.pairings)>1:
            self.transmit()
            self.truncate()
        return self.pairings[0].label




class Triangle(object):
    def __init__(self, left, bottom, right):
        if not self._is_consistent(left,bottom, right):
            raise Exception("Intersection numbers inconsistent")
        self.left = left
        self.bottom = bottom
        self.right = right
        self.intersection_numbers = [left, bottom, right]
        self._find_transition_numbers()

    def _is_consistent(self, left, bottom, right):
         if left > right+bottom or bottom > left+right or right > left+bottom:
             return False
         if (left+bottom+right)%2 == 1:
             return False
         return True
     
    def _find_transition_numbers(self):
         x = (self.left+self.right-self.bottom)/2
         y = (self.left+self.bottom-self.right)/2
         z = (self.bottom+self.right-self.left)/2

         assert x >= 0 and y >= 0 and z >= 0
         self.transition_numbers = [x, y, z]

    def _cross_triangle(self,side, position):
        if position < self.transition_numbers[side]:
            return (side-1)%3, self.intersection_numbers[side-1] - position - 1
        else:
            return (side+1)%3, self.intersection_numbers[side] - position - 1

e = WordMonoid.identity
a = MonoidElement('a',WordMonoid)
b = MonoidElement('b',WordMonoid)
c = MonoidElement('c',WordMonoid)
d = MonoidElement('d',WordMonoid)
def genus_2_curve_pairings(wl, wm, wr, twl, twm, twr):
    
    t1, t2, t3 = Triangle(wl, wl, wm).transition_numbers
    t4, t5, t6 = Triangle(wm, wr, wr).transition_numbers

    print('t1, t2, t2: {}, {}, {}'.format(t1,t2,t3))
    print('t4, t5, t6: {}, {}, {}'.format(t4,t5,t6))
    #start of each interval
    s1 = 1
    e1 = wl
    
    s2 = s1+wl
    e2 = e1+wl
    
    s3 = s2+wl
    e3 = e2+wm
    
    s4 = s3+wm
    e4 = e3+wm
    
    s5 = s4+wm
    e5 = e4+wr
    
    s6 = s5+wr
    e6 = e5+wr
    print('s1, s2, s3, s4, s5, s6: {}, {}, {}, {}, {}, {}'.format(s1,s2,s3, s4,s5,s6))
    print('e1, e2, e3, e4, e5, e6: {}, {}, {}, {}, {}, {}'.format(e1,e2,e3, e4,e5,e6))
    pairings = []
    
    pairings.append(Flip([s1,s1+t1-1],[s3+t3,e3],e))
#    print(pairings)
    
    pairings.append(Flip([s2+t2,e2],[s3,s3+t3-1],e))
#    print(pairings)
    
    if t2>0:        
        pairings.append(Flip([s1+t1,e1],[s2,s2+t2-1],e))
#    print(pairings)
    
    pairings.append(Flip([s4,s4+t4-1],[s6+t6,e6],e))
#    print(pairings)

    if t6>0:
        pairings.append(Flip([s5+t5,e5],[s6,s6+t6-1],e))
#    print(pairings)
    
    pairings.append(Flip([s4+t4,e4],[s5,s5+t5-1],e))
#    print(pairings)
    
    pairings.append(Flip([s1, e1-twl],[s2,e2-twl],a))
#    print(pairings)
    
    if twl>0:
        pairings.append(Flip([e1-twl+1, e1],[e2-twl+1,e2],a))
#    print(pairings)

    pairings.append(Flip([s3, e3-twm],[s4,e4-twm],e))
#    print(pairings)
    if twm>0:
        pairings.append(Flip([e3-twm+1, e3],[e4-twm+1,e4],e))
#    print(pairings)

    pairings.append(Flip([s5, e5-twr],[s6,e6-twr],b))
#    print(pairings)
    if twr>0:
        pairings.append(Flip([e5-twr+1, e5],[e6-twr+1,e6],b))
#    print(pairings)
        
    return pairings

#H = Pseudogroup([ Flip([11,16],[13,18], a), Shift([1,4],[3,6], b), Shift([10,14],[14,18], c), Shift([8,10],[14,16], d), Flip([6,7],[8,9], e) ], StringMonoid, Interval(1,18))

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
    while True:
        alpha = randint(1,n-1)
        #    beta = randint(1,n)
        beta = randint(1,n-1)
        gamma = randint(1, min(2*alpha,2*beta)/2)*2
        alpha_twist = randint(0,alpha-1)
        beta_twist = randint(0,beta-1)
        gamma_twist = randint(0, gamma-1)
        G = Pseudogroup(genus_2_curve_pairings(alpha,gamma,beta,alpha_twist,gamma_twist,beta_twist),WordMonoid)
        if G.reduce() == 1:
            print(alpha, gamma, beta, alpha_twist, gamma_twist, beta_twist)
            G = Pseudogroup(genus_2_curve_pairings(alpha,gamma,beta,alpha_twist,gamma_twist,beta_twist),WordMonoid)
            print(G)
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
    ax.axis('off')
    walk = relator_walk(relator)
    for i in range(len(walk)-1):
        x1, y1 = walk[i]
        x2, y2 = walk[i+1]
        ax.plot([x1,x2],[y1,y2], c='black')
    ax.set_axis_bgcolor('white')
    with open(filename,'w') as outfile:
        fig.canvas.print_png(outfile)
