from types import IntType, LongType
from monoids import *
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

    def ordered_labels(self):
        self.pairings.sort()
        return [pairing.label.data for pairing in self.pairings]

    def ordered_top_markings(self):
        self.pairings.sort()
        return [int(pairing.label.top_marked()) for pairing in self.pairings]

    def ordered_bottom_markings(self):
        self.pairings.sort()
        return [int(pairing.label.bottom_marked()) for pairing in self.pairings]
    
    def reduce_to_single_pairing(self, show_labels=False):
        labels = []
        top_markings = []
        bottom_markings = []
        while len(self.pairings)>1:
            if show_labels:
                labels.append(self.ordered_labels())
                top_markings.append(self.ordered_top_markings())
                bottom_markings.append(self.ordered_bottom_markings())
            self.transmit()
            self.truncate()
        if show_labels:
            return self.pairings[0].label, labels, top_markings, bottom_markings
        else:
            return self.pairings[0].label
        

    def reducible_to_marked_labels(self):
        while len(self.pairings)>1:
#            print('pairings')
#            print(self.pairings)
            self.transmit()
            self.truncate()
            all_tops_marked = True
            all_bottoms_marked = True
            for pairing in self.pairings:
#                print(pairing)
#                print('top marked: '+str(pairing.label.top_marked()))
#                print('bottom marked: '+str(pairing.label.bottom_marked()))
                if not pairing.label.top_marked():
                    all_tops_marked = False
                if not pairing.label.bottom_marked():
                    all_bottoms_marked = False
                if not all_tops_marked and not all_bottoms_marked:
                    break
            if all_tops_marked or all_bottoms_marked:
                return True
        final_pairing = self.pairings[0]
        return final_pairing.label.top_marked() or final_pairing.label.bottom_marked()
        


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

def genus_2_curve_pairings(wl, wm, wr, twl, twm, twr, e, a, b):
    
    t1, t2, t3 = Triangle(wl, wl, wm).transition_numbers
    t4, t5, t6 = Triangle(wm, wr, wr).transition_numbers

#    print('t1, t2, t2: {}, {}, {}'.format(t1,t2,t3))
#    print('t4, t5, t6: {}, {}, {}'.format(t4,t5,t6))
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
#    print('s1, s2, s3, s4, s5, s6: {}, {}, {}, {}, {}, {}'.format(s1,s2,s3, s4,s5,s6))
#    print('e1, e2, e3, e4, e5, e6: {}, {}, {}, {}, {}, {}'.format(e1,e2,e3, e4,e5,e6))
    pairings = []
    
    pairings.append(Flip([s1,s1+t1-1],[s3+t3,e3],e))

    
    pairings.append(Flip([s2+t2,e2],[s3,s3+t3-1],e))

    
    if t2>0:        
        pairings.append(Flip([s1+t1,e1],[s2,s2+t2-1],e))

    
    pairings.append(Flip([s4,s4+t4-1],[s6+t6,e6],e))


    if t6>0:
        pairings.append(Flip([s5+t5,e5],[s6,s6+t6-1],e))

    
    pairings.append(Flip([s4+t4,e4],[s5,s5+t5-1],e))

    
    pairings.append(Flip([s1, e1-twl],[s2,e2-twl],a))

    
    if twl>0:
        pairings.append(Flip([e1-twl+1, e1],[e2-twl+1,e2],a))


    pairings.append(Flip([s3, e3-twm],[s4,e4-twm],e))

    if twm>0:
        pairings.append(Flip([e3-twm+1, e3],[e4-twm+1,e4],e))


    pairings.append(Flip([s5, e5-twr],[s6,e6-twr],b))

    if twr>0:
        pairings.append(Flip([e5-twr+1, e5],[e6-twr+1,e6],b))

        
    return pairings

def genus_2_pseudogroup(wl, wm, wr, twl, twm, twr, e, a, b, Monoid):
    return Pseudogroup(genus_2_curve_pairings(wl, wm, wr, twl, twm, twr, e, a, b), Monoid)


def relator(wl, wm, wr, twl, twm, twr):
    return genus_2_pseudogroup(wl, wm, wr, twl, twm, twr, e, a, b, WordMonoid).reduce_to_single_pairing().data

def homology(wl, wm, wr, twl, twm, twr):
    return genus_2_pseudogroup(wl, wm, wr, twl, twm, twr, eh, ah, bh, Z2Monoid).reduce_to_single_pairing().data


def starting_fibering_boxes(phi_a, phi_b):
    if phi_a == 0 or phi_b == 0:
        if phi_a == 0 and phi_b == 0:
            raise Exception("phi is zero")
        if phi_a == 0:
            box_a = MarkedMonoidElement((0,0,0,2,2),BoxMonoid)
            box_b = MarkedMonoidElement((1,1,0,0,0),BoxMonoid)        
        if phi_b == 0:
            box_a = MarkedMonoidElement((1,1,0,0,0),BoxMonoid)
            box_b = MarkedMonoidElement((0,0,0,2,2),BoxMonoid)        
    else:
        if phi_a > 0:
            box_a = MarkedMonoidElement((phi_a,phi_a,0,1,1),BoxMonoid)
        if phi_a < 0:
            box_a = MarkedMonoidElement((phi_a,0,phi_a,1,1),BoxMonoid)
        if phi_b > 0:
            box_b = MarkedMonoidElement((phi_b,phi_b,0,1,1),BoxMonoid)
        if phi_b < 0:
            box_b = MarkedMonoidElement((phi_b,0,phi_b,1,1),BoxMonoid)
    return box_a, box_b

def fibers(wl, wm, wr, twl, twm, twr):
    x, y = homology(wl, wm, wr, twl, twm, twr)
    phi_a, phi_b = -y, x
    box_a, box_b = starting_fibering_boxes(phi_a, phi_b)
    box_id = BoxMonoid.identity
    return not genus_2_pseudogroup(wl, wm, wr, twl, twm, twr, box_id, box_a, box_b, BoxMonoid).reducible_to_marked_labels()

def fibering_box(wl, wm, wr, twl, twm, twr, show_labels=False):
    x, y = homology(wl, wm, wr, twl, twm, twr)
    phi_a, phi_b = -y, x
    box_a, box_b = starting_fibering_boxes(phi_a, phi_b)
    box_id = BoxMonoid.identity
    return genus_2_pseudogroup(wl, wm, wr, twl, twm, twr, box_id, box_a, box_b, BoxMonoid).reduce_to_single_pairing(show_labels=show_labels)


def starting_signed_fibering_boxes(phi_a, phi_b):
    if phi_a == 0 or phi_b == 0:
        if phi_a == 0 and phi_b == 0:
            raise Exception("phi is zero")
        if phi_a == 0:
            box_a = MarkedMonoidElement((0,0,0,2,0,2,0),SignedBoxMonoid)
            box_b = MarkedMonoidElement((1,1,0,0,0,0,0),SignedBoxMonoid)        
        if phi_b == 0:
            box_a = MarkedMonoidElement((1,1,0,0,0,0,0),SignedBoxMonoid)
            box_b = MarkedMonoidElement((0,0,0,0,2,0,2),SignedBoxMonoid)        
    else:
        if phi_a > 0:
            box_a = MarkedMonoidElement((phi_a,phi_a,0,1,0,0,1),SignedBoxMonoid)
        if phi_a < 0:
            box_a = MarkedMonoidElement((phi_a,0,phi_a,0,1,1,0),SignedBoxMonoid)
        if phi_b > 0:
            box_b = MarkedMonoidElement((phi_b,phi_b,0,0,1,1,0),SignedBoxMonoid)
        if phi_b < 0:
            box_b = MarkedMonoidElement((phi_b,0,phi_b,1,0,0,1),SignedBoxMonoid)
    return box_a, box_b

def alexander_poly_top_coeff_cancels(wl, wm, wr, twl, twm, twr):
    x, y = homology(wl, wm, wr, twl, twm, twr)
    phi_a, phi_b = -y, x
    box_a, box_b = starting_signed_fibering_boxes(phi_a, phi_b)
    box_id = SignedBoxMonoid.identity
    return genus_2_pseudogroup(wl, wm, wr, twl, twm, twr, box_id, box_a, box_b, SignedBoxMonoid).reducible_to_marked_labels()


def signed_fibering_box(wl, wm, wr, twl, twm, twr, show_labels=False):
    x, y = homology(wl, wm, wr, twl, twm, twr)
    phi_a, phi_b = -y, x
    box_a, box_b = starting_signed_fibering_boxes(phi_a, phi_b)
    box_id = SignedBoxMonoid.identity
    return genus_2_pseudogroup(wl, wm, wr, twl, twm, twr, box_id, box_a, box_b, SignedBoxMonoid).reduce_to_single_pairing(show_labels=show_labels)


#H = Pseudogroup([ Flip([11,16],[13,18], a), Shift([1,4],[3,6], b), Shift([10,14],[14,18], c), Shift([8,10],[14,16], d), Flip([6,7],[8,9], e) ], StringMonoid, Interval(1,18))
