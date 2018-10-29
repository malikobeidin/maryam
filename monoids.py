#from orbits import *

class Monoid(object):
    def __init__(self, identity, product_map, involution):
        self.product_map = product_map
        self.involution = involution
        self.identity = MonoidElement(identity, self)

class MarkedMonoid(object):
    def __init__(self, identity, product_map, involution, top_marking, bottom_marking):
        self.product_map = product_map
        self.involution = involution
        self.identity = MarkedMonoidElement(identity, self)
        self.top_marking = top_marking
        self.bottom_marking = bottom_marking

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

class MarkedMonoidElement(MonoidElement):
    def __init__(self, data, monoid):
        self.data = data
        self.monoid = monoid

    def __mul__(self, other):
        if self.monoid != other.monoid:
            raise TypeError("Elements are not in the same monoid")
        else:
            new_data = self.monoid.product_map(self.data, other.data)
            return MarkedMonoidElement(new_data, self.monoid)

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
        return MarkedMonoidElement(self.monoid.involution(self.data), self.monoid)
    
    def top_marked(self):
        return self.monoid.top_marking(self.data)

    def bottom_marked(self):
        return self.monoid.bottom_marking(self.data)


                        
StringMonoid = Monoid( '' , lambda x,y: x+y, lambda t: ''.join(reversed(t)))
WordMonoid = Monoid( '' , lambda x,y: x+y, lambda t: ''.join(reversed(t)).swapcase())
IntMonoid = Monoid( 0, lambda x,y: x+y, lambda x: -x)
PositiveIntMonoid = Monoid( 0, lambda x,y: x+y, lambda x: x)
#IsometryMonoid = Monoid(Isometry(0), lambda x,y: x*y, lambda x: x**(-1))
TrivialMonoid = Monoid(0, lambda x,y: 0, lambda x: 0)

def pair_sum(x,y):
    return (x[0]+y[0], x[1]+y[1])

def flip(x):
    return (-x[0], -x[1])


Z2Monoid = Monoid((0,0), pair_sum, flip)


e = WordMonoid.identity
a = MonoidElement('a',WordMonoid)
b = MonoidElement('b',WordMonoid)
c = MonoidElement('c',WordMonoid)
d = MonoidElement('d',WordMonoid)

eh = Z2Monoid.identity
ah = MonoidElement((1,0),Z2Monoid)
bh = MonoidElement((0,1),Z2Monoid)


def box_product(box1, box2):
    s1, t1, b1, nt1, nb1 = box1
    s2, t2, b2, nt2, nb2 = box2
    s_new = s1+s2
    
    if t1>s1+t2:
        t_new = t1
        nt_new = nt1
    elif t1<s1+t2:
        t_new = s1+t2
        nt_new = nt2
    else:
        t_new = t1
        nt_new = nt1+nt2

    if b1<s1+b2:
        b_new = b1
        nb_new = nb1
    elif b1>s1+b2:
        b_new = s1+b2
        nb_new = nb2
    else:
        b_new = b1
        nb_new = nb1+nb2

    return (s_new, t_new, b_new, nt_new, nb_new)

def box_reverse(box):
    s, t, b, nt, nb = box
    return (-s,t-s,b-s,nt,nb)

def top_marked(box):
    return box[3]>2
    
def bottom_marked(box):
    return box[4]>2
    
BoxMonoid = MarkedMonoid((0,0,0,0,0), box_product, box_reverse, top_marked, bottom_marked)
