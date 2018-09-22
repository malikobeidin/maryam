# This is interval-exchange-o-matic by Marc Culler and Nathan Dunfield
# (mostly all Marc's work).   Written for Maryam Mizakhani 2005/2/4
#
# To use, in Terminal go to the directory containing this file, and type
# "python maryam.py"  (w/o the quotes).
#

import os, sys, re, random, math, time

# return a parition where each segement has randomly choosen length
# between 1 and 2 N/n

def square_random_partition(N, n):
    lengths = [ random.randrange(1, 2 * N//n)  for i in range(n) ]
    spaces = [  sum(lengths[:i]) for i in range(n+1) ]
    return [ [spaces[i] + 1, spaces[i+1]] for i in range(n) ]
    

# a random partition of [1..N] into n pieces

def basic_random_partition(N, n):
    spaces = random.sample(xrange(1, N), n-1) + [0, N]
    spaces.sort()
    return [ [spaces[i] + 1, spaces[i+1]] for i in range(n) ]

# choose a partition into n piece among all such with total length < N.
# compute choice of length using volumes of simplices in R^n

def random_partition(N, n):
    M = int( math.floor( random.random()**(1.0/n) * (N - n)) + n)
    return basic_random_partition(M, n)

def create_pseudogroup(perm, partition):
    n = len(perm)
    lengths = [ x[1] - x[0] + 1 for x in partition]
    permlengths = [ lengths[perm.index(i)] for i in range(n)]
    spaces = [  sum(permlengths[:i]) for i in range(n+1) ]
    permimages =  [ [spaces[i] + 1, spaces[i+1]] for i in range(n) ]

    return Pseudogroup(
        [Shift( partition[i], permimages[perm[i]]) for i in range(len(perm))],
        Interval(partition[0][0], partition[-1][-1]))

def random_interval_exchange( perm, N ):
    return create_pseudogroup(perm, random_partition(N, len(perm)))

def run_trials(perm, N, num_trials, outfile_name):
    f = open(outfile_name, "a")
    for i in range(num_trials):
        p = random_interval_exchange(perm, N)
        f.write( "%s\t%s\t%d\n" % (perm, N, p.reduce2()))

# this command assumes all guys in the file correspond to the same
# permutation

def examine_trials(file_name):
    ans = load_trials(file_name)
    keys = ans.keys()
    keys.sort()
    for N in keys:
        print_stats(ans[N], N)

def print_stats(data, N):
    num_trials = sum(data.values())
    probone = (data[1]*1.0)/num_trials
    ave =(1.0/num_trials) *  sum( [i * data[i] for i in data.keys()] )
    sigma = math.sqrt(probone* (1 - probone)/num_trials)
    print "N = %d  numtrials = %d  Prob conn: %f  Ave num components: %f  Std. Dev %f" % (N, num_trials, probone, ave, sigma)

    probs = []
    for n in data.keys():
        p = data[n]/(num_trials*1.0)
        if p > 3 * math.sqrt(p* (1 - p)/num_trials):
            probs.append( (n, p) )

    probs.sort()
    for n, p in probs:
        print "%d: %f " % (n,p), 
    print
    
def load_trials(file_name):
    ans = {}
    for line in open(file_name).xreadlines():
        perm, N, components = map(eval, line.split("\t"))
        if not N in ans.keys():
            ans[N] = {components:1}
        else:
            if components not in ans[N].keys():
                ans[N][components] = 1
            else:
                ans[N][components] += 1

    return ans

#  Here's the interactive code
#

def get_permutation():
    found = 0
    while not found:
        try:
            ans = raw_input( "Enter permuation as list of images, e.g. 3 1 2 : ")
            ans = [ int(a) - 1 for a in ans.split()]
            check = ans[:]
            check.sort()
            if not check == range(len(ans)):
                print "Not a bijection of [1, 2, ..., n]"
            elif len(ans) > 0:
                found = 1

        except:
            print "Invalid input"

    return ans

def run_trials(perm, N, trials):
    completed_trials = 0
    data = {}
    last_print_time = time.time()

    while completed_trials < trials:
        exch = random_interval_exchange(perm, N)
        components = exch.reduce()
        if components not in data.keys():
            data[components] = 1
        else:
            data[components] += 1

        if time.time() - last_print_time > 30:
            print_stats(data, N)
            last_print_time = time.time()

        completed_trials += 1

    print_stats(data, N)
        

def main():
    print "Welcome to interval-exchange-o-matic\n    by Marc Culler and Nathan Dunfield\n "
    while 1:
        perm = get_permutation()
        N = int(raw_input("Enter the max points for the interval: "))
        trials = int(raw_input("Enter number of trials: "))
        print "\nStaring computation, will output intermediate results every 30 seconds\n  Hit ctrl-C to stop at any time"
        try:
            run_trials(perm, N, trials)
        except KeyboardInterrupt:
            print "You asked me to stop\n"

        ans = raw_input("\nDo you want to do another compuation? (y/n)")
        if ans[0] != "y":
            break
        
    print "Bye!"


#--------------
#
#  Below code is orbits.py by Marc Culler, version of 2005/2/10
#
#-------------------

from types import IntType, LongType
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
    def __init__(self, domain, isometry):
        self.domain, self.isometry = domain, isometry
        self.range = self(self.domain)

    def __repr__(self):
        if self.isometry.flip:
            op = ' ~> '
        else:
            op = ' -> '
        return str(self.domain) + op + str(self.range)

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
            return Pairing(shift(self.domain), shift * self.isometry * ~shift)
        elif I.end < self.range.start:
            return Pairing(self.domain, shift * self.isometry)
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
                return Pairing(domain, self.isometry)
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
        if other.range not in trim.range:
            return other
        domain = other.domain
        if not trim.is_preserving():
            isometry = trim.isometry * other.isometry
            if domain in trim.range:
                isometry = isometry * trim.isometry
                domain = trim.isometry(domain)
        else:
            shift = trim.isometry.shift
            post = -(1 + (other.range.start - trim.range.start)/shift)
            isometry = (trim.isometry**post) * other.isometry
            if domain in trim.range:
                pre = 1 + (other.domain.start - trim.range.start)/shift
                isometry = isometry * (trim.isometry**pre)
                domain = (trim.isometry**(-pre))(domain)
        range = isometry(domain)
        if range.start < domain.start:
            isometry = isometry**(-1)
            domain = range
        return Pairing(domain, isometry)
            
def Shift(domain, range):
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
    return Pairing(domain, isometry)

def Flip(domain, range):
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
    return Pairing(domain, isometry)

class Pseudogroup:
    """
    Pseudogroup(P,U) is the pseudogroup of maps of the interval U which
    is generated by the Pairings in the list P.  
    """
    def __init__(self, pairings, universe=None):
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
        self.pairings = [Pairing(domain, g.isometry)] + self.pairings

    def simplify(self):
        """
        Do one cycle of the orbit counting reduction algorithm due
        to Agol, Hass and Thurston.
        """
        self.clean()
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

    def reduce(self):
        """
        Reduce the pseudogroup to nothing.  Return the number of orbits.
        """
        count = 0
        while self.pairings != None:
            count += self.simplify()
        return count


#--------- end Marcs code ------


if __name__ == "__main__":
    main()
        
    
    
        
            
