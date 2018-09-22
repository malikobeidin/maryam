from orbits import *
G = Pseudogroup([ Flip([11,16],[13,18]), Shift([1,4],[3,6]), Shift([10,14],[14,18]), Shift([8,10],[14,16]) ], Interval(1,18))
print 'G=', G
print 'Reducing G:', G.reduce(), 'orbits.'
print

H = Pseudogroup([ Flip([11,16],[13,18]), Shift([1,4],[3,6]), Shift([10,14],[14,18]), Shift([8,10],[14,16]), Flip([6,7],[8,9]) ], Interval(1,18))
print 'H=', H
print 'Reducing H:', H.reduce(), 'orbits.'

