#!/usr/bin/env python
import math
import cmath
from collections import Counter

#start = input("start value here:")

start = 100
start = int(start)
limit = start+100
limit = int(limit)
#list of TRUE to depth of search space
amplitude_map = [0]*int(limit-start)
#get RANGE for number of iterations through functions
#constants for QUADRATIC
a = 90
b = -300
c = 250 - limit 
# calculate the discriminant
d = (b**2) - (4*a*c)
# find two solutions
sol1 = (-b-cmath.sqrt(d))/(2*a)
sol2 = (-b+cmath.sqrt(d))/(2*a)
print('The solution are {0} and {1}'.format(sol1,sol2))
#the integer REAL part of this value is the limit for RANGE
new_limit = sol2
#here comes a function for generating composites
def drLD(x, l, m, z):   
  y = 90*(x*x) - l*x + m #y is the "location operator"
  if start <= y < limit: 
    amplitude_map[y-start] = amplitude_map[y-start]+1
 #cancellation operators p,q:
  p = z+(90*(x-1))  
  newp_start = int(((start-y)/p))
  newp_lim = int(((limit-y)/p)+1)
  if int(new_limit.real)/x < 1:
    newp_start = newp_start-1
  for n in range(newp_start, newp_lim):
    new_y = y+(p*n)
    if start <= new_y < limit: 
      amplitude_map[new_y - start] = amplitude_map[new_y-start]+1
 

for x in range(1, int(new_limit.real)): 
#11 A201804:
    drLD(x, 78, -1, 11)  #11,91 @11, 203 5    564
    drLD(x, 78, -1, 91)  #11,91 @11, 203 5    564
    
    drLD(x, 120, 34, 7)  #7,53  @4,  154 1    831
    drLD(x, 120, 34, 53)  #7,53  @4,  154 1    831   
    
    drLD(x, 132, 48, 19) #19,29 @6,  144 2    371
    drLD(x, 132, 48, 29) #19,29 @6,  144 2    371
    
    drLD(x, 120, 38, 17) #17,43 @8,  158 3    403 
    drLD(x, 120, 38, 43) #17,43 @8,  158 3    403     
    
    drLD(x, 90, 11, 13)  #13,77 @11, 191 4    493
    drLD(x, 90, 11, 77)  #13,77 @11, 191 4    493    
    
    drLD(x, 108, 32, 31) #31,41 @14, 176 6    262
    drLD(x, 108, 32, 41) #31,41 @14, 176 6    262    
    
    drLD(x, 90, 17, 23)  #23,67 @17, 197 7    320
    drLD(x, 90, 17, 67)  #23,67 @17, 197 7    320    
    
    drLD(x, 72, 14, 49)  #49,59 @32, 230 8    192
    drLD(x, 72, 14, 59)  #49,59 @32, 230 8    192    
    
    drLD(x, 60, 4, 37)   #37,83 @34, 244 9    229
    drLD(x, 60, 4, 83)   #37,83 @34, 244 9    229    
    
    drLD(x, 60, 8, 47)   #47,73 @38, 248 10   195   
    drLD(x, 60, 8, 73)   #47,73 @38, 248 10   195     
    
    drLD(x, 48, 6, 61)   #61,71 @48, 270 11   165
    drLD(x, 48, 6, 71)   #61,71 @48, 270 11   165    
    
    drLD(x, 12, 0, 79)   #79,89 @78, 336 12   139
    drLD(x, 12, 0, 89)   #79,89 @78, 336 12   139
    
    
print("The total energy to produce this list:", sum(amplitude_map))
A201804c = [i for i,x in enumerate(amplitude_map) if x == 0]
newA = [(i+(start)) for i in A201804c]
newA = newA#[1:]
print("These are terms in OEIS sequence A201804", newA)
print("This is the number of A201804 amplitude_map zeroes", len(newA), "between", start, "and", limit)
A201804b = [(i*90)+11 for i in newA]
A201804b = A201804b#[:1] #defect of inplementation 
amplitude_map = amplitude_map#[1:]
print("This is the base-10 expression of A201804 (see A142317)", A201804b)#[-50:])
base10lim = (limit*90)+11
operators = int(new_limit.real)*12
print("This is the number of LOCATION OPERATORS :", operators)
print("This is the number of operations:", sum(amplitude_map))
print("This is the limit in base-10:", base10lim)
print("This is the ratio of base10limit to OPERATORS:", base10lim/operators)
print("This is the square root of the base10limit:", math.sqrt(base10lim))
print("This is the ratio of the sq.rt. to the OPERATORS:", math.sqrt(base10lim)/operators)
print("This is the ratio sqrt(lim) to ACTUAL list ops:", math.sqrt(base10lim)/sum(amplitude_map))
newdict = Counter(amplitude_map)
bigOmega = dict(newdict)
print("This is the distributino of composite types:", bigOmega)



#where it BREAKS: 100000000000000000

# testing @ 1000000000000000
#4333 ops on list based on changed algo


#13
#    drLD(x, 76, -1, 13, 91)   #13,91
#    drLD(x, 94, 18, 19, 67)  #19,67
#    drLD(x, 94, 24, 37, 49)  #37,49
#    drLD(x, 76, 11, 31, 73)  #31,73
#    drLD(x, 86, 6, 11, 83)   #11,83
#    drLD(x, 104, 29, 29, 47) #29,47
#    drLD(x, 86, 14, 23, 71)  #23,71
#    drLD(x, 86, 20, 41, 53)  #41,53
#    drLD(x, 104, 25, 17, 59) #17,59
#    drLD(x, 14, 0, 77, 89)   #77,89
#    drLD(x, 94, 10, 7, 79)   #7,79
#    drLD(x, 76, 15, 43, 61)  #43,61