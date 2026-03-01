import cmath
import math

def compute_holes(h):
    limit = h
    epoch = 90*(h*h) - 12*h + 1
    limit = epoch
    a = 90
    b = -300
    c = 250 - limit 
    d = (b**2) - (4*a*c)
    sol1 = (-b-cmath.sqrt(d))/(2*a)
    sol2 = (-b+cmath.sqrt(d))/(2*a)
    new_limit = sol2
    list17 = [0]*int(limit+100) 
    def drLD(x, l, m, z, listvar): 
      y = 90*(x*x) - l*x + m
      if 0 <= y < len(listvar):
          listvar[y] = listvar[y]+1   
      p = z+(90*(x-1))
      for n in range (1, int(((limit-y)/p)+1)):  
        yy = y+(p*n)
        if 0 <= yy < len(listvar):
            listvar[yy] = listvar[yy]+1

    for x in range(1, int(new_limit.real)+1):
        drLD(x, 72, -1, 17, list17)   
        drLD(x, 72, -1, 91, list17)   
        
        drLD(x, 108, 29, 19, list17) 
        drLD(x, 108, 29, 53, list17) 
        
        drLD(x, 72, 11, 37, list17)   
        drLD(x, 72, 11, 71, list17)   
        
        drLD(x, 18, 0, 73, list17)   
        drLD(x, 18, 0, 89, list17)   
        
        drLD(x, 102, 20, 11, list17)  
        drLD(x, 102, 20, 67, list17)  
        
        drLD(x, 138, 52, 13, list17) 
        drLD(x, 138, 52, 29, list17) 
        
        drLD(x, 102, 28, 31, list17)  
        drLD(x, 102, 28, 47, list17)  
        
        drLD(x, 48, 3, 49, list17)  
        drLD(x, 48, 3, 83, list17)  
        
        drLD(x, 78, 8, 23, list17)  
        drLD(x, 78, 8, 79, list17)  
        
        drLD(x, 132, 45, 7, list17) 
        drLD(x, 132, 45, 41, list17) 
        
        drLD(x, 78, 16, 43, list17)   
        drLD(x, 78, 16, 59, list17)   
        
        drLD(x, 42, 4, 61, list17) 
        drLD(x, 42, 4, 77, list17) 

    list17 = list17[:int(limit)] 
    list17a = [i for i,x in enumerate(list17) if x == 0] 
    quantity = len(list17a)
    density = quantity / epoch if epoch > 0 else 0
    return h, epoch, quantity, density

hs = [5, 10, 20, 50, 100, 150, 200, 210, 220, 230, 240, 250]
results = []
for h in hs:
    results.append(compute_holes(h))
print(results)

import math
def formula_density(h):
    return 3.75 / (9 + 2 * math.log(h)) if h > 1 else 0

for r in results:
    print(f"h={r[0]}, empirical density={r[3]}, formula={formula_density(r[0])}, error={abs(r[3] - formula_density(r[0])) / formula_density(r[0])}")