#!/usr/bin/env python
#use timidity output.mid to play audio in Terminal
#use timidity output.mid -Ow -o primes.wav in Terminal to convert *.mid to *.wav

from midiutil.MidiFile import MIDIFile
#import py_midicsv #--not working, library seems to have errors



trunc_var = "False" #whether to trim the list by amount = new_limit


#when new_limit > 0, the first enumerated element of the returned list is not to be included as a part of the prime sequence. It is currently included in the output. Later versions will likely remove this datapoint. 
                 #3333336
new_test  = 1000
old_test =  1000

new_limit =  900


test_list = list(range(new_test, old_test+1))

target_var = "Yes"
if target_var == "Yes":
  new_test = int(new_test - new_limit)
  print new_test





A181732 = [91]*int(new_test) #A181732 Numbers n such that 90n + 1 are prime CHECKD
A202116 = [89]*int(new_test) #A202116 Numbers n such that 90n + 89 is prime. CHECKD
A196007 = [83]*int(new_test) #A196007 Numbers n such that 90n + 83 is prime. CHECKD
A202112 = [79]*int(new_test) #A202112 Numbers n such that 90n + 79 is prime. CHECKD
A201822 = [77]*int(new_test) #A201822 Numbers k such that 90*k + 77 is prime. CHECKD
A195993 = [73]*int(new_test) #A195993 Numbers n such that 90n + 73 is prime. CHECKD
A202129 = [71]*int(new_test) #A202129 Numbers n such that 90n + 71 is prime. CHECKD
A201817 = [67]*int(new_test) #A201817 Numbers k such that 90*k + 67 is prime. CHECKD
A202113 = [61]*int(new_test) #A202113 Numbers n such that 90n + 61 is prime. CHECKD
A202101 = [59]*int(new_test) #A202101 Numbers n such that 90*n + 59 is prime. CHECKD
A202114 = [53]*int(new_test) #A202114 Numbers n such that 90n + 53 is prime. CHECKD see:A142316
A201818 = [49]*int(new_test) #A201818 Numbers k such that 90*k + 49 is prime.CHECKD
A201734 = [47]*int(new_test) #A201734 Numbers n such that 90*n + 47 is prime. CHECKD 
A202105 = [43]*int(new_test) #A202105 Numbers n such that 90*n + 43 is prime. CHECKD
A202104 = [41]*int(new_test) #A202104 Numbers n such that 90*n + 41 is prime. CHECKD
A198382 = [37]*int(new_test) #A198382 Numbers n such that 90n + 37 is prime. CHECKD
A201819 = [31]*int(new_test) #A201819 Numbers n such that 90*n + 31 is prime. CHECKD
A201739 = [29]*int(new_test) #A201739 Numbers n such that 90*n + 29 is prime. CHECKD
A201820 = [23]*int(new_test) #A201820 Numbers k such that 90*k + 23 is prime. CHECKD
A196000 = [19]*int(new_test) #A196000 Numbers n such that 90n + 19 is prime. CHECKD
A202115 = [17]*int(new_test) #A202115 Numbers n such that 90n + 17 is prime. CHECKD
A201816 = [13]*int(new_test) #A201816 Numbers k such that 90*k + 13 is prime. CHECKD
A201804 = [11]*int(new_test) #A201804 Numbers k such that 90*k + 11 are prime. CHECKD
A202110 = [7]*int(new_test) #A202110 Numbers n such that 90*n + 7 is prime. CHECKD see:A142315




#A181732 = [91]*int(new_test) #A181732 Numbers n such that 90n + 1 are prime CHECKD
#A202116 = [89]*int(new_test) #A202116 Numbers n such that 90n + 89 is prime. CHECKD
#A196007 = [83]*int(new_test) #A196007 Numbers n such that 90n + 83 is prime. CHECKD
#A202112 = [79]*int(new_test) #A202112 Numbers n such that 90n + 79 is prime. CHECKD
#A201822 = [77]*int(new_test) #A201822 Numbers k such that 90*k + 77 is prime. CHECKD
#A195993 = [73]*int(new_test) #A195993 Numbers n such that 90n + 73 is prime. CHECKD
#A202129 = [71]*int(new_test) #A202129 Numbers n such that 90n + 71 is prime. CHECKD
#A201817 = [67]*int(new_test) #A201817 Numbers k such that 90*k + 67 is prime. CHECKD
#A202113 = [61]*int(new_test) #A202113 Numbers n such that 90n + 61 is prime. CHECKD
#A202101 = [59]*int(new_test) #A202101 Numbers n such that 90*n + 59 is prime. CHECKD
#A202114 = [53]*int(new_test) #A202114 Numbers n such that 90n + 53 is prime. CHECKD see:A142316
#A201818 = [49]*int(new_test) #A201818 Numbers k such that 90*k + 49 is prime.CHECKD
#A201734 = [47]*int(new_test) #A201734 Numbers n such that 90*n + 47 is prime. CHECKD 
#A202105 = [43]*int(new_test) #A202105 Numbers n such that 90*n + 43 is prime. CHECKD
#A202104 = [41]*int(new_test) #A202104 Numbers n such that 90*n + 41 is prime. CHECKD
#A198382 = [37]*int(new_test) #A198382 Numbers n such that 90n + 37 is prime. CHECKD
#A201819 = [31]*int(new_test) #A201819 Numbers n such that 90*n + 31 is prime. CHECKD
#A201739 = [29]*int(new_test) #A201739 Numbers n such that 90*n + 29 is prime. CHECKD
#A201820 = [23]*int(new_test) #A201820 Numbers k such that 90*k + 23 is prime. CHECKD
#A196000 = [19]*int(new_test) #A196000 Numbers n such that 90n + 19 is prime. CHECKD
#A202115 = [17]*int(new_test) #A202115 Numbers n such that 90n + 17 is prime. CHECKD
#A201816 = [13]*int(new_test) #A201816 Numbers k such that 90*k + 13 is prime. CHECKD
#A201804 = [11]*int(new_test) #A201804 Numbers k such that 90*k + 11 are prime. CHECKD
#A202110 = [7]*int(new_test) #A202110 Numbers n such that 90*n + 7 is prime. CHECKD see:A142315




#[81053, 81233, 81773, 81953, 82223, 82493, 82763, 83663, 83843, 83933, 84653, 85103, 85193, 85643, 85733, 86183, 86453, 86813, 86993, 87083, 87443, 87623, 87803, 88523, 88793, 88883, 89153, 89513, 89603, 89783, 89963]




def drLD(x, l, m, z, o, list_var, primitive_var):   
  """This is where we create a family tree / player piano for primes."""   
  y = 90*(x*x) - l*x + m 
  #print x
  if y > new_limit:
    #print x
    #step_var = raw_input("Stop here to recoed limit")
    try:
      list_var[y-(new_limit)] = 0
      pass
    except:
      pass
  if y > old_test: #should be value of original new_test input
    print("This is where the function overshoots", x, primitive_var)
    return
  else:
    pass



  for n in xrange (1, new_test+1):
    w1 = (n+new_limit)-y 
    if w1 <= 0:
      print("What x is this hitting at:", x)
      return
    new_y = ((z+(90*(x-1))))
    if w1%new_y == 0:
      try:
        list_var[n] = 0
        if primitive_var == 59:
          print("This is the composite x", n, new_y, x, z, primitive_var, w1, y) 
        pass
      except:
        if primitive_var == 59:
          print("This failed, too large", x, n, new_y)
          pass
        else:
          pass
        #pass

    new2_y = ((o+(90*(x-1))))  
    if w1%new2_y == 0:
      try:
        list_var[n] = 0
        if primitive_var == 59:
          print("This is the composite x", n, new2_y, x, o, primitive_var, w1, y) 
         
      except:
        if primitive_var == 59:
          print("This failed, too large", x, n, new2_y)
          pass
        #return #trying something new here
      else:
        pass

    if new_y > old_test:
      print("new_y overflows at this n, z, o", n, z, o)
      return

    else:
      pass

# 14 digit range = 33333336
# 13 digit range =  1056666
# 12 digit range =   33336

#for limit 100000000000000 #16 digit base 10
#for 59 @ x_limit = 1000008:
#stops factoring at x=30815 [verify this]



#for limit 1000000000000000 #17 digit base 10
#for 59 @ x_limit = 3333336


# 89999999999999153 = 17 digit [real 102m11.220s] (last 9 positions in seive)
# 8999999999999153 = 16 digit [real 30m18.274s] (last 9 positions in seive)
# 899999999999153 = 333336 = 15 digit [real 10m11.786s]
# 89999999999333 = 105666 = 14-digit

#for the function below
#1000 = 4
#10000 = 11
#100000 = 36
#1000000 = 106
#10000000 = 333
#100000000 = 1055
#1000000000 = 3334
#10000000000 = 10541 [899,999,999,153] real	0m20.076s
#100000000000 = 33334 [8,999,999,999,153] real	1m1.018s
#1000000000000 = 105410 [89,999,999,999,153] real 3m23.014s
#10000000000000 = 333334 [899,999,999,999,153]

for x in xrange(1, 4): # 10000 = 12; 100000 = 36; 1,000,000 = 108; 10,000,000 = 336; 100,000,000 = 1056; 1,000,000,000 = "Hey I need to checksum these tables, all hand-entered
 

#11 = a
    drLD(x, 120, 34, 7, 53, A201804, 11)  #7,53  @4,  154 1
    drLD(x, 132, 48, 19, 29, A201804, 11) #19,29 @6,  144 2
    drLD(x, 120, 38, 17, 43, A201804, 11) #17,43 @8,  158 3
    drLD(x, 90, 11, 13, 77, A201804, 11)  #13,77 @11, 191 4
    drLD(x, 78, -1, 11, 91, A201804, 11)  #11,91 @11, 203 5
    drLD(x, 108, 32, 31, 41, A201804, 11) #31,41 @14, 176 6
    drLD(x, 90, 17, 23, 67, A201804, 11)  #23,67 @17, 197 7
    drLD(x, 72, 14, 49, 59, A201804, 11)  #49,59 @32, 230 8
    drLD(x, 60, 4, 37, 83, A201804, 11)   #37,83 @34, 244 9
    drLD(x, 60, 8, 47, 73, A201804, 11)   #47,73 @38, 248 10
    drLD(x, 48, 6, 61, 71, A201804, 11)   #61,71 @48, 270 11
    drLD(x, 12, 0, 79, 89, A201804, 11)   #79,89 @78, 336 12


#13 = b
    drLD(x, 76, -1, 13, 91, A201816, 13)   #13,91
    drLD(x, 94, 18, 19, 67, A201816, 13)  #19,67
    drLD(x, 94, 24, 37, 49, A201816, 13)  #37,49
    drLD(x, 76, 11, 31, 73, A201816, 13)  #31,73
    drLD(x, 86, 6, 11, 83, A201816, 13)   #11,83
    drLD(x, 104, 29, 29, 47, A201816, 13) #29,47
    drLD(x, 86, 14, 23, 71, A201816, 13)  #23,71
    drLD(x, 86, 20, 41, 53, A201816, 13)  #41,53
    drLD(x, 104, 25, 17, 59, A201816, 13) #17,59
    drLD(x, 14, 0, 77, 89, A201816, 13)   #77,89
    drLD(x, 94, 10, 7, 79, A201816, 13)   #7,79
    drLD(x, 76, 15, 43, 61, A201816, 13) #43,61

#17 = c
    drLD(x, 72, -1, 17, 91, A202115, 17)   #17,91
    drLD(x, 108, 29, 19, 53, A202115, 17) #19,53
    drLD(x, 72, 11, 37, 71, A202115, 17)   #37,71
    drLD(x, 18, 0, 73, 89, A202115, 17)   #73,89
    drLD(x, 102, 20, 11, 67, A202115, 17)  #11,67
    drLD(x, 138, 52, 13, 29, A202115, 17) #13,29
    drLD(x, 102, 28, 31, 47, A202115, 17)  #31,47
    drLD(x, 48, 3, 49, 83, A202115, 17)  #49,83
    drLD(x, 78, 8, 23, 79, A202115, 17)  #23,79
    drLD(x, 132, 45, 7, 41, A202115, 17) #7,41
    drLD(x, 78, 16, 43, 59, A202115, 17)   #43,59
    drLD(x, 42, 4, 61, 77, A202115, 17) #61,77   

#91 = A181732
    drLD(x, -2, 0, 91, 91, A181732, 91) #91,91
    drLD(x, 142, 56, 19, 19, A181732, 91) #19,19
    drLD(x, 70, 10, 37, 73, A181732, 91) #37, 73
    drLD(x, 128, 43, 11, 41, A181732, 91) #11, 41
    drLD(x, 92, 21, 29, 59, A181732, 91) #29,59
    drLD(x, 110, 32, 23, 47, A181732, 91) #23,47
    drLD(x, 20, 1, 77, 83, A181732, 91) #77,83
    drLD(x, 160, 71, 7, 13, A181732, 91) #7,13
    drLD(x, 88, 19, 31, 61, A181732, 91) #31,61
    drLD(x, 52, 5, 49, 79, A181732, 91) #49,79
    drLD(x, 70, 12, 43, 67, A181732, 91) #43,67
    drLD(x, 110, 30, 17, 53, A181732, 91) #17,53
    drLD(x, 38, 4, 71, 71, A181732, 91) #71,71
    drLD(x, 2, 0, 89, 89, A181732, 91) #89,89

# 89 = A202116
    drLD(x, 1, 0, 89, 91, A202116, 89) #89,91
    drLD(x, 90, 14, 19, 71, A202116, 89) #19,71
    drLD(x, 126, 42, 17, 37, A202116, 89) #17,37
    drLD(x, 54, 6, 53, 73, A202116, 89) #53,73
    drLD(x, 120, 35, 11, 49, A202116, 89) #11,49
    drLD(x, 120, 39, 29, 31, A202116, 89) #29,31
    drLD(x, 66, 10, 47, 67, A202116, 89) #47,67
    drLD(x, 84, 5, 13, 83, A202116, 89) #13,83
    drLD(x, 114, 34, 23, 43, A202116, 89) #23,43
    drLD(x, 60, 5, 41, 79, A202116, 89) #41,79
    drLD(x, 60, 9, 59, 61, A202116, 89) #59,61
    drLD(x, 96, 11, 7, 77, A202116, 89) #7,77

# 83 = f
    drLD(x, 6, -1, 83, 91, A196007, 83) #83,91
    drLD(x, 114, 33, 19, 47, A196007, 83) #19,47
    drLD(x, 114, 35, 29, 37, A196007, 83) #37,29
    drLD(x, 96, 14, 11, 73, A196007, 83) #11,73
    drLD(x, 126, 41, 13, 41, A196007, 83) #13,41
    drLD(x, 126, 43, 23, 31, A196007, 83) #23,31
    drLD(x, 54, 5, 49, 77, A196007, 83) #49,77
    drLD(x, 54, 7, 59, 67, A196007, 83) #59,67
    drLD(x, 84, 0, 7, 89, A196007, 83) #7,89
    drLD(x, 66, 9, 43, 71, A196007, 83) #43,71
    drLD(x, 66, 11, 53, 61, A196007, 83) #53,61
    drLD(x, 84, 8, 17, 79, A196007, 83) #17,79



# 79 = A202112
    drLD(x, 10, -1, 79, 91, A202112, 79) #79,91
    drLD(x, 100, 22, 19, 61, A202112, 79) #19,61
    drLD(x, 136, 48, 7, 37, A202112, 79) #7,37
    drLD(x, 64, 8, 43, 73, A202112, 79) #43,73
    drLD(x, 80, 0, 11, 89, A202112, 79) #11,89
    drLD(x, 80, 12, 29, 71, A202112, 79) #29,71
    drLD(x, 116, 34, 17, 47, A202112, 79) #17,47
    drLD(x, 44, 2, 53, 83, A202112, 79) #53,83
    drLD(x, 154, 65, 13, 13, A202112, 79) #13,13
    drLD(x, 100, 26, 31, 49, A202112, 79) #31,49
    drLD(x, 46, 5, 67, 67, A202112, 79) #67,67
    drLD(x, 134, 49, 23, 23, A202112, 79) #23,23
    drLD(x, 80, 16, 41, 59, A202112, 79) #41,59
    drLD(x, 26, 1, 77, 77, A202112, 79) #77,77


# 77 = h
    drLD(x, 12, -1, 77, 91, A201822, 77) #77,91
    drLD(x, 138, 52, 19, 23, A201822, 77) #19,23
    drLD(x, 102, 28, 37, 41, A201822, 77) #37,41
    drLD(x, 48, 5, 59, 73, A201822, 77) #59,73
    drLD(x, 162, 72, 7, 11, A201822, 77) #7,11
    drLD(x, 108, 31, 29, 43, A201822, 77) #29,43
    drLD(x, 72, 13, 47, 61, A201822, 77) #47,61
    drLD(x, 18, 0, 79, 83, A201822, 77) #79,83
    drLD(x, 78, 0, 13, 89, A201822, 77) #13,89
    drLD(x, 132, 47, 17, 31, A201822, 77) #17,31
    drLD(x, 78, 16, 49, 53, A201822, 77) #49,53
    drLD(x, 42, 4, 67, 71, A201822, 77) #67,71

# 73 = A195993
    drLD(x, 16, -1, 73, 91, A195993, 73) #73,91
    drLD(x, 124, 41, 19, 37, A195993, 73) #19,37
    drLD(x, 146, 58, 11, 23, A195993, 73) #11,23
    drLD(x, 74, 8, 29, 77, A195993, 73) #29,77
    drLD(x, 74, 14, 47, 59, A195993, 73) #47,59
    drLD(x, 56, 3, 41, 83, A195993, 73) #41,83
    drLD(x, 106, 24, 13, 61, A195993, 73) #13,61
    drLD(x, 106, 30, 31, 43, A195993, 73) #31,43
    drLD(x, 124, 37, 7, 49, A195993, 73) #7,49
    drLD(x, 34, 2, 67, 79, A195993, 73) #67,79
    drLD(x, 74, 0, 17, 89, A195993, 73) #17,89
    drLD(x, 56, 7, 53, 71, A195993, 73) #53,71


# 71 = j
    drLD(x, 18, -1, 71, 91, A202129, 71) #71,91
    drLD(x, 72, 0, 19, 89, A202129, 71) #19,89
    drLD(x, 90, 21, 37, 53, A202129, 71) #37,53
    drLD(x, 90, 13, 17, 73, A202129, 71) #17,73
    drLD(x, 138, 51, 11, 31, A202129, 71) #11,31
    drLD(x, 102, 27, 29, 49, A202129, 71) #29,49
    drLD(x, 120, 36, 13, 47, A202129, 71) #13,47
    drLD(x, 30, 1, 67, 83, A202129, 71) #67,83
    drLD(x, 150, 61, 7, 23, A202129, 71) #7,23
    drLD(x, 78, 15, 41, 61, A202129, 71) #41,61
    drLD(x, 42, 3, 59, 79, A202129, 71) #59,79
    drLD(x, 60, 6, 43, 77, A202129, 71) #43,77

# 67 = A201817
    drLD(x, 22, -1, 67, 91, A201817, 67) #67,91
    drLD(x, 148, 60, 13, 19, A201817, 67) #13,19
    drLD(x, 112, 34, 31, 37, A201817, 67) #31,37
    drLD(x, 58, 7, 49, 73, A201817, 67) #49,73
    drLD(x, 122, 37, 11, 47, A201817, 67) #11,47
    drLD(x, 68, 4, 29, 83, A201817, 67) #29,83
    drLD(x, 122, 39, 17, 41, A201817, 67) #17,41
    drLD(x, 68, 12, 53, 59, A201817, 67) #53,59
    drLD(x, 32, 2, 71, 77, A201817, 67) #71,77
    drLD(x, 112, 26, 7, 61, A201817, 67) #7,61
    drLD(x, 58, 5, 43, 79, A201817, 67) #43,79
    drLD(x, 68, 0, 23, 89, A201817, 67) #23,89

# 61 = q
    drLD(x, 28, -1, 61, 91, A202113, 61) #61,91
    drLD(x, 82, 8, 19, 79, A202113, 61) #19,79
    drLD(x, 100, 27, 37, 43, A202113, 61) #37,43)
    drLD(x, 100, 15, 7, 73, A202113, 61) #7,73
    drLD(x, 98, 16, 11, 71, A202113, 61) #11,71
    drLD(x, 62, 0, 29, 89, A202113, 61) #29,89
    drLD(x, 80, 17, 47, 53, A202113, 61) #47,53
    drLD(x, 80, 5, 17, 83, A202113, 61) #17,83
    drLD(x, 100, 19, 13, 67, A202113, 61) #13,67
    drLD(x, 118, 38, 31, 31, A202113, 61) #31,31
    drLD(x, 82, 18, 49, 49, A202113, 61) #49,49
    drLD(x, 80, 9, 23, 77, A202113, 61) #23,77
    drLD(x, 98, 26, 41, 41, A202113, 61) #41,41
    drLD(x, 62, 10, 59, 59, A202113, 61) #59,59


# 59 = p
    drLD(x, 30, -1, 59, 91, A202101, 59) #59,91
    drLD(x, 120, 38, 19, 41, A202101, 59) #19,41
    drLD(x, 66, 7, 37, 77, A202101, 59) #37,77
    drLD(x, 84, 12, 23, 73, A202101, 59) #23,73
    drLD(x, 90, 9, 11, 79, A202101, 59) #11,79
    drLD(x, 90, 19, 29, 61, A202101, 59) #29,61
    drLD(x, 126, 39, 7, 47, A202101, 59) #7,47
    drLD(x, 54, 3, 43, 83, A202101, 59) #43,83
    drLD(x, 114, 31, 13, 53, A202101, 59) #13,53
    drLD(x, 60, 0, 31, 89, A202101, 59) #31,89
    drLD(x, 60, 8, 49, 71, A202101, 59) #49,71
    drLD(x, 96, 18, 17, 67, A202101, 59) #17,67

# 53 = r
    drLD(x, 36, -1, 53, 91, A202114, 53) #53,91
    drLD(x, 144, 57, 17, 19, A202114, 53) #17,19
    drLD(x, 54, 0, 37, 89, A202114, 53) #37,89
    drLD(x, 36, 3, 71, 73, A202114, 53) #71,73
    drLD(x, 156, 67, 11, 13, A202114, 53) #11,13
    drLD(x, 84, 15, 29, 67, A202114, 53) #29,67
    drLD(x, 84, 19, 47, 49, A202114, 53) #47,49
    drLD(x, 66, 4, 31, 83, A202114, 53) #31,83
    drLD(x, 96, 21, 23, 61, A202114, 53) #23,61
    drLD(x, 96, 25, 41, 43, A202114, 53) #41,43
    drLD(x, 114, 28, 7, 59, A202114, 53) #7,59
    drLD(x, 24, 1, 77, 79, A202114, 53) #77,79

# 49 = s
    drLD(x, 40, -1, 49, 91, A201818, 49) #49,91
    drLD(x, 130, 46, 19, 31, A201818, 49) #19,31
    drLD(x, 76, 13, 37, 67, A201818, 49) #37,67
    drLD(x, 94, 14, 13, 73, A201818, 49) #13,73
    drLD(x, 140, 53, 11, 29, A201818, 49) #11,29
    drLD(x, 86, 20, 47, 47, A201818, 49) #47,47
    drLD(x, 14, 0, 83, 83, A201818, 49) #83,83
    drLD(x, 104, 27, 23, 53, A201818, 49) #23,53
    drLD(x, 50, 0, 41, 89, A201818, 49) #41,89
    drLD(x, 50, 6, 59, 71, A201818, 49) #59,71
    drLD(x, 86, 10, 17, 77, A201818, 49) #17,77
    drLD(x, 166, 76, 7, 7, A201818, 49) #7,7
    drLD(x, 94, 24, 43, 43, A201818, 49) #43,43
    drLD(x, 40, 3, 61, 79, A201818, 49) #61,79

# 47 = t
    drLD(x, 42, -1, 47, 91, A201734, 47) #47,91
    drLD(x, 78, 5, 19, 83, A201734, 47) #19,83
    drLD(x, 132, 46, 11, 37, A201734, 47) #11,37
    drLD(x, 78, 11, 29, 73, A201734, 47) #29,73
    drLD(x, 108, 26, 13, 59, A201734, 47) #13,59
    drLD(x, 72, 8, 31, 77, A201734, 47) #31,77
    drLD(x, 108, 30, 23, 49, A201734, 47) #23,49
    drLD(x, 102, 17, 7, 71, A201734, 47) #7,71
    drLD(x, 48, 0, 43, 89, A201734, 47) #43,89
    drLD(x, 102, 23, 17, 61, A201734, 47) #17,61
    drLD(x, 48, 4, 53, 79, A201734, 47) #53,79
    drLD(x, 72, 12, 41, 67, A201734, 47) #41,67


# 43 = u
    drLD(x, 46, -1, 43, 91, A202105, 43) #43,91
    drLD(x, 154, 65, 7, 19, A202105, 43) #7,19
    drLD(x, 64, 6, 37, 79, A202105, 43) #37,79
    drLD(x, 46, 5, 61, 73, A202105, 43) #61,73
    drLD(x, 116, 32, 11, 53, A202105, 43) #11,53
    drLD(x, 134, 49, 17, 29, A202105, 43) #17,29
    drLD(x, 44, 0, 47, 89, A202105, 43) #47,89
    drLD(x, 26, 1, 71, 83, A202105, 43) #71,83
    drLD(x, 136, 50, 13, 31, A202105, 43) #13,31
    drLD(x, 64, 10, 49, 67, A202105, 43) #49,67
    drLD(x, 116, 36, 23, 41, A202105, 43) #23,41
    drLD(x, 44, 4, 59, 77, A202105, 43) #59,77

# 41 = v
    drLD(x, 48, -1, 41, 91, A202104, 41) #41,91
    drLD(x, 42, 0, 49, 89, A202104, 41) #49,89
    drLD(x, 102, 24, 19, 59, A202104, 41) #19,59
    drLD(x, 120, 39, 23, 37, A202104, 41) #23,37
    drLD(x, 108, 25, 11, 61, A202104, 41) #11,61
    drLD(x, 72, 7, 29, 79, A202104, 41) #29,79
    drLD(x, 90, 22, 43, 47, A202104, 41) #43,47
    drLD(x, 150, 62, 13, 17, A202104, 41) #13,17
    drLD(x, 78, 12, 31, 71, A202104, 41) #31,71
    drLD(x, 30, 2, 73, 77, A202104, 41) #73, 77
    drLD(x, 60, 9, 53, 67, A202104, 41) #53,67
    drLD(x, 90, 6, 7, 83, A202104, 41) #7,83

# 37 = w
    drLD(x, 52, -1, 37, 91, A198382, 37) #37,91
    drLD(x, 88, 13, 19, 73, A198382, 37) #19,73
    drLD(x, 92, 11, 11, 77, A198382, 37) #11,77
    drLD(x, 128, 45, 23, 29, A198382, 37) #23,29
    drLD(x, 92, 23, 41, 47, A198382, 37) #41,47
    drLD(x, 38, 2, 59, 83, A198382, 37) #59,83
    drLD(x, 88, 9, 13, 79, A198382, 37) #13,79
    drLD(x, 142, 54, 7, 31, A198382, 37) #7,31
    drLD(x, 88, 21, 43, 49, A198382, 37) #43,49
    drLD(x, 52, 7, 61, 67, A198382, 37) #61,67
    drLD(x, 92, 15, 17, 71, A198382, 37) #17,71
    drLD(x, 38, 0, 53, 89, A198382, 37) #53,89

# 31 = ww
    drLD(x, 58, -1, 31, 91, A201819, 31) #31,91
    drLD(x, 112, 32, 19, 49, A201819, 31) #19,49
    drLD(x, 130, 45, 13, 37, A201819, 31) #13,37
    drLD(x, 40, 4, 67, 73, A201819, 31) #67,73
    drLD(x, 158, 69, 11, 11, A201819, 31) #11,11
    drLD(x, 122, 41, 29, 29, A201819, 31) #29,29
    drLD(x, 50, 3, 47, 83, A201819, 31) #47,83
    drLD(x, 140, 54, 17, 23, A201819, 31) #17,23
    drLD(x, 68, 10, 41, 71, A201819, 31) #41,71
    drLD(x, 32, 0, 59, 89, A201819, 31) #59,89
    drLD(x, 50, 5, 53, 77, A201819, 31) #53,77
    drLD(x, 130, 43, 7, 43, A201819, 31) #7,43
    drLD(x, 58, 9, 61, 61, A201819, 31) #61,61
    drLD(x, 22, 1, 79, 79, A201819, 31) #79,79

# 29 = www
    drLD(x, 60, -1, 29, 91, A201739, 29) #29,91
    drLD(x, 150, 62, 11, 19, A201739, 29) #11,19
    drLD(x, 96, 25, 37, 47, A201739, 29) #37,47
    drLD(x, 24, 1, 73, 83, A201739, 29) #73,83
    drLD(x, 144, 57, 13, 23, A201739, 29) #13,23
    drLD(x, 90, 20, 31, 59, A201739, 29) #31,59
    drLD(x, 90, 22, 41, 49, A201739, 29) #41,49
    drLD(x, 36, 3, 67, 77, A201739, 29) #67,77
    drLD(x, 156, 67, 7, 17, A201739, 29) #7,17
    drLD(x, 84, 19, 43, 53, A201739, 29) #43,53
    drLD(x, 30, 0, 61, 89, A201739, 29) #61,89
    drLD(x, 30, 2, 71, 79, A201739, 29) #71,79

# 23 = wx
    drLD(x, 66, -1, 23, 91, A201820, 23) #23,91
    drLD(x, 84, 10, 19, 77, A201820, 23) #19,77
    drLD(x, 84, 18, 37, 59, A201820, 23) #37,59
    drLD(x, 66, 9, 41, 73, A201820, 23) #41,73
    drLD(x, 126, 41, 11, 43, A201820, 23) #11,43
    drLD(x, 144, 56, 7, 29, A201820, 23) #7,29
    drLD(x, 54, 5, 47, 79, A201820, 23) #47,79
    drLD(x, 36, 2, 61, 83, A201820, 23) #61,83
    drLD(x, 96, 16, 13, 71, A201820, 23) #13,71
    drLD(x, 96, 24, 31, 53, A201820, 23) #31,53
    drLD(x, 114, 33, 17, 49, A201820, 23) #17,49
    drLD(x, 24, 0, 67, 89, A201820, 23) #67,89

# 19 = gg
    drLD(x, 70, -1, 19, 91, A196000, 19) #19,91
    drLD(x, 106, 31, 37, 37, A196000, 19) #37,73
    drLD(x, 34, 3, 73, 73, A196000, 19) #73,73
    drLD(x, 110, 27, 11, 59, A196000, 19) #11,59
    drLD(x, 110, 33, 29, 41, A196000, 19) #29,41
    drLD(x, 56, 6, 47, 77, A196000, 19) #47,77
    drLD(x, 74, 5, 23, 83, A196000, 19) #23,83
    drLD(x, 124, 40, 13, 43, A196000, 19) #13,43
    drLD(x, 70, 7, 31, 79, A196000, 19) #31,79
    drLD(x, 70, 13, 49, 61, A196000, 19) #49,61
    drLD(x, 106, 21, 7, 67, A196000, 19) #7,67
    drLD(x, 20, 0, 71, 89, A196000, 19) #71,89
    drLD(x, 74, 15, 53, 53, A196000, 19) #53,53
    drLD(x, 146, 59, 17, 17, A196000, 19) #17,17


# 7 = xy
    drLD(x, 82, -1, 7, 91, A202110, 7) #7,91
    drLD(x, 118, 37, 19, 43, A202110, 7) #19,43
    drLD(x, 82, 17, 37, 61, A202110, 7) #37,61
    drLD(x, 28, 2, 73, 79, A202110, 7) #73,79
    drLD(x, 152, 64, 11, 17, A202110, 7) #11,17
    drLD(x, 98, 25, 29, 53, A202110, 7) #29,53
    drLD(x, 62, 9, 47, 71, A202110, 7) #47,71
    drLD(x, 8, 0, 83, 89, A202110, 7) #83,89
    drLD(x, 118, 35, 13, 49, A202110, 7) #13,49
    drLD(x, 82, 15, 31, 67, A202110, 7) #31,67
    drLD(x, 98, 23, 23, 59, A202110, 7) #23,59
    drLD(x, 62, 7, 41, 77, A202110, 7) #41,77



print("Let's make some music!")

##### LET'S MAKE SOME MIDI MUSIC #####

# Create the MIDIFile Object with 24 tracks
MyMIDI = MIDIFile(23)

#Create Instrument Variable Declarations
#example variable:
instrument_var = 105 #a sitar

#piano
instrument_var_acousticgrand = 1 #acoustic grand piano
instrument_var_brightacousticpiano = 2
instrument_var_electricgrandpiano = 3
instrument_var_honkytonkpiano = 4
instrument_var_electricpiano1 = 5
instrument_var_electricpiano2 = 6
instrument_var_harpsichord = 7
instrument_var_clavinet = 8 #Clavinet

#chromatic percussion
instrument_var_celestra = 9
instrument_var_glockenspiel = 10 #Glockenspiel
instrument_var_musicbox = 11
instrument_var_vibraphone = 12
instrument_var_marima = 13
instrument_var_xylophone = 14 #xylophone
instrument_var_tubularbells = 15
instrument_var_dulcimer = 16


#organ
instrument_var_drawbarorgan = 17
instrument_var_percussiveorgan = 18
instrument_var_rockorgan = 19
instrument_var_churchorgan = 20
instrument_var_reedorgan = 21
instrument_var_accordion = 22
instrument_var_harmonica = 23
instrument_var_tangoaccordion = 34

#guitar
instrument_var_acousticnylonguitar = 25
instrument_var_acousticsteelguitar = 26
instrument_var_electricjazzguitar = 27
instrument_var_electriccleanguitar = 28
instrument_var_electricmutedguitar = 29
instrument_var_overdrivenguitar = 30
instrument_var_distortionguitar = 31
instrument_var_harmonicsguitar = 32

#bass
instrument_var_acousticbass = 33
instrument_var_fingerelectricbass = 34
instrument_var_pickelectricbass = 35
instrument_var_fretlessbass = 36
instrument_var_slap1bass = 37
instrument_var_slap2bass = 38
instrument_var_synth1bass = 38
instrument_var_synth2bass = 40

#strings
instrument_var_violin = 41 #violin
instrument_var_viola = 42 #viola
instrument_var_cello = 43 #Cello
instrument_var_contrabass = 44 #Contrabase
instrument_var_tremulostrings = 45
instrument_var_pizzicatostrings = 46
instrument_var_orchestralharp = 47 #orchestral harp
instrument_var_timpani = 48
instrument_var_stringensemble1 = 49 #string ensemble
instrument_var_stringensemble2 = 50
instrument_var_synthstrings1 = 51
instrument_var_synthstrings2 = 52
instrument_var_choralahh = 53
instrument_var_voiceooh = 54
instrument_var_synthvoice = 55
instrument_var_orchestrahit = 56

#brass
instrument_var_trumpet = 57 #trumpet
instrument_var_trombone = 58 #trombone
instrument_var_tuba = 59
instrument_var_mutedtrumpet = 60
instrument_var_frenchhorn = 61 #french horn
instrument_var_brasssection = 62
instrument_var_synthbrass1 = 63
instrument_var_synthbrass2 = 64

#reed
instrument_var_sopranosax = 65 #sopranosax
instrument_var_altosax = 66
instrument_var_tenorsax = 67
instrument_var_baritonesax = 68
instrument_var_oboe = 69 #Oboe
instrument_var_englishhorn = 70 #English Horn
instrument_var_bassoon = 71 #bassoon
instrument_var_clarinet = 72

#pipe
instrument_var_piccolo = 73
instrument_var_flute = 74
instrument_var_recorder = 75
instrument_var_panflute = 76
instrument_var_blownbottle = 77
instrument_var_shakuhachi = 78
instrument_var_whistle = 79
instrument_var_ocarina = 80

#synthlead
instrument_var_square = 81
instrument_var_sawtooth = 82
instrument_var_calliope = 83
instrument_var_chiff = 84
instrument_var_charang = 85
instrument_var_voice = 86
instrument_var_fifths = 87
instrument_var_baselead = 88

#synthpad
instrument_var_newagepad = 89 
instrument_var_warmpad = 90
instrument_var_polysynthpad = 91
instrument_var_choirpad = 92
instrument_var_bowedpad = 93
instrument_var_metallicpad = 94 
instrument_var_halopad = 95
instrument_var_sweeppad = 96

#syntheffects
instrument_var_rain = 97
instrument_var_soundtrack = 98
instrument_var_crystal = 99
instrument_var_atmosphere =100 
instrument_var_brightness = 101
instrument_var_goblins = 102
instrument_var_echoes = 103
instrument_var_scifi = 104

#ethnic
instrument_var_sitar = 105
instrument_var_banjo = 106
instrument_var_shamisen = 107
instrument_var_koto = 108
instrument_var_kalimba = 109
instrument_var_bagpipe = 110
instrument_var_fiddle = 111
instrument_var_shanai = 112

#percussive
instrument_var_tinklebell = 113
instrument_var_agogo = 114
instrument_var_steeldrums = 115
instrument_var_woodblock = 116
instrument_var_taikodrum = 117
instrument_var_melodictom = 118
instrument_var_synthdrum = 119

#soundeffects
instrument_var_reversesymbal = 120 
instrument_var_guitarfretnoise = 121
instrument_var_breathnoise = 122
instrument_var_seashore = 123
instrument_var_birdtweet = 124
instrument_var_telephonering = 125 
instrument_var_helicopter = 126
instrument_var_applause = 127
instrument_var_gunshot = 128




MyMIDI.addProgramChange(0, 0, 0, int(instrument_var_electricpiano1))
MyMIDI.addProgramChange(0, 1, 0, int(instrument_var_electricpiano1))
MyMIDI.addProgramChange(0, 2, 0, int(instrument_var_electricpiano1))
MyMIDI.addProgramChange(0, 3, 0, int(instrument_var_electricpiano1))
MyMIDI.addProgramChange(0, 4, 0, int(instrument_var_electricpiano1))
MyMIDI.addProgramChange(0, 5, 0, int(instrument_var_electricpiano1))
MyMIDI.addProgramChange(0, 6, 0, int(instrument_var_electricpiano1))
MyMIDI.addProgramChange(0, 7, 0, int(instrument_var_electricpiano1))
MyMIDI.addProgramChange(0, 8, 0, int(instrument_var_electricpiano1))
MyMIDI.addProgramChange(0, 9, 0, int(instrument_var_electricpiano1))
MyMIDI.addProgramChange(0, 10, 0, int(instrument_var_electricpiano1))
MyMIDI.addProgramChange(0, 11, 0, int(instrument_var_electricpiano1))
MyMIDI.addProgramChange(0, 12, 0, int(instrument_var_electricpiano1))
MyMIDI.addProgramChange(0, 13, 0, int(instrument_var_electricpiano1))
MyMIDI.addProgramChange(0, 14, 0, int(instrument_var_electricpiano1))
MyMIDI.addProgramChange(0, 15, 0, int(instrument_var_electricpiano1))
MyMIDI.addProgramChange(0, 16, 0, int(instrument_var_electricpiano1))
MyMIDI.addProgramChange(0, 17, 0, int(instrument_var_electricpiano1))
MyMIDI.addProgramChange(0, 18, 0, int(instrument_var_electricpiano1))
MyMIDI.addProgramChange(0, 19, 0, int(instrument_var_electricpiano1))
MyMIDI.addProgramChange(0, 20, 0, int(instrument_var_electricpiano1))
MyMIDI.addProgramChange(0, 21, 0, int(instrument_var_electricpiano1))
MyMIDI.addProgramChange(0, 22, 0, int(instrument_var_electricpiano1))
MyMIDI.addProgramChange(0, 23, 0, int(instrument_var_electricpiano1))




track1 = 0
channel1 = 0
track2 = 0
channel2 = 0
track3 = 0
channel3 = 0
track4 = 0
channel4 = 0
track5 = 0
channel5 = 0
track6 = 0
channel6 = 0
track7 = 0
channel7 = 0
track8 = 0
channel8 = 0
track9 = 0
channel9 = 0
track10 = 0
channel10 = 0
track11 = 0
channel11 = 0
track12 = 0
channel12 = 0
track13 = 0
channel13 = 0
track14 = 0
channel14 = 0
track15 = 0
channel15 = 0
track16 = 0
channel16 = 0
track17 = 0
channel17 = 0
track18 = 0
channel18 = 0
track19 = 0
channel19 = 0
track20 = 0
channel20 = 0
track21 = 0
channel21 = 0
track22 = 0
channel22 = 0
track23 = 0
channel23 = 0
track24 = 0
channel24 = 0

#set your variable global time variable. The beat the write output starts at
time_var = 0
#time_var = 0
#time_var_delta = 111 #example

#duration variable setup. This is how long you hold the note
duration_var = 1

#volume variable setup. This is how loud it is "recorded"
volume_var = 60

#pitch variable setup #leave this at 0 to ignore 0 on midi. Above 0, this gets 0 as noise
pitch_var = 0

#example time variable
#time2_var = 10000

pitch1 = int(pitch_var)
time1 = int(time_var) # OR int(time2_var)
duration1 = int(duration_var) 
volume1 = int(volume_var)

pitch2 = int(pitch_var)
time2 = int(time_var)
duration2 = int(duration_var)
volume2 = int(volume_var)

pitch3 = int(pitch_var)
time3 = int(time_var)
duration3 = int(duration_var)
volume3 = int(volume_var)

pitch4 = int(pitch_var)
time4 = int(time_var)
duration4 = int(duration_var)
volume4 = int(volume_var)

pitch5 = int(pitch_var)
time5 = int(time_var)
duration5 = int(duration_var)
volume5 = int(volume_var)

pitch6 = int(pitch_var)
time6 = int(time_var)
duration6 = int(duration_var)
volume6 = int(volume_var)

pitch7 = int(pitch_var)
time7 = int(time_var)
duration7 = int(duration_var)
volume7 = int(volume_var)

pitch8 = int(pitch_var)
time8 = int(time_var)
duration8 = int(duration_var)
volume8 = int(volume_var)

pitch9 = int(pitch_var)
time9 = int(time_var)
duration9 = int(duration_var)
volume9 = int(volume_var)

pitch10 = int(pitch_var)
time10 = int(time_var)
duration10 = int(duration_var)
volume10 = int(volume_var)

pitch11 = int(pitch_var)
time11 = int(time_var)
duration11 = int(duration_var)
volume11 = int(volume_var)

pitch12 = int(pitch_var)
time12 = int(time_var)
duration12 = int(duration_var)
volume12 = int(volume_var)

pitch13 = int(pitch_var)
time13 = int(time_var)
duration13 = int(duration_var)
volume13 = int(volume_var)

pitch14 = int(pitch_var)
time14 = int(time_var)
duration14 = int(duration_var)
volume14 = int(volume_var)

pitch15 = int(pitch_var)
time15 = int(time_var)
duration15 = int(duration_var)
volume15 = int(volume_var)

pitch16 = int(pitch_var)
time16 = int(time_var)
duration16 = int(duration_var)
volume16 = int(volume_var)

pitch17 = int(pitch_var)
time17 = int(time_var)
duration17 = int(duration_var)
volume17 = int(volume_var)

pitch18 = int(pitch_var)
time18 = int(time_var)
duration18 = int(duration_var)
volume18 = int(volume_var)

pitch19 = int(pitch_var)
time19 = int(time_var)
duration19 = int(duration_var)
volume19 = int(volume_var)

pitch20 = int(pitch_var)
time20 = int(time_var)
duration20 = int(duration_var)
volume20 = int(volume_var)

pitch21 = int(pitch_var)
time21 = int(time_var)
duration21 = int(duration_var)
volume21 = int(volume_var)

pitch22 = int(pitch_var)
time22 = int(time_var)
duration22 = int(duration_var)
volume22 = int(volume_var)

pitch23 = int(pitch_var)
time23 = int(time_var)
duration23 = int(duration_var)
volume23 = int(volume_var)

pitch24 = int(pitch_var)
time24 = int(time_var)
duration24 = int(duration_var)
volume24 = int(volume_var)




#set the maximum number of beats to process
time_limit_var = 10000

#setup track inputs from our lists
print("The truncation variable is equal to:", trunc_var)
if trunc_var == "Yes":
  print("Let's delete some list elements using del()")
  var = raw_input("testing here at trunc_var")
  del A181732[:new_limit]
  del A202116[:new_limit]
  del A196007[:new_limit]
  del A202112[:new_limit]
  del A201822[:new_limit]
  del A195993[:new_limit]
  del A202129[:new_limit]
  del A201817[:new_limit]
  del A202113[:new_limit]
  del A202101[:new_limit]
  del A202114[:new_limit]
  del A201818[:new_limit]
  del A201734[:new_limit]
  del A202105[:new_limit]
  del A202104[:new_limit]
  del A198382[:new_limit]
  del A201819[:new_limit]
  del A201739[:new_limit]
  del A201820[:new_limit]
  del A196000[:new_limit]
  del A202115[:new_limit]
  del A201816[:new_limit]
  del A201804[:new_limit]
  del A202110[:new_limit]
  pass















for digit in A201804:
    #if digit == '.':
    #    continue
    MyMIDI.addNote(track1, channel1, pitch1 + int(digit), time1, duration1, volume1)
    time1 += 1
    if time1 == int(time_limit_var):
        pass


for digit in A201816:
    #if digit == '.':
    #    continue
    #print digit
    MyMIDI.addNote(track2, channel2, pitch2 + int(digit), time2, duration2, volume2)
    time2 += 1
    if time2 == int(time_limit_var):
        pass



for digit in A202115:
    #if digit == '.':
    #    continue
    MyMIDI.addNote(track3, channel3, pitch3 + int(digit), time3, duration3, volume3)
    time3 += 1
    if time3 == int(time_limit_var):
        pass



for digit in A181732:
    #if digit == '.':
    #    continue
    MyMIDI.addNote(track4, channel4, pitch4 + int(digit), time4, duration4, volume4)
    time4 += 1
    if time4 == int(time_limit_var):
        pass


for digit in A202116:
    #if digit == '.':
    #    continue
    MyMIDI.addNote(track5, channel5, pitch5 + int(digit), time5, duration5, volume5)
    time5 += 1
    if time5 == int(time_limit_var):
        pass



for digit in A196007:
    #if digit == '.':
    #    continue
    MyMIDI.addNote(track6, channel6, pitch6 + int(digit), time6, duration6, volume6)
    time6 += 1
    if time6 == int(time_limit_var):
        pass



for digit in A202112:
    #if digit == '.':
    #    continue
    MyMIDI.addNote(track7, channel7, pitch7 + int(digit), time7, duration7, volume7)
    time7 += 1
    if time7 == int(time_limit_var):
        pass

for digit in A201822:
    #if digit == '.':
    #    continue
    MyMIDI.addNote(track8, channel8, pitch8 + int(digit), time8, duration8, volume8)
    time8 += 1
    if time8 == int(time_limit_var):
        pass


for digit in A195993:
    #if digit == '.':
    #    continue
    MyMIDI.addNote(track9, channel9, pitch9 + int(digit), time9, duration9, volume9)
    time9 += 1
    if time9 == int(time_limit_var):
        pass



for digit in A202129:
    #if digit == '.':
    #    continue
    MyMIDI.addNote(track10, channel10, pitch10 + int(digit), time10, duration10, volume10)
    time10 += 1
    if time10 == int(time_limit_var):
        pass


for digit in A201817:
    #if digit == '.':
    #    continue
    MyMIDI.addNote(track11, channel11, pitch11 + int(digit), time11, duration11, volume11)
    time11 += 1
    if time11 == int(time_limit_var):
        pass



for digit in A202113:
    #if digit == '.':
    #    continue
    MyMIDI.addNote(track12, channel12, pitch12 + int(digit), time12, duration12, volume12)
    time12 += 1
    if time12 == int(time_limit_var):
        pass


for digit in A202101:
    #if digit == '.':
    #    continue
    MyMIDI.addNote(track13, channel13, pitch13 + int(digit), time13, duration13, volume13)
    time13 += 1
    if time13 == int(time_limit_var):
        pass


for digit in A202114:
    #if digit == '.':
    #    continue
    MyMIDI.addNote(track14, channel14, pitch14 + int(digit), time14, duration14, volume14)
    time14 += 1
    if time14 == int(time_limit_var):
        pass


for digit in A201818:
    #if digit == '.':
    #    continue
    MyMIDI.addNote(track15, channel15, pitch15 + int(digit), time15, duration15, volume15)
    time15 += 1
    if time15 == int(time_limit_var):
        pass



for digit in A201734:
    #if digit == '.':
    #    continue
    MyMIDI.addNote(track16, channel16, pitch16 + int(digit), time16, duration16, volume16)
    time16 += 1
    if time16 == int(time_limit_var):
        pass


for digit in A202105:
    #if digit == '.':
    #    continue
    MyMIDI.addNote(track17, channel17, pitch17 + int(digit), time17, duration17, volume17)
    time17 += 1
    if time17 == int(time_limit_var):
        pass


for digit in A202104:
    #if digit == '.':
    #    continue
    MyMIDI.addNote(track18, channel18, pitch18 + int(digit), time18, duration18, volume18)
    time18 += 1
    if time18 == int(time_limit_var):
        pass


for digit in A198382:
    #if digit == '.':
    #    continue
    MyMIDI.addNote(track19, channel19, pitch19 + int(digit), time19, duration19, volume19)
    time19 += 1
    if time19 == int(time_limit_var):
        pass


for digit in A201819:
    #if digit == '.':
    #    continue
    MyMIDI.addNote(track20, channel20, pitch20 + int(digit), time20, duration20, volume20)
    time20 += 1
    if time20 == int(time_limit_var):
        pass


for digit in A201739:
    #if digit == '.':
    #    continue
    MyMIDI.addNote(track21, channel21, pitch21 + int(digit), time21, duration21, volume21)
    time21 += 1
    if time21 == int(time_limit_var):
        pass


for digit in A201820:
    #if digit == '.':
    #    continue
    MyMIDI.addNote(track22, channel22, pitch22 + int(digit), time22, duration22, volume22)
    time22 += 1
    if time22 == int(time_limit_var):
        pass


for digit in A196000:
    #if digit == '.':
    #    continue
    MyMIDI.addNote(track23, channel23, pitch23 + int(digit), time23, duration23, volume23)
    time23 += 1
    if time23 == int(time_limit_var):
        pass


for digit in A202110:
    #if digit == '.':
    #    continue
    MyMIDI.addNote(track24, channel24, pitch24 + int(digit), time24, duration24, volume24)
    time24 += 1
    if time24 == int(time_limit_var):
        break




# And write it to disk.
binfile = open("output.mid", 'wb')
MyMIDI.writeFile(binfile)
binfile.close()



print(A202114)
list_A = [i for i,b in enumerate(A202114) if b == 53]
#list_B = [((x*90)+7) for x in list_A] #compare to http://oeis.org/A142317
list_B = [((x1+(new_limit))*90)+53 for x1 in list_A] #compare to http://oeis.org/A142317/b142317.txt
print list_A
print list_B
print len(list_B)
#print len(list_B)



#composites_base10.sort(reverse = True)
#print(composites_base10)
#zed = [x for x in composite_base10]



#print(gg)
##list_A = [i for i,b in enumerate(gg) if b == 0]
#list_B = [((x+(new_limit))*90)+19 for x in list_A] #compare to http://oeis.org/A142317/b142317.txt
#print list_A
#print list_B
#print len(list_A)
#print len(list_B)
#print sorted(composites_base10)
#print a  #uncomment this line to see the full list (None=primes and 0=composite)
#print list_A201804
#print "%d is the %dth term of Sloane's A201804." % (list_A201804[-1], len(list_A201804))




# Load the MIDI file and parse it into CSV format
#csv_string = py_midicsv.midi_to_csv("output.mid")
#print csv_string
# Parse the CSV output of the previous command back into a MIDI file
#midi_object = py_midicsv.csv_to_midi(csv_string)

# Save the parsed MIDI file to disk
#with open("example_converted.mid", "wb") as output_file:
#    midi_writer = py_midicsv.FileWriter(output_file)
#    midi_writer.write(midi_object)







