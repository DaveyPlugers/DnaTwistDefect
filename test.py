import random
import numpy as np
Z = 1000
List = []
String = ''
#String += '\''
N = 148
for g in range(Z):
   for k in range(N):
      Index = random.randint(0,3)
      if Index == 0:
         String += 'A'
      elif Index == 1:
         String += 'T'
      elif Index == 2:
         String += 'C'
      elif Index == 3:
         String += 'G'
   #String += '\''
   #print(String)
   List.append(String)
   String = ""
   #String += '\''
#print(List)


MatlabString = "["
for k in range(len(List)):
   MatlabString += '\'' +  List[k] + '\''
   if k < len(List)-1:
      MatlabString += ';'
   else:
      MatlabString += ']'
print(MatlabString)