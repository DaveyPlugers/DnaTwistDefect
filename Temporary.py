import numpy as np
import random
import matplotlib.pyplot as plt
#Taken from position 51 to 60, repeat at modulo 10


AA_Array = np.array([0.0159,    0.0311,   0.0646,    0.1152,    0.2071 ,   0.3038  ,  0.2520 ,   0.1153  ,  0.0417 ,   0.0170])
TA_Array = np.array([0.0146,    0.0177,    0.0177,    0.0146,    0.0212,    0.0657,    0.1706,    0.1706,    0.0657,    0.0212])
CA_Array = np.array([0.0449,    0.0354,    0.0339,    0.0359 ,   0.0392 ,   0.0454 ,   0.0601 ,   0.0821   , 0.0845  ,  0.0630])
GA_Array = np.array([0.0086,    0.0322,    0.0875,    0.1173,    0.0853,    0.0409,    0.0169 ,   0.0079,    0.0046 ,   0.0040])

AT_Array = np.array([0.0238,    0.0507 ,   0.0507 ,   0.0238  ,  0.0161 ,   0.0231  ,  0.0361   , 0.0361 ,   0.0231  ,  0.0161])
TT_Array = np.array([0.1152,    0.0646,    0.0311 ,   0.0159 ,   0.0170  ,  0.0417  ,  0.1153,   0.2520  ,  0.3038  ,  0.2071])
CT_Array = np.array([0.1603,    0.0944,    0.0433 ,   0.0174  ,  0.0094 ,   0.0095 ,   0.0181 ,   0.0542   , 0.1304 ,   0.1878])
GT_Array = np.array([0.0534,    0.0733 ,   0.0786 ,   0.0592 ,   0.0415 ,   0.0309 ,   0.0269 ,   0.0335 ,   0.0423 ,   0.0447])

AC_Array = np.array([0.0592 ,   0.0786 ,   0.0733   , 0.0534   , 0.0447 ,   0.0423 ,   0.0335  ,  0.0269  ,  0.0309  ,  0.0415])
TC_Array = np.array([0.1173 ,   0.0875 ,   0.0322 ,   0.0086 ,   0.0040   , 0.0046  ,  0.0079  ,  0.0169  ,  0.0409  ,  0.0853])
CC_Array = np.array([0.1229 ,   0.1143 ,   0.0758  ,  0.0381  , 0.0202  ,  0.0159  ,  0.0182  ,  0.0282  ,  0.0502   , 0.0894])
GC_Array = np.array([0.1050 ,   0.1151  ,  0.1151 ,   0.1050   , 0.0900   , 0.0702  , 0.0522  ,  0.0522  ,  0.0702   , 0.0900])

AG_Array = np.array([0.0174  ,  0.0433   , 0.0944  ,  0.1603  ,  0.1878  ,  0.1304  ,  0.0542  ,  0.0181 ,   0.0095   , 0.0094])
TG_Array = np.array([0.0359 ,   0.0339   , 0.0354  ,  0.0449  ,  0.0630  ,  0.0845 ,   0.0821   , 0.0601 ,   0.0454 ,   0.0392])
CG_Array = np.array([0.0674  ,  0.0521  ,  0.0521  ,  0.0674  ,  0.0642  ,  0.0410 ,   0.0278  ,  0.0278  ,  0.0410  ,  0.0642])
GG_Array = np.array([0.0381  ,  0.0758  ,  0.1143    ,0.1229  ,  0.0894  ,  0.0502  ,  0.0282  ,  0.0182   , 0.0159  ,  0.0202])

#plt.plot(AA_Array)
#plt.plot(TT_Array)
#plt.plot(TA_Array)
#plt.show()

ASummed_Array = AA_Array + AT_Array + AC_Array + AG_Array
TSummed_Array = TA_Array + TT_Array + TC_Array + TG_Array
CSummed_Array = CA_Array + CT_Array + CC_Array + CG_Array
GSummed_Array = GA_Array + GT_Array + GC_Array + GG_Array


AllSummed_Array = ASummed_Array + TSummed_Array + CSummed_Array + GSummed_Array
print(AllSummed_Array)

Startindexation = -1
N = 151
Z = 100000

#print(ASummed_Array)
#print(ASummed_Array[0])
#print(TSummed_Array)
#print(CSummed_Array)
#print(GSummed_Array)

DNA_Arrays = [0]*Z

for z in range(Z):
    Random_Value = random.random()

    if Random_Value < ASummed_Array[Startindexation%10]: #It starts with an A
        if Random_Value < AA_Array[Startindexation%10]:
            DNAString = 'AA'
            FirstIndex = 0
        elif Random_Value < AA_Array[Startindexation%10] + AT_Array[Startindexation%10]:
            DNAString = 'AT'
            FirstIndex = 1
        elif Random_Value < AA_Array[Startindexation%10] + AT_Array[Startindexation%10] + AC_Array[Startindexation%10]:
            DNAString = 'AC'
            FirstIndex = 2
        else:
            DNAString = 'AG'
            FirstIndex = 3
    elif Random_Value < ASummed_Array[Startindexation%10] + TSummed_Array[Startindexation%10]: #It starts with a T
        if Random_Value < ASummed_Array[Startindexation%10] + TA_Array[Startindexation%10]:
            DNAString = 'TA'
            FirstIndex = 0
        elif Random_Value < ASummed_Array[Startindexation%10] + TA_Array[Startindexation%10] + TT_Array[Startindexation%10]:
            DNAString = 'TT'
            FirstIndex = 1
        elif Random_Value < ASummed_Array[Startindexation%10] + TA_Array[Startindexation%10] + TT_Array[Startindexation%10] + TC_Array[Startindexation%10]:
            DNAString = 'TC'
            FirstIndex = 2
        else:
            DNAString = 'TG'
            FirstIndex = 3
    elif Random_Value < ASummed_Array[Startindexation%10] + TSummed_Array[Startindexation%10] + CSummed_Array[Startindexation%10]: #It starts with a C
        if Random_Value < ASummed_Array[Startindexation%10] + TSummed_Array[Startindexation%10] +CA_Array[Startindexation%10]:
            DNAString = 'CA'
            FirstIndex = 0
        elif Random_Value < ASummed_Array[Startindexation%10] + TSummed_Array[Startindexation%10] +CA_Array[Startindexation%10] + CT_Array[Startindexation%10]:
            DNAString = 'CT'
            FirstIndex = 1
        elif Random_Value < ASummed_Array[Startindexation%10] + TSummed_Array[Startindexation%10] + CA_Array[Startindexation%10] + CT_Array[Startindexation%10] + CC_Array[Startindexation%10]:
            DNAString = 'CC'
            FirstIndex = 2
        else:
            DNAString = 'CG'
            FirstIndex = 3
    else:
        if Random_Value < ASummed_Array[Startindexation%10] + TSummed_Array[Startindexation%10] + CSummed_Array[Startindexation%10] + GA_Array[Startindexation%10]:
            DNAString = 'GA'
            FirstIndex = 0
        elif Random_Value < ASummed_Array[Startindexation%10] + TSummed_Array[Startindexation%10] + CSummed_Array[Startindexation%10] + GA_Array[Startindexation%10] + GT_Array[Startindexation%10]:
            DNAString = 'GT'
            FirstIndex = 1
        elif Random_Value < ASummed_Array[Startindexation%10] + TSummed_Array[Startindexation%10] + CSummed_Array[Startindexation%10] +  GA_Array[Startindexation%10] + GT_Array[Startindexation%10] + GC_Array[Startindexation%10]:
            DNAString = 'GC'
            FirstIndex = 2
        else:
            DNAString = 'GG'
            FirstIndex = 3

    #print(DNAString)

    for k in range(2+Startindexation,N+Startindexation):

        Random_Value = random.random()

        if FirstIndex==0:  # It starts with an A
            Random_Value = Random_Value*ASummed_Array[k%10]

            if Random_Value < AA_Array[k%10]:
                DNAString += 'A'
                FirstIndex = 0
            elif Random_Value < AA_Array[k%10] + AT_Array[k%10]:
                DNAString += 'T'
                FirstIndex = 1
            elif Random_Value < AA_Array[k%10] + AT_Array[k%10] + AC_Array[k%10]:
                DNAString += 'C'
                FirstIndex = 2
            else:
                DNAString += 'G'
                FirstIndex = 3
        elif FirstIndex==1:  # It starts with a T
            Random_Value = Random_Value * TSummed_Array[k % 10]

            if Random_Value < TA_Array[k % 10]:
                DNAString += 'A'
                FirstIndex = 0
            elif Random_Value < TA_Array[k % 10] + TT_Array[k % 10]:
                DNAString += 'T'
                FirstIndex = 1
            elif Random_Value < TA_Array[k % 10] + TT_Array[k % 10] + TC_Array[k % 10]:
                DNAString += 'C'
                FirstIndex = 2
            else:
                DNAString += 'G'
                FirstIndex = 3
        elif FirstIndex ==2:  # It starts with a C
            Random_Value = Random_Value * CSummed_Array[k % 10]

            if Random_Value < CA_Array[k % 10]:
                DNAString += 'A'
                FirstIndex = 0
            elif Random_Value < CA_Array[k % 10] + CT_Array[k % 10]:
                DNAString += 'T'
                FirstIndex = 1
            elif Random_Value < CA_Array[k % 10] + CT_Array[k % 10] + CC_Array[k % 10]:
                DNAString += 'C'
                FirstIndex = 2
            else:
                DNAString += 'G'
                FirstIndex = 3
        else:
            Random_Value = Random_Value * GSummed_Array[k % 10]

            if Random_Value < GA_Array[k % 10]:
                DNAString += 'A'
                FirstIndex = 0
            elif Random_Value < GA_Array[k % 10] + GT_Array[k % 10]:
                DNAString += 'T'
                FirstIndex = 1
            elif Random_Value < GA_Array[k % 10] + GT_Array[k % 10] + GC_Array[k % 10]:
                DNAString += 'C'
                FirstIndex = 2
            else:
                DNAString += 'G'
                FirstIndex = 3
    DNA_Arrays[z] = DNAString

#print(DNA_Arrays)

A_Check_Array = np.array([[0]*10]*4,dtype=float)
T_Check_Array = np.array([[0]*10]*4,dtype=float)
C_Check_Array = np.array([[0]*10]*4,dtype=float)
G_Check_Array = np.array([[0]*10]*4,dtype=float)


for i in range(Z):
    for j in range(N-1):
        String = DNA_Arrays[i][j:j+2]
        if String == 'AA':
            A_Check_Array[0][(j-1)%10] += 10/((N-1)*Z)
        elif String == 'AT':
            A_Check_Array[1][(j-1)%10] += 10 / ((N - 1) * Z)
        elif String == 'AC':
            A_Check_Array[2][(j-1)%10] += 10 / ((N - 1) * Z)
        elif String == 'AG':
            A_Check_Array[3][(j-1)%10] += 10 / ((N - 1) * Z)
        elif String == 'TA':
            T_Check_Array[0][(j-1)%10] += 10 / ((N - 1) * Z)
        elif String == 'TT':
            T_Check_Array[1][(j-1)%10] += 10 / ((N - 1) * Z)
        elif String == 'TC':
            T_Check_Array[2][(j-1)%10] += 10 / ((N - 1) * Z)
        elif String == 'TG':
            T_Check_Array[3][(j-1)%10] += 10 / ((N - 1) * Z)
        elif String == 'CA':
            C_Check_Array[0][(j-1)%10] += 10 / ((N - 1) * Z)
        elif String == 'CT':
            C_Check_Array[1][(j-1)%10] += 10 / ((N - 1) * Z)
        elif String == 'CC':
            C_Check_Array[2][(j-1)%10] += 10 / ((N - 1) * Z)
        elif String == 'CG':
            C_Check_Array[3][(j-1)%10] += 10 / ((N - 1) * Z)
        elif String == 'GA':
            G_Check_Array[0][(j-1)%10] += 10 / ((N - 1) * Z)
        elif String == 'GT':
            G_Check_Array[1][(j-1)%10] += 10 / ((N - 1) * Z)
        elif String == 'GC':
            G_Check_Array[2][(j-1)%10] += 10 / ((N - 1) * Z)
        elif String == 'GG':
            G_Check_Array[3][(j-1)%10] += 10 / ((N - 1) * Z)

print(A_Check_Array)

print("AA check =>   " + str(A_Check_Array[0]-AA_Array) + "\nWhich is normalised as:    " + str(np.linalg.norm(A_Check_Array[0]-AA_Array)))
print("AT check =>   " + str(A_Check_Array[1]-AT_Array) + "\nWhich is normalised as:    " + str(np.linalg.norm(A_Check_Array[1]-AT_Array)))
print("AC check =>   " + str(A_Check_Array[2]-AC_Array) + "\nWhich is normalised as:    " + str(np.linalg.norm(A_Check_Array[2]-AC_Array)))
print("AG check =>   " + str(A_Check_Array[3]-AG_Array) + "\nWhich is normalised as:    " + str(np.linalg.norm(A_Check_Array[3]-AG_Array)))
print("TA check =>   " + str(T_Check_Array[0]-TA_Array) + "\nWhich is normalised as:    " + str(np.linalg.norm(T_Check_Array[0]-TA_Array)))
print("TT check =>   " + str(T_Check_Array[1]-TT_Array) + "\nWhich is normalised as:    " + str(np.linalg.norm(T_Check_Array[1]-TT_Array)))
print("TC check =>   " + str(T_Check_Array[2]-TC_Array) + "\nWhich is normalised as:    " + str(np.linalg.norm(T_Check_Array[2]-TC_Array)))
print("TG check =>   " + str(T_Check_Array[3]-TG_Array) + "\nWhich is normalised as:    " + str(np.linalg.norm(T_Check_Array[3]-TG_Array)))
print("CA check =>   " + str(C_Check_Array[0]-CA_Array) + "\nWhich is normalised as:    " + str(np.linalg.norm(C_Check_Array[0]-CA_Array)))
print("CT check =>   " + str(C_Check_Array[1]-CT_Array) + "\nWhich is normalised as:    " + str(np.linalg.norm(C_Check_Array[1]-CT_Array)))
print("CC check =>   " + str(C_Check_Array[2]-CC_Array) + "\nWhich is normalised as:    " + str(np.linalg.norm(C_Check_Array[2]-CC_Array)))
print("CG check =>   " + str(C_Check_Array[3]-CG_Array) + "\nWhich is normalised as:    " + str(np.linalg.norm(C_Check_Array[3]-CG_Array)))
print("GA check =>   " + str(G_Check_Array[0]-GA_Array) + "\nWhich is normalised as:    " + str(np.linalg.norm(G_Check_Array[0]-GA_Array)))
print("GT check =>   " + str(G_Check_Array[1]-GT_Array) + "\nWhich is normalised as:    " + str(np.linalg.norm(G_Check_Array[1]-GT_Array)))
print("GC check =>   " + str(G_Check_Array[2]-GC_Array) + "\nWhich is normalised as:    " + str(np.linalg.norm(G_Check_Array[2]-GC_Array)))
print("GG check =>   " + str(G_Check_Array[3]-GG_Array) + "\nWhich is normalised as:    " + str(np.linalg.norm(G_Check_Array[3]-GG_Array)))
#print(A_Check_Array)
#print(T_Check_Array)
#print(C_Check_Array)
#print(G_Check_Array)


#print(DNA_Arrays)


MatlabString = "["
for k in range(len(DNA_Arrays)):
   MatlabString += '\'' +  DNA_Arrays[k] + '\''
   if k < len(DNA_Arrays)-1:
      MatlabString += ';'
   else:
      MatlabString += ']'
#print(MatlabString)