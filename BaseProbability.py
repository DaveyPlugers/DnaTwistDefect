import numpy as np
import math
import matplotlib.pyplot as plt
# {A,T,C,G} = [0,1,2,3]

First_Enzyme = 0
Second_Enzyme = 0

phi=0
N=147
Gamma = 4.46
Theta_Twist = 35.575


Theta_Tilt_Array = [[-1.4,0.0,-0.1,-1.7],[0.0,1.4,1.5,-0.5],[0.5,1.7,0.1,0.0],[-1.5,0.1,0.0,-0.1]]
K_Tilt_Array = [[0.100,0.166,0.111,0.149],[0.148,0.100,0.087,0.082],[0.082,0.149,0.119,0.068],[0.087,0.111,0.082,0.119]]
Theta_Roll_Array = [[0.7,1.1,0.7,4.5],[3.3,0.7,1.9,4.7],[4.7,4.5,3.6,5.4],[1.9,0.7,0.3,3.6]]
K_Roll_Array = [[0.049,0.055,0.080,0.096],[0.029,0.049,0.046,0.048],[0.048,0.096,0.064,0.050],[0.046,0.080,0.082,0.064]]


#Geometric_Properties[a][b][c]   %c is first index, b is second (cb = 00 =AA, 23 = CG,
#a is [0,1,2,3] = [theta_tilt,K_tilt,Theta_Roll,K_Roll]

Geometric_Properties = [Theta_Tilt_Array,K_Tilt_Array,Theta_Roll_Array,K_Roll_Array]
Print1 = False
if Print1: #Is test, verwijder
    print(Theta_Roll_Array)
    print(Theta_Roll_Array[0])
    print(Theta_Roll_Array[0][1])
    print(Geometric_Properties)
    print(Geometric_Properties[2][0])
    print(Geometric_Properties[2][0][1])


#Theta_Tilt = Gamma*sin(2*pi*n/10 - phi)
#Theta_roll = Gamma*cos(2*pi*n/10 - phi)
#Energy^n(N^n N^(n+1)) = E_tilt^n(N^n N^(n+1)) + E_roll^n(N^n N^(n+1))
#E_tilt^n(N^n N^(n+1)) = 0.5* K_tilt(N^nN^(n+1) [Theta_tilt^n - %Theta_tilt(N^nN^(n+1))]^2
#T^n_GC = exp^(-Beta E^n(GC))


#Numpy version
Theta_Tilt_Array = np.array([[-1.4,0.0,-0.1,-1.7],[0.0,1.4,1.5,-0.5],[0.5,1.7,0.1,0.0],[-1.5,0.1,0.0,-0.1]])
K_Tilt_Array = np.array([[0.100,0.166,0.111,0.149],[0.148,0.100,0.087,0.082],[0.082,0.149,0.119,0.068],[0.087,0.111,0.082,0.119]])
Theta_Roll_Array = np.array([[0.7,1.1,0.7,4.5],[3.3,0.7,1.9,4.7],[4.7,4.5,3.6,5.4],[1.9,0.7,0.3,3.6]])
K_Roll_Array = np.array([[0.049,0.055,0.080,0.096],[0.029,0.049,0.046,0.048],[0.048,0.096,0.064,0.050],[0.046,0.080,0.082,0.064]])

Numpy_Geometric_Properties = np.array([Theta_Tilt_Array,K_Tilt_Array,Theta_Roll_Array,K_Roll_Array])

#print(Numpy_Geometric_Properties)

Transfer_Matrices = np.zeros(N).astype(list)



def Energy_Tilt(Angle,First_Index,Second_Index):
    E_Tilt = 0.5*Numpy_Geometric_Properties[1,Second_Index,First_Index]*(
                Gamma*math.sin(Angle - phi) - Numpy_Geometric_Properties[0,Second_Index,First_Index])**2
    return E_Tilt

def Energy_Roll(Angle,First_Index,Second_Index):
    E_Roll = 0.5*Numpy_Geometric_Properties[3,Second_Index,First_Index]*(
                Gamma*math.cos(Angle - phi) - Numpy_Geometric_Properties[2,Second_Index,First_Index])**2
    return E_Roll



for i in range(N):
    Angle = 2*math.pi*i/10
    Transfer_Values_At_i = np.array([[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]]).astype(float)
    for j in range(4):
        for k in range(4):
            Transfer_Values_At_i[j,k] = math.e**(-(Energy_Tilt(Angle,k,j)+Energy_Roll(Angle,k,j)))
    Transfer_Matrices[i] = Transfer_Values_At_i

Z_function = Transfer_Matrices[0]
for i in range(N-2):
    Z_function = np.matmul(Z_function,Transfer_Matrices[i+1])

print(Z_function)
Z = np.sum(np.sum(Z_function))
print(Z)

Probabilities = np.zeros(N-2)
for First_Enzyme in range(4):
    for Second_Enzyme in range(4):
        for i in range(N-2):
            Lower_Matrix = np.zeros(4)
            Upper_Matrix = np.zeros(4)

            if i==0: #I skip these edges for now
                1+1
            elif i==N-2: #I skip these edges for now
                1+2
            else:
                Lower_Matrix = Transfer_Matrices[0]
                for k in range(1,i):
                    Lower_Matrix = np.matmul(Lower_Matrix,Transfer_Matrices[k])
                Higher_Matrix = Transfer_Matrices[i+1]
                for k in range(i+2,N-1):
                    Higher_Matrix = np.matmul(Higher_Matrix,Transfer_Matrices[k])

                SumValues = 0
                for m in range(4):
                    for n in range(4):
                        SumValues += Lower_Matrix[m,First_Enzyme]*Transfer_Matrices[i][First_Enzyme,Second_Enzyme]*Higher_Matrix[Second_Enzyme,n]
                Probabilities[i] +=SumValues/Z

print(Probabilities[1:N-1])
plt.plot(Probabilities[1:N-1])
plt.ylim([0,1.1])
plt.show()