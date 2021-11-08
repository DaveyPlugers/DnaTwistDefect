%
% {A,T,C,G} = {1,2,3,4]

First_Enzyme = 3
Second_Enzyme = 4

phi=0
N=147;
Gamma = 4.46;
Theta_Twist = 35.575;


Theta_Tilt_Array = [-1.4,0.0,-0.1,-1.7;0.0,1.4,1.5,-0.5;0.5,1.7,0.1,0.0;-1.5,0.1,0.0,-0.1];
K_Tilt_Array = [0.100,0.166,0.111,0.149;0.148,0.100,0.087,0.082;0.082,0.149,0.119,0.068;0.087,0.111,0.082,0.119];
Theta_Roll_Array = [0.7,1.1,0.7,4.5;3.3,0.7,1.9,4.7;4.7,4.5,3.6,5.4;1.9,0.7,0.3,3.6];
K_Roll_Array = [0.049,0.055,0.080,0.096;0.029,0.049,0.046,0.048;0.048,0.096,0.064,0.050;0.046,0.080,0.082,0.064];

Geometric_Properties = Theta_Tilt_Array;
Geometric_Properties(:,:,2) = K_Tilt_Array;
Geometric_Properties(:,:,3) = Theta_Roll_Array;
Geometric_Properties(:,:,4) = K_Roll_Array;

Geometric_Properties(2,3,2)
%Geometric_Properties(a,b,c)   %a is first index, b is second (ab = 11 =AA, 34 = CG, 
%c is [1,2,3,4] = [theta_tilt,K_tilt,Theta_Roll,K_Roll]



%Theta_Tilt = Gamma*sin(2*pi*n/10 - phi)
%Theta_roll = Gamma*cos(2*pi*n/10 - phi)

%Energy^n(N^n N^(n+1)) = E_tilt^n(N^n N^(n+1)) + E_roll^n(N^n N^(n+1))

%E_tilt^n(N^n N^(n+1)) = 0.5* K_tilt(N^nN^(n+1) [Theta_tilt^n - %Theta_tilt(N^nN^(n+1))]^2

%T^n_GC = exp^(-Beta E^n(GC))

%Making the transfer matrices
Transfer_Matrices = zeros(4);
for i=0:N-2
    Transfer_Matrices(:,:,i+1) = [0,0,0,0;0,0,0,0;0,0,0,0;0,0,0,0];
    for j=1:4 %Loop over first index
        for k=1:4 %Second index
            E_Tilt = 0.5*Geometric_Properties(j,k,2)*(Gamma*sin(2*pi*i/10-phi)-Geometric_Properties(j,k,1))^2;
            E_Roll = 0.5*Geometric_Properties(j,k,4)*(Gamma*cos(2*pi*i/10-phi)-Geometric_Properties(j,k,3))^2;
            Transfer_Matrices(j,k,i+1) = exp(-(E_Tilt + E_Roll));
        end
    end   
end

Z_function = Transfer_Matrices(:,:,1);
for i=1:N-2
    Z_function
    Transfer_Matrices(:,:,i+1)
    Z_function = Z_function*Transfer_Matrices(:,:,i+1);
    
end

Z = sum(sum(Z_function))


Sum_Z = 0
for i=1:4
   for j=1:4
       Z_Twee = Transfer_Matrices(:,:,2);
       Low = Transfer_Matrices(i,:,1);
       High = Transfer_Matrices(:,j,N-1);
       for k=3:N-2
           Z_Twee = Z_Twee*Transfer_Matrices(:,:,k);
       end
       Sum_Z = Sum_Z + Low * Z_Twee * High;
   end
end


Probabilities = zeros(1,N-1);


for i=1:N-1
    Lower_Matrix = zeros(4);
    Higher_Matrix = zeros(4);

    if i==1
        1+1;
    elseif i ==N-1
        1+2;
    else 
        
        Lower_Matrix = Transfer_Matrices(:,:,1);
        for k = 2:i-1
            Lower_Matrix = Lower_Matrix*Transfer_Matrices(:,:,k);
        end
        Higher_Matrix = Transfer_Matrices(:,:,i+1);
        for k =i+2:N-1

            Higher_Matrix = Higher_Matrix*Transfer_Matrices(:,:,k);
        end
        Transfer_Matrices(:,:,N-1)
        SumValues = 0;
        for m=1:4
            for n=1:4
                if i==3
                    1+1;
                end
                Lower_Matrix(m,First_Enzyme)
                Transfer_Matrices(First_Enzyme,Second_Enzyme,i)
                Higher_Matrix(Second_Enzyme,n)
                SumValues = SumValues + Lower_Matrix(m,First_Enzyme)*Transfer_Matrices(First_Enzyme,Second_Enzyme,i)*Higher_Matrix(Second_Enzyme,n);
            end
        end
        Probabilities(i) = SumValues/Z;
    end     
end

plot(Probabilities(2:N-2))

