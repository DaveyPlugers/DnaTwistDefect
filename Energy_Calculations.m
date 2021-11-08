
DNAString = 'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA'
DNA = char(DNAString)
DNAIndex = []

for i=1:N+1
    if DNA(i) == 'A'    
        DNAIndex = [DNAIndex,1];
    elseif DNA(i) == 'T'
        DNAIndex = [DNAIndex,2];
    elseif DNA(i) == 'C'
        DNAIndex = [DNAIndex,3];
    else
        DNAIndex = [DNAIndex,4];
    end
end

% {A,T,C,G} = {1,2,3,4]

%Different functions (Use solid geometric shape, only include some energies
%etc.)
Probability_Mode = true;

Use_Tilt_Energy = true;
Use_Roll_Energy = true;
Use_Twist_Energy = true;
Use_Raise_Energy = true;

Plot_Energy = true;
Plot_Individual_Energy = true;
Plot_DNA = true;

Use_Static_Geometrics = true;
Use_Variable_Raise_Twist = false;

Twist_Defect_Type1 = true; %5/5 split into 4/1/4
Twist_Defect_Type2 = false;  %10/10 split into 9/1/9
Defect_Location = 41;
%Variables
N=147;

%These are for probability mode
First_Enzyme = 3;
Second_Enzyme = 4;

Gamma = 4.46;
Theta_Twist = 35.575;
Raise = 3.4; %Angstrom
Phase_Shift = 0;
%Initialisation


Tilt_Array = [[-1.4,0.0,-0.1,-1.7];[0.0,1.4,1.5,-0.5];[0.5,1.7,0.1,0.0];[-1.5,0.1,0.0,-0.1]];
K_Tilt_Array = [[0.100,0.166,0.111,0.149];[0.148,0.100,0.087,0.082];[0.082,0.149,0.119,0.068];[0.087,0.111,0.082,0.119]];
Roll_Array = [[0.7,1.1,0.7,4.5];[3.3,0.7,1.9,4.7];[4.7,4.5,3.6,5.4];[1.9,0.7,0.3,3.6]];
K_Roll_Array = [[0.049,0.055,0.080,0.096];[0.029,0.049,0.046,0.048];[0.048,0.096,0.064,0.050];[0.046,0.080,0.082,0.064]];
Twist_Array = [[35.1,29.3,31.5,31.9];[37.8,35.1,36.3,37.3];[37.3,31.9,32.9,36.1];[36.3,31.5,33.6,32.9]];
K_Twist_Array = [[0.092,0.070,0.073,0.064];[0.052,0.092,0.071,0.043];[0.043,0.064,0.041,0.047];[0.071,0.073,0.055,0.041]];
Raise_Array = [[3.27,3.31,3.36,3.34];[3.42,3.27,3.37,3.33];[3.33,3.34,3.42,3.39];[3.37,3.36,3.40,3.42]];
K_Raise_Array = [[21.748,25.547,23.860,29.496];[21.914,21.748,22.820,18.235];[18.235,29.496,30.312,14.164];[22.820,23.860,25.860,30.312]];


Geometric_Properties = Tilt_Array;
Geometric_Properties(:,:,2) = K_Tilt_Array;
Geometric_Properties(:,:,3) = Roll_Array;
Geometric_Properties(:,:,4) = K_Roll_Array;
Geometric_Properties(:,:,5) = Twist_Array;
Geometric_Properties(:,:,6) = K_Twist_Array;
Geometric_Properties(:,:,7) = Raise_Array;
Geometric_Properties(:,:,8) = K_Raise_Array;




if Use_Static_Geometrics==1
    %[slide,shift,raise, twist,tilt,roll]
    DNA_Geometry = zeros(N,6);
    for i=1:N
        DNA_Geometry(i,:) = [0,0,Raise,Theta_Twist,Gamma*cos(2*pi*i/10-Phase_Shift),Gamma*sin(2*pi*i/10-Phase_Shift)];
    end
elseif Use_Variable_Raise_Twist==1
    %[slide,shift,raise, twist,tilt,roll]
    DNA_Geometry = zeros(N,6);
    for i=1:N
        DNA_Geometry(i,:) = [0,0,Geometric_Properties(DNAIndex(i),DNAIndex(i+1),7),Geometric_Properties(DNAIndex(i),DNAIndex(i+1),5),Gamma*cos(2*pi*i/10-Phase_Shift),Gamma*sin(2*pi*i/10-Phase_Shift)];
    end
end


if Plot_DNA
    figure()
    PlotterFunc(DNA_Geometry,2)
end


if Twist_Defect_Type1;
    for i=Defect_Location-4:Defect_Location+4
        
        DNA_Geometry(i,:) = [0,0,Raise*10/9,Theta_Twist*10/9,Gamma*cos(2*pi*((i-Defect_Location+5)*10/9+Defect_Location-5)/10 -Phase_Shift),Gamma*sin(2*pi*((i-Defect_Location+5)*10/9+Defect_Location-5)/10-Phase_Shift)];
    end
    %for i=Defect_Location+1:Defect_Location+10
       
        %DNA_Geometry(i,:) = [0,0,Raise*10/9,Theta_Twist*10/9,Gamma*cos(2*pi*i/10-Phase_Shift),Gamma*sin(2*pi*i/10-Phase_Shift)];
    %end
    for i=Defect_Location+5:N-1
        DNA_Geometry(i,:) = DNA_Geometry(i+1,:);
    end
    DNA_Geometry(N,:) = [0,0,Raise,Theta_Twist,Gamma*cos(2*pi*N/10-Phase_Shift),Gamma*sin(2*pi*N/10-Phase_Shift)];
end

if Twist_Defect_Type2;
    for i=Defect_Location-9:Defect_Location+9
        
        DNA_Geometry(i,:) = [0,0,Raise*20/19,Theta_Twist*20/19,Gamma*cos(2*pi*((i-Defect_Location+10)*20/19+Defect_Location-10)/10 -Phase_Shift),Gamma*sin(2*pi*((i-Defect_Location+10)*20/19+Defect_Location-10)/10-Phase_Shift)];
    end
    %for i=Defect_Location+1:Defect_Location+10
       
        %DNA_Geometry(i,:) = [0,0,Raise*10/9,Theta_Twist*10/9,Gamma*cos(2*pi*i/10-Phase_Shift),Gamma*sin(2*pi*i/10-Phase_Shift)];
    %end
    for i=Defect_Location+10:N-1
        DNA_Geometry(i,:) = DNA_Geometry(i+1,:);
    end
    DNA_Geometry(N,:) = [0,0,Raise,Theta_Twist,Gamma*cos(2*pi*N/10-Phase_Shift),Gamma*sin(2*pi*N/10-Phase_Shift)];
end

if Plot_DNA
    PlotterFunc(DNA_Geometry,1)
end
figure()
%Geometric_Properties(a,b,c)   %a is first index, b is second (ab = 11 =AA, 34 = CG, 
%c is [1,2,3,4,5,6,7,8] = [Angle_tilt,K_tilt,Angle_Roll,K_Roll,Angle_Twist,K_Twist,Angstrom_Raise,K_Raise]

%Theta_Tilt = Gamma*sin(2*pi*n/10 - phi)
%Theta_roll = Gamma*cos(2*pi*n/10 - phi)

%Energy^n(N^n N^(n+1)) = E_tilt^n(N^n N^(n+1)) + E_roll^n(N^n N^(n+1))

%E_tilt^n(N^n N^(n+1)) = 0.5* K_tilt(N^nN^(n+1) [Theta_tilt^n - %Theta_tilt(N^nN^(n+1))]^2

%T^n_GC = exp^(-Beta E^n(GC))



%Energy Calculation
Energy = zeros(N-1,6);

if Use_Tilt_Energy
    for i=1:N-1
        Energy(i,5) = 0.5*Geometric_Properties(DNAIndex(i),DNAIndex(i+1),2)*(DNA_Geometry(i,5)-Geometric_Properties(DNAIndex(i),DNAIndex(i+1),1))^2;
    end       
end

if Use_Roll_Energy
    for i=1:N-1
        Energy(i,6) = 0.5*Geometric_Properties(DNAIndex(i),DNAIndex(i+1),4)*(DNA_Geometry(i,6)-Geometric_Properties(DNAIndex(i),DNAIndex(i+1),3))^2;
    end       
end

if Use_Twist_Energy
    for i=1:N-1
        Energy(i,4) = 0.5*Geometric_Properties(DNAIndex(i),DNAIndex(i+1),6)*(DNA_Geometry(i,4)-Geometric_Properties(DNAIndex(i),DNAIndex(i+1),5))^2;
    end       
end

if Use_Raise_Energy
    for i=1:N-1
        Energy(i,3) = 0.5*Geometric_Properties(DNAIndex(i),DNAIndex(i+1),8)*(DNA_Geometry(i,3)-Geometric_Properties(DNAIndex(i),DNAIndex(i+1),7))^2;
    end       
end

if Plot_Energy;
    plot(sum(Energy,2))
    if Plot_Individual_Energy;
        hold on
        plot(Energy(:,3))
        plot(Energy(:,4))
        plot(Energy(:,6)) 
        plot(Energy(:,5))
        Energy_String_Legend = strcat('Total Energy (k_BT) = ' , num2str(sum(sum(Energy))))
        legend('Total Energy','Raise Energy','Twist Energy','Roll Energy','Tilt Energy')
        title(Energy_String_Legend)
    end
end



%Making the transfer matrices
if Probability_Mode
    Transfer_Matrices = zeros(4);
    %Making Transfer Matrices
    for i=1:N-1
        Transfer_Matrices(:,:,i) = [0,0,0,0;0,0,0,0;0,0,0,0;0,0,0,0];
        for j=1:4 %Loop over first index
            for k=1:4 %Second index
                Energy_Individual = 0;
                if Use_Tilt_Energy
                    E_Tilt = 0.5*Geometric_Properties(j,k,2)*(DNA_Geometry(i,5)-Geometric_Properties(j,k,1))^2;
                    Energy_Individual = Energy_Individual + E_Tilt;
                end
                if Use_Roll_Energy
                    E_Roll = 0.5*Geometric_Properties(j,k,4)*(DNA_Geometry(i,6)-Geometric_Properties(j,k,3))^2;;
                    Energy_Individual = Energy_Individual + E_Roll;
                end
                if Use_Twist_Energy
                    E_Twist = 0.5*Geometric_Properties(j,k,6)*(DNA_Geometry(i,4)-Geometric_Properties(j,k,5))^2;;
                    Energy_Individual = Energy_Individual + E_Twist;
                end
                if Use_Raise_Energy
                    E_Raise = 0.5*Geometric_Properties(j,k,8)*(DNA_Geometry(i,3)-Geometric_Properties(j,k,7))^2;
                    Energy_Individual = Energy_Individual + E_Raise;
                end
                Transfer_Matrices(j,k,i) = exp(-(Energy_Individual));
            end
        end   
    end

    %Normalisation
    Z_function = Transfer_Matrices(:,:,1);
    for i=1:N-2
        Transfer_Matrices(:,:,i+1);
        Z_function = Z_function*Transfer_Matrices(:,:,i+1);

    end
    
    Z = sum(sum(Z_function));


    Sum_Z = 0;
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

    %Calculating probabilities
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
            Transfer_Matrices(:,:,N-1);
            SumValues = 0;
            for m=1:4
                for n=1:4
                    if i==3
                        1+1;
                    end
                    Lower_Matrix(m,First_Enzyme);
                    Transfer_Matrices(First_Enzyme,Second_Enzyme,i);
                    Higher_Matrix(Second_Enzyme,n);
                    SumValues = SumValues + Lower_Matrix(m,First_Enzyme)*Transfer_Matrices(First_Enzyme,Second_Enzyme,i)*Higher_Matrix(Second_Enzyme,n);
                end
            end
            Probabilities(i) = SumValues/Z;
        end     
    end
       figure()
    plot(Probabilities(2:N-2))
end



