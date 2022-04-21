function [GeometricArray,DNAIndexation] = GeomArrayMaker(DNAString,Mode,Geometric_Properties)
%GEOMARRAYMAKER 
%   Input: 
%   DNAString = 148 character string of DNA
%   Mode = What type of geometry (Mode=1 static raise and twist, 2 variable raise and twists)   
%   Geometric_Properties = The geometric properties array so it doesn't
%   have to be remade every run
%   Output: 
%   GeometricArray = Corresponding 147 length array with the geometric properties
%   (raise,twist,tilt,roll)
%   DNAIndexation = Corresponding 147 character string with numbers

%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%Initial values, remove later? Matlab is being stupid and doesn't recognize
%these unless I parse them through or redefine, check computing cost later
Gamma = 4.46;
Theta_Twist = 35.575;
Raise = 3.4; %Angstrom
Phase_Shift = -6*pi/20;


%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

N=147;
DNAIndexation = []; %Speed up possible but negligeble
for i=1:N+1
        if DNAString(i) == 'A'    
            DNAIndexation = [DNAIndexation,1];
        elseif DNAString(i) == 'T'
            DNAIndexation = [DNAIndexation,2];
        elseif DNAString(i) == 'C'
            DNAIndexation = [DNAIndexation,3];
        else
            DNAIndexation = [DNAIndexation,4];
        end
end


if Mode==1
    %[slide,shift,raise, twist,tilt,roll]
    GeometricArray = zeros(N,6);
    for i=1:N
        GeometricArray(i,:) = [0,0,Raise,Theta_Twist,Gamma*cos(2*pi*i/10-Phase_Shift),Gamma*sin(2*pi*i/10-Phase_Shift)];
    end
elseif Mode==2
    %[slide,shift,raise, twist,tilt,roll]
    GeometricArray = zeros(N,6);
    for i=1:N
        GeometricArray(i,:) = [0,0,Geometric_Properties(DNAIndexation(i),DNAIndexation(i+1),7),Geometric_Properties(DNAIndexation(i),DNAIndexation(i+1),5),Gamma*cos(2*pi*i/10-Phase_Shift),Gamma*sin(2*pi*i/10-Phase_Shift)];
    end
end



end

