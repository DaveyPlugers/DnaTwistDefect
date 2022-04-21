function [Energy_Difference,Defected_Geom] = Introduction_Energy_Cost(SHL,Direction,DNA_Geometry,DNAIndexation,Geometric_Properties)
%INTRODUCTION_ENERGY_COST Summary of this function goes here
%   Detailed explanation goes here

%We will calculate the energy cost to introduce an overtwist and an
%undertwist on both sides of the remodeller, once we have this we will use
%this energy value to calculate the introduction rate

%We give it the DNA_geometry, we return the energy cost
Gamma = 4.46;
Theta_Twist = 35.575;
Raise = 3.4; %Angstrom
Phase_Shift = -6*pi/20;
N=147;
AminoBP = [6,16,26,36,46,56,66,76,86,96,106,116,126,136,146];
%SHL where our overtwist is introduced, the one we actually use for the dynamics, 
%Position = SHL - 1 is then the location of the undertwist

%6 to 16 is -7
%56 to 66 is -2
%For SHL = -2 we want 46 to 55 undertwist and 55 to 66 overtwist




BeginIndex = 10*(SHL+7)+6; 
Defected_Geom = DNA_Geometry;
if Direction==0
    k=0; %This is used as a correction of the indexation to fit the geometry (see plotterfunc for visualisation)
    for i=BeginIndex-9:BeginIndex-1
       Defected_Geom(i,:) = [0,0,Raise*10/9,Theta_Twist*10/9,5/4.5*Gamma*cos(2*pi*(i-BeginIndex+9)*10/(9*10) + 2*pi*(BeginIndex-9)/10 -Phase_Shift),5/4.5*Gamma*sin(2*pi*(i-BeginIndex+9)*10/(9*10)+ 2*pi*(BeginIndex-9)/10-Phase_Shift)];
    end

    for i=BeginIndex+k:BeginIndex+10+k
       Defected_Geom(i,:) = [0,0,Raise*10/11,Theta_Twist*10/11,5/5.5*Gamma*cos(2*pi*(i-BeginIndex)*10/(11*10) + 2*pi*(BeginIndex)/10 -Phase_Shift),5/5.5*Gamma*sin(2*pi*(i-BeginIndex)*10/(11*10)+ 2*pi*(BeginIndex)/10-Phase_Shift)];
    end

    AminoBP(SHL+7+1) = AminoBP(SHL+7+1)-1;
elseif Direction==1
    k=1;
    for i=BeginIndex-9+k:BeginIndex+1+k
       Defected_Geom(i-1+k,:) = [0,0,Raise*10/11,Theta_Twist*10/11,5/5.5*Gamma*cos(2*pi*(i-BeginIndex-1)*10/(11*10) + 2*pi*(BeginIndex-1)/10 -Phase_Shift),5/5.5*Gamma*sin(2*pi*(i-BeginIndex-1)*10/(11*10)+ 2*pi*(BeginIndex-1)/10-Phase_Shift)];
    end

    for i=BeginIndex+2+k:BeginIndex+10+k
       Defected_Geom(i-1+k,:) = [0,0,Raise*10/9,Theta_Twist*10/9,5/4.5*Gamma*cos(2*pi*(i-BeginIndex-1)*10/(9*10) + 2*pi*(BeginIndex-1)/10 -Phase_Shift),5/4.5*Gamma*sin(2*pi*(i-BeginIndex-1)*10/(9*10)+ 2*pi*(BeginIndex-1)/10-Phase_Shift)];
    end

    AminoBP(SHL+7+1) = AminoBP(SHL+7+1)-1;
end


Undefected_Energy = sum(sum(EnergyValuesCalculator(DNA_Geometry,DNAIndexation,Geometric_Properties,true,true,true,true)));
Defected_Energy = sum(sum(EnergyValuesCalculator(Defected_Geom,DNAIndexation,Geometric_Properties,true,true,true,true)));

Energy_Difference = Defected_Energy - Undefected_Energy;




end

