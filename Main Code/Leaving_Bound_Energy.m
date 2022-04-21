function [Energy] = Leaving_Bound_Energy(Direction,DNAString,Geometric_Properties)
%LEAVING_BOUND_ENERGY Summary of this function goes here
%   Detailed explanation goes here

%Going to be using this to update the energy leaving bound later

%DNAString = 'ACGCCGATCGCGCTACGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG';
[DNA_Geometry,DNAIndexation] = GeomArrayMaker(DNAString,1,Geometric_Properties);



Gamma = 4.46;
Theta_Twist = 35.575;
Raise = 3.4; %Angstrom
Phase_Shift = -6*pi/20;
N=147;

%Direction=1;




if Direction==0
    Defect_Location = 14;
    AminoBP = [6,16,26,36,46,56,66,76,86,96,106,116,126,136,146];
    BeginIndexType6 = 10*(Defect_Location-1)+6;
    Defected_Geom = DNA_Geometry;
    TempGeomType6 = DNA_Geometry(BeginIndexType6+10,:);
    for i=BeginIndexType6:BeginIndexType6+10
        Defected_Geom(i,:) = [0,0,Raise*20/21,Theta_Twist*20/21,20/21*Gamma*cos(2*pi*(i-BeginIndexType6)*20/(21*10) + 2*pi*BeginIndexType6/10 -Phase_Shift),20/21*Gamma*sin(2*pi*(i-BeginIndexType6)*20/(21*10)+ 2*pi*BeginIndexType6/10-Phase_Shift)];
    end

    for i=0:N-BeginIndexType6-11
        Defected_Geom(N-i,:) = Defected_Geom(N-i-1,:);
    end
    
    Defected_Geom(BeginIndexType6+11,:) = TempGeomType6;
    AminoBP(Defect_Location+1:15)=AminoBP(Defect_Location+1:15)+1;
    AminoBP = [AminoBP(1:Defect_Location), 10*Defect_Location+6 , AminoBP(Defect_Location+1:15)];
    Defect_Temp = EnergyValuesCalculator(Defected_Geom,DNAIndexation,Geometric_Properties,true,true,true,true);
    Normal_Temp = EnergyValuesCalculator(DNA_Geometry,DNAIndexation,Geometric_Properties,true,true,true,true);
    Energy = sum(sum(Defect_Temp(136:147,:) - Normal_Temp(136:147,:)));
else
    AminoBP = [7,17,27,37,47,57,67,77,87,97,107,117,127,137,147];
    Defect_Location = 0;
    BeginIndexType6 = 7;
    Defected_Geom = DNA_Geometry;
    TempGeomType6 = DNA_Geometry(BeginIndexType6-1,:);
    for i=BeginIndexType6+1:BeginIndexType6+11
        Defected_Geom(i-2,:) = [0,0,Raise*20/21,Theta_Twist*20/21,20/21*Gamma*cos(2*pi*(i-BeginIndexType6-1)*20/(21*10) + 2*pi*(BeginIndexType6-1)/10 -Phase_Shift),20/21*Gamma*sin(2*pi*(i-BeginIndexType6-1)*20/(21*10)+ 2*pi*(BeginIndexType6-1)/10-Phase_Shift)];
    end
    
    for i=1:BeginIndexType6-3
        Defected_Geom(i,:) = Defected_Geom(i+1,:);
    end
    Defected_Geom(BeginIndexType6-2,:) = TempGeomType6;
    AminoBP(1:Defect_Location)=AminoBP(1:Defect_Location)-1;
    AminoBP = [AminoBP(1:Defect_Location), 10*Defect_Location+6 , AminoBP(Defect_Location+1:15)];
    Defect_Temp = EnergyValuesCalculator(Defected_Geom,DNAIndexation,Geometric_Properties,true,true,true,true);
    Normal_Temp = EnergyValuesCalculator(DNA_Geometry,DNAIndexation,Geometric_Properties,true,true,true,true);
    Energy = sum(sum(Defect_Temp(6:17,:) - Normal_Temp(6:17,:)));
end



end

