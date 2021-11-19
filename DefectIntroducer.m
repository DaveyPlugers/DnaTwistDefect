function [Defected_Geom,AminoBP] = DefectIntroducer(DNA_Geometry,Defect_Location,Defect_Type)
%DEFECTINTRODUCER 
%   Input:
%   DNA_Geometry = The geometry of the DNA that we will change
%   Defect_Location = Between which BP's will we add this defect (1 => 6-16
%   for K10, 1 => 6-26 for K20 then +10*n for location 1 + n)
%   Defect_type = What type of defect are we adding (1 => K10 undertwist,
%   2=> K20 undertwist) Overtwist not added in this function yet.
%   Output:
%   Defected_Geom = The geometry with the defect integrated into it
%   AminoBP = Array saving all the positions of binding sites

Gamma = 4.46;
Theta_Twist = 35.575;
Raise = 3.4; %Angstrom
Phase_Shift = 0;
N=147;
AminoBP = [6,16,26,36,46,56,66,76,86,96,106,116,126,136,146];


if Defect_Type == 1 %Undertwist 9 to 10
    BeginIndexType4 = 10*(Defect_Location-1)+6;
    Defected_Geom = DNA_Geometry;
    TempGeomType4 = DNA_Geometry(BeginIndexType4+10,:);
    for i=BeginIndexType4+1:BeginIndexType4+10
        Defected_Geom(i,:) = [0,0,Raise*9/10,Theta_Twist*9/10,5/5.5*Gamma*cos(2*pi*(i-BeginIndexType4)*9/(10*10) + 2*pi*BeginIndexType4/10 -Phase_Shift),5/5.5*Gamma*sin(2*pi*(i-BeginIndexType4)*9/(10*10)+ 2*pi*BeginIndexType4/10-Phase_Shift)];
    end
    
    for i=0:N-BeginIndexType4-12
        Defected_Geom(N-i,:) = Defected_Geom(N-i-1,:);
    end
    Defected_Geom(BeginIndexType4+11,:) = TempGeomType4;
    AminoBP(Defect_Location+1:15)=AminoBP(Defect_Location+1:15)+1;

elseif Defect_Type == 2
    BeginIndexType6 = 10*(Defect_Location-1)+6;
    Defected_Geom = DNA_Geometry;
    TempGeomType6 = DNA_Geometry(BeginIndexType6+20,:);
    for i=BeginIndexType6+1:BeginIndexType6+20
        Defected_Geom(i,:) = [0,0,Raise*19/20,Theta_Twist*19/20,20/21*Gamma*cos(2*pi*(i-BeginIndexType6)*19/(20*10) + 2*pi*BeginIndexType6/10 -Phase_Shift),20/21*Gamma*sin(2*pi*(i-BeginIndexType6)*19/(20*10)+ 2*pi*BeginIndexType6/10-Phase_Shift)];
    end
    
    for i=0:N-BeginIndexType6-22
        Defected_Geom(N-i,:) = Defected_Geom(N-i-1,:);
    end
    Defected_Geom(BeginIndexType6+21,:) = TempGeomType6;
    AminoBP(Defect_Location+1:15)=AminoBP(Defect_Location+1:15)+1;
    AminoBP = [AminoBP(1:Defect_Location), 10*Defect_Location+6 , AminoBP(Defect_Location+1:15)];
    
    
end

end

