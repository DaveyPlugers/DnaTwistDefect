function [Defected_Geom,AminoBP] = ReversedDefectIntroducer(DNA_Geometry,Defect_Location,Defect_Type)
%REVERSEDDEFECTINTRODUCER 
%   Input:
%   DNA_Geometry = The geometry of the DNA that we will change
%   Defect_Location = Between which BP's will we add this defect (1 => 7-17
%   for K10, 1 => 7-27 for K20 then +10*n for location 1 + n)
%   Defect_type = What type of defect are we adding (1 => K10 undertwist,
%   2=> K20 undertwist) Overtwist not yet added in this function.
%   Output:
%   Defected_Geom = The geometry with the defect integrated into it (but
%   now we do it in a way our defect "moves" the other way
%   So it takes the bp from the left, where the original takes the bp from
%   the right
%   AminoBP = Array saving all the positions of binding sites

Gamma = 4.46;
Theta_Twist = 35.575;
Raise = 3.4; %Angstrom
Phase_Shift = -6*pi/20;
N=147;
AminoBP = [7,17,27,37,47,57,67,77,87,97,107,117,127,137,147];


if Defect_Type == 1 %Undertwist 9 to 10
    BeginIndexType4 = 10*(Defect_Location-1)+7;
    Defected_Geom = DNA_Geometry;
    TempGeomType4 = DNA_Geometry(BeginIndexType4,:);
    for i=BeginIndexType4:BeginIndexType4+9
        Defected_Geom(i,:) = [0,0,Raise*9/10,Theta_Twist*9/10,5/5.5*Gamma*cos(2*pi*(i-1-BeginIndexType4)*9/(10*10) + 2*pi*BeginIndexType4/10 -Phase_Shift),5/5.5*Gamma*sin(2*pi*(i-1-BeginIndexType4)*9/(10*10)+ 2*pi*BeginIndexType4/10-Phase_Shift)];
    end
    
    for i=1:BeginIndexType4-2
        Defected_Geom(i,:) = Defected_Geom(i+1,:);
    end
    Defected_Geom(BeginIndexType4-1,:) = TempGeomType4;
    AminoBP(1:Defect_Location)=AminoBP(1:Defect_Location)-1;

elseif Defect_Type == 2
    BeginIndexType6 = 10*(Defect_Location-1)+7;
    Defected_Geom = DNA_Geometry;
    TempGeomType6 = DNA_Geometry(BeginIndexType6,:);
    for i=BeginIndexType6:BeginIndexType6+19
        Defected_Geom(i,:) = [0,0,Raise*19/20,Theta_Twist*19/20,20/21*Gamma*cos(2*pi*(i-1-BeginIndexType6)*19/(20*10) + 2*pi*BeginIndexType6/10 -Phase_Shift),20/21*Gamma*sin(2*pi*(i-1-BeginIndexType6)*19/(20*10)+ 2*pi*BeginIndexType6/10-Phase_Shift)];
    end
    
    for i=1:BeginIndexType6-2
        Defected_Geom(i,:) = Defected_Geom(i+1,:);
    end
    Defected_Geom(BeginIndexType6-1,:) = TempGeomType6;
    AminoBP(1:Defect_Location)=AminoBP(1:Defect_Location)-1;
    AminoBP = [AminoBP(1:Defect_Location), 10*Defect_Location+6 , AminoBP(Defect_Location+1:15)];


end

