function [Defected_Geom,AminoBP] = ReversedDefectIntroducer(DNA_Geometry,Defect_Location,Defect_Type)
%REVERSEDDEFECTINTRODUCER 
%   Input:
%   DNA_Geometry = The geometry of the DNA that we will change
%   Defect_Location = Between which BP's will we add this defect (1 => 6-16
%   for K10, 1 => 6-26 for K20 then +10*n for location 1 + n)
%   Defect_type = What type of defect are we adding (1 => K10 undertwist,
%   2=> K20 undertwist) Overtwist not yet added in this function.
%   Output:
%   Defected_Geom = The geometry with the defect integrated into it (but
%   now we do it in a way our defect "moves" the other way
%   So it takes the bp from the left, where the original takes the bp from
%   the right
%   AminoBP = Array saving all the positions of binding sites
Defect_Type = Defect_Type-4; %We adjust since in other code we do +4, it is sloppy but we can fix it if we wish to continue using the code
Gamma = 4.46;
Theta_Twist = 35.575;
Raise = 3.4; %Angstrom
Phase_Shift = -6*pi/20;
N=147;
AminoBP = [6,16,26,36,46,56,66,76,86,96,106,116,126,136,146];


if Defect_Type == 1 %Undertwist 10 to 11
    BeginIndexType5 = 10*(Defect_Location-1)+6;
    Defected_Geom = DNA_Geometry;
    TempGeomType5 = DNA_Geometry(BeginIndexType5-1,:);
    for i=BeginIndexType5:BeginIndexType5+10
        Defected_Geom(i-1,:) = [0,0,Raise*10/11,Theta_Twist*10/11,5/5.5*Gamma*cos(2*pi*(i-BeginIndexType5)*10/(11*10) + 2*pi*(BeginIndexType5)/10 -Phase_Shift),5/5.5*Gamma*sin(2*pi*(i-BeginIndexType5)*10/(11*10)+ 2*pi*(BeginIndexType5)/10-Phase_Shift)];
        %Defected_Geom(i,:) = [0,0,Raise*9/10,Theta_Twist*9/10,5/5.5*Gamma*cos(2*pi*(i-1-BeginIndexType4)*9/(10*10) + 2*pi*BeginIndexType4/10 -Phase_Shift),5/5.5*Gamma*sin(2*pi*(i-1-BeginIndexType4)*9/(10*10)+ 2*pi*BeginIndexType4/10-Phase_Shift)];
    end
    
    for i=1:BeginIndexType5-3
        Defected_Geom(i,:) = Defected_Geom(i+1,:);
    end
    Defected_Geom(BeginIndexType5-2,:) = TempGeomType5;
    AminoBP(1:Defect_Location)=AminoBP(1:Defect_Location)-1;

elseif Defect_Type == 2
    BeginIndexType6 = 10*(Defect_Location-1)+6;
    Defected_Geom = DNA_Geometry;
    TempGeomType6 = DNA_Geometry(BeginIndexType6-1,:);
    for i=BeginIndexType6:BeginIndexType6+20
        Defected_Geom(i-1,:) = [0,0,Raise*20/21,Theta_Twist*20/21,20/21*Gamma*cos(2*pi*(i-BeginIndexType6)*20/(21*10) + 2*pi*(BeginIndexType6)/10 -Phase_Shift),20/21*Gamma*sin(2*pi*(i-BeginIndexType6)*20/(21*10)+ 2*pi*(BeginIndexType6)/10-Phase_Shift)];
        %Defected_Geom(i,:) = [0,0,Raise*19/20,Theta_Twist*19/20,20/21*Gamma*cos(2*pi*(i-1-BeginIndexType6)*19/(20*10) + 2*pi*BeginIndexType6/10 -Phase_Shift),20/21*Gamma*sin(2*pi*(i-1-BeginIndexType6)*19/(20*10)+ 2*pi*BeginIndexType6/10-Phase_Shift)];
    end
    
    for i=1:BeginIndexType6-3
        Defected_Geom(i,:) = Defected_Geom(i+1,:);
    end
    Defected_Geom(BeginIndexType6-2,:) = TempGeomType6;
    AminoBP(1:Defect_Location)=AminoBP(1:Defect_Location)-1;
    AminoBP = [AminoBP(1:Defect_Location), 10*Defect_Location+6 , AminoBP(Defect_Location+1:15)];


end

