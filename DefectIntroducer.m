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
Phase_Shift = -6*pi/20;
N=147;
AminoBP = [6,16,26,36,46,56,66,76,86,96,106,116,126,136,146];


if Defect_Type == 1 %Undertwist 9 to 10, legacy code
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
elseif Defect_Type==3  %Overtwist 9 to 8 
    
    BeginIndexType5 = 10*(Defect_Location-1)+6;
    Defected_Geom = DNA_Geometry;
    for i=BeginIndexType5+1:BeginIndexType5+8
        Defected_Geom(i,:) = [0,0,Raise*9/8,Theta_Twist*9/8,5/4.5*Gamma*cos(2*pi*(i-BeginIndexType5)*9/(8*10) + 2*pi*BeginIndexType5/10 -Phase_Shift),5/4.5*Gamma*sin(2*pi*(i-BeginIndexType5)*9/(8*10)+ 2*pi*BeginIndexType5/10-Phase_Shift)];
    end
    for i=BeginIndexType5+9:N-1
        Defected_Geom(i,:) = Defected_Geom(i+1,:);
    end
    Defected_Geom(N,:) = [0,0,Raise,Theta_Twist,Gamma*cos(2*pi*(N+1)/10-Phase_Shift),Gamma*sin(2*pi*(N+1)/10-Phase_Shift)];
    AminoBP(Defect_Location+1:15)=AminoBP(Defect_Location+1:15)-1;
    
elseif Defect_Type ==4 %19 to 18
    Correctiefactor = 1.02;
    BeginIndexType7 = 10*(Defect_Location-1)+6;
    Defected_Geom = DNA_Geometry;
    for i=BeginIndexType7+1:BeginIndexType7+18
        Defected_Geom(i,:) = [0,0,Raise*19/18,Theta_Twist*19/18,Correctiefactor*19/18*Gamma*cos(2*pi*(i-BeginIndexType7)*19/(18*10) + 2*pi*BeginIndexType7/10 -Phase_Shift),Correctiefactor*19/18*Gamma*sin(2*pi*(i-BeginIndexType7)*19/(18*10)+ 2*pi*BeginIndexType7/10-Phase_Shift)];
    end
    for i=BeginIndexType7+19:N-1
        Defected_Geom(i,:) = Defected_Geom(i+1,:);
    end
    Defected_Geom(N,:) = [0,0,Raise,Theta_Twist,Gamma*cos(2*pi*(N+1)/10-Phase_Shift),Gamma*sin(2*pi*(N+1)/10-Phase_Shift)];
    AminoBP(Defect_Location+1:15)=AminoBP(Defect_Location+1:15)-1;
    AminoBP = [AminoBP(1:Defect_Location), AminoBP(Defect_Location+2:15)];
elseif Defect_Type==7  %Adding new version of the defects, these are correct, ADD MORE IN ALL DEFECT SCRIPTS
    %Correct overtwist (10 to 9 bp's)
    BeginIndexType7 = 10*(Defect_Location-1)+6;
    Defected_Geom = DNA_Geometry;
    for i=BeginIndexType7:BeginIndexType7+8
        Defected_Geom(i,:) = [0,0,Raise*10/9,Theta_Twist*10/9,5/4.5*Gamma*cos(2*pi*(i-BeginIndexType7)*10/(9*10) + 2*pi*BeginIndexType7/10 -Phase_Shift),5/4.5*Gamma*sin(2*pi*(i-BeginIndexType7)*10/(9*10)+ 2*pi*BeginIndexType7/10-Phase_Shift)];
    end
    for i=BeginIndexType7+9:N-1
        Defected_Geom(i,:) = Defected_Geom(i+1,:);
    end
    Defected_Geom(N,:) = [0,0,Raise,Theta_Twist,Gamma*cos(2*pi*(N+1)/10-Phase_Shift),Gamma*sin(2*pi*(N+1)/10-Phase_Shift)];
    AminoBP(Defect_Location+1:15)=AminoBP(Defect_Location+1:15)-1;

elseif Defect_Type ==8 %K20 overtwist (20 to 19 bp's)
    BeginIndexType8 = 10*(Defect_Location-1)+6;
    Defected_Geom = DNA_Geometry;
    
    for i=BeginIndexType8:BeginIndexType8+18
        Defected_Geom(i,:) = [0,0,Raise*20/19,Theta_Twist*20/19,20/19*Gamma*cos(2*pi*(i-BeginIndexType8)*20/(19*10) + 2*pi*BeginIndexType8/10 -Phase_Shift),20/19*Gamma*sin(2*pi*(i-BeginIndexType8)*20/(19*10)+ 2*pi*BeginIndexType8/10-Phase_Shift)];
    end
    
    for i=BeginIndexType8+19:N-1
        Defected_Geom(i,:) = Defected_Geom(i+1,:);
    end
    Defected_Geom(N,:) = [0,0,Raise,Theta_Twist,Gamma*cos(2*pi*(N+1)/10-Phase_Shift),Gamma*sin(2*pi*(N+1)/10-Phase_Shift)];
    AminoBP(Defect_Location+1:15)=AminoBP(Defect_Location+1:15)-1;
    AminoBP = [AminoBP(1:Defect_Location), AminoBP(Defect_Location+2:15)];    
    
elseif Defect_Type ==5 %10 to 11
    BeginIndexType5 = 10*(Defect_Location-1)+6;
    Defected_Geom = DNA_Geometry;
    TempGeomType4 = DNA_Geometry(BeginIndexType5+10,:);
    for i=BeginIndexType5:BeginIndexType5+10
        Defected_Geom(i,:) = [0,0,Raise*10/11,Theta_Twist*10/11,5/5.5*Gamma*cos(2*pi*(i-BeginIndexType5)*10/(11*10) + 2*pi*BeginIndexType5/10 -Phase_Shift),5/5.5*Gamma*sin(2*pi*(i-BeginIndexType5)*10/(11*10)+ 2*pi*BeginIndexType5/10-Phase_Shift)];
    end
    
    for i=0:N-BeginIndexType5-12
        Defected_Geom(N-i,:) = Defected_Geom(N-i-1,:);
    end
    Defected_Geom(BeginIndexType5+11,:) = TempGeomType4;
    AminoBP(Defect_Location+1:15)=AminoBP(Defect_Location+1:15)+1;
elseif Defect_Type ==6 %K20 20 to 21    
    BeginIndexType6 = 10*(Defect_Location-1)+6;
    Defected_Geom = DNA_Geometry;
    TempGeomType6 = DNA_Geometry(BeginIndexType6+20,:);
    for i=BeginIndexType6:BeginIndexType6+20
        Defected_Geom(i,:) = [0,0,Raise*20/21,Theta_Twist*20/21,20/21*Gamma*cos(2*pi*(i-BeginIndexType6)*20/(21*10) + 2*pi*BeginIndexType6/10 -Phase_Shift),20/21*Gamma*sin(2*pi*(i-BeginIndexType6)*20/(21*10)+ 2*pi*BeginIndexType6/10-Phase_Shift)];
    end
    
    for i=0:N-BeginIndexType6-21
        Defected_Geom(N-i,:) = Defected_Geom(N-i-1,:);
    end
    Defected_Geom(BeginIndexType6+21,:) = TempGeomType6;
    AminoBP(Defect_Location+1:15)=AminoBP(Defect_Location+1:15)+1;
    AminoBP = [AminoBP(1:Defect_Location), 10*Defect_Location+6 , AminoBP(Defect_Location+1:15)];
end

end

