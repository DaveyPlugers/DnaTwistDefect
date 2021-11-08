function [outputArg1,outputArg2] = PlotterFunc(inputArg12,ColourType,AminoBP)
%PLOTTERFUNC Summary of this function goes here
%   Detailed explanation goes here
if ColourType == 1
    Colour = [0,1,0];
elseif ColourType == 2
    Colour = [1,0,0]
end
Geometric_Values=inputArg12


Length = length(Geometric_Values(:,1));
HalfWidth = 2;
HalfHeight = 0.5;
HalfDepth=5;


Start = [[0,0,0]];

q1 = (Start + [-HalfWidth,-HalfDepth,-HalfHeight]).';
q2 = (Start + [HalfWidth,-HalfDepth,-HalfHeight]).';
q3 = (Start + [HalfWidth,HalfDepth,-HalfHeight]).';
q4 = (Start + [-HalfWidth,HalfDepth,-HalfHeight]).';
q5 = (Start + [-HalfWidth,-HalfDepth,HalfHeight]).';
q6 = (Start + [HalfWidth,-HalfDepth,HalfHeight]).';
q7 = (Start + [HalfWidth,HalfDepth,HalfHeight]).';
q8 = (Start + [-HalfWidth,HalfDepth,HalfHeight]).';

p1=q1;
p2=q2;
p3=q3;
p4=q4;
p5=q5;
p6=q6;
p7=q7;
p8=q8;

%fill3([p1(1) p2(1) p3(1) p4(1)], [p1(2) p2(2) p3(2) p4(2)], [p1(3) p2(3) p3(3) p4(3)], 1);
%hold on
%fill3([p5(1) p6(1) p7(1) p8(1)], [p5(2) p6(2) p7(2) p8(2)], [p5(3) p6(3) p7(3) p8(3)], 1);
%fill3([p2(1) p6(1) p7(1) p3(1)], [p2(2) p6(2) p7(2) p3(2)], [p2(3) p6(3) p7(3) p3(3)], 1);
%fill3([p2(1) p6(1) p5(1) p1(1)], [p2(2) p6(2) p5(2) p1(2)], [p2(3) p6(3) p5(3) p1(3)], 1);
%fill3([p4(1) p8(1) p5(1) p1(1)], [p4(2) p8(2) p5(2) p1(2)], [p4(3) p8(3) p5(3) p1(3)], 1);
%fill3([p4(1) p8(1) p7(1) p3(1)], [p4(2) p8(2) p7(2) p3(2)], [p4(3) p8(3) p7(3) p3(3)], 1);
Position = [0,0,0].';
Rectangle(p1,p2,p3,p4,p5,p6,p7,p8,Colour,false)
hold on

for i=1:Length
    if ColourType == 1
        Colour = Colour + [0,-1/Length,1/Length];
        if i==Length
            Colour = [0,0,1];
        end
    end
    Axis_Rotation = [[-(p8-p7)/(2*HalfWidth)],[(p8-p5)/(2*HalfDepth)],[(p8-p4)/(2*HalfHeight)]];
    Position = Position +  Axis_Rotation*Geometric_Values(i,1:3).';
    Twistangle = -deg2rad(Geometric_Values(i,4));
    Tiltangle = -deg2rad(Geometric_Values(i,5));
    Rollangle = deg2rad(Geometric_Values(i,6));
    TwistMatrix = [[cos(Twistangle),-sin(Twistangle),0];[sin(Twistangle),cos(Twistangle),0];[0,0,1]];
    RollMatrix = [[1,0,0];[0,cos(Rollangle),-sin(Rollangle)];[0,sin(Rollangle),cos(Rollangle)]];
    TiltMatrix = [[cos(Tiltangle),0,sin(Tiltangle)];[0,1,0];[-sin(Tiltangle),0,cos(Tiltangle)]];


    p1 = Axis_Rotation*TiltMatrix*RollMatrix*TwistMatrix*q1 + Position;
    p2 = Axis_Rotation*TiltMatrix*RollMatrix*TwistMatrix*q2 + Position;
    p3 = Axis_Rotation*TiltMatrix*RollMatrix*TwistMatrix*q3 + Position;
    p4 = Axis_Rotation*TiltMatrix*RollMatrix*TwistMatrix*q4 + Position;
    p5 = Axis_Rotation*TiltMatrix*RollMatrix*TwistMatrix*q5 + Position;
    p6 = Axis_Rotation*TiltMatrix*RollMatrix*TwistMatrix*q6 + Position;
    p7 = Axis_Rotation*TiltMatrix*RollMatrix*TwistMatrix*q7 + Position;
    p8 = Axis_Rotation*TiltMatrix*RollMatrix*TwistMatrix*q8 + Position;
    if ismember(i,AminoBP)
        Rectangle(p1,p2,p3,p4,p5,p6,p7,p8,Colour,true);
    else
        Rectangle(p1,p2,p3,p4,p5,p6,p7,p8,Colour,false);
    end

%Vanaf hier, bepaal normaal vector van het nieuwe vlak, laat de
%verplaatsing van de raise, shift en slide hierop verder gaan.
%Slide zal in dezelfde richting als een p1,p3 paar zijn en shift p1 p2 dan
%Raise zal loodrecht staan ofwel p1 p5 paar richting

%Normaliseer deze vectoren en vermenigvuldig ze dan met de
%raise/slide/shift parameters

%p8-p7 width, p8-p5 depth, p8-p4 height

end



end

