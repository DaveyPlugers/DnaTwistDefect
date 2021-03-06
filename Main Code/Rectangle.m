function [p1,p2,p3,p4,p5,p6,p7,p8,Colour] = TransferMatrix(p1,p2,p3,p4,p5,p6,p7,p8,Colour,Side)
%TRANSFERMATRIX Summary of this function goes here
%   Detailed explanation goes here

fill3([p1(1) p2(1) p3(1) p4(1)], [p1(2) p2(2) p3(2) p4(2)], [p1(3) p2(3) p3(3) p4(3)], Colour);
hold on
fill3([p5(1) p6(1) p7(1) p8(1)], [p5(2) p6(2) p7(2) p8(2)], [p5(3) p6(3) p7(3) p8(3)], Colour);
fill3([p2(1) p6(1) p7(1) p3(1)], [p2(2) p6(2) p7(2) p3(2)], [p2(3) p6(3) p7(3) p3(3)], Colour);
fill3([p2(1) p6(1) p5(1) p1(1)], [p2(2) p6(2) p5(2) p1(2)], [p2(3) p6(3) p5(3) p1(3)], Colour);
fill3([p4(1) p8(1) p5(1) p1(1)], [p4(2) p8(2) p5(2) p1(2)], [p4(3) p8(3) p5(3) p1(3)], Colour);
if Side
    fill3([p4(1) p8(1) p7(1) p3(1)], [p4(2) p8(2) p7(2) p3(2)], [p4(3) p8(3) p7(3) p3(3)], [1,0,1]);
else
    fill3([p4(1) p8(1) p7(1) p3(1)], [p4(2) p8(2) p7(2) p3(2)], [p4(3) p8(3) p7(3) p3(3)], Colour)
end
end

