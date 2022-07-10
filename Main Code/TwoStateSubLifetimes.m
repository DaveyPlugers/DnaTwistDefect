function [Expected_Lifetime] = TwoStateSubLifetimes(Side,Occupancy,knmin,knplus,Previous_Lifetime)
%TWOSTATESUBLIFETIMES Summary of this function goes here
%   Detailed explanation goes here

%knmin = 2;
%knplus = 2;
%Previous_Lifetime = 10;
%Occupancy = 0.1;
%Side = 1;  
%Side 1 means we have a particle to the left, Side 2 means to the right
%For the second most right position we don't give the actual occupancy but
%an approximation

if Side==1
    %Expected_Lifetime = 1/(knplus + knmin) + knmin/knplus*(1/(knplus+knmin) + (1-Occupancy)*(Previous_Lifetime ));
    Expected_Lifetime = 1/knplus + knmin/knplus*(1-Occupancy)*(Previous_Lifetime);
    
elseif Side==2
    %Expected_Lifetime = 1/(1-Occupancy)*(Occupancy/(knplus + knmin) + knmin/knplus*(Previous_Lifetime + 1/(knplus + knmin)) + (1-Occupancy)/(knplus + knmin));
    Expected_Lifetime = 1/(1-Occupancy)*(1/knplus + knmin/knplus*Previous_Lifetime);
end





end

