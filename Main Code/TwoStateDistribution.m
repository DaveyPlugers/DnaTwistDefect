function [Occupancy] = TwoStateDistribution(Barriers,Side,Positions)
%TWOSTATEDISTRIBUTION Summary of this function goes here
%   Detailed explanation goes here

Filler = zeros(2,9); %Should replace 9 with n later
Filler(:,:) = Barriers(:,:,:);
Barriers = Filler;

if Side==1
    Barriers = Barriers(:,1:Positions);
elseif Side==2
    Barriers = Barriers(:,end-Positions:end);
end


if Positions==1
    'This is an error from the TwoStateDistribution code, Tried to calculate for one position, double check values'  
end    
%Barriers = zeros(2,Positions);
%Barriers(1,1) = 2;
%Barriers(1,2) = 2;
%Barriers(2,1) = 2;
%Barriers(2,2) = 2;
%Barriers(1,3) = 8;
%Barriers(2,3) = 4;
%Barriers(1,4) = 2;
%Barriers(2,4) = 2;

Barrier_Parameters = zeros(Positions,Positions);

Barrier_Parameters(Positions,:) = 1;

Solutions = zeros(Positions,1);
Solutions(Positions) = 1;
for i=1:Positions-1
    if i==1
        Barrier_Parameters(1,1) = Barriers(2,1);% + k+_1
        Barrier_Parameters(1,2) = -Barriers(1,2);% - k-_2
    elseif i==Positions-1
        Barrier_Parameters(Positions-1,Positions-1) = Barriers(2,Positions-1);% + k+_(n-1)
        Barrier_Parameters(Positions-1,Positions) = -Barriers(1,Positions);% - k-_n
    else
        Barrier_Parameters(i,i) = Barriers(2,i) + Barriers(1,i);% + (k+_i + k-_i)
        Barrier_Parameters(i,i-1) = - Barriers(2,i-1);% - k+_(i-1)
        Barrier_Parameters(i,i+1) = - Barriers(1,i+1);% - k-_(i+1)
    end
    
end

X = linsolve(Barrier_Parameters,Solutions);

%Side=0, we are to the left so we need the right most value
%Side=1, we are to the right so we need the left most value
if Side==1
    Occupancy = X(Positions);
elseif Side==2
    Occupancy = X(1);
end




%end
