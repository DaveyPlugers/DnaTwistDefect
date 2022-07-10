function [ExpectedLifetime] = TwoStateTheoreticalTime(Side,knplus,knmin,krplus,krmin,krplus2,krmin2,ExpPrevious,ExpN,ExpNPrevious)
%TWOSTATETHEORETICALTIME Summary of this function goes here
%   Detailed explanation goes here


%Side=0;


if Side==1
    %knplus = 2; %Rate to go to the right for our particle in position n
    %knmin = 2; %Rate to go to the left for our particle in position n
    %krplus = 2; %Rate to go to the right for boundary particle in position n-1
    %krmin = 2; %Rate to go to the left for boundary particle in position n-1
    %krplus2 = 2; %Rate to go to the right for boundary particle in position n-2
    %krmin2 = 2; %Rate to go to the left for boundary particle in position n-2
    %ExpN = 3; %Expected value E[n^*]
    %ExpPrevious = 3; %E[(n-2,n-1)]

    
     

    ExpectedLifetime = 1/(1-krmin*krplus2/((knplus + krmin)*(knplus + knmin + krplus2 + krmin2)))*( ...   %Pre-factor in front
    knplus/(knplus + krmin)^2 + ...   %Just step to the right
    krmin/(knplus + krmin)*(... %Pre factor to the 4 terms
    knplus/(knplus + knmin + krplus2 + krmin2)*(1/(knplus + krmin) + 1/(knplus + knmin + krplus2 + krmin2))  + ... %left then n to right
    knmin/(knplus + knmin + krplus2 + krmin2)*(1/(knplus + krmin) + 1/(knplus + knmin + krplus2 + krmin2) + ExpPrevious + ExpN) + ...   %left then n to left
    krplus2/(knplus + knmin + krplus2 + krmin2)*(1/(knplus + krmin) + 1/(knplus + knmin + krplus2 + krmin2)) + ... %left then n-2 to right
    krmin2/(knplus + knmin + krplus2 + krmin2)*(1/(knplus + krmin) + 1/(knplus + knmin + krplus2 + krmin2) + ExpN)));    %left then n-2 to left

elseif Side==2
    
    %knplus = 2; %Rate to go to the right for our particle in position n
    %knmin = 2; %Rate to go to the left for our particle in position n
    %krplus = 2; %Rate to go to the right for boundary particle in position n+1
    %krmin = 2; %Rate to go to the left for boundary particle in position n+1
    %krplus2 = 2; %Rate to go to the right for boundary particle in position n+2
    %krmin2 = 2; %Rate to go to the left for boundary particle in position n+2
    %ExpN = 3; %Expected value E[n^#]
    %ExpNPrevious %Expected value E[n-1^#]
    
    ExpectedLifetime = 1/(1-krplus*krmin2/((knmin + krplus)*(knplus + knmin + krplus2 + krmin2)))*( ...   %Pre-factor in front
    knmin/(knmin + krplus)*(1/(knmin + krplus) + ExpNPrevious + ExpN) + ... %We just take step to the left
    krplus/(knmin + krplus)*( ... %pre-factor to the 4 terms
    knplus/(knplus + knmin + krplus2 + krmin2)*(1/(knmin + krplus) + 1/(knplus + knmin + krplus2 + krmin2)) + ... %Other steps to the right, then main particle steps to the right
    knmin/(knplus + knmin + krplus2 + krmin2)*(1/(knmin + krplus) + 1/(knplus + knmin + krplus2 + krmin2) + ExpNPrevious + ExpN) + ...%Other steps to the right, then main particle steps to the left
    krplus2/(knplus + knmin + krplus2 + krmin2)*(1/(knmin + krplus) + 1/(knplus + knmin + krplus2 + krmin2) + ExpN) + ... %Other steps to the right, then other particle steps to the right
    krmin2/(knplus + knmin + krplus2 + krmin2)*(1/(knmin + krplus) + 1/(knplus + knmin + krplus2 + krmin2)))); %Other steps to the right, then other particle steps to the left


end



end

