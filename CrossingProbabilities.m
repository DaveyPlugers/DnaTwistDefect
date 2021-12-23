function [Crossing_Probabilities] = CrossingProbabilities(Peak_Defect_Energies)
%CROSSINGPROBABILITIES 
%   Input:
%   Peak_Defect_Energies = 14 wide array with the energies of the energy peaks that need to be crossed
%   Output:
%   Crossing_Probabilities = 2 by 13 wide array that gives the chance to go
%   left (First index= 1) or right (First index= 2) in the different energy
%   values between the peaks.
%   Crossing_Probabilities(1,i) Gives the chance to go left over peak i and
%   Crossing_Probabilities(2,i) gives the probability to go right over peak i+1

Crossing_Probabilities = zeros(2,length(Peak_Defect_Energies)-1);
for i=1:length(Peak_Defect_Energies)-1
    Crossing_Probabilities(1,i) = 1/(1+exp(-Peak_Defect_Energies(i+1)+Peak_Defect_Energies(i)));
    Crossing_Probabilities(2,i) = exp(-Peak_Defect_Energies(i+1)+Peak_Defect_Energies(i))/(1+exp(-Peak_Defect_Energies(i+1)+Peak_Defect_Energies(i)));
end

end

