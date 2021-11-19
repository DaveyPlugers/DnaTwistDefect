function [Energy] = EnergyValuesCalculator(DNA_Geometry,DNAIndexation,Geometric_Properties,UseRaise,UseTwist,UseTilt,UseRoll)
%ENERGYVALUESCALCULATOR 
%   Input:
%   DNA_Geometry = The Dna's geometrics which we use to calculate the
%   energy cost of forming
%   DNAIndexation = N+1 = 148 Numbers representing the DNA basepairs (ATCG => 1234)
%   Geometric_Properties = The geometric properties array so it doesn't
%   have to be remade every run
%   UseXXXX = Boolean to see if we want to use these values in our energy
%   calculation
%   Output:
%   Energy: A (N=147 x 6) array that gives the individual energy components
%   (shift and slide are added too but no option to calculate cost yet so always 0)
N = 147;
Energy = zeros(N,6);

if UseTilt
    for i=1:N
        Energy(i,5) = 0.5*Geometric_Properties(DNAIndexation(i),DNAIndexation(i+1),2)*(DNA_Geometry(i,5)-Geometric_Properties(DNAIndexation(i),DNAIndexation(i+1),1))^2;
    end
end

if UseRoll
    for i=1:N
        Energy(i,6) = 0.5*Geometric_Properties(DNAIndexation(i),DNAIndexation(i+1),4)*(DNA_Geometry(i,6)-Geometric_Properties(DNAIndexation(i),DNAIndexation(i+1),3))^2;
    end
end

if UseTwist
    for i=1:N
        Energy(i,4) = 0.5*Geometric_Properties(DNAIndexation(i),DNAIndexation(i+1),6)*(DNA_Geometry(i,4)-Geometric_Properties(DNAIndexation(i),DNAIndexation(i+1),5))^2;
    end
end

if UseRaise
    for i=1:N
        Energy(i,3) = 0.5*Geometric_Properties(DNAIndexation(i),DNAIndexation(i+1),8)*(DNA_Geometry(i,3)-Geometric_Properties(DNAIndexation(i),DNAIndexation(i+1),7))^2;
    end
end


end

