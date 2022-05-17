


%We take DNA string, we calculate the energy of the wrapping around the
%nucleosomes in all positions and use it to get a statistical distribution
%of the position where we would find it

DNA = 'GGTGGGCTCCGAAAAATTCGCAGGGCGACCGGCGAAATGCTCAAAAAATCAAAAAATATTCCCTGAAACAAAAGCAACTCTTGCGAAACGGGCGCATTTTTAAAAAATTCGGAGAAAATTTTTGAAATCGTGACGCATCTCTTGCCGTTCGCCAAACAATCCAGAAATTTTATTGTTCAGGCAAAATTCACTAGGTTTTATGGTGAAACGCAAAAAATTCTGACGTTTTCATGAACGTCTTTTAGGATTTTCAGGTTAATGCGGTTGGGCTCCAAAATCTCAAGCAGTCTCGTAGAAATTTCAAAAATTTCGGCATTCTTGGAGAATCACGGGAAGATACTCGCCGGTTT'; 
M = length(DNA);
DNA = [DNA DNA];

Energies = zeros(1,M);
N=147;
for i=1:M
    DNAString = DNA(i:i+M);
    [DNA_Geometry,DNAIndexation] = GeomArrayMaker(DNAString,1,Geometric_Properties);%1 is to denote the mode, see function description
    Undefected_Energy = EnergyValuesCalculator(DNA_Geometry,DNAIndexation,Geometric_Properties,true,true,true,true);
    Energies(i) = sum(sum(Undefected_Energy));   
end

%Now we have the raw energies but we want to damp them in a smart way so
%that we preserve the ratio of lowest to highest energy similar to in the
%literature
x = linspace(1,350,350)
myfittype = fittype('b + a*sin(2*pi*x/10+g)*(1+c*x+d*x^2)','dependent',{'Energies'},'independent',{'x'},'coefficients',{'a','b','c','d','g'});
startPoints = [200 40 0 0 0];
myfit = fit(x',Energies',myfittype,'Start',startPoints)
%6 = a*Damp, Damp = 6/a
%New_Energies = Damp*(Energies - b) + b
plot(Energies)
hold on
NewEnergies = 0.18*(Energies - 206) + 206;
plot(NewEnergies)
figure
%Now that we have done this we calculate the actual statistical
%distribution

Occupancy = zeros(1,M);
Z = sum(exp(-NewEnergies));
Occupancy = exp(-NewEnergies)/Z;
figure
plot(Occupancy)