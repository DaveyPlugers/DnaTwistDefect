%function [outputArg1,outputArg2] = Full_DNA_Energy_Landscape(inputArg1,inputArg2)
%FULL_DNA_ENERGY_LANDSCAPE Summary of this function goes here
%   Takes a DNA string as input and calculates the energy landscape for
%   every possible position

DNA = 'ACATACATTCTTTTAAATTTCGAAGACCTCGTGATTCTTCCCGCGTTTACTGGCAAAATCTTAGAAATCGCGGCAAAGTCAAAAAAATTTTACAAATTGTCGTCGGGTCCTTAGAATTCGGCAGTCGTTTCTGGGTTCACGCACGACTTCTTTGAAAGGTCAAAAGCCGCCGAGTTTTTTCGGCGTAACTCCCGATCAGACCGAAACATCTGAAAACTTTCCTGTGACTACGCAGACGTCCGGGAAGTTCAAATCGTTTCGGGAGTTATCTTGAGATTTCTTCAAAACACACCAGCGTTTTCGAAAAATCTCCGTAATCACGCAAAACACGCTACGAATCAGAAAACGAAC';
M = length(DNA);
DNA_Extended = [DNA DNA]; %Make longer since we want every position and we take it to be periodic

N=147;
Gamma = 4.46;
Theta_Twist = 35.575;
Raise = 3.4; %Angstrom
Phase_Shift = 0;
Defect_Type = 1; %1 is for overtwist, 3 is for undertwist (Defect_Type +1 is for their respective K20 versions)
Negative_Barrier_Replacement = 0.1;
%Sometimes there is no actual barrier and it seems to be negative, to
%resolve this we replace this by a really small barrier value given here


Tilt_Array = [[-1.4,0.0,-0.1,-1.7];[0.0,1.4,1.5,-0.5];[0.5,1.7,0.1,0.0];[-1.5,0.1,0.0,-0.1]];
K_Tilt_Array = [[0.100,0.166,0.111,0.149];[0.148,0.100,0.087,0.082];[0.082,0.149,0.119,0.068];[0.087,0.111,0.082,0.119]];
Roll_Array = [[0.7,1.1,0.7,4.5];[3.3,0.7,1.9,4.7];[4.7,4.5,3.6,5.4];[1.9,0.7,0.3,3.6]];
K_Roll_Array = [[0.049,0.055,0.080,0.096];[0.029,0.049,0.046,0.048];[0.048,0.096,0.064,0.050];[0.046,0.080,0.082,0.064]];
Twist_Array = [[35.1,29.3,31.5,31.9];[37.8,35.1,36.3,37.3];[37.3,31.9,32.9,36.1];[36.3,31.5,33.6,32.9]];
K_Twist_Array = [[0.092,0.070,0.073,0.064];[0.052,0.092,0.071,0.043];[0.043,0.064,0.041,0.047];[0.071,0.073,0.055,0.041]];
Raise_Array = [[3.27,3.31,3.36,3.34];[3.42,3.27,3.37,3.33];[3.33,3.34,3.42,3.39];[3.37,3.36,3.40,3.42]];
K_Raise_Array = [[21.748,25.547,23.860,29.496];[21.914,21.748,22.820,18.235];[18.235,29.496,30.312,14.164];[22.820,23.860,25.860,30.312]];


Geometric_Properties = Tilt_Array;
Geometric_Properties(:,:,2) = K_Tilt_Array;
Geometric_Properties(:,:,3) = Roll_Array;
Geometric_Properties(:,:,4) = K_Roll_Array;
Geometric_Properties(:,:,5) = Twist_Array;
Geometric_Properties(:,:,6) = K_Twist_Array;
Geometric_Properties(:,:,7) = Raise_Array;
Geometric_Properties(:,:,8) = K_Raise_Array;

Amino_Unbinding_Cost = [0,repmat([12,0],1,13)];
Full_Energy_Landscape = zeros(2*M,27); %there is redundancy in here but makes it easier to follow
%First M is movement to the left, second M is movement to the right of the
%defect, making it a 3d array makes matlab act stupid


for w=1:M %We go over every position, this is only 1 of the directions though (defect moving to the left)
    
   DNAString = DNA_Extended(w:N+w);  %WE NEED TO CHANGE SOMETHING HERE
   %IT SEEMS WE START AT THE FIRST INDEX OF THE DNA, IN REALITY THE DEFECT
   %STARTS FROM INDEX 6 UP TO 146!!!!!!!!!!!!!
   %!!!!!!!!!!!!!!!
   %!!!!!!!!!!!!!!!
   %!!!!!!!!!!!!!!!
   
   [DNA_Geometry,DNAIndexation] = GeomArrayMaker(DNAString,1,Geometric_Properties);%1 is to denote the mode, see function description
   Undefected_Energy = EnergyValuesCalculator(DNA_Geometry,DNAIndexation,Geometric_Properties,true,true,true,true);

   for k=1:14 %K10 loop
        [Defected_DNA_Geom,AminoBP] = DefectIntroducer(DNA_Geometry,k,Defect_Type);
        Temporary_Code = EnergyValuesCalculator(Defected_DNA_Geom,DNAIndexation,Geometric_Properties,true,true,true,true);
        
        
        %Here we only take a small part of the energy landscape which
        %we use for the later movement
        Full_Energy_Landscape(w,(2*k-1)) = sum(sum(Temporary_Code(:)-Undefected_Energy(:)));
   end
   
   for k=1:13 %K10 loop
        [Defected_DNA_Geom,AminoBP] = DefectIntroducer(DNA_Geometry,k,Defect_Type+1);
        Temporary_Code = EnergyValuesCalculator(Defected_DNA_Geom,DNAIndexation,Geometric_Properties,true,true,true,true);
        
        
        %Here we only take a small part of the energy landscape which
        %we use for the later movement
        Full_Energy_Landscape(w,(2*k)) = sum(sum(Temporary_Code(:)-Undefected_Energy(:)));
   end
   Full_Energy_Landscape(w,:) = Full_Energy_Landscape(w,:) + Amino_Unbinding_Cost;
    
end



DNA_Extended_Moved = [DNA(M) DNA_Extended]; %Here we add the amino at the end to the front as well
%Because when it moves the other way we want to make sure we are looking at
%the same base DNA, ReversedDefectIntroducer also used that everything has
%already been moved by one (amino binding sites at 10*k + 7 instead of 10*k + 6)

%What used to be on position 6 is now on position 7 which is what we want

for w=1:M
    DNAString = DNA_Extended_Moved(w:N+w);
    [DNA_Geometry,DNAIndexation] = GeomArrayMaker(DNAString,3,Geometric_Properties);%3 is to denote the mode, see function description
    Undefected_Energy = EnergyValuesCalculator(DNA_Geometry,DNAIndexation,Geometric_Properties,true,true,true,true);

    
    for k=1:14 %K10 loop
        [Defected_DNA_Geom,AminoBP] = ReversedDefectIntroducer(DNA_Geometry,k,Defect_Type);
        Temporary_Code = EnergyValuesCalculator(Defected_DNA_Geom,DNAIndexation,Geometric_Properties,true,true,true,true);
        
        
        %Here we only take a small part of the energy landscape which
        %we use for the later movement
        Full_Energy_Landscape(M+w,(2*k-1)) = sum(sum(Temporary_Code(:)-Undefected_Energy(:)));
   end
   
   for k=1:13 %K10 loop
        [Defected_DNA_Geom,AminoBP] = ReversedDefectIntroducer(DNA_Geometry,k,Defect_Type+1);
        Temporary_Code = EnergyValuesCalculator(Defected_DNA_Geom,DNAIndexation,Geometric_Properties,true,true,true,true);
        
        
        %Here we only take a small part of the energy landscape which
        %we use for the later movement
        Full_Energy_Landscape(M+w,(2*k)) = sum(sum(Temporary_Code(:)-Undefected_Energy(:)));
   end
   Full_Energy_Landscape(M+w,:) = Full_Energy_Landscape(M+w,:) + Amino_Unbinding_Cost;
    
end

%Energy barrier calculations + checking how often it is negative, we need
%to fix this issue

Barrier_Right = zeros(2*M,14);
Barrier_Left = zeros(2*M,14);
Negatives = 0;
Positives = 0;
for w=1:2*M
    
    Barrier_Right(w,1) = Full_Energy_Landscape(w,2) - Full_Energy_Landscape(w,1);
    if Barrier_Right(w,1)<0
        Negatives = Negatives + 1;
         Barrier_Right(w,1) = Negative_Barrier_Replacement;
    else
        Positives = Positives + 1;
    end
    
    for k=2:13
        Barrier_Right(w,k) = Full_Energy_Landscape(w,2*k) - Full_Energy_Landscape(w,2*k-1);
        Barrier_Left(w,k) = Full_Energy_Landscape(w,2*(k-1)) - Full_Energy_Landscape(w,2*k-1);
        
        if Barrier_Right(w,k)<0
            Negatives = Negatives + 1;
            Barrier_Right(w,k) = Negative_Barrier_Replacement;
        else
            Positives = Positives + 1;
        end
        
        if Barrier_Left(w,k)<0
            Negatives = Negatives + 1;
            Barrier_Left(w,k) = Negative_Barrier_Replacement;
        else
            Positives = Positives + 1;
        end
        
    end
    
    Barrier_Left(w,14) = Full_Energy_Landscape(w,26) - Full_Energy_Landscape(w,27);
    if Barrier_Left(w,14)<0
        Barrier_Left(w,14) = Negative_Barrier_Replacement;
        Negatives = Negatives + 1;
    else
        Positives = Positives + 1;
    end
    
    
end

%end

