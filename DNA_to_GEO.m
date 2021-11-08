function [outputArg1] = DNA_to_GEO(DNAString)
%DNA_TO_GEO Summary of this function goes here
%{A,T,C,G} = {1,2,3,4]

%[slide,shift,raise, twist,roll,tilt]
% [[aa,at,ac,ag];[ta,tt,tc,tg];[ca,ct,cc,cg];[ga,gt,gc,gg]]
Geometrics = [[-0.08,-0.59,-0.58,-0.25];[0.05,-0.08,0.09,0.53];[0.53,-0.25,-0.22,0.41];[0.09,-0.58,-0.38,-0.22]];
Geometrics(:,:,2) = [[-0.03,0,0.13,0.09];[0,0.03,0.28,-0.09];[0.09,-0.09,-0.05,0];[-0.28,-0.13,0,0.05]];
Geometrics(:,:,3) = [[3.27,3.31,3.36,3.34];[3.42,3.27,3.37,3.33];[3.33,3.34,3.42,3.39];[3.37,3.36,3.40,3.42]];
Geometrics(:,:,4) = [[35.1,29.3,31.5,31.9];[37.8,35.1,36.3,37.3];[37.3,31.9,32.9,36.1];[36.3,31.5,33.6,32.9]];
Geometrics(:,:,5) = [[0.7,1.1,0.7,4.5];[3.3,0.7,1.9,4.7];[4.7,4.5,3.6,5.4];[1.9,0.7,0.3,3.6]];
Geometrics(:,:,6) = [[-1.4,0,-0.1,-1.7];[0,1.4,1.5,-0.5];[0.5,1.7,0.1,0];[-1.5,0.1,0,-0.1]];

DNA = char(DNAString);

%First 2 strings
if DNA(1) == 'A'
    i=1;
    if DNA(2) == 'A'    
        j=1;
    elseif DNA(2) == 'T'
        j=2;
    elseif DNA(3) == 'C'
        j=3;
    else
        j=4;
    end
elseif DNA(1) == 'T'
    i=2;
    if DNA(2) == 'A'    
        j=1;
    elseif DNA(2) == 'T'
        j=2;
    elseif DNA(3) == 'C'
        j=3;
    else
        j=4;
    end
elseif DNA(1) == 'C'
    i=3;
    if DNA(2) == 'A'    
        j=1;
    elseif DNA(2) == 'T'
        j=2;
    elseif DNA(3) == 'C'
        j=3;
    else
        j=4;
    end
else 
    i=4;
    if DNA(2) == 'A'    
        j=1;
    elseif DNA(2) == 'T'
        j=2;
    elseif DNA(3) == 'C'
        j=3;
    else
        j=4;
    end
end

Array = [0,0,0,0,0,0];
for z=1:6
    Array(z) = (Geometrics(i,j,z));
end
Geometric_Values = Array;

for k=3:length(DNAString)
    i=j;
    if DNA(k) == 'A'
        j=1;
    elseif DNA(k) == 'T'
        j=2;
    elseif DNA(k) == 'C'
        j=3;
    else
        j=4;
    end
        
    Array = [0,0,0,0,0,0];
    for z=1:6
        Array(z) = (Geometrics(i,j,z));
    end
    Geometric_Values = [Geometric_Values;Array];
end
    

outputArg1 =  Geometric_Values;
end

