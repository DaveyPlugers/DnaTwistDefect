%function [outputArg1,outputArg2] = Chromatin_Transition_Script(inputArg1,inputArg2)
%CHROMATIN_TRANSITION_SCRIPT Summary of this function goes here
%   Detailed explanation goes here

M=351; %Length of DNA string
n = 9;
p = 14-n;
Defect_Array = zeros(1,n);
Time = 0;
IntroductionRate = 10^(-5); %we change this later

Iterations = 10000000;
Introduction_Time = zeros(1,floor(Iterations/n^2));
Leave_Time = zeros(1,floor(Iterations/n^2));

Chr_Dir = 0; %Short for Chromatin_Direction, 0 is defect moves to the left, 1 is defect moves to the right
Passable_Chromatin = true;
Chromatin_Barrier = exp(-12); %Make this a better value later


k = 328; %what position does our simulation start at

Failed_Defect=zeros(2,1);
g=1;
h=1;
r=1;
for q=1:Iterations
    [Left_Movement,Right_Movement,Introduction_Possible] = LeftRightMovement(Defect_Array);
    
    %Here we can change some stuff to have our introduction rate be
    %sequence dependent, for now we have it constant, implement next time
    TransitionRate = Introduction_Possible*IntroductionRate;
    %!!!!!
    %!!!!!
    %!!!!!
    
    for y=1:n
        if y==1    %Chr_Dir*M gives 0 or M which denotes left or right, mod(k-1)+1 makes our position between 1 and M + our previous from Chr_Dir
            Right_Movement(y) = Right_Movement(y)*exp(-(Barrier_Right(Chr_Dir*M+mod(k-1,M)+1,p+1)));
            %p+1, p is because we have our chrom. at position 14-n=p
            if Passable_Chromatin %If we can go through we have to add this
                Left_Movement(y) = Defect_Array(y)*Chromatin_Barrier;
            end
            
            TransitionRate = TransitionRate + Right_Movement(y) + Left_Movement(y);
        elseif y==n
            l = mod(k - sum(Defect_Array(1:n-1))-1,M)+1; %Double check, something maybe wrong with index since it doesn't correspond with theoretical peak
            Left_Movement(y) = Left_Movement(y)*exp(-(Barrier_Left(Chr_Dir*M+l,14)));
            Right_Movement(y) = Defect_Array(y)*exp(-(3)); %Note this is the leaving one, change later
            TransitionRate = TransitionRate + Left_Movement(y) + Right_Movement(y);
        else
            
            l = mod(k - sum(Defect_Array(1:y-1))-1,M)+1; %Double check, something maybe wrong with index since it doesn't correspond
            Left_Movement(y) = Left_Movement(y)*exp(-(Barrier_Left(Chr_Dir*M+l,p+y)));
            Right_Movement(y) = Right_Movement(y)*exp(-(Barrier_Right(Chr_Dir*M+mod(l-1,M)+1,p+y)));
            TransitionRate = TransitionRate + Right_Movement(y) + Left_Movement(y);
        end
    end   



    TransitieCoefficient = rand();
    %TransitionRate
    TimeStep = -log(1-TransitieCoefficient)/TransitionRate;
    Time = Time + TimeStep;
    TransitieCoefficient = TransitieCoefficient*TransitionRate;

    

    if TransitieCoefficient < Introduction_Possible*IntroductionRate
        Defect_Array(1) = Defect_Array(1) + 1; %Here we introduce a twist, doing like this just to test if it works, replace with =1 later
        Introduction_Time(g) = Time;
        g = g+1;
        k = mod(k,M) +1;
    else
        TransitieCoefficient = TransitieCoefficient - Introduction_Possible*IntroductionRate;
        GillespieLoop = true;
        y=1;
        while GillespieLoop
            if TransitieCoefficient < Left_Movement(y)
                %Defect moves to left for this y
                GillespieLoop = false;
                Defect_Array(y) = Defect_Array(y)-1;
                if y==1
                    %We save when a defect failed, might be interesting
                    %to check later when it happens more often
                    Failed_Defect(:,r) = [k,Introduction_Time(g-1)]';
                    g=g-1;
                    k = mod(k-2,M) +1;
                    r = r+1;
                else
                    Defect_Array(y-1) = Defect_Array(y-1) + 1;
                end
            elseif TransitieCoefficient < Left_Movement(y) + Right_Movement(y)
                %Defect moves to right for this y
                GillespieLoop = false;
                Defect_Array(y) = Defect_Array(y)-1;
                if y==n
                    Leave_Time(h) = Time;
                    h = h+1;
                    %if Time > 10^10 
                    %   Time = Time - 10^10;
                    %end
                else
                    Defect_Array(y+1) = Defect_Array(y+1) + 1;
                end
            else
                TransitieCoefficient = TransitieCoefficient - (Left_Movement(y) + Right_Movement(y));
                y = y+1;
            end

        end
    end


end









%end

