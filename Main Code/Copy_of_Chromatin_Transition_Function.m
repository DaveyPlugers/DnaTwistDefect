function [Introduction_Time,Leave_Time,k,Time] = Copy_of_Chromatin_Transition_Function(Duration,DNA,Position_Remodeller,Start_Position,Chr_Dir,IntroductionRate,Variable_Introrate,Passable_Chromatin,Chromatin_Barrier,Barrier_Left,Barrier_Right)
%CHROMATIN_TRANSITION_FUNCTION Creates a subsimulation of the twist defects
%motion. Takes begin position, duration and energy barriers, returns the
%times at which it was at a certain location once done
%   Input:
%   Duration = How long do we want to keep running
%   DNA = String that we are working on
%   Position_Remodeller = What SHL is our remodeller located
%   Start_Position = Where is the first bp around the nucleosome w.r.t. the
%   DNA
%   Chr_Dir = Chromatin_Direction, 0 is nucleosome moves to the right, 1 to
%   the left
%   IntroductionRate = How often do defects get introduced
%   Variable_Introrate = true/false Do we want it to be sequent dependent
%   Passable_Chromatin = true/false Do we want to allow defects to pass through
%   remodeller?
%   Chromatin_Barrier = What is the rate they can pass through the barrier
%   Barrier_Left = Energy barriers to move to the left
%   Barrier_Right = Energy barriers to move to the right
%   Output:
%   Introduction_Time = 2 by N array with timestamps and locations of when
%   we introduced defects
%   Leave_Time = 2 by N array with timestamps and locations of when a
%   defect got out
%   k = Position of the nucleosome once the simulation is done (can also
%   get it from previous arrays)
%   Time: What was the last time value before simulation ended
M = length(DNA);
DNA = [DNA DNA];

p = Position_Remodeller + 7;
n = 14 - p; %How many positions are there for the defect to move in
Defect_Array = zeros(1,n);
Time = 0;

% I just have a big number for this, if we end up having overflows fix this
Introduction_Time = zeros(2,50000);
Leave_Time = zeros(2,50000);


if Variable_Introrate
    IntroductionRate = zeros(1,M);
    SHL = -2;
    [DNA_Geometry,~] = GeomArrayMaker(DNAString,1,Geometric_Properties);
    
    DNAIndexation = []; %We do a different DNAIndexation here so we don't have to do it again constantly
    for i=1:2*M %Can pre-allocate DNAIndexation to be more efficient
            if DNA(i) == 'A'    
                DNAIndexation = [DNAIndexation,1];
            elseif DNA(i) == 'T'
                DNAIndexation = [DNAIndexation,2];
            elseif DNA(i) == 'C'
                DNAIndexation = [DNAIndexation,3];
            else
                DNAIndexation = [DNAIndexation,4];
            end
    end
    for k=1:M
        IntroductionRate(k) = exp(-Introduction_Energy_Cost(SHL,DNA_Geometry,DNAIndexation(k:k+N),Geometric_Properties)/2);
    end
    
    Average_Introduction_Rates = zeros(1,10);

    for i=1:floor(M/10)*10

        Average_Introduction_Rates(mod(i-1,10)+1) = Average_Introduction_Rates(mod(i-1,10)+1) + IntroductionRate(i)/floor(M/10);
    end

end
%Chromatin_Barrier = exp(-12); %Make this a better value later


k = Start_Position + (-1)^(Chr_Dir+1); %what position does our simulation start at

Failed_Defect=zeros(2,1);
g=1; %Index for introduction_time
h=1; %Index for the leave_time
r=1; %Index for failed defects
Looping = true;
while Looping
    [Left_Movement,Right_Movement,Introduction_Possible] = LeftRightMovement(Defect_Array);
    
    %Here we can change some stuff to have our introduction rate be
    %sequence dependent, for now we have it constant, implement next time
    TransitionRate = Introduction_Possible*IntroductionRate;
    %!!!!!
    %!!!!!
    %!!!!!
    
    for y=1:n
        if y==1    %Chr_Dir*M gives 0 or M which denotes left or right, mod(k-1)+1 makes our position between 1 and M + our previous from Chr_Dir
            Right_Movement(y) = Right_Movement(y)*exp(-(Barrier_Right(Chr_Dir*M+mod(k-1,M)+1,p)));
            %p+1, p is because we have our chrom. at position 14-n=p
            if Passable_Chromatin %If we can go through we have to add this
                Left_Movement(y) = Defect_Array(y)*Chromatin_Barrier;
            end
            
            TransitionRate = TransitionRate + Right_Movement(y) + Left_Movement(y);
        elseif y==n
            if Chr_Dir == 0%2 directions work differently for our new index
                l = mod(k + sum(Defect_Array(1:n-1))-1,M)+1; %Double check, something maybe wrong with index since it doesn't correspond with theoretical peak
            else
                l = mod(k - sum(Defect_Array(1:n-1))-1,M)+1;
            end %l is used to take into account the different defects can be present and thus we shouldn't
            %use k but the actual position the defect is living in
            Left_Movement(y) = Left_Movement(y)*exp(-(Barrier_Left(Chr_Dir*M+l,14)));
            Right_Movement(y) = Defect_Array(y)*exp(-(3)); %Note this is the leaving one, change later
            TransitionRate = TransitionRate + Left_Movement(y) + Right_Movement(y);
        else
            if Chr_Dir == 0 %2 directions work differently for our new index
                l = mod(k + sum(Defect_Array(1:y-1))-1,M)+1; %Double check, something maybe wrong with index since it doesn't correspond
            else
                l = mod(k - sum(Defect_Array(1:y-1))-1,M)+1;
            end
            Left_Movement(y) = Left_Movement(y)*exp(-(Barrier_Left(Chr_Dir*M+l,p+y-1)));
            Right_Movement(y) = Right_Movement(y)*exp(-(Barrier_Right(Chr_Dir*M+mod(l-1,M)+1,p+y-1)));
            TransitionRate = TransitionRate + Right_Movement(y) + Left_Movement(y);
        end
    end   



    TransitieCoefficient = rand();
    %TransitionRate
    TimeStep = -log(1-TransitieCoefficient)/TransitionRate;
    Time = Time + TimeStep;
    TransitieCoefficient = TransitieCoefficient*TransitionRate;

    if Time > Duration
       Looping = false; %Stop the loop
       k = mod(k - sum(Defect_Array)*(-1)^(Chr_Dir+1)-1,M)+1; %We undo the still present defects, we can make this more accurate later
       %By seeing if they would still complete it or not, we do - for
       %Chr_Dir= 0 and + for Chr_Dir = 1
       
       
       %Removing trailing zero's
       Length_Intro = find(Introduction_Time,1,'last')/2;
       Length_Leave = find(Leave_Time,1,'last')/2; %Finds last non-zero element
       
       Introduction_Time = Introduction_Time(:,1:Length_Intro);
       
       %This next line is optional if we want to remove them already, I use
       %Leave_Time length in Simulation_Nucleosome instead so I don't have
       %to do it here and can still see how many were present
       %Introduction_Time = Introduction_Time(:,1:Length(Introduction_Time)-sum(Defect_Array)); %We have to remove the ones that were introduced
       Leave_Time = Leave_Time(:,1:Length_Leave);
       
    elseif TransitieCoefficient < Introduction_Possible*IntroductionRate
        Defect_Array(1) = Defect_Array(1) + 1; %Here we introduce a twist, doing like this just to test if it works, replace with =1 later
        if Chr_Dir ==0
            k = mod(k-2,M) +1; %k -> k-1
        else
            k = mod(k,M) +1; %k -> k+1
        end
        Introduction_Time(:,g) = [Time,k];
        g = g+1;
        
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
                    if Chr_Dir == 0
                        k = mod(k,M) +1;
                    else
                        k = mod(k-2,M) +1;
                    end
                    r = r+1;
                else
                    Defect_Array(y-1) = Defect_Array(y-1) + 1;
                end
            elseif TransitieCoefficient < Left_Movement(y) + Right_Movement(y)
                %Defect moves to right for this y
                GillespieLoop = false;
                Defect_Array(y) = Defect_Array(y)-1;
                if y==n
                    
                    if Chr_Dir == 0%2 directions work differently for our new index
                        l = mod(k + sum(Defect_Array(1:n-1))-1,M)+1; %Double check, something maybe wrong with index since it doesn't correspond with theoretical peak
                    else
                        l = mod(k - sum(Defect_Array(1:n-1))-1,M)+1;
                    end
                    
                    
                    Leave_Time(:,h) = [Time;l];
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









end

