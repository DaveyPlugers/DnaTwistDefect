%rng(15071544)

%What is the idea:
%Here we generate a n-length array which represents where the defects are
%located 0 is empty, 1 is defect
%Then we make 2 more arrays of the same length, 1 asks if there is a defect
%there and if it can move to the left (1) or none or can't move (0)
%The other does the same for the right

%Once these arrays are made we also check if there is a defect in the left
%most to ask if our defect generator works
%Then we combine all these effects for our Gillespie algorithm, calculating
%the expected time and checking which one actually happens


%Can uncomment some of these if I use basic version
Energy_Landscape = Chromatin_Energy_Landscape(:,1:18);
%Energy_Landscape = [repmat([2,5],1,n)];
%Defect_Array = [0,0,0,0,0,0,0,0,0,0];
%Use a filler landscape for now
n=length(Energy_Landscape(1,:))/2;
Defect_Array = [repmat([0],n)];


Time = 0;
IntroductionRate = 0.1; %Small value so it doesn't happen often

Iterations = 1;
Introduction_Time = [];
Leave_Time = [];

for q=1:Iterations
    [Left_Movement,Right_Movement,Introduction_Possible] = LeftRightMovement(Defect_Array);

    TransitionRate = Introduction_Possible*IntroductionRate;
    for y=1:n
        if y==1
            Right_Movement(y) = Right_Movement(y)*(Energy_Landscape(2*y) - Energy_Landscape(2*y-1));
            TransitionRate = TransitionRate + Right_Movement(y);
        elseif y==n
            Left_Movement(y) = Left_Movement(y)*(Energy_Landscape(2*(y-1)) - Energy_Landscape(2*(y-1)+1));
            Right_Movement(y) = Defect_Array(y)*(3); %Note this is the leaving one, change later
            TransitionRate = TransitionRate + Left_Movement(y) + Right_Movement(y);
        else
            Left_Movement(y) = Left_Movement(y)*(Energy_Landscape(2*(y-1)) - Energy_Landscape(2*(y-1)+1));
            Right_Movement(y) = Right_Movement(y)*(Energy_Landscape(2*y) - Energy_Landscape(2*y-1));
            TransitionRate = TransitionRate + Right_Movement(y) + Left_Movement(y);
        end
    end   



    TransitieCoefficient = rand();
    TimeStep = -log(1-TransitieCoefficient)/TransitionRate;
    Time = Time + TimeStep;
    TransitieCoefficient = TransitieCoefficient*TransitionRate;

    

    if TransitieCoefficient < Introduction_Possible*IntroductionRate
        Defect_Array(1) = Defect_Array(1) + 1; %Here we introduce a twist, doing like this just to test if it works, replace with =1 later
        Introduction_Time = [Introduction_Time,Time];
    else
        TransitieCoefficient = TransitieCoefficient - Introduction_Possible*IntroductionRate;
        GillespieLoop = true;
        y=1;
        while GillespieLoop
            if TransitieCoefficient < Left_Movement(y)
                %Defect moves to left for this y
                GillespieLoop = false;
                Defect_Array(y) = Defect_Array(y)-1;
                Defect_Array(y-1) = Defect_Array(y-1) + 1;
            elseif TransitieCoefficient < Left_Movement(y) + Right_Movement(y)
                %Defect moves to right for this y
                GillespieLoop = false;
                Defect_Array(y) = Defect_Array(y)-1;
                if y==n
                    Leave_Time = [Leave_Time, Time];
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


Introduction_Time;
Leave_Time;
Defect_Array

