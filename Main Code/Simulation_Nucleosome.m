%This is the main script, we will be running the simulation through this
%main code. Let me give a quick description of what the code does. For
%further details it is recommended to look into the function descriptions

%% Code explanation:

%We start out by getting the full energy landscape for twists, This is done
%for every single position in the given sequence and both directions. We do
%it for overtwist now, further development could be to do undertwist or the
%average energy landscape by having multiple sequences or adding the damp
%factor.

%After this we need to prepare the loop defect values, I will do this later
%I don't have this code yet so I will just use some filler values for now


%Then a random position gets chosen followed by a random direction for the
%twists, sub-simulation is run which does twist propagation, this will
%return the position w.r.t. time and once it is over we save these values
%here. Completing this sub simulation a new sub simulation will commence,
%now doing loop propagation, this will also return the positions w.r.t.
%time to the main simulation. This will be completed and a new random
%direction is chosen a sub simulation of twist defects starts again and we
%repeat this progress until a sub-simulation causes us to reach the final
%time (once the sub-simulation is over, not while running).

%When the simulations have been completed we end up with an array of the
%position w.r.t. time and a note in what mode our mode the simulation was
%running. We use this to calculate correlation functions, we can check if
%the mode has a big effect on the correlation time too. We also calculate
%an approximate probability distribution function of where we are more
%likely to find the nucleosome. This is done by summing the time duration
%during which they were in that position for all the positions and dividing
%by the total duration, This is only valid however for sufficiently long
%simulations. BUT STILL CHECK EPS(x) WE MIGHT HAVE BIG ERRORS 
%(eps(10^20) = 16000!!!



%Parameters which we can adjust:
%Passable_Chromatin = true/false to adjust the simulation type
%Variable_Introduction = true/false to adjust the simulation type
%DNA = ... To change the sequence
%Duration_Twists = ... To change time of Twist propagation (Exp distr.)
%Duration_Loops = ... To change time of loop propagation
%Chromatin_Barrier = ... what rate do we cross the remodeller if allowed
%Intro_Rate = ... (if we don't vary it)
%Damp_Factor = ... (Between 0 to 1) How big effect displacement energy
%Position_Remodeller = ... we can move it a bit but very likely SHL=-2
%Loop_Velocity = ... Change timescale on which the loops move w.r.t. twists
%Total_Duration = ... Time we want to run our simulation for
%Note that high timesteps will cause errors, USE EPS(x) TO CHECK!!!

%Parameters which we can add later if interested
%Assymetric_Remodeller = true/false Different chance for the directions
%Overtwist = true/false we could write the code to allow this too
%Freeze_Out = true/false we could freeze the simulation dynamics with time
%% Parameter initialisation

Passable_Chromatin = false;
Variable_Introduction = false;
Duration_Twists = 10^8;
Duration_Loops = 10^5;
Chromatin_Barrier = 10^(-12);
Intro_Rate = 10^(-6);
Damp_Factor = 0.9;
Position_Remodeller = -2;
Loop_Velocity = 1;
Total_Duration = 10^(9);


DNA = 'GGTGGGCTCCGAAAAATTCGCAGGGCGACCGGCGAAATGCTCAAAAAATCAAAAAATATTCCCTGAAACAAAAGCAACTCTTGCGAAACGGGCGCATTTTTAAAAAATTCGGAGAAAATTTTTGAAATCGTGACGCATCTCTTGCCGTTCGCCAAACAATCCAGAAATTTTATTGTTCAGGCAAAATTCACTAGGTTTTATGGTGAAACGCAAAAAATTCTGACGTTTTCATGAACGTCTTTTAGGATTTTCAGGTTAATGCGGTTGGGCTCCAAAATCTCAAGCAGTCTCGTAGAAATTTCAAAAATTTCGGCATTCTTGGAGAATCACGGGAAGATACTCGCCGGTTT'; 
%% Calculating EnergyLandscape
%Later write some more code here to verify validity

[Barrier_Right,Barrier_Left,Geometric_Properties,Positives,Negatives,Full_Energy_Landscape] = Full_DNA_Energy_Landscape(DNA,Damp_Factor);


%%

%% Simulation
Length_DNA = length(DNA);
Position = floor(length(DNA)*rand)+1;
Looping = true;
Time_Passed = 0;

Size = 100000;

Position_History= zeros(1,Size,'int16');
Duration_History = zeros(1,Size);
Iterations = 0;

while Looping
    %Direction = floor(2*rand);
    Direction = 1;
    Duration = exprnd(Duration_Twists);
    
    
    [Introduction_Time,Leave_Time,k,Temp_Time] = Chromatin_Transition_Function(Duration,DNA,Position_Remodeller,Position,Direction,Intro_Rate,Variable_Introduction,Passable_Chromatin,Chromatin_Barrier,Barrier_Left,Barrier_Right);
    
    %Saving duration it was at some location and at what time
    Position_History(Iterations+1) = Position;
    if length(Leave_Time(1,:)) > 1
        Duration_History(Iterations+1) = Leave_Time(1,1);
    else
        Duration_History(Iterations+1) = Temp_Time;
    end
    if length(Leave_Time(1,:)) > 1
        for i=1:length(Leave_Time)-1
            Position_History(Iterations+i+1) = Leave_Time(2,i);
            Duration_History(Iterations+i+1) = Leave_Time(1,i+1) - Leave_Time(1,i);
        end
        
    end
    if length(Leave_Time(1,:)) > 1 
        Position_History(Iterations+length(Leave_Time(1,:))+1) = Leave_Time(2,length(Leave_Time(1,:)));
        Duration_History(Iterations+length(Leave_Time(1,:))+1) = Temp_Time - Leave_Time(1,length(Leave_Time(1,:)));
    end
    Position = k;
    Iterations = Iterations + length(Leave_Time(1,:))+1;
    
    Time_Passed = Time_Passed + Duration
    if Time_Passed > Total_Duration
        Looping = false;
    else
        1+1; %Write loop twist here later
        
        if Time_Passed > Total_Duration
            Looping = false;
        end
    end
    
end    

%%



Position_History = Position_History(1:find(Position_History,1,'last'));
Duration_History = Duration_History(1:find(Duration_History,1,'last'));

Fractional_Duration = zeros(1,Length_DNA);

for j=1:length(Position_History)
    if not(Position_History(j) ==0)
        Fractional_Duration(Position_History(j)) = Fractional_Duration(Position_History(j)) + Duration_History(j);
    end
end
Fractional_Duration = Fractional_Duration/Time_Passed;
plot(Fractional_Duration)