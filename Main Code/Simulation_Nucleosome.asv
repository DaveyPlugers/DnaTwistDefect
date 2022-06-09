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
%Amount_Of_Simulations = ... (integer) how many simulations will we have,
%can do this to compare changing variables

%Parameters which we can add later if interested
%Assymetric_Remodeller = true/false Different chance for the directions
%Overtwist = true/false we could write the code to allow this too
%Freeze_Out = true/false we could freeze the simulation dynamics with time
%% Parameter initialisation
rng(10071584) %Seed
tic %Timing

Passable_Chromatin = false;
Variable_Introduction = false;
Duration_Twists = 10^8;
Duration_Loops = 10^5;
Chromatin_Barrier = 10^(-12);
Intro_Rate = 10^(-6);

Damp_Factor = 1;
Position_Remodeller = -2;
Loop_Velocity = 1;
Total_Duration = 10^(11);
Amount_Of_Simulations = 1;
Sequence_Difference = [10^(-6);10^(-8);10^(-10);10^(-12)];%Here we can change some parameters for the different simulations
Sequence_Difference2 = [10^7; 10^12;10^14;10^16];
Sequence_Difference3 = [10^7;0.5*10^10;0.5*10^12;0.5*10^14];
DNA = 'GGTGGGCTCCGAAAAATTCGCAGGGCGACCGGCGAAATGCTCAAAAAATCAAAAAATATTCCCTGAAACAAAAGCAACTCTTGCGAAACGGGCGCATTTTTAAAAAATTCGGAGAAAATTTTTGAAATCGTGACGCATCTCTTGCCGTTCGCCAAACAATCCAGAAATTTTATTGTTCAGGCAAAATTCACTAGGTTTTATGGTGAAACGCAAAAAATTCTGACGTTTTCATGAACGTCTTTTAGGATTTTCAGGTTAATGCGGTTGGGCTCCAAAATCTCAAGCAGTCTCGTAGAAATTTCAAAAATTTCGGCATTCTTGGAGAATCACGGGAAGATACTCGCCGGTTT'; 
%% Calculating EnergyLandscape
%Later write some more code here to verify validity

[Barrier_Right,Barrier_Left,Geometric_Properties,Positives,Negatives,Full_Energy_Landscape] = Full_DNA_Energy_Landscape(DNA,Damp_Factor);
"We had " + Negatives + " energy barrier corrections"

%%

%% Simulation



Fractional_Duration_History = zeros(Amount_Of_Simulations,length(DNA));
for z=1:Amount_Of_Simulations
Intro_Rate = Sequence_Difference(z);
Total_Duration = Sequence_Difference2(z);
Duration_Twists = Sequence_Difference3(z);
Length_DNA = length(DNA);
Position = floor(length(DNA)*rand)+1
Looping = true;
Time_Passed = 0;

Size = 100000;

Position_History= zeros(1,Size,'int16');
Duration_History = zeros(1,Size);
Iterations = 0;

WaitingBar = waitbar(0,'Simulating DNA, please wait');
s = clock;
Steps =0;
while Looping
    Direction = floor(2*rand);
    %Direction = 0; %If we want to try one direction %0 is to the left 1 to
    %the right
    Duration = exprnd(Duration_Twists);
    
    
    [Introduction_Time,Leave_Time,k,Temp_Time] = Copy_of_Nucleosome_Transition_Function(Duration,DNA,Position_Remodeller,Position,Direction,Intro_Rate,Variable_Introduction,Passable_Chromatin,Chromatin_Barrier,Barrier_Left,Barrier_Right);
    Steps = Steps + abs(Position-k);
    %Saving duration it was at some location and at what time
    Position_History(Iterations+1) = Position;
    if length(Leave_Time(1,:)) > 0
        Duration_History(Iterations+1) = Leave_Time(1,1);
        
        for i=1:length(Leave_Time(1,:))-1
            Position_History(Iterations+i+1) = Leave_Time(2,i);
            Duration_History(Iterations+i+1) = Leave_Time(1,i+1) - Leave_Time(1,i);
        end
        
        Position_History(Iterations+length(Leave_Time(1,:))+1) = Leave_Time(2,length(Leave_Time(1,:)));
        Duration_History(Iterations+length(Leave_Time(1,:))+1) = Temp_Time - Leave_Time(1,length(Leave_Time(1,:)));
    else
        Duration_History(Iterations+1) = Temp_Time;
    end
    
  
    Position = k;
    Iterations = Iterations + length(Leave_Time(1,:))+1;
    
    Time_Passed = Time_Passed + Duration;
    
    is = etime(clock,s);
    esttime = is * Total_Duration / Time_Passed;
    
    WaitingBar = waitbar(Time_Passed / Total_Duration,WaitingBar,['Simulating DNA, time remaining: ', num2str(esttime/60 - is/60), ' min']);
    
    if Time_Passed > Total_Duration
        Looping = false;
    else
        1+1; %Write loop twist here later
        
        if Time_Passed > Total_Duration
            Looping = false;
        end
    end
    
end   
close(WaitingBar)
WaitingBar = waitbar(1,'Finishing Data Analysis');
%%
close(WaitingBar)


Position_History = Position_History(1:find(Position_History,1,'last'));
Duration_History = Duration_History(1:find(Duration_History,1,'last'));

Fractional_Duration = zeros(1,Length_DNA);

for j=1:length(Position_History)
    if not(Position_History(j) ==0)
        Fractional_Duration(Position_History(j)) = Fractional_Duration(Position_History(j)) + Duration_History(j);
    end
end
Fractional_Duration = Fractional_Duration/Time_Passed;

%Avoid red and green for colour blindness, red in line so no green here
%Uncomment some of these if no package available
%Colours = [[0  0.75 0.75]; [0 0 1]; [0.25 0.25 0.25]; [0.6350, 0.0780, 0.1840]; [0.8500, 0.3250, 0.0980]]; %Cyan, Dark blue, Grey, Bourdeaux, Orange, extra colour if needed: magenta [0.75 0 0.75]
%Colours = summer(Amount_Of_Simulations);  %This is more a sequential way of doing it with 5 colours, should also be good for colour blind? Either uncomment this and the line below or the one above
%Colours(5,:) = [0.75, 0.75, 0]; %Change the last yellow to a darker one
Colours = inferno(Amount_Of_Simulations+1); %+1 since the colormap becomes too white in the end
Lbl_text = ["10^{-4}";"10{-6}"; "10^{-8}";"10^{-10}";"10^{-12}"];

txt = "Introduction rate = " + Lbl_text(z);

plot(Fractional_Duration,'color',Colours(z,:), 'linewidth', 1.2, 'DisplayName', txt)
hold on
Fractional_Duration_History(z,:) = Fractional_Duration;
end
 
title("Occupation distribution remodeller for different passable remodeller rates")
xlabel('Begin Index of 147 length DNA string')
ylabel('Fractional probability')
Uniform = 1/length(DNA);
line([0 350],[Uniform Uniform],'Color','red','LineStyle','--','DisplayName','Uniform probability' )
legend show
grid minor

toc
Steps