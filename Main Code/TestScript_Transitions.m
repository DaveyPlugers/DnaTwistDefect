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
n = 9;
Energy_Landscape = fliplr(Full_Energy_Landscape(1:147,1:2*n-1));

%Use a filler landscape for now
%Energy_Landscape = [repmat([2,5],1,n)];
%Defect_Array = [0,0,0,0,0,0,0,0,0,0];


Defect_Array = zeros(1,n);


Time = 0;
IntroductionRate = 10^(-7); %Small value so it doesn't happen often

Iterations = 10000000;
Introduction_Time = zeros(1,floor(Iterations/n^2));
Leave_Time = zeros(1,floor(Iterations/n^2));
k = 10;
g=1;
h=1;
for q=1:Iterations
    [Left_Movement,Right_Movement,Introduction_Possible] = LeftRightMovement(Defect_Array);
    
    TransitionRate = Introduction_Possible*IntroductionRate;
    for y=1:n
        if y==1
            Right_Movement(y) = Right_Movement(y)*exp(-(Energy_Landscape(k,2*y) - Energy_Landscape(k,2*y-1)));
            TransitionRate = TransitionRate + Right_Movement(y);
        elseif y==n
            l = mod(k - sum(Defect_Array(1:n-1))-1,10)+1; %Double check, something maybe wrong with index since it doesn't correspond with theoretical peak
            Left_Movement(y) = Left_Movement(y)*exp(-(Energy_Landscape(l,2*(y-1)) - Energy_Landscape(l,2*(y-1)+1)));
            Right_Movement(y) = Defect_Array(y)*exp(-(3)); %Note this is the leaving one, change later
            TransitionRate = TransitionRate + Left_Movement(y) + Right_Movement(y);
        else
            
            l = mod(k - sum(Defect_Array(1:y-1))-1,10)+1; %Double check, something maybe wrong with index since it doesn't correspond
            Left_Movement(y) = Left_Movement(y)*exp(-(Energy_Landscape(l,2*(y-1)) - Energy_Landscape(l,2*(y-1)+1)));
            Right_Movement(y) = Right_Movement(y)*exp(-(Energy_Landscape(l,2*y) - Energy_Landscape(l,2*y-1)));
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
        k = mod(k,10) +1;
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


Introduction_Time = Introduction_Time(Introduction_Time ~= 0);;
Leave_Time = Leave_Time(Leave_Time ~= 0);
Defect_Array






%Here statistical analysis of the time values

TenSteps = floor(length(Leave_Time)/10);

Average_Lifetime = zeros(1,10);
Lifetime = zeros(TenSteps,10);
for i=1:TenSteps
    for j=1:10
        z = 10*(i-1) + j;
        Average_Lifetime(j) = Average_Lifetime(j) + (Leave_Time(z) - Introduction_Time(z))/TenSteps;
        Lifetime(i,j) = Leave_Time(z) - Introduction_Time(z);
    end
end
x = 1:10;
errorbar(x,Average_Lifetime,std(Lifetime)/sqrt(TenSteps))
set(gca,'YScale','log')

%Calculating theoretical value
TTR = zeros(11,2,n); %TTR = theoretical transition rate
for y=1:10
    TTR(y,1,1) = exp(-(Energy_Landscape(y,2)-Energy_Landscape(y,1)));
    for z=2:n-1
        TTR(y,1,z) = exp(-(Energy_Landscape(y,2*z)-Energy_Landscape(y,2*z-1)));
        TTR(y,2,z) = exp(-(Energy_Landscape(y,2*(z-1))-Energy_Landscape(y,2*z-1)));
    end
    z=9;
    TTR(y,2,9) = exp(-(Energy_Landscape(y,2*(z-1))-Energy_Landscape(y,2*z-1)));
    TTR(y,1,9) = exp(-3); %We just use an energy barrier of 3 for all of these, change later
end
TTR(11,:,:) = exp(-3);
Theoretical_Lifetimes = zeros(11,n);

for y=1:11
    Theoretical_Lifetimes(y,1) = 1/TTR(y,1,1);
    for z = 2:n
        %Theoretical_Lifetimes(y,z) = Theoretical_Lifetimes(y,z-1) + 1/(TTR(y,1,z) + TTR(y,2,z))*(1+ TTR(y,2,z)/TTR(y,1,z) + TTR(y,2,z)*(TTR(y,1,z)+TTR(y,2,z))/TTR(y,1,z)*(1/TTR(y,2,z) + Theoretical_Lifetimes(y,z-1)));
        Theoretical_Lifetimes(y,z) = TTR(y,1,z)/(TTR(y,1,z)+TTR(y,2,z))^2*(1+TTR(y,2,z)/TTR(y,1,z) + (TTR(y,2,z)*(TTR(y,1,z)+TTR(y,2,z)))/(TTR(y,1,z))^2* ( 1+(TTR(y,1,z)+TTR(y,2,z))*Theoretical_Lifetimes(y,z-1)));
    end
end
hold on
plot(sum(Theoretical_Lifetimes(1:10,:),2))

legend('Single Sided Simulation Lifetime','Theoretical single particle Lifetime')



if floor(length(Leave_Time(1:end-2))/10) == TenSteps
    Time_Delays_Double = zeros(1,10);
    Time_Delays_Triple = zeros(1,10);

    for i=1:TenSteps
        for j=1:10
            Time_Delays_Double(j) = (Leave_Time((i-1)*10+j+1) - Leave_Time((i-1)*10+j))/TenSteps;
            Time_Delays_Triple(j) = (Leave_Time((i-1)*10+j+2) - Leave_Time((i-1)*10+j))/TenSteps;
        end
    end
    
else
    Time_Delays_Double = zeros(1,10);
    Time_Delays_Triple = zeros(1,10);

    for i=1:TenSteps-1
        for j=1:10
            Time_Delays_Double(j) = (Leave_Time((i-1)*10+j+1) - Leave_Time((i-1)*10+j))/(TenSteps-1);
            Time_Delays_Triple(j) = (Leave_Time((i-1)*10+j+2) - Leave_Time((i-1)*10+j))/(TenSteps-1);
        end
    end
end



figure
plot(x,Time_Delays_Double)

hold on
plot(x,Time_Delays_Triple)
set(gca,'YScale','log')