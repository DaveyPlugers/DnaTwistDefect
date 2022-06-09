



n = 9;
m = 14-n;
Energy_Landscape = fliplr(Chromatin_Energy_Landscape(:,1:2*n-1));
Other_Side_Landscape = Chromatin_Energy_Landscape_Opp(:,1:2*(m)-1);
Defect_Array = zeros(1,n);
Defect_Array_Opp = zeros(1,m);

Time = 0;
IntroductionRate = 10^(-1); %Small value so it doesn't happen often

Iterations = 10000000;
Introduction_Time = zeros(1,Iterations/10);
Leave_Time = zeros(1,Iterations/10);
k = 10;
p = 10;
g=1;
h=1;
for q=1:Iterations
    [Left_Movement,Right_Movement,Introduction_Possible,Left_Movement_Opp,Right_Movement_Opp] = LeftRightMovement_Two_Sided(Defect_Array,Defect_Array_Opp);
    
    TransitionRate = Introduction_Possible*IntroductionRate;
    for y=1:n %Loop over right one
        if y==1
            Right_Movement(y) = Right_Movement(y)*exp(-(Energy_Landscape(k,2*y) - Energy_Landscape(k,2*y-1)));
            TransitionRate = TransitionRate + Right_Movement(y);
        elseif y==n
            l = mod(k + sum(Defect_Array(1:n-1))-1,10)+1; %Double check, something maybe wrong with index since it doesn't correspond with theoretical peak
            
            Left_Movement(y) = Left_Movement(y)*exp(-(Energy_Landscape(l,2*(y-1)) - Energy_Landscape(l,2*(y-1)+1)));
            Right_Movement(y) = Defect_Array(y)*exp(-(3)); %Note this is the leaving one, change later
            TransitionRate = TransitionRate + Left_Movement(y) + Right_Movement(y);
        else
            l = mod(k + sum(Defect_Array(1:y-1))-1,10)+1; %Double check, something maybe wrong with index since it doesn't correspond
            Left_Movement(y) = Left_Movement(y)*exp(-(Energy_Landscape(l,2*(y-1)) - Energy_Landscape(l,2*(y-1)+1)));
            Right_Movement(y) = Right_Movement(y)*exp(-(Energy_Landscape(l,2*y) - Energy_Landscape(l,2*y-1)));
            TransitionRate = TransitionRate + Right_Movement(y) + Left_Movement(y);
        end
    end   

    for y=1:m %Loop over left one
        if y==1
            l = mod(p + sum(Defect_Array_Opp(2:m))-1,10)+1;  %Maybe we have to swtich sign on the sum?
            Left_Movement_Opp(y) = Defect_Array_Opp(y)*exp(-3);
            Right_Movement_Opp(y) = Right_Movement_Opp(y)*exp(-(Other_Side_Landscape(l,2*y)-Other_Side_Landscape(l,2*y-1)));
            TransitionRate = TransitionRate + Left_Movement_Opp(y) + Right_Movement_Opp(y);
        elseif y==m
            Left_Movement_Opp(y) = Left_Movement_Opp(y)*exp(-(Other_Side_Landscape(p,2*(y-1))-Other_Side_Landscape(p,2*y-1)));
            TransitionRate = TransitionRate + Left_Movement_Opp(y);
        else
            l = mod(p + sum(Defect_Array_Opp(y+1:m))-1,10)+1; %Maybe we have to swtich sign on the sum?
            Right_Movement_Opp(y) = Right_Movement_Opp(y)*exp(-(Other_Side_Landscape(l,2*y)-Other_Side_Landscape(l,2*y-1)));
            Left_Movement_Opp(y) = Left_Movement_Opp(y)*exp(-(Other_Side_Landscape(l,2*(y-1))-Other_Side_Landscape(l,2*y-1)));
            TransitionRate = TransitionRate + Left_Movement_Opp(y) + Right_Movement_Opp(y);
        
        end
    end
    

    TransitieCoefficient = rand();
    %TransitionRate
    TimeStep = -log(1-TransitieCoefficient)/TransitionRate;
    Time = Time + TimeStep;
    TransitieCoefficient = TransitieCoefficient*TransitionRate;

    

    if TransitieCoefficient < Introduction_Possible*IntroductionRate
        Defect_Array(1) = Defect_Array(1) + 1; %Here we introduce a twist, doing like this just to test if it works, replace with =1 later
        Defect_Array_Opp(m) = Defect_Array_Opp(m) + 1;
        Introduction_Time(g) = Time;
        g = g+1;
        k = mod(k-2,10) +1;
        p = mod(p-2,10) + 1; %Don't need p, we can use k but for clarity
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
                    if Time > 10^10
                       Time = Time - 10^10;
                    end
                else
                    Defect_Array(y+1) = Defect_Array(y+1) + 1;
                end
            elseif y==n
                
                TransitieCoefficient = TransitieCoefficient - (Left_Movement(y) + Right_Movement(y));
                z=1;
                while GillespieLoop
                    if TransitieCoefficient < Right_Movement_Opp(z)
                        %Defect moves to right for this z
                        GillespieLoop = false;
                        Defect_Array_Opp(z) = Defect_Array_Opp(z)-1;
                        Defect_Array_Opp(z+1) = Defect_Array_Opp(z+1) + 1;
                    elseif TransitieCoefficient < Right_Movement_Opp(z) + Left_Movement_Opp(z)
                        %Defect moves to left for this z
                        GillespieLoop = false;
                        Defect_Array_Opp(z) = Defect_Array_Opp(z)-1;
                        if not(z==1)
                            Defect_Array_Opp(z-1) = Defect_Array_Opp(z-1) + 1;
                        end
                    else
                       TransitieCoefficient = TransitieCoefficient - (Left_Movement_Opp(z) + Right_Movement_Opp(z));
                       z = z+1;
                    end
                end
            else
                TransitieCoefficient = TransitieCoefficient - (Left_Movement(y) + Right_Movement(y));
                y = y+1;
            end

        end
    end


end


Introduction_Time = Introduction_Time(Introduction_Time ~= 0);
Leave_Time = Leave_Time(Leave_Time ~= 0);
Defect_Array

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



