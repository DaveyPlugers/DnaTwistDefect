Intro_Rate = 10^(-6);
Initial_Lifetimes = sum(Theoretical_Lifetimes(1:10,:),2);
n=9;   %Some places I still use 9 instead of n, correct this later if you want to use this code

Transition_Rates_Two_State = zeros(10,2,n);
%10 is the amount of periodicity in our landscape
%1 is left (minus) 2 is right (plus) movements, note this is opposite as
%what TTR does
%n is for the different positions it can be in

Transition_Rates_Two_State(1:10,1,:) = TTR(1:10,2,:);
Transition_Rates_Two_State(1:10,2,:) = TTR(1:10,1,:);

Transition_Rates_Two_State(1:10,1,1) = NaN; %We should never use this value, we use NaN to make debugging easier


Occupancies=zeros(10,2,9); %If there are 2 particles and we know it's left or right, we calculate the probability we have one to the left of ours or the the right of ours at these positions
for k=1:10
for i=1:n
    if i==1
        Occupancies(k,2,i) = TwoStateDistribution(Transition_Rates_Two_State(k,:,:),2,n-1);
        Occupancies(k,1,i) = NaN;
    elseif i==9
        Occupancies(k,1,i) = TwoStateDistribution(Transition_Rates_Two_State(k,:,:),1,n-1);
        Occupancies(k,2,i) = NaN;
    else
        Occupancies(k,1,i) = TwoStateDistribution(Transition_Rates_Two_State(k,:,:),1,i-1);
        Occupancies(k,2,i) = TwoStateDistribution(Transition_Rates_Two_State(k,:,:),2,n-i);
    end
end
end

%We manually change an occupancy since otherwise we get an infinite lifetime in our TwoStateSubLifetimes calculator for right movement
Occupancies(1:10,2,8) = 1-0.5*(1-Occupancies(1:10,2,7));
%In reality occupancy is 1 since we say there is already one present, however our system can evolve and the particle can leave. Our analytical
%solution isn't that smart though and lives in the past. We approximate it by using the previous occupancy value so it is still position dependent

%We should correct the Occupancies later in a better way


SubLifeTimes = zeros(10,2,n);
Occupancies(:,2,:) = Occupancies(:,2,:)/2;


for k=1:10 
    for i=1:n
        if i==1
            SubLifeTimes(k,1,i) = Theoretical_Lifetimes(k,i);
            SubLifeTimes(k,2,i) = TwoStateSubLifetimes(2,Occupancies(k,2,i),0,Transition_Rates_Two_State(k,2,i),0);
        elseif i==9
            SubLifeTimes(k,1,i) = TwoStateSubLifetimes(1,Occupancies(k,1,i),Transition_Rates_Two_State(k,1,i),Transition_Rates_Two_State(k,2,i),SubLifeTimes(k,1,i-1));
            SubLifeTimes(k,2,i) = NaN;
        else
            SubLifeTimes(k,1,i) = TwoStateSubLifetimes(1,Occupancies(k,1,i),Transition_Rates_Two_State(k,1,i),Transition_Rates_Two_State(k,2,i),SubLifeTimes(k,1,i-1));
            SubLifeTimes(k,2,i) = TwoStateSubLifetimes(2,Occupancies(k,2,i),Transition_Rates_Two_State(k,1,i),Transition_Rates_Two_State(k,2,i),SubLifeTimes(k,2,i-1));
        end
    end
end





Transition_Rates_Two_State(1:10,1,1) = 0; %We set it to 0 now because of the TwoStateLifetimeArray(k,2,1) calculation
ETR = zeros(10,2,10);  %Expanded Transition Rates
ETR(1:10,:,1:9) = Transition_Rates_Two_State;
%we expand 9 to 10 for when we calculate having another particle to the right when we have i=8, we can do the same
%thing for i=2 and particle to the left but this shifts all indices so we don't


TwoStateLifetimeArray = zeros(10,2,n);
for k=1:10
    %i=1
    TwoStateLifetimeArray(k,1,1) = NaN;   %mod(k-2,10)+1 to go from 1 to 10, 10 to 9 etc. previous position
    TwoStateLifetimeArray(k,2,1) = TwoStateTheoreticalTime(2,ETR(k,2,1),ETR(k,1,1),ETR(mod(k-2,10)+1,2,2),ETR(mod(k-2,10)+1,1,2),ETR(mod(k-2,10)+1,2,3),ETR(mod(k-2,10)+1,1,3),NaN,SubLifeTimes(k,2,1),0);
    %i=2 %mod(k,10) + 1 to go from 10 to 1 9 to 10 etc.
    TwoStateLifetimeArray(k,1,2) = TwoStateTheoreticalTime(1,ETR(k,2,2),ETR(k,1,2),ETR(mod(k,10)+1,2,1),ETR(mod(k,10)+1,1,1),0,0,0,SubLifeTimes(k,1,2),NaN); %We insert 0 for expprevious since we cannot get here
    TwoStateLifetimeArray(k,2,2) = TwoStateTheoreticalTime(2,ETR(k,2,2),ETR(k,1,2),ETR(mod(k-2,10)+1,2,3),ETR(mod(k-2,10)+1,1,3),ETR(mod(k-2,10)+1,2,4),ETR(mod(k-2,10)+1,1,4),NaN,SubLifeTimes(k,2,2),SubLifeTimes(k,2,1));

    for i=3:n-1
        TwoStateLifetimeArray(k,1,i) = TwoStateTheoreticalTime(1,ETR(k,2,i),ETR(k,1,i),ETR(mod(k,10)+1,2,i-1),ETR(mod(k,10)+1,1,i-1),ETR(mod(k,10)+1,2,i-2),ETR(mod(k,10)+1,1,i-2),TwoStateLifetimeArray(k,1,i-1),SubLifeTimes(k,1,i),NaN);
        TwoStateLifetimeArray(k,2,i) = TwoStateTheoreticalTime(2,ETR(k,2,i),ETR(k,1,i),ETR(mod(k-2,10)+1,2,i+1),ETR(mod(k-2,10)+1,1,i+1),ETR(mod(k-2,10)+1,2,i+2),ETR(mod(k-2,10)+1,1,i+2),NaN,SubLifeTimes(k,2,i),SubLifeTimes(k,2,i-1));
    end
    %i=n
    TwoStateLifetimeArray(k,1,n) = TwoStateTheoreticalTime(1,ETR(k,2,n),ETR(k,1,n),ETR(mod(k,10)+1,2,n-1),ETR(mod(k,10)+1,1,n-1),ETR(mod(k,10)+1,2,n-2),ETR(mod(k,10)+1,1,n-2),TwoStateLifetimeArray(k,1,n-1),SubLifeTimes(k,1,n),NaN);
    TwoStateLifetimeArray(k,2,n) = NaN;
end    





Probability_Left = zeros(10,n);
Probability_Right = zeros(10,n);

for k=1:10
    %i=1
    Probability_Left(k,1) = 0;
    Probability_Right(k,1) = (1-((1-exp(-Intro_Rate*Initial_Lifetimes(mod(k-2,10)+1)))/(Intro_Rate*Initial_Lifetimes(mod(k-2,10)+1))))*Occupancies(k,2,1);
    for i=2:n-1
        Probability_Left(k,i) = (1-(1-exp(-Intro_Rate*Initial_Lifetimes(k)))/(Intro_Rate*Initial_Lifetimes(k)))*Occupancies(k,1,i);
        %Probability_Right(k,i) = (Initial_Lifetimes(mod(k-2,10)+1)/Initial_Lifetimes(k)*(exp(-1/(Intro_Rate*Initial_Lifetimes(mod(k-2,10)+1))) - exp(-(1/Intro_Rate + Initial_Lifetimes(k))/Initial_Lifetimes(mod(k-2,10)+1))))*Occupancies(k,2,i);
        Probability_Right(k,i) = (1-((1-exp(-Intro_Rate*Initial_Lifetimes(mod(k-2,10)+1)))/(Intro_Rate*Initial_Lifetimes(mod(k-2,10)+1))))*Occupancies(k,2,i);
    end
    %i=n
    Probability_Left(k,n) = (1-(1-exp(-Intro_Rate*Initial_Lifetimes(k)))/(Intro_Rate*Initial_Lifetimes(k)))*Occupancies(k,1,n);
    Probability_Right(k,n) = 0;
end

Probability_Right(:,7:9) = Probability_Right(:,7:9)*0.5;

Step_LifeTimes = zeros(10,9);

for k=1:10
   i=1;
   Step_LifeTimes(k,1) = Probability_Right(k,i)*TwoStateLifetimeArray(k,2,i) + (1-Probability_Left(k,i)-Probability_Right(k,i))*Theoretical_Lifetimes(k,i);

   for i=2:8 
        Step_LifeTimes(k,i) = Probability_Left(k,i)*TwoStateLifetimeArray(k,1,i) + Probability_Right(k,i)*TwoStateLifetimeArray(k,2,i) + (1-Probability_Left(k,i)-Probability_Right(k,i))*Theoretical_Lifetimes(k,i);
       
   end  
   i=9;
   Step_LifeTimes(k,i) = Probability_Left(k,i)*TwoStateLifetimeArray(k,1,i) + (1-Probability_Left(k,i)-Probability_Right(k,i))*Theoretical_Lifetimes(k,i);

end

DifferenceLeft = zeros(10,9);
for k=1:10
    for j=1:9
        DifferenceLeft(k,j) = Theoretical_Lifetimes(k,j) - TwoStateLifetimeArray(k,1,j);
    end
end
DifferenceRight = zeros(10,9);
for k=1:10
    for j=1:9
        DifferenceRight(k,j) = Theoretical_Lifetimes(k,j) - TwoStateLifetimeArray(k,2,j);
    end
end

errorbar(x,Average_Lifetime,std(Lifetime)/sqrt(TenSteps))
set(gca,'YScale','log')
hold on
plot(sum(Theoretical_Lifetimes(1:10,:),2))
plot(sum(Step_LifeTimes(1:10,:),2))
legend('Single Sided Simulation Lifetime','Theoretical single particle Lifetime','Theoretical Two-side')