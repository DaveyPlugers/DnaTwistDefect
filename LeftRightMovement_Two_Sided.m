function [Left_Movement,Right_Movement,Introduction_Possible,Left_Movement_Opp,Right_Movement_Opp] = LeftRightMovement_Two_Sided(Defect_Array,Defect_Array_Opp)
%LEFTRIGHTMOVEMENT 
%   Input:
%   Defect_Array = The right-side array that shows where we have defects
%   Defect_Array_Opp = The left-side array that shows where we have defects
%   Output:
%   Left_Movement: On which index can defect move to the left for the
%   right-side array
%   Right_Movement: On which index can defect move to the right
%   right-side array
%   Left_Movement_opp: On which index can defect move to the left for the
%   left-side array
%   Right_Movement_opp: On which index can defect move to the right for the
%   left-side array
%   Introduction_Possible: Can the defect introducer do it's job


Amount_Of_Positions = length(Defect_Array);
Amount_Of_Positions_Opp = length(Defect_Array_Opp);
Left_Movement = zeros(1,Amount_Of_Positions);
Right_Movement = zeros(1,Amount_Of_Positions);
Left_Movement_Opp = zeros(1,Amount_Of_Positions_Opp);
Right_Movement_Opp = zeros(1,Amount_Of_Positions_Opp);


for z=1:Amount_Of_Positions %~0 =1 and ~1 = 0, 
    if z==1 %There already is defect in one of them so we can't move them
        if or(Defect_Array(z)==1, Defect_Array_Opp(Amount_Of_Positions_Opp)==1)
           Introduction_Possible = false;
        else
            Introduction_Possible = true;
        end
        Left_Movement(z) = 0;
        Right_Movement(z) = ~Defect_Array(z+1)*Defect_Array(z);
    elseif z==Amount_Of_Positions
        Left_Movement(z) = ~Defect_Array(z-1)*Defect_Array(z);
        Right_Movement(z) = 0;
    else
        Left_Movement(z) = ~Defect_Array(z-1)*Defect_Array(z);
        Right_Movement(z) = ~Defect_Array(z+1)*Defect_Array(z);
    end 
end


for z=1:Amount_Of_Positions_Opp
   if z==1
       Left_Movement_Opp(z) = Defect_Array_Opp(z);
       Right_Movement_Opp(z) = ~Defect_Array_Opp(z+1)*Defect_Array_Opp(z);
   elseif z==Amount_Of_Positions_Opp
       Right_Movement_Opp(z) =0;
       Left_Movement_Opp(z) = ~Defect_Array_Opp(z-1)*Defect_Array_Opp(z);
   else
       Left_Movement_Opp(z) = ~Defect_Array_Opp(z-1)*Defect_Array_Opp(z);
       Right_Movement_Opp(z) = ~Defect_Array_Opp(z+1)*Defect_Array_Opp(z);  
   end
end


end


