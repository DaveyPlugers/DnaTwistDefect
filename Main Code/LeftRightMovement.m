function [Left_Movement,Right_Movement,Introduction_Possible] = LeftRightMovement(Defect_Array,Left_Movement,Right_Movement)
%LEFTRIGHTMOVEMENT 
%   Input:
%   Defect_Array = The array that shows where we have defects
%   Left_Movement: Pre-made array, see output
%   Right_Movement: Pre-made array, see output
%   Output:
%   Left_Movement: On which index can defect move to the left
%   Right_Movement: On which index can defect move to the right
%   Introduction_Possible: Can the defect introducer do it's job


Amount_Of_Positions = length(Defect_Array);
%Left_Movement = zeros(1,Amount_Of_Positions); 
%we premake them in our
%Simulation function now to save computation time
%Right_Movement = zeros(1,Amount_Of_Positions);


for z=1:Amount_Of_Positions %~0 =1 and ~1 = 0
    if z==1
        if Defect_Array(z)==1
           Introduction_Possible = false;
        else
            Introduction_Possible = true;
        end
        Left_Movement(z) = 0;
        Right_Movement(z) = ~Defect_Array(z+1)*Defect_Array(z);
    elseif z==Amount_Of_Positions
        Left_Movement(z) = ~Defect_Array(z-1)*Defect_Array(z);
        Right_Movement(z) = Defect_Array(z);
    else
        Left_Movement(z) = ~Defect_Array(z-1)*Defect_Array(z);
        Right_Movement(z) = ~Defect_Array(z+1)*Defect_Array(z);
    end 
end


end

