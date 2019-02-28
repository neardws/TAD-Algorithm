function output = isTaskDone(taskFinish)
%ISTASKDONE Summary of this function goes here
%   Detailed explanation goes here
    for i = 1 : numel(taskFinish)
        if taskFinish(i) == 0
            output = 0;
            break;
        end
    end
    output = 1;
end

