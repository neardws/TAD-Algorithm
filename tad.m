% import frome TXT

VehicleSum = 110;
TimeSum = 300;

origindata = importdata("data.txt");
[m,n] = size(origindata);
VehicleTime = zeros(VehicleSum,TimeSum);
VehicleTrace = zeros(VehicleSum*TimeSum,2);
for i = 1 : m
    trace = origindata(i,:);
    vehicleID = trace(1);
    time = trace(2);
    x = trace(3);
    y = trace(4);
    if vehicleID ~= 0 && time ~= 0
        VehicleTime(vehicleID,time) = 1;
        VehicleTrace((vehicleID-1)*TimeSum+time,1) = x;
        VehicleTrace((vehicleID-1)*TimeSum+time,2) = y;
    end
end

%{
    Init some date
%}

% Vehicle
TaskEtimes = 1e9;
MobileFogNum = 20;
VehicleNum = VehicleSum - MobileFogNum;
VehicleTask = randi([1,1],VehicleNum,1);

% Task
TaskNum = sum(VehicleTask(:));
TaskVehicle = zeros(TaskNum,VehicleNum);
TaskSize = randi([10,25],TaskNum,1);
TaskCpu = randi([10,25],TaskNum,1)*TaskEtimes;
TaskEndTime = randi([10,300],TaskNum,1);

% Fog
FogEtimes = 1e-10;

FixedFogNum = 4;
FixedFogLocal = [500 500; 500 1000; 1000 500; 1000 1000];
FixedFogSize = randi([80,150],FixedFogNum,1);
FixedFogCompu = randi([5,10],FixedFogNum,1)*FogEtimes;
FixedFogTrans = randi([20,40],FixedFogNum,1);

MobileFogID = randperm(110, MobileFogNum);
MobileFogSize = randi([50,100],MobileFogNum,1);
MobileFogCompu = randi([5,10],MobileFogNum,1)*FogEtimes;
MobileFogTrans = randi([10,20],MobileFogNum,1);

FogNum = MobileFogNum + FixedFogNum;
FogSize = cat(1,MobileFogSize,FixedFogSize);
FogCompu = cat(1,MobileFogCompu,FixedFogCompu);
FogTrans = cat(1,MobileFogTrans,FixedFogTrans);

% Task & Vehicle
for i = 1 : TaskNum
    for n = 1 : VehicleNum
        if i <= sum(VehicleTask(1:n,:),1)
            TaskVehicle(i,n)=1;
            break;
        end
    end
end

MF = []; % Mobile Fog ID
V = []; % Vehicles ID
for id = 1 : VehicleSum
    if ismember(id,MobileFogID) 
        MF = [MF,id];
    else
        V = [V,id];
    end
end


%{ 
    NowTime Update Every Arragement
%}
% nowTime = 1;

% % Task & Fog
% TaskFog = zeros(TaskNum,FogNum);
% 
% for i = 1: TaskNum
%     for n = 1 : VehicleNum
%         % Task i is in Vehicle n
%         if TaskVehicle(i,n) == 1
%             vehicleId = V(n);
%             if VehicleTime(vehicleId, nowTime)
%                 vehicleLoc = VehicleTrace(((vehicleId-1)*300+nowTime),:);
%                 if MobileFogNum ~= 0
%                     for m = 1 : MobileFogNum
%                         fogId = MF(m);
%                         if VehicleTime(fogId, nowTime)
%                             mobileFogLoc = VehicleTrace((fogId-1)*300+nowTime,:);
%                             if isIn(vehicleLoc, mobileFogLoc)
%                                 TaskFog(i,m)=1;
%                             end
%                         end
%                     end
%                     for m = 1 : FixedFogNum
%                         if isIn(vehicleLoc, FixedFogLocal(m,:))
%                             TaskFog(i,MobileFogNum+m) = 1;
%                         end
%                     end
%                 end
%             end
%         end
%     end
% end
% 
% % Profits of task in fog 
% TaskFogProfit = zeros(TaskNum,FogNum);
% 
% if MobileFogNum ~= 0
%     for m = 1 : MobileFogNum
%         taskInFog = TaskFog(:,m);
%         sumSize = sum(TaskSize.*taskInFog);
%         sumCom = sum(TaskCpu.*taskInFog);
%         for i = 1 : TaskNum
%             if TaskFog(i,m) == 1
%                 TaskFogProfit(i,m) = TaskSize(i)/sumSize + TaskCpu(i)/sumCom;
%             end
%         end
%     end
% end
% 
% for m = MobileFogNum +1 : FogNum
%     taskInFog = TaskFog(:,m);
%     sumSize = sum(TaskSize.*taskInFog);
%     sumCom = sum(TaskCpu.*taskInFog);
%     for i = 1 : TaskNum
%         if TaskFog(i,m) == 1
%             TaskFogProfit(i,m) = TaskSize(i)/sumSize + TaskCpu(i)/sumCom;
%         end
%     end
% end
% 
% % Mini time of task
% TaskFogMiniTime = zeros(TaskNum, FogNum);
% 
% for i = 1 : TaskNum
%     for m = 1 : FogNum
%         if TaskFog(i,m) == 1
%             TaskFogMiniTime(i,m) = TaskEndTime(i);
%         end
%     end
% end
% 
% for i = 1 : TaskNum
%     % find vehicleID
%     for n = 1 : VehicleNum
%         if TaskVehicle(i,n) == 1
%             vehicleId = V(n);
%         end
%     end
%     if MobileFogNum ~= 0
%         for m = 1 : MobileFogNum
%             if TaskFog(i,m) == 1
%                 for t = nowTime : TimeSum
%                     vehicleLoc = VehicleTrace(((vehicleId-1)*TimeSum + t),:);
%                     mobileFogLoc = VehicleTrace(((MF(m)-1)*TimeSum + t),:);
%                     if isIn(vehicleLoc,mobileFogLoc) == 0
%                         leaveTime = t;
%                         TaskFogMiniTime(i,m) = min(leaveTime,taskFogMiniTime(i));
%                         break;
%                     end
%                 end
%             end
%         end
%     end
%     for m = MobileFogNum + 1 : FogNum
%         if TaskFog(i,m) == 1
%             for t = nowTime : TimeSum
%                 vehicleLoc = VehicleTrace(((vehicleId-1)*TimeSum + t),:);
%                 fixedFogLoc = FixedFogLocal(m-MobileFogNum,:);
%                 if isIn(vehicleLoc,fixedFogLoc) == 0
%                     leaveTime = t;
%                     TaskFogMiniTime(i,m) = min(leaveTime,taskEndTime(i));
%                     break;
%                 end
%             end
%         end
%     end
% end

% taskSumInFog = zeros(FogNum,1);
% for i = 1 : FogNum
%     taskInFog = TaskFog(:,i);
%     taskSumInFog(i) = sum(taskInFog);
% end
% 
% maxTaskSumInFog = max(taskSumInFog);

%{
    Write File
%}
% fileName = 'datas.txt';
% writeNum(fileName,taskNum);
% writeNum(fileName,fogNum);
% writeNum(fileName,maxTaskSumInFog);
% writeMatrix(fileName,taskSize');
% writeMatrix(fileName,taskCpu');
% writeMatrix(fileName,fogSize');
% writeMatrix(fileName,fogCompu');
% writeMatrix(fileName,fogTrans');
% writeMatrix(fileName,taskSumInFog');
% writeMatrix(fileName,taskFog');
% writeMatrix(fileName,taskFogProfit');
% writeMatrix(fileName,taskFogMiniTime');


%{
    TAD-Algorithm
%}

TaskFinish = zeros(TaskNum,1);
TaskChoosed = zeros(TaskNum,FogNum);
startTime = 1;
lastComplete = -1;
profitsSum = zeros(100,1);
endTime = zeros(100,1);
arragementTime = 1;
taskFogFinishTime = zeros(TaskNum, FogNum);

while isTaskDone(TaskFinish)
    %{ 
        NowTime Update Every Arragement
    %} 
    % Task & Fog
    nowTime = startTime;
    TaskFog = zeros(TaskNum,FogNum);
    
    for i = 1: TaskNum
        for n = 1 : VehicleNum
            % Task i is in Vehicle n
            if TaskVehicle(i,n) == 1
                vehicleId = V(n);
                if VehicleTime(vehicleId, nowTime) == 1
                    vehicleLoc = VehicleTrace(((vehicleId-1)*300+nowTime),:);
                    if MobileFogNum ~= 0
                        for m = 1 : MobileFogNum
                            fogId = MF(m);
                            if VehicleTime(fogId, nowTime) == 1
                                mobileFogLoc = VehicleTrace((fogId-1)*300+nowTime,:);
                                if isIn(vehicleLoc, mobileFogLoc)
                                    TaskFog(i,m)=1;
                                end
                            end
                        end
                        for m = 1 : FixedFogNum
                            if isIn(vehicleLoc, FixedFogLocal(m,:))
                                TaskFog(i,MobileFogNum+m) = 1;
                            end
                        end
                    end
                end
            end
        end
    end

    % Profits of task in fog 
    TaskFogProfit = zeros(TaskNum,FogNum);

    if MobileFogNum ~= 0
        for m = 1 : MobileFogNum
            taskInFog = TaskFog(:,m);
            sumSize = sum(TaskSize.*taskInFog);
            sumCom = sum(TaskCpu.*taskInFog);
            for i = 1 : TaskNum
                if TaskFog(i,m) == 1
                    TaskFogProfit(i,m) = TaskSize(i)/sumSize + TaskCpu(i)/sumCom;
                end
            end
        end
    end

    for m = MobileFogNum +1 : FogNum
        taskInFog = TaskFog(:,m);
        sumSize = sum(TaskSize.*taskInFog);
        sumCom = sum(TaskCpu.*taskInFog);
        for i = 1 : TaskNum
            if TaskFog(i,m) == 1
                TaskFogProfit(i,m) = TaskSize(i)/sumSize + TaskCpu(i)/sumCom;
            end
        end
    end

    % Mini time of task
    TaskFogMiniTime = zeros(TaskNum, FogNum);

    for i = 1 : TaskNum
        for m = 1 : FogNum
            if TaskFog(i,m) == 1
                TaskFogMiniTime(i,m) = TaskEndTime(i);
            end
        end
    end

    for i = 1 : TaskNum
        % find vehicleID
        for n = 1 : VehicleNum
            if TaskVehicle(i,n) == 1
                vehicleId = V(n);
            end
        end
        if MobileFogNum ~= 0
            for m = 1 : MobileFogNum
                if TaskFog(i,m) == 1
                    for t = nowTime : TimeSum
                        vehicleLoc = VehicleTrace(((vehicleId-1)*TimeSum + t),:);
                        mobileFogLoc = VehicleTrace(((MF(m)-1)*TimeSum + t),:);
                        if isIn(vehicleLoc,mobileFogLoc) == 0
                            leaveTime = t;
                            TaskFogMiniTime(i,m) = min(leaveTime,TaskFogMiniTime(i,m));
                            break;
                        end
                    end
                end
            end
        end
        for m = MobileFogNum + 1 : FogNum
            if TaskFog(i,m) == 1
                for t = nowTime : TimeSum
                    vehicleLoc = VehicleTrace(((vehicleId-1)*TimeSum + t),:);
                    fixedFogLoc = FixedFogLocal(m-MobileFogNum,:);
                    if isIn(vehicleLoc,fixedFogLoc) == 0
                        leaveTime = t;
                        TaskFogMiniTime(i,m) = min(leaveTime,TaskFogMiniTime(i,m));
                        break;
                    end
                end
            end
        end
    end
    % once arragement
    
    for i = 1 : FogNum
        disp('FogNum is ');
        disp(i);
        
        theTaskFogMiniTime = TaskFogMiniTime(:,i);
        timeMax = max(theTaskFogMiniTime);
        stopTime = timeMax;
        if timeMax+1 < startTime
            disp('timeMax+1 < startTime');
            break;
        end
        if startTime == 0
            disp('startTime == 0');
            break;
        end
        if stopTime == 0
            disp('stopTime == 0');
            continue;
        end
        theFogSize = FogSize(i);
        theTaskInFog = TaskFog(:,i);
        [sortedTask, marked] = sort(theTaskFogMiniTime);
        choosedIdSet = cell(TaskNum,stopTime-startTime,theFogSize);
        profits = zeros(TaskNum,timeMax-startTime-1,theFogSize); 
        for b = 1 : TaskNum
            if sortedTask(b,1) ~= 0
                theTaskInFog = marked([b,TaskNum],:);
                break;
            end
        end
        for b = 1 : numel(theTaskInFog)
            for t = startTime : stopTime
                for s = 1 : theFogSize
                    taskId = theTaskInFog(b);
                    if s < TaskSize(taskId)+1
                        if b == 1
                            profits(b,t,s) = 0;
                            choosedIdSet{b,t,s}=[];
                        else
                            profits(b,t,s) = profits(b-1,t,s);
                            choosedIdSet(b,t,s) = choosedIdSet(b-1,t,s);
                        end
                    else
                        timeTrans = round(TaskSize(taskId) / FogTrans(i));
                        timeProcess = round(TaskCpu(taskId) * FogCompu(i));
                        if numel(choosedIdSet(b,t,s)) == 0
                            if t >= startTime+timeTrans+timeProcess && t <= theTaskFogMiniTime(taskId)
                                profits(b,t,s) = TaskFogProfit(taskId,i);
                                choosedIdSet(b,t,s) = [choosedIdSet{b-1,t,s},taskId];
                                taskFogFinishTime(taskId, i) = t;
                            end
                        else
                            if t >= startTime + timeProcess && t <= theTaskFogMiniTime(taskId)
                                if b > 1
                                    if profits(b-1,t-timeProcess,s-TaskSize(taskId)) + TaskFogProfit(taskId,i) > profits(b,t,s)
                                        profits(b,t,s) = profits(b-1,t-timeProcess,s-TaskSize(taskId)) + TaskFogProfit(taskId,i);
                                        choosedIdSet{b,t,s} = [choosedIdSet{b,t,s},taskId];
                                        taskFogFinishTime(taskId, i) = t;
                                    end
                                end
                            end
                        end
                    end
                    
                end
            end
        end
        disp('profits is');
        disp(profits(b,t,s));
        profitsSum(arragementTime) = profitsSum(arragementTime) + profits(b,t,s);
    end
    
    % vehicle choose the mini time
    for i = 1 : TaskNum
        theTaskFinishTime = taskFogFinishTime(i,:);
        [sortedFinishTime, fogMarked] = sort(theTaskFinishTime);
        for j = 1 : FogNum
            if sortedFinishTime(1,j) ~= 0
                TaskFinish(i) = 1;
                TaskChoosed(i,j) = 1;
                break;
            end
        end
    end
    % when the task is not down and the time is not over
    % renew the startTime
    taskTimeSumInFog = zeros(FogNum,1);
    for i = 1 : FogNum
        for j = 1 : TaskNum
            if TaskChoosed(j,i) == 1
                taskTimeSumInFog(i) = taskTimeSumInFog(i) + round(TaskCpu(j) * FogCompu(i));
            end
        end
    end
    startTime = max(taskTimeSumInFog);
    
    disp('StartTime = ');
    disp(startTime);
    
    disp('Result = ');
    disp(TaskFinish);
    
    complete = sum(TaskFinish) / TaskNum;
    
    disp('Complete rate is');
    disp(complete);
    
    if startTime >= 300
        disp('startTime >= 300');
        break;
    end
    
    if lastComplete == complete
        disp('lastComplete = complete');
        break;
    end
    
    endTime(arragementTime) = startTime;
    lastComplete = complete;
    arragementTime = arragementTime + 1;
    disp('new arrage');
end

complete = sum(TaskFinish) / TaskNum;
disp('Final complete rate is');
disp(complete);





