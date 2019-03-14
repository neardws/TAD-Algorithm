% % % % import frome TXT

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
    if vehicleID ~= 0 && time ~= 0
        VehicleTime(vehicleID,time) = 1;
        VehicleTrace((vehicleID-1)*TimeSum+time,1) = trace(3);
        VehicleTrace((vehicleID-1)*TimeSum+time,2) = trace(4);
    end
end

%{
    Init some date
%}

% Vehicle
TaskEtimes = 1e8;
MobileFogNum = 20;
VehicleNum = VehicleSum - MobileFogNum;


% Fog
FogEtimes = 1e-12;

FixedFogNum = 4;
FixedFogLocal = [500 500; 500 1000; 1000 500; 1000 1000];
FixedFogSize = randi([800,1200],FixedFogNum,1);
FixedFogCompu = randi([100,200],FixedFogNum,1)*FogEtimes;
FixedFogTrans = randi([20,40],FixedFogNum,1);

MobileFogID = randperm(110, MobileFogNum);
MobileFogSize = randi([150,300],MobileFogNum,1);
MobileFogCompu = randi([5,10],MobileFogNum,1)*FogEtimes;
MobileFogTrans = randi([10,20],MobileFogNum,1);

FogNum = MobileFogNum + FixedFogNum;
FogSize = cat(1,MobileFogSize,FixedFogSize);
FogCompu = cat(1,MobileFogCompu,FixedFogCompu);
FogTrans = cat(1,MobileFogTrans,FixedFogTrans);

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
    TAD-Algorithm
%}

VehicleTask = randi([0,1],VehicleNum,1);

% Task
TaskNum = sum(VehicleTask(:));
TaskVehicle = zeros(TaskNum,VehicleNum);
TaskSize = randi([10,25],TaskNum,1);
TaskCpu = randi([10,25],TaskNum,1)*TaskEtimes;
TaskEndTime = randi([100,300],TaskNum,1);

% Task & Vehicle
for i = 1 : TaskNum
    for n = 1 : VehicleNum
        if i <= sum(VehicleTask(1:n,:),1)
            TaskVehicle(i,n)=1;
            break;
        end
    end
end

TaskFinish = zeros(TaskNum,1);
TaskChoosed = zeros(TaskNum,FogNum);
startTime = 1;
lastComplete = -1;
profitsSum = zeros(500,1);
endTime = zeros(500,1);
arragementTime = 0;
taskFogFinishTime = zeros(TaskNum, FogNum);

while isTaskDone(TaskFinish)
    %{ 
        NowTime Update Every Arragement
    %} 
    if startTime >= 300
        disp('startTime >= 300');
        break;
    end
    % Task & Fog
    nowTime = startTime;
    arragementTime = arragementTime + 1;

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
    
    
    for i = 1 : TaskNum
        if TaskFinish(i) == 1
            for m = 1 : FogNum
                TaskFog(i,m) = 0;
                TaskFogMiniTime(i,m) = 0;
            end
        end
    end
    

    % Profits of task in fog 
    TaskFogProfit = zeros(TaskNum,FogNum);

    if MobileFogNum ~= 0
        for m = 1 : MobileFogNum
            taskInFog = TaskFog(:,m);
            taskInFogNum = sum(taskInFog(:)==1);
            sumSize = sum(TaskSize.*taskInFog)/taskInFogNum;
            sumCom = sum(TaskCpu.*taskInFog)/taskInFogNum;
            for i = 1 : TaskNum
                if TaskFog(i,m) == 1
                    TaskFogProfit(i,m) = TaskSize(i)/sumSize + TaskCpu(i)/sumCom;
                end
            end
        end
    end

    for m = MobileFogNum +1 : FogNum
        taskInFog = TaskFog(:,m);
        taskInFogNum = sum(taskInFog(:)==1);
        sumSize = sum(TaskSize.*taskInFog)/taskInFogNum;
        sumCom = sum(TaskCpu.*taskInFog)/taskInFogNum;
        for i = 1 : TaskNum
            if TaskFog(i,m) == 1
                TaskFogProfit(i,m) = TaskSize(i)/sumSize + TaskCpu(i)/sumCom;
            end
        end
    end
    
    TaskSumInFog = zeros(FogNum,1);
    for i = 1 : FogNum
        TaskInFog = TaskFog(:,i);
        TaskSumInFog(i) = sum(TaskInFog);
    end

    MaxTaskSumInFog = max(TaskSumInFog);
    
    %{
        Write File
    %}
    if arragementTime == 1
        fileName = 'datas.txt';
        writeNum(fileName,TaskNum);
        writeNum(fileName,FogNum);
        writeNum(fileName,MaxTaskSumInFog);
        writeMatrix(fileName,TaskSize);
        writeMatrix(fileName,TaskCpu);
        writeMatrix(fileName,FogSize);
        writeMatrix(fileName,FogCompu);
        writeMatrix(fileName,FogTrans);
        writeMatrix(fileName,TaskSumInFog);
        writeMatrix(fileName,TaskFog');
        writeMatrix(fileName,TaskFogProfit');
        writeMatrix(fileName,TaskFogMiniTime');
    end
    
    % once arragement
    
    for i = 1 : FogNum
%         disp('FogNum is ');
%         disp(i);
        
        theTaskFogMiniTime = TaskFogMiniTime(:,i);
        timeMax = max(theTaskFogMiniTime);
        stopTime = timeMax;
        if timeMax+1 <= startTime
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
%         disp('profits is');
        if s > theFogSize
%             disp(profits(b,t,theFogSize));
        	profitsSum(arragementTime) = profitsSum(arragementTime) + profits(b,t,theFogSize);
        else
%             disp(profits(b,t,s));
        	profitsSum(arragementTime) = profitsSum(arragementTime) + profits(b,t,s);
        end
        
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
                taskTimeSumInFog(i) = taskTimeSumInFog(i) + TaskCpu(j) * FogCompu(i);
            end
        end
    end
    
    
%     disp('StartTime = ');
%     disp(startTime);
    
%     disp('Result = ');
%     disp(TaskFinish);
    
    complete = sum(TaskFinish) / TaskNum;
    
    disp('Complete rate is');
    disp(complete);
    
%     disp('Task Finish is');
%     disp(numel(TaskFinish));
    
    disp('StartTime = ');
    disp(startTime);
    disp('======================================')
    
    if lastComplete == complete
        startTime = startTime + 1;
        disp('lastComplete = complete');
    else
        if max(taskTimeSumInFog) == 0
            startTime = startTime + 1;
            endTime(arragementTime) = 0;
        else
            startTime = startTime + round(max(taskTimeSumInFog));
            endTime(arragementTime) =  max(taskTimeSumInFog);
        end
        disp('lastComplete ~= complete');
    end
%     disp(startTime);
    
    lastComplete = complete;
%     disp('new arrage');
end

complete = sum(TaskFinish) / TaskNum;
disp('Final complete rate is');
disp(complete);





