% init

origindata = importdata("data.txt");
[m,n] = size(origindata);
vehicleTime = zeros(110,300);
vehicleTrace = zeros(110*300,2);
for i = 1 : m
    trace = origindata(i,:);
    vehicleID = trace(1);
    time = trace(2);
    x = trace(3);
    y = trace(4);
    if vehicleID ~= 0 && time ~= 0
        vehicleTime(vehicleID,time) = 1;
        vehicleTrace((vehicleID-1)*300+time,1) = x;
        vehicleTrace((vehicleID-1)*300+time,2) = y;
    end
end

mobileFogID = randi([1,110],10,1);

nowTime = 1;
MF = []; % Mobile Fog
V = []; % Vehicles
for i = 1 : nowTime
    inMapID = vehicleTime(:,i);
    for j = 1 : vehicleID
        if inMapID(j) == 1
            % j is a mobile Fog
            if ismember(j,mobileFogID) 
                MF = [MF,j];
            % j is a client vehicle
            else
                V = [V,j];
            end
        end  
    end
end

vehicleNum = numel(V);
vehicleTask = randi([1,1l],vehicleNum,1);
taskNum = sum(vehicleTask(:));
taskVehicle = zeros(taskNum,vehicleNum);
taskSize = randi([10,25],taskNum,1);
taskCpu = randi([10,25],taskNum,1);
taskCpu = taskCpu.*1e9;

etimes = 1e-10;
fixedFogNum = 4;
mobileFogNum = numel(MF);
fogNum = fixedFogNum + mobileFogNum;
fixedFogSize = randi([80,150],fixedFogNum,1);
mobileFogSize = randi([50,100],mobileFogNum,1);
fixedFogCompu = randi([5,10],fixedFogNum,1);
fixedFogCompu = fixedFogCompu.*etimes;
mobileFogCompu = randi([1,5],mobileFogNum,1);
mobileFogCompu = mobileFogCompu.*etimes;
fixedFogTrans = randi([20,40],fixedFogNum,1);
mobileFogTrans = randi([10,20],mobileFogNum,1);

taskFog = zeros(taskNum,mobileFogNum+fixedFogNum);
fixedFogLocal = [500 500; 500 1000; 1000 500; 1000 1000];

for i = 1 : taskNum
    for n = 1 : vehicleNum
        if i <= sum(vehicleTask(1:n,:),1)
            taskVehicle(i,n)=1;
            break;
        end
    end
end

for i = 1: taskNum
    for n = 1 : vehicleNum
        if taskVehicle(i,n) == 1
            vehicleId = V(n);
            vehicleLoc = vehicleTrace(((vehicleId-1)*300+nowTime),:);
            if mobileFogNum ~= 0
                for j = 1 : mobileFogNum
                    fogId = MF(j);
                    mobileFogLoc = vehicleTrace((fogId-1)*300+nowTime,:);
                    if isIn(vehicleLoc, mobileFogLoc)
                        taskFog(i,j)=1;
                    end
                end
            end
            for m = 1 : fixedFogNum
                if isIn(vehicleLoc, fixedFogLocal(m,:))
                    taskFog(i,mobileFogNum+m) = 1;
                end
            end
        end
    end
end

taskFogProfit = zeros(taskNum,mobileFogNum+4);

if mobileFogNum ~= 0
    for i = 1 : mobileFogNum
        taskInFog = taskFog(:,i);
        sumSize = sum(taskSize.*taskInFog);
        sumCom = sum(taskCpu.*taskInFog);
        for j = 1 : taskNum
            if taskFog(j,i) == 1
                taskFogProfit(j,i) = taskSize(j)/sumSize + taskCpu(j)/sumCom;
            end
        end
    end
end

for i = mobileFogNum+1 : mobileFogNum+fixedFogNum
    taskInFog = taskFog(:,i);
    sumSize = sum(taskSize.*taskInFog);
    sumCom = sum(taskCpu.*taskInFog);
    for j = 1 : taskNum
        if taskFog(j,i) == 1
            taskFogProfit(j,i) = taskSize(j)/sumSize + taskCpu(j)/sumCom;
        end
    end
end

taskEndTime = randi([1,300],taskNum,1);
taskFogMiniTime = zeros(taskNum, mobileFogNum+fixedFogNum);

for i = 1 : taskNum
    for j = 1 : mobileFogNum + fixedFogNum
        if taskFog(i,j) == 1
            taskFogMiniTime(i,j) = taskEndTime(i);
        end
    end
end

for i = 1 : taskNum
    % find vehicleID
    for n = 1 : vehicleNum
        if taskVehicle(i,n) == 1
            vehicleId = V(n);
        end
    end
    if mobileFogNum ~= 0
        for j = 1 : mobileFogNum
            if taskFog(i,j) == 1
                for t = nowTime : 300
                    vehicleLoc = vehicleTrace(((vehicleId-1)*300+t),:);
                    mobileFogLoc = vehicleTrace(((MF(j)-1)*300+t),:);
                    if isIn(vehicleLoc,mobileFogLoc) == 0
                        leaveTime = t;
                        taskFogMiniTime(i,j) = min(leaveTime,taskFogMiniTime(i,j));
                        break;
                    end
                end
            end
        end
    end
    for j = mobileFogNum + 1 : mobileFogNum + fixedFogNum
        if taskFog(i,j) == 1
            for t = nowTime : 300
                vehicleLoc = vehicleTrace(((vehicleId-1)*300+t),:);
                fixedFogLoc = fixedFogLocal(j-mobileFogNum,:);
                if isIn(vehicleLoc,fixedFogLoc) == 0
                    leaveTime = t;
                    taskFogMiniTime(i,j) = min(leaveTime,taskEndTime(i));
                    break;
                end
            end
        end
    end
end

taskSumInFog = zeros(fogNum,1);
for i = 1 : fogNum
    taskInFog = taskFog(:,i);
    taskSumInFog(i) = sum(taskInFog);
end

maxTaskSumInFog = max(taskSumInFog);

fogSize = cat(1,mobileFogSize,fixedFogSize);
fogCompu = cat(1,mobileFogCompu,fixedFogCompu);
fogTrans = cat(1,mobileFogTrans,fixedFogTrans);

fileName = 'datas.txt';
writeNum(fileName,taskNum);
writeNum(fileName,fogNum);
writeNum(fileName,maxTaskSumInFog);
writeMatrix(fileName,taskSize');
writeMatrix(fileName,taskCpu');
writeMatrix(fileName,fogSize');
writeMatrix(fileName,fogCompu');
writeMatrix(fileName,fogTrans');
writeMatrix(fileName,taskSumInFog');
writeMatrix(fileName,taskFog');
writeMatrix(fileName,taskFogProfit');
writeMatrix(fileName,taskFogMiniTime');
            
% begin algorithm
% for one fog node

taskFinish = zeros(taskNum,1);
taskChoosed = zeros(taskNum,fogNum);
startTime = 1;
lastComplete = -1;

disp('Start\n');

while isTaskDone(taskFinish)
     % once arragement
    taskFogFinishTime = zeros(taskNum, fogNum);
    for i = 1 : fogNum
        
        disp('FogNum is ');
        disp(i);
        
        theTaskFogMiniTime = taskFogMiniTime(:,i);
        timeMax = max(theTaskFogMiniTime);
        stopTime = timeMax-startTime+1;
        theFogSize = fogSize(i);
        theTaskInFog = taskFog(:,i);
        [sortedTask, marked] = sort(theTaskFogMiniTime);
        choosedIdSet = cell(taskNum,stopTime-startTime,theFogSize);
        profits = zeros(taskNum,timeMax-startTime-1,theFogSize); 
        for b = 1 : taskNum
            if sortedTask(b,1) ~= 0
                theTaskInFog = marked([b,taskNum],:);
                break;
            end
        end
        for b = 1 : numel(theTaskInFog)
            for t = startTime : stopTime
                for s = 1 : theFogSize
                    taskId = theTaskInFog(b);
                    if s < taskSize(taskId)+1
                        if b == 1
                            profits(b,t,s) = 0;
                            choosedIdSet{b,t,s}=[];
                        else
                            profits(b,t,s) = profits(b-1,t,s);
                            choosedIdSet(b,t,s) = choosedIdSet(b-1,t,s);
                        end
                    else
                        timeTrans = round(taskSize(taskId) / fogTrans(i));
                        timeProcess = round(taskCpu(taskId) * fogCompu(i));
                        if numel(choosedIdSet(b,t,s)) == 0
                            if t >= startTime+timeTrans+timeProcess && t <= theTaskFogMiniTime(taskId)
                                profits(b,t,s) = taskFogProfit(taskId,i);
                                choosedIdSet(b,t,s) = [choosedIdSet{b-1,t,s},taskId];
                                taskFogFinishTime(taskId, i) = t;
                            end
                        else
                            if t >= startTime + timeProcess && t <= theTaskFogMiniTime(taskId)
                                if b > 1
                                    if profits(b-1,t-timeProcess,s-taskSize(taskId)) + taskFogProfit(taskId,i) > profits(b,t,s)
                                        profits(b,t,s) = profits(b-1,t-timeProcess,s-taskSize(taskId)) + taskFogProfit(taskId,i);
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
    end
    % vehicle choose the mini time
    for i = 1 : taskNum
        theTaskFinishTime = taskFogFinishTime(i,:);
        [sortedFinishTime, fogMarked] = sort(theTaskFinishTime);
        for j = 1 : fogNum
            if sortedFinishTime(1,j) ~= 0
                taskFinish(i) = 1;
                taskChoosed(i,j) = 1;
                break;
            end
        end
    end
    % when the task is not down and the time is not over
    % renew the startTime
    taskTimeSumInFog = zeros(fogNum,1);
    for i = 1 : fogNum
        for j = 1 : taskNum
            if taskChoosed(j,i) == 1
                taskTimeSumInFog(i) = taskTimeSumInFog(i) + round(taskCpu(j) * fogCompu(i));
            end
        end
    end
    startTime = max(taskTimeSumInFog);
    
    disp('StartTime = ');
    disp(startTime);
    
    disp('Result = ');
    disp(taskFinish);
    
    complete = sum(taskFinish) / taskNum;
    
    disp('Complete rate is');
    disp(complete);
    
    if startTime >= 300
        break;
    end
    
    if lastComplete == complete
        break;
    end
    
    % renew dates
    nowTime = startTime;
    lastComplete = complete;
    for i = 1: taskNum
        for n = 1 : vehicleNum
            if taskVehicle(i,n) == 1
                vehicleId = V(n);
                vehicleLoc = vehicleTrace(((vehicleId-1)*300+nowTime),:);
                if mobileFogNum ~= 0
                    for j = 1 : mobileFogNum
                        fogId = MF(j);
                        mobileFogLoc = vehicleTrace((fogId-1)*300+nowTime,:);
                        if isIn(vehicleLoc, mobileFogLoc)
                            taskFog(i,j)=1;
                        end
                    end
                end
                for m = 1 : fixedFogNum
                    if isIn(vehicleLoc, fixedFogLocal(m,:))
                        taskFog(i,mobileFogNum+m) = 1;
                    end
                end
            end
        end
    end
    
    
    for i = 1 : taskNum
        if taskFinish(i) == 1
            for j = 1 : fogNum
                taskFog(i,j) = 0;
                taskFogMiniTime(i,j) = 0;
            end
        end
    end
    
    if mobileFogNum ~= 0
        for i = 1 : mobileFogNum
            taskInFog = taskFog(:,i);
            sumSize = sum(taskSize.*taskInFog);
            sumCom = sum(taskCpu.*taskInFog);
            for j = 1 : taskNum
                if taskFog(j,i) == 1
                    taskFogProfit(j,i) = taskSize(j)/sumSize + taskCpu(j)/sumCom;
                end
            end
        end
    end

    for i = mobileFogNum+1 : mobileFogNum+fixedFogNum
        taskInFog = taskFog(:,i);
        sumSize = sum(taskSize.*taskInFog);
        sumCom = sum(taskCpu.*taskInFog);
        for j = 1 : taskNum
            if taskFog(j,i) == 1
                taskFogProfit(j,i) = taskSize(j)/sumSize + taskCpu(j)/sumCom;
            end
        end
    end
    
    for i = 1 : taskNum
    % find vehicleID
        for n = 1 : vehicleNum
            if taskVehicle(i,n) == 1
                vehicleId = V(n);
            end
        end
        if mobileFogNum ~= 0
            for j = 1 : mobileFogNum
                if taskFog(i,j) == 1
                    for t = nowTime : 300
                        vehicleLoc = vehicleTrace(((vehicleId-1)*300+t),:);
                        mobileFogLoc = vehicleTrace(((MF(j)-1)*300+t),:);
                        if isIn(vehicleLoc,mobileFogLoc) == 0
                            leaveTime = t;
                            taskFogMiniTime(i,j) = min(leaveTime,taskFogMiniTime(i,j));
                            break;
                        end
                    end
                end
            end
        end
        for j = mobileFogNum + 1 : mobileFogNum + fixedFogNum
            if taskFog(i,j) == 1
                for t = nowTime : 300
                    vehicleLoc = vehicleTrace(((vehicleId-1)*300+t),:);
                    fixedFogLoc = fixedFogLocal(j-mobileFogNum,:);
                    if isIn(vehicleLoc,fixedFogLoc) == 0
                        leaveTime = t;
                        taskFogMiniTime(i,j) = min(leaveTime,taskEndTime(i));
                        break;
                    end
                end
            end
        end
    end
end

complete = sum(taskFinish) / taskNum;
disp('Final complete rate is');
disp(complete);

                        