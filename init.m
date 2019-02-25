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
        vehicleTrace(vehicleID*time,1) = x;
        vehicleTrace(vehicleID*time,2) = y;
    end
end


    
    

