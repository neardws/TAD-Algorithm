function writeNum(fileName, num)
%WRITE Summary of this function goes here
%   Detailed explanation goes here
    fid=fopen(fileName,'a');
    fprintf(fid,'%d',num);
    fprintf(fid,"~\n");
    fclose(fid);
end
