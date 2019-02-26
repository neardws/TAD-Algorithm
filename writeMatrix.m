function writeMatrix(fileName, matrix)
%WRITE Summary of this function goes here
%   Detailed explanation goes here
    [r,c]=size(matrix);
    fid=fopen(fileName,'a');
    for i=1:r
        for j=1:c
            fprintf(fid,'%d\t',matrix(i,j));
        end
        if i == r
            fprintf(fid,"~\n");
        else
            fprintf(fid,"\n");
        end
    end
    
    fclose(fid);
end

