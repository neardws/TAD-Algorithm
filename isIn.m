function [in] = isIn(a,b)
    x1 = a(1);
    y1 = a(2);
    x2 = b(1);
    y2 = b(2);
    x = abs(x1-x2);
    y = abs(y1-y2);
    if x*x + y*y >= 500*500
        in = 0;
    else
        in = 1;
    end
end

