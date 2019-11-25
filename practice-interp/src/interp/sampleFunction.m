function [sx, sy] = sampleFunction(x, coeff, density)

[~,n] = size(x);
n=n-1;

[~,m] = size(coeff);

sx = [];
sy = [];

for i = 1:n
    for j = 0:density
        cur_x = x(i)+(x(i+1)-x(i))/density*j;
        cur_y = 0;
        for k = 1:m
            cur_y = cur_y*cur_x+coeff(i,m-k+1);
        end
        sx = [sx cur_x];
        sy = [sy cur_y];
    end
end
    
end

