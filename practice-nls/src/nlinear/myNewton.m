function [x, x_series] = myNewton(x0, f, df, eps, kmax)
% Newton法 解 f(x)=0
% df: f的导函数
% kmax: 最大迭代次数
% eps: 相邻两次迭代结果小于eps，则终止

x = x0;
x_series = [x0];
for k = 1:kmax
    last_x = x;
    x = x-f(x)/df(x);
    x_series = [x_series x];
    if abs(x-last_x) < eps
        break
    end
end

end

