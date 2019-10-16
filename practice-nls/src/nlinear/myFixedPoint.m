function [x, x_series] = myFixedPoint(x0, phi, eps, kmax)
% 不动点迭代法 解 x=phi(x)
% kmax: 最大迭代次数
% eps: 相邻两次迭代结果小于eps，则终止

x = x0;
x_series = [x0];
for k = 1:kmax
    last_x = x;
    x = phi(last_x);
    x_series = [x_series x];
    if abs(x-last_x) < eps
        break
    end
end

end

