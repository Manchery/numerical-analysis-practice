function [lambda, v] = inversePowerMethod(A, q, eps)
% 原点位移的反幂法求解特征值和特征向量
% eps: 相邻两次迭代结果小于eps，则终止

[n,~] = size(A);
A = A-q*eye(n);
v = ones(n, 1);

k=0;

while (true)
    k = k+1;
    z = A\v;
    [~,pos] = max(abs(z));
    m = z(pos); % 最大模元素
    v = z/m;
    if (k>1 && abs(1/m-1/last_m)<eps) 
        break;
    end
    last_m = m;
end

lambda = 1/m+q;

end