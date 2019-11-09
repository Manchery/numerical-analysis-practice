function [lambda, Q] = myJacobiClassic(A, eps)
% 经典Jacobi方法
% eps: 非对角元素绝对值小于eps，则终止
% lambda: 特征值
% Q: 对应特征值的特征向量
% 满足 inv(Q)*A*Q = diag(lambda)

[n,~] = size(A);

Q = eye(n);

while (true)
    C = abs(A - diag(diag(A)));
    val = max(max(C)); % 非对角元素的最大模
    if (val<eps)
        break
    end
    [k, l] = find(val==C); % 最大模的位置
    J = givensForJiacobi(A, k(1), l(1));
    A = J*A*J';
    Q = Q*J';
end

lambda = diag(A);    

end

