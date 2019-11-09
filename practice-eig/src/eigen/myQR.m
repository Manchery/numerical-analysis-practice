function [lambda, X] = myQR(A, eps)
% QR算法求解特征值
% eps: 下对角元素绝对值小于eps，则终止
% lambda: 特征值
% X: 对应特征值的特征向量
% 满足 inv(X)*A*X = diag(lambda)

[n,~] = size(A);

% 为了减少运算量，可将A先转化为上Hessenberg阵
% 在本次作业中 A 已经为该形式
% H = toHessenberg(A);
H = A;

while (max(max( tril( abs( H-diag(diag(H)) )) ))>eps)
    [U,R] = QRForHessenberg(H);
    H = R*U;
end

lambda = diag(H);

% 特征向量可用反幂法求得
X = zeros(n,n);
for i=1:n
    [lambda(i), X(:,i)] = inversePowerMethod(A, lambda(i)+(1e-5*randn()), eps);
end

end

