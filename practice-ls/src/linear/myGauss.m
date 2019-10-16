function [x] = myGauss(A, b)
% 选取列主元的高斯消元法 Ax=b
% 假定参数的size满足 A: [n,n], b: [n,1]

[n,~] = size(A);
A = [A, b]; % 增广矩阵

% 选取列主元的高斯消元
for i = 1:n
    [~, pivot] = max(abs(A(i:end, i)));
    pivot = pivot + i - 1; % 列主元
    A([i,pivot],:) = A([pivot,i],:);
    scale = A(i+1:end, i) / A(i,i);
    A(i+1:end, :) = A(i+1:end, :) - scale * A(i, :);    
end

% 解 Ux=b1

x = zeros(n,1);
x(n) = A(n,end) / A(n,n);
for i = n-1:-1:1
    x(i) = (A(i, end) - A(i, i+1:end-1) * x(i+1:end)) / A(i,i);
end

end
