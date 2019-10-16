function [x] = myCholesky(A, b)
% cholesky分解方法求解 Ax=LL'x=b
% 假定参数的size满足 A: [n,n], b: [n,1]

[n,~] = size(A);

% Cholesky 分解
L = zeros(n,n);
for j = 1:n
    L(j,j) = sqrt(A(j,j) - sum(L(j, 1:j-1).^2));
    L(j+1:n, j) = ( A(j+1:n,j) -  L(j+1:n,1:j-1)*L(j,1:j-1)') / L(j,j); 
end

% Ly=b
y = zeros(n,1);
y(1) = b(1) / L(1,1);
for i = 2:n
    y(i) = (b(i) - L(i,1:i-1)*y(1:i-1)) / L(i,i);
end

% L'x=y
x = zeros(n,1);
L_t = L';
x(n) = y(n) / L_t(n,n);
for i = n-1:-1:1
    x(i) = (y(i) - L_t(i, i+1:end) * x(i+1:end)) / L_t(i,i);
end

end