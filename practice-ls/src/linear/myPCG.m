function [x, iters] = myPCG(A, b, kmax, eps)
% 预处理的共轭梯度法 Ax=b
% 假定参数的size满足 A: [n,n], b: [n,1]
% 返回值 iter 表示收敛时迭代次数

[n,~] = size(A);

M = diag(diag(A)); % Jacobi迭代的分裂矩阵

x = zeros(n,1);
r = b-A*x; z = M\r; p = z;

for k = 1:kmax
    alpha = (z'*r)/(p'*(A*p));
    x = x+alpha*p;
    last_r = r; r = r-alpha*A*p;
    last_z = z; z = M\r;
    beta = (z'*r)/(last_z'*last_r);
    p = z+beta*p;
    
    if norm(r)/norm(b) < eps
        iters = k;
        return;
    end
end

iters = kmax;

end