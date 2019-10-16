function [x, iters] = myCG(A, b, kmax, eps)
% 共轭梯度法 Ax=b
% 假定参数的size满足 A: [n,n], b: [n,1]
% 返回值 iter 表示收敛时迭代次数

[n,~] = size(A);

x = zeros(n,1);
r = b-A*x;  p = r;

for k = 1:kmax
    alpha = (r'*r)/(p'*(A*p));
    x = x+alpha*p;
    last_r = r; r = r-alpha*A*p;
    beta = (r'*r)/(last_r'*last_r);
    p = r*beta+p;
    
    if norm(r)/norm(b) < eps
        iters = k;
        return;
    end
end

iters = kmax;

end