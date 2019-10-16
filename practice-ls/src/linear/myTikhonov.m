function [x] = myTikhonov(A, b, delta, directMethod)
% Tikhonov正则化方法改进条件数求解 Ax=b
% 假定参数的size满足 A: [n,n], b: [n,1]
% directMethod为直接解法的函数指针
%  如 @myGauss @myCholesky

[n,~] = size(A);
eigen = eig(A);
miu1 = eigen(end);
alpha = miu1^2*delta;

% 调用直接解法求解改进后的方程组
x = directMethod(alpha*eye(n)+A'*A, A'*b);

end
