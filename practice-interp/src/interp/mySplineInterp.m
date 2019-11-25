function coeff = mySplineInterp(x,y)
% 三次自然样条插值
% 返回值 coeff 为 (n,4) 的矩阵，分别表示 n 段的常数项到三次项系数

[~,n] = size(x);
n = n-1;

coeff = zeros(n,4);

h = x(2:end) - x(1:end-1);
mu = h(1:end-1)./(h(1:end-1)+h(2:end));
lambda = 1-mu;
df = (y(2:end)-y(1:end-1))./(x(2:end)-x(1:end-1));
d = 6*(df(2:end)-df(1:end-1))./(x(3:end)-x(1:end-2));

% 转角方程
A = 2*eye(n-1)+diag(mu(2:end), 1)+diag(lambda(1:end-1), -1);
b = d;

M = [ 0; A\(b'); 0];

for j = 1:n
    coeff(j,:) = coeff(j,:) + [x(j+1)^3, -3*x(j+1)^2, 3*x(j+1), -1]*M(j)/6/h(j);
    coeff(j,:) = coeff(j,:) + [-x(j)^3, 3*x(j)^2, -3*x(j), 1]*M(j+1)/6/h(j);
    coeff(j,:) = coeff(j,:) + [x(j+1), -1, 0, 0]*(y(j)-M(j)*h(j)^2 / 6) / h(j);
    coeff(j,:) = coeff(j,:) + [-x(j), 1, 0, 0]*(y(j+1)-M(j+1)*h(j)^2 / 6) / h(j);
end

end

