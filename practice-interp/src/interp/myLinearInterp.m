function coeff = myLinearInterp(x,y)
% 分段线性插值
% 返回值 coeff 为 (n,2) 的矩阵，分别表示 n 段的常数项和一次项系数

[~,n] = size(x);
n = n-1;

coeff = zeros(n,2);

for i = 1:n
    h = x(i+1)-x(i);
    coeff(i,1) = (y(i)*x(i+1)-y(i+1)*x(i))/h;
    coeff(i,2) = (y(i+1)-y(i))/h;
    
end

