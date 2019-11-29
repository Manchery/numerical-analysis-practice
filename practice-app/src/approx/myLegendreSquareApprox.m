function coeff = myLegendreSquareApprox(f)
% 用Legendre多项式求最佳平方三次逼近
% 返回值 coeff 为 (1,4) 的向量，分别表示从常数项到三次项的系数

% 归一化的 Legendre 多项式
P0 = @(x)(sqrt(1/2));
P1 = @(x)(sqrt(3/2).*x);
P2 = @(x)(sqrt(5/8)*3.*(x.^2-1/3));
P3 = @(x)(sqrt(7/8).*x.*(5*x.^2-3));

a0 = integral(@(x)(f(x).*P0(x)),-1,1);
a1 = integral(@(x)(f(x).*P1(x)),-1,1);
a2 = integral(@(x)(f(x).*P2(x)),-1,1);
a3 = integral(@(x)(f(x).*P3(x)),-1,1);

coeff = [sqrt(1/2)*a0-sqrt(5/8)*a2, sqrt(3/2)*a1-3*sqrt(7/8)*a3, 3*sqrt(5/8)*a2, 5*sqrt(7/8)*a3];

end

