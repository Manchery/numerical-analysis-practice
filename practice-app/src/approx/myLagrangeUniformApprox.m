function coeff = myLagrangeUniformApprox(f)
% 插值余项极小化 近似最佳一致三次逼近
% 返回值 coeff 为 (1,4) 的向量，分别表示从常数项到三次项的系数

n = 4;
% Tchebychev 多项式零点
x = cos((2*(1:n)-1)/2/n*pi);
coeff = myLagrangeInterp(x, f(x));

end

