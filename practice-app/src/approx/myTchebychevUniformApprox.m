function coeff = myTchebychevUniformApprox(f)
% 截断 Tchebychev 级数近似最佳一致三次逼近
% 返回值 coeff 为 (1,4) 的向量，分别表示从常数项到三次项的系数

% Tchebychev 多项式
T0 = @(x)(1);
T1 = @(x)(x);
T2 = @(x)(2*x.^2-1);
T3 = @(x)(4*x.^3-3.*x);

c0 = 2/pi*integral(@(x)(f(x).*T0(x)./sqrt(1-x.^2)),-1,1);
c1 = 2/pi*integral(@(x)(f(x).*T1(x)./sqrt(1-x.^2)),-1,1);
c2 = 2/pi*integral(@(x)(f(x).*T2(x)./sqrt(1-x.^2)),-1,1);
c3 = 2/pi*integral(@(x)(f(x).*T3(x)./sqrt(1-x.^2)),-1,1);

coeff = [c0/2 - c2, c1-3*c3, 2*c2, 4*c3];

end

