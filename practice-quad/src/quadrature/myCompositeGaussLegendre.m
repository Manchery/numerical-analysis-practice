function I = myCompositeGaussLegendre(f,a,b,n)
% 五点 Gauss-Legendre 求积公式复合求 f 在 [a,b] 上积分
% 将 [a,b] 等分成 n 段

I = 0;
for i=0:n-1
    ai = a+(b-a)/n*i;
    bi = a+(b-a)/n*(i+1);
    I = I + myGaussLegendre(f,ai,bi);
end

end
