function I = myRomberg(f,a,b,eps)
% Romberg 求积方法, f 在 [a,b] 上积分
% eps 为指定的停止时的误差控制

T = (b-a)/2*(f(a)+f(b));
l = 2;
while true
    T = [T, zeros(l-1,1); zeros(1,l)];
    n = 2^(l-1);
    T(l,1) = (f(a)+f(b)+2*sum(f(a+(b-a)/n*(1:n-1))))*(b-a)/n/2;
    for j = 2:l
        T(l,j) = (4^(j-1)*T(l,j-1)-T(l-1,j-1)) / (4^(j-1)-1);
    end
    if abs(T(l,l)-T(l-1,l-1))<eps
        break
    end
    l = l+1;
end

% disp(T);
I = T(l,l);

end

