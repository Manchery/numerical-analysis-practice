function y = myRungeKutta24(f, y0, a, b, h, eps)
% 隐式二级四阶 Runge-Kutta 方法
% 求解微分方程 dy/dx = f(x, y), x in [a,b]; y(a) = y0
% h 为步长

c1 = 1/2-sqrt(3)/6; c2 = 1/2+sqrt(3)/6; 
b1 = 1/2; b2 = 1/2;
a11 = 1/4; a12 = 1/4-sqrt(3)/6;
a21 = 1/4+sqrt(3)/6; a22 = 1/4;

n = floor((b-a)/h);
x = a + h*(0:n); 
y = y0; 
yn = y0;
for i = 1:n
    xn = x(i);
    k1 = f(xn, yn); k2 = k1;
    while true % 迭代
        new_k1 = f(xn+c1*h, yn + a11*h*k1 + a12*h*k2);
        new_k2 = f(xn+c2*h, yn + a21*h*k1 + a22*h*k2);
        if (max_error(k1, new_k1)< eps && max_error(k2, new_k2)< eps)
            break
        end
        k1 = new_k1; k2 = new_k2;
    end
    yn1 = yn + h*(b1*k1+b2*k2);
    y = [y yn1]; yn = yn1;
end

end

