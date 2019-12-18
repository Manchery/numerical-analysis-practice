function y = myRungeKutta44(f, y0, a, b, h)
% 古典四级四阶 Runge-Kutta 方法
% 求解微分方程 dy/dx = f(x, y), x in [a,b]; y(a) = y0
% h 为步长

n = floor((b-a)/h);
x = a + h*(0:n); 
y = y0; 
yn = y0;
for i = 1:n
    xn = x(i);
    k1 = f(xn, yn);
    k2 = f(xn + h/2, yn + h/2*k1);
    k3 = f(xn + h/2, yn + h/2*k2);
    k4 = f(xn + h, yn + h*k3);
    yn1 = yn + h/6*(k1+2*k2+2*k3+k4);
    y = [y yn1]; yn = yn1;
end

end