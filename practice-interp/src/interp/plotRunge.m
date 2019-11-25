hold on;

sx = (0:100)*2/100-1;
sy = f(sx);
plot(sx, sy);

n = 6;
x = (0:n)*2/n-1; y = f(x);
coeff_lagrange = myLagrangeInterp(x, y);
[sx, sy] = sampleFunction([-1, 1], coeff_lagrange, 100);
plot(sx, sy);

n = 10;
x = (0:n)*2/n-1; y = f(x);
coeff_lagrange = myLagrangeInterp(x, y);
[sx, sy] = sampleFunction([-1, 1], coeff_lagrange, 100);
plot(sx, sy);

n = 14;
x = (0:n)*2/n-1; y = f(x);
coeff_lagrange = myLagrangeInterp(x, y);
[sx, sy] = sampleFunction([-1, 1], coeff_lagrange, 100);
plot(sx, sy);

xlabel('x');
ylabel('y');
lgd = legend('f(x)', 'L_{6}(x)', 'L_{10}(x)', 'L_{14}(x)');
lgd.Location = 'north';

function y=f(x)
    y = 1./(1+25*x.^2);
end