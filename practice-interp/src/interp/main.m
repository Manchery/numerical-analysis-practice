n = 10;
fprintf("n = %d\n", n);
x = (0:n)*2/n-1;
y = f(x);

hold on;

sx = (0:100)*2/100-1;
sy = f(sx);
plot(sx, sy);

coeff_lagrange = myLagrangeInterp(x, y);
% disp(coeff_lagrange);
[sx, sy] = sampleFunction([-1, 1], coeff_lagrange, 100);
plot(sx, sy);
fprintf("Lagrange error: %f\n", compute_error(sx,sy))

coeff_linear = myLinearInterp(x, y);
% disp(coeff_linear);
[sx, sy] = sampleFunction(x, coeff_linear, 10);
plot(sx, sy);
fprintf("Linear error: %f\n", compute_error(sx,sy))

coeff_spline = mySplineInterp(x, y);
% disp(coeff_spline);
[sx, sy] = sampleFunction(x, coeff_spline, 10);
plot(sx, sy);
fprintf("Spline error: %f\n", compute_error(sx,sy))

xlabel('x');
ylabel('y');
lgd = legend('f(x)', 'L_{n}(x)', 'I_h(x)', 'S_3(x)');
lgd.Location = 'southeast';

function y=f(x)
    y = 1./(1+25*x.^2);
end

function d = compute_error(sx, sy)
y = f(sx);
d = max(y-sy);
end