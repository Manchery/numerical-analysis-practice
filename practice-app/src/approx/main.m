hold on;

density = 100;
sx = (0:density)*2/density-1;
sy = f(sx);
plot(sx,sy);

coeff_leg = myLegendreSquareApprox(@f);
fprintf("Legendre: %fx^3+%fx^2+%fx+%f\n", coeff_leg(4), coeff_leg(3), coeff_leg(2), coeff_leg(1));
[sx, sy] = sampleFunction([-1,1], coeff_leg, density);
plot(sx,sy);
% plot(sx,f(sx)-sy);
fprintf("Legendre error: %f\n", compute_error(sx,sy))
fprintf("Legendre l2 error: %f\n", compute_l2_error(@f, @(x)(coeff_leg(1)+coeff_leg(2).*x+coeff_leg(3).*x.^2+coeff_leg(4).*x.^3), @(x)(1)));

coeff_tche = myTchebychevUniformApprox(@f);
fprintf("Tchebychev: %fx^3+%fx^2+%fx+%f\n", coeff_tche(4), coeff_tche(3), coeff_tche(2), coeff_tche(1));
[sx, sy] = sampleFunction([-1,1], coeff_tche, density);
plot(sx,sy);
% plot(sx,f(sx)-sy);
fprintf("Tchebychev error: %f\n", compute_error(sx,sy))

coeff_lag = myLagrangeUniformApprox(@f);
fprintf("Lagrange: %fx^3+%fx^2+%fx+%f\n", coeff_lag(4), coeff_lag(3), coeff_lag(2), coeff_lag(1));
[sx, sy] = sampleFunction([-1,1], coeff_lag, density);
plot(sx,sy);
% plot(sx,f(sx)-sy);
fprintf("Lagrange error: %f\n", compute_error(sx,sy))

lgd = legend('f','Legendre', 'Tchebychev', 'Lagrange');
lgd.Location = 'northwest';

% lgd = legend('Legendre', 'Tchebychev', 'Lagrange');
% lgd.Location = 'south';

function y=f(x)
    y = (x.^2).*log(2+x);
end

function d = compute_error(sx, sy)
y = f(sx);
d = max(abs(y-sy));
end

function d = compute_l2_error(f,g,rho)
d = integral(@(x)(rho(x).*(f(x)-g(x)).^2), -1, 1)^0.5;
end

