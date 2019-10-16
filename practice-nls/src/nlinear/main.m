% main

phi1 = @(x)(20-2*x^2-x^3)/10.0;
phi2 = @(x)nthroot(20-10*x-2*x^2,3);
f = @(x)x^3+2*x^2+10*x-20;
df = @(x)3*x^2+4*x+10;

[x_fp1, xs_fp1] = myFixedPoint(1, phi1, 1e-8, 50);
[x_fp2, xs_fp2] = myFixedPoint(1, phi2, 1e-8, 50);

[x_steff1, xs_steff1] = mySteffensen(1, phi1, 1e-8, 50);
[x_steff2, xs_steff2] = mySteffensen(1, phi2, 1e-8, 50);
[x_newton, xs_newton] = myNewton(1, f, df, 1e-8, 50);

disp('x series fp1:'); disp(vpa(xs_fp1,12));
disp('x series fp2:'); disp(vpa(xs_fp2,12));
disp('x series steff1:'); disp(vpa(xs_steff1,12));
disp('x series steff2:'); disp(vpa(xs_steff2,12));
disp('x series newton:'); disp(vpa(xs_newton,12));

hold on;
% plot(1:length(xs_fp1), xs_fp1, 'LineWidth',1.5);
% plot(1:length(xs_fp2), xs_fp2, 'LineWidth',1.5);
plot(1:length(xs_steff1), xs_steff1, 'LineWidth',1.5);
plot(1:length(xs_steff2), xs_steff2, 'LineWidth',1.5);
plot(1:length(xs_newton), xs_newton, 'LineWidth',1.5);

legend('steff1', 'steff2', 'newton');
% legend('fp1', 'fp2', 'steff1', 'steff2', 'newton');

xlabel('k');
ylabel('x_k');