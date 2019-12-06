I_true = -0.238732414637843;
I_gauss = myCompositeGaussLegendre(@f,1,3,4);
I_romberg = myRomberg(@f,1,3,1e-7);

fprintf("gauss: %.15f\n", I_gauss);
fprintf("romberg: %.15f\n", I_romberg);

Is_gauss = [];
for n = 1:20
    Is_gauss = [Is_gauss myCompositeGaussLegendre(@f,1,3,n)];
end
Is_romberg = [];
for n = 3:13
    Is_romberg = [Is_romberg myRomberg(@f,1,3,10^(-n))];
end

plot((1:20), log10(abs(Is_gauss-I_true)/abs(I_true)), 'LineWidth',1.5);
xlabel('n');
ylabel('log_{10}(relative error)');

% plot((3:13), log10(abs(Is_romberg-I_true)/abs(I_true)), 'LineWidth',1.5);
% xlabel('eps = 10^{-x}');
% ylabel('log_{10}(relative error)');

function y = f(x)
    y = sin(2*pi./x)./(x.^2);
end