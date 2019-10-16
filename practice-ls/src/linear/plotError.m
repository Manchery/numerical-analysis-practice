% 比较结果相对误差

t = 2:20;

relative_error = zeros(6, length(t));

for i=1:length(t)
    n=t(i);
    A = hilb(n);
    x_true = ones(n,1);
    b = A*x_true;
    
    relative_error(:,i) = getRelativeErrors(A, b, x_true);
end

plot(t,log10(relative_error),'LineWidth',1.5);
xlabel('n');
ylabel('log_{10}(relative error)');
lgd = legend('gauss', 'cholesky', 'tikhonov+gauss', 'tikhonov+cholesky', 'pcg', 'gmres');
lgd.Location = 'southeast';

n = 100;
A = hilb(n);
x_true = ones(n,1);
b = A*x_true;
res = getRelativeErrors(A, b, x_true);
for i=1:6
    disp(res(i));
end

function re = relativeError(x, x_true)
    re = norm(x_true - x) / norm(x_true);
end

function re = getRelativeErrors(A, b, x_true)
    x_gauss = myGauss(A, b);
    x_choles = myCholesky(A, b);
    x_tik_gauss = myTikhonov(A, b, 1e-13, @myGauss);
    x_tik_chole = myTikhonov(A, b, 1e-13, @myCholesky);
    [x_pcg,~] = myPCG(A, b, 10000, 1e-15);
    [x_gmres,~] = myGMRESm(A, b, 5, 10000, 1e-15);
%     x_gmres = gmres(A, b);
    
    re = [
        relativeError(x_gauss, x_true);
        relativeError(x_choles, x_true);
        relativeError(x_tik_gauss, x_true);
        relativeError(x_tik_chole, x_true);
        relativeError(x_pcg, x_true);
        relativeError(x_gmres, x_true);
    ];
end