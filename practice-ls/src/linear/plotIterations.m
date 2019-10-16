% 比较迭代次数

t = 2:20;

iterations = zeros(5, length(t));

for i=1:length(t)
    n=t(i);
    A = hilb(n);
    x_true = ones(n,1);
    b = A*x_true;
    
    iterations(:,i) = getIterations(A, b);
end

plot(t,log10(iterations),'LineWidth',1.5);
xlabel('n');
ylabel('log_{10}(iterations)');
lgd = legend('cg', 'pcg', 'gmres m=2', 'gmres m=5', 'gmres m=10');
lgd.Location = 'southeast';

function iters = getIterations(A, b)
    [~,iter_cg] = myCG(A, b, 10000, 1e-8);
    [~,iter_pcg] = myPCG(A, b, 10000, 1e-8);
    [~,iter_gmres2] = myGMRESm(A, b, 2, 10000, 1e-8);
    [~,iter_gmres5] = myGMRESm(A, b, 5, 10000, 1e-8);
    [~,iter_gmres10] = myGMRESm(A, b, 10, 10000, 1e-8);
    
    iters = [
        iter_cg;
        iter_pcg;
        iter_gmres2;
        iter_gmres5;
        iter_gmres10;
    ];
end