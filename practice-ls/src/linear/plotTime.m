% 比较运行时间

t = 2:200;

times = zeros(2, length(t));

for i=1:length(t)
    n=t(i);
    A = hilb(n);
    x_true = ones(n,1);
    b = A*x_true;
    
    times(:,i) = getTimes(A, b, n);
end

plot(t,times,'LineWidth',1.5);
xlabel('n');
ylabel('time (s)');
lgd = legend('pcg', 'gmres');
lgd.Location = 'southeast';

function times = getTimes(A, b, n)
    tic;
    myPCG(A, b, 10000, 1e-8);
    time_pcg = toc;
    tic;
    myGMRESm(A, b, max(int32(n/10),10), 10000, 1e-8);
    time_gmres = toc;
    
    times = [
        time_pcg;
        time_gmres;
    ];
end