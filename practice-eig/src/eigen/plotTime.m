% 比较运行时间

t = 3:100;

times = zeros(3, length(t));

for i=1:length(t)
    n = i;
    A = diag(ones(1,n)*2) +diag(-ones(1,n-1),-1) +diag(-ones(1,n-1),1);  
    times(:,i) = getTimes(A);
end

plot(t,times,'LineWidth',1.5);
xlabel('n');
ylabel('time (s)');
lgd = legend('Jacobi classic', 'Jacobi threshold', 'QR');
lgd.Location = 'northwest';

function times = getTimes(A)
    tic;
    myJacobiClassic(A, 1e-7);
    time_jacobi_cls = toc;
    
    tic;
    myJacobiThreshold(A, 1e-7);
    time_jacobi_thr = toc;
   
    tic;
    myQR(A, 1e-7);
    time_qr = toc;
    
    times = [
        time_jacobi_cls;
        time_jacobi_thr;
        time_qr
    ];
end