n=20;
A = hilb(n);
x_true = ones(n,1);
b = A*x_true;

t = 2:n;
iterations = zeros(1, length(t));

for i = 1:length(t)
    m = t(i);
    [~,iterations(i)] = myGMRESm(A, b, m, 10000, 1e-8);
end

plot(t,log10(iterations),'LineWidth',1.5);
xlabel('m');
ylabel('log_{10}(iterations)');
