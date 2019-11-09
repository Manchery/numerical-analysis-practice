n = 15;
A = diag(ones(1,n)*2) +diag(-ones(1,n-1),-1) +diag(-ones(1,n-1),1);
eig_true = eig(A);
% disp(eig_true);

disp("# Jacobi classic");
tic;
[lambda_cls, Q_cls] = myJacobiClassic(A, 1e-7);
t = toc;
% disp(sort(lambda_cls));
fprintf("* sf  : %d\n", min(sigFigures(eig_true, sort(lambda_cls))));
fprintf("* time: %f\n", t);
% disp(Q_cls);
% disp(Q_cls\(A*Q_cls));

disp("# Jacobi threshold");
tic;
[lambda_thr, Q_thr] = myJacobiThreshold(A, 1e-7);
t = toc;
% disp(sort(lambda_thr));
fprintf("* sf  : %d\n", min(sigFigures(eig_true, sort(lambda_thr))));
fprintf("* time: %f\n", t);
% disp(Q_thr);
% disp(Q_thr\(A*Q_thr));

disp("# QR")
tic;
[lambda_qr, Q_qr] = myQR(A, 1e-7);
t = toc;
% disp(sort(lambda_qr));
fprintf("* sf  : %d\n", min(sigFigures(eig_true, sort(lambda_qr))));
fprintf("* time: %f\n", t);
% disp(Q_qr);
% disp(Q_qr\(A*Q_qr));

function sf = sigFigures(a,b)
sf = floor(log10(a))-floor(log10(abs(a-b)));
end