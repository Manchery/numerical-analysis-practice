function [x, iters] = myGMRESm(A, b, m, kmax, eps)
% GMRES算法求解 Ax=b
% 假定参数的size满足 A: [n,n], b: [n,1]

[n,~] = size(A);
m = min(m,n);
x0 = zeros(n,1);
% x0 = randn(n,1)*0.0001;

for k = 1:kmax
    r0 = b-A*x0;
    v1 = r0/norm(r0);

    % Arnoldi过程
    Vm = zeros(n,m+1); Vm(:,1) = v1;
    Hm = zeros(m+1,m);
    success = true;

    for i = 1:m
        for j = 1:i
            Hm(j,i) = dot(A*Vm(:,i), Vm(:,j)) / dot(Vm(:,j),Vm(:,j));
        end
        ri = A*Vm(:,i) - Vm(:,1:i)*Hm(1:i,i);
        if ~any(ri(:)) 
            success = false;
            break
        end % ri==0，中断
        if i<n
            Hm(i+1,i) = norm(ri);
            Vm(:,i+1) = ri/norm(ri);
        end
    end
    
    if ~success 
        x0 = randn(n,1)*0.0001;
        continue
    end % 若Arnoldi过程中断，重新选取x0

    Vm = Vm(:,1:m);
    % Arnoldi过程结束

    ym = Hm\(norm(r0)*speye(m+1,1));
    xm = x0+Vm*ym;
    rm = b-A*xm;
    
    if norm(rm)/norm(b)<eps
        x = xm;
        iters = k;
        return
    else
        x0 = xm;
    end
    
end

iters = kmax;
x = xm;

end

