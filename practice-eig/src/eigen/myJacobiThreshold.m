function [lambda, Q] = myJacobiThreshold(A, eps)
% Jacobi过关算法
% eps: 过关法阈值小于eps，则终止
% lambda: 特征值
% Q: 对应特征值的特征向量
% 满足 inv(Q)*A*Q = diag(lambda)

[n,~] = size(A);

Q = eye(n);
% 初始阈值 sqrt(N(A))/2
delta = sqrt(sum(sum(abs(A)))-sum(sum(diag(abs(A)))))/n;

while (delta>eps)
    while (true)
        mapped = false;
        % 过关扫描
        for i = 1:n
            for j = i+1:n
                if (abs(A(i,j))>delta)
                    % givens 变换
                    J = givensForJiacobi(A, i, j);
                    A = J*A*J';
                    Q = Q*J';
                    mapped = true;
                end
            end
        end
        % 所有元素都绝对值小于阈值
        if (~mapped)
            break
        end
    end
    delta = delta/n;
end

lambda = diag(A);    

end

