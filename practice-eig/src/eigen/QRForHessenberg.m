function [U, R] = QRForHessenberg(H)
% 上Hessenberg矩阵的QR分解
% 结果满足 H = UR

[n,~] = size(H);
Ut = eye(n);

for i=1:n-1
    t = H(i+1,i)/H(i,i);
    c = 1/sqrt(1+t^2);
    s = t*c;
    
    % givens 变换
    H(i,:) = H(i,:)*c + H(i+1,:)*s;
    H(i+1,:) = (c+s^2/c)*H(i+1,:) - s/c*H(i,:);
    Ut(i,:) = Ut(i,:)*c + Ut(i+1,:)*s;
    Ut(i+1,:) = (c+s^2/c)*Ut(i+1,:) - s/c*Ut(i,:);
    
    % J: givens 变换矩阵
    % J = eye(n);
    % J(i,i) = c; J(i+1,i+1) = c;
    % J(i+1,i) = -s; J(i,i+1) = s;
    % H = J*H;
    % Ut = J*Ut;
end

U = Ut';
R = H;

end

