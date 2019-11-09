function [J] = givensForJiacobi(A,k_,l_)
% Jiacobi方法中的givens变换矩阵

k = min(k_,l_); l = max(k_,l_);

[n,~] = size(A);

if (A(k,k) ~= A(l,l))
    ct2 = (A(k,k)-A(l,l))/2/A(k,l);
    t = sign(ct2)/(abs(ct2) + sqrt(1+ct2^2));
    c = 1/sqrt(1+t^2);
    s = c*t;
else
    c = cos(pi/4);
    s = sin(pi/4);
end

J = eye(n);
J(k,k) = c; J(l,l) = c;
J(k,l) = s; J(l,k) = -s;

end

