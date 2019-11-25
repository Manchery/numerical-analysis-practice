function coeff = myLagrangeInterp(x,y)
% lagrange 插值法
% 返回格式：coeff 为 (1,n+1) 大小的向量，分别表示从常数项到最高次项的系数

[~,n] = size(x);
n = n-1;

coeff = zeros(1,n+1);

l = zeros(n+1,n+1); % 基函数

for i = 1:n+1
    l(i,1) = 1; prod = 1;
    for j = 1:n+1
        if (i==j)
            continue
        end
        l(i,:) = [0, l(i, 1:end-1)] - x(j)*l(i,:);  
        prod = prod*(x(i)-x(j));
    end
    coeff = coeff + y(i)*l(i,:)/prod;
end

end

