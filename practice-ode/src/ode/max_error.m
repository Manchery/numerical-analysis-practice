function r = max_error(y1, y2)
% 最大绝对误差
A = abs(y1-y2);
r = max(reshape(A,numel(A),1));
end
