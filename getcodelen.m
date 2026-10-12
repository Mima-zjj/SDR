
function L = getcodelen(a, b, t)
%getcodelen - 计算a到b的十进制数的最大二进制编码长度
%
% Syntax: L = getcodelen(a, b, t)
%
% Long description
    n = length(a);
    L = zeros(n, 1);
    for i = 1:n
       while ((b(i) - a(i))*10.^t) > ((2^L(i)) - 1)
          L(i) = L(i) + 1;
       end
    end
end
