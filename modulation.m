function [tx] = modulation(x,M,com_order, complex_value)
%MODULATION 调制函数
%   M 调制阶数  tx 调制信号 默认 tx 按照一维向量形式默认输入
% com_order  星座点顺序  complex_value  对应的星座点值
l = length(x);
m = log2(M);
ord = zeros(1,l/m);
for i =1:1:l/log2(M)
    ord(i) = bi2de(x((i-1)*m+1:i*m),"left-msb");
end
[~,index] = sort(com_order);
sort_complex_value = complex_value(index);
tx = sort_complex_value(ord+1);
end

