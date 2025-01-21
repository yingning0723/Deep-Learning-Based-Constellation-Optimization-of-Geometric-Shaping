function [tx] = modulation(x,M,com_order, complex_value)
%MODULATION ���ƺ���
%   M ���ƽ���  tx �����ź� Ĭ�� tx ����һά������ʽĬ������
% com_order  ������˳��  complex_value  ��Ӧ��������ֵ
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

