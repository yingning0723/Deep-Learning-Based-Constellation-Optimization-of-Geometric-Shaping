function [llr] = df_qamdemod(rx,M,bits_seq, complex_value)
%与matlab qamdemod  approxllr 一模一样
%M 阶数  rx 接收复信号  bits_seq  M个log2（M）长度比特序列  complex_value：对应的星座点 mode：max-log-map or log-map
%每个bit输出的llr 等价于 1星座点集最小距离 - 0星座点集最小距离
signal_len = length(rx);
m = log2(M);
bit1_location = zeros(m,M/2);
bit0_location = zeros(m,M/2);

for i =  1:1:m
    bit1_location(i,:) = find(bits_seq(:,i));
    bit0_location(i,:) =  setdiff(1:M, bit1_location(i,:));
end

llr = zeros(m,signal_len);
for i= 1:1:signal_len
    distance = abs(complex_value - rx(i)).^2;
    % 近似计算概率正比于e^(-d^2)
    for j = 1:1:m
       llr(j,i) = -min(distance(bit0_location(j,:))) +  min(distance(bit1_location(j,:)));
    end
end


