function [llr] = df_qamdemod(rx,M,bits_seq, complex_value)
%��matlab qamdemod  approxllr һģһ��
%M ����  rx ���ո��ź�  bits_seq  M��log2��M�����ȱ�������  complex_value����Ӧ�������� mode��max-log-map or log-map
%ÿ��bit�����llr �ȼ��� 1�����㼯��С���� - 0�����㼯��С����
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
    % ���Ƽ������������e^(-d^2)
    for j = 1:1:m
       llr(j,i) = -min(distance(bit0_location(j,:))) +  min(distance(bit1_location(j,:)));
    end
end


