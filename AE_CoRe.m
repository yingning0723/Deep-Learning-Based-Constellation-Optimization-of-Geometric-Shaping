clc
clear all;
num_conste_points = 256;
com_order = 0:num_conste_points-1;
bit_seq = de2bi(0:1:num_conste_points-1,"left-msb");
bit_seq = circshift(bit_seq,[0,-2]);
bit_seq(:,1:6) = ~bit_seq(:,1:6);
order = bi2de(bit_seq,"left-msb");
load('qam256Train_SNR4complex_values.mat');
complex_value_x = complex_values(:,1);
complex_value_y = complex_values(:,2);

m1 = zeros(256,2);
figure(1);
scatter(complex_value_x,complex_value_y,"filled");
axis([-2 2 -2 2]);
grid on;
for i= 1:num_conste_points
%     m1(i,1) = complex_value_x(order(i));
%     m1(i,2) = complex_value_y(order(i));
    text (complex_value_x(i),complex_value_y(i),num2str(order(i)),'FontSize',12,'FontWeight','normal');
%    text (complex_value_x(i),complex_value_y(i),num2str(order(i)));
end
%save 256AE_CoRe.mat m1;
% figure(2);
% scatter(m1(:,1),m1(:,2),"filled");
% axis([-2 2 -2 2]);
% grid on;
% for i= 1:num_conste_points
%     text (m1(i,1),m1(i,2),num2str(com_order(i)),'FontSize',14,'FontWeight','bold');
% %    text (complex_value_x(i),complex_value_y(i),num2str(order(i)));
% end