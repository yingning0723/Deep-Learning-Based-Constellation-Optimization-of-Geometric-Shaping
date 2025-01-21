% clc;clear all;close all;
% load("qam64_SNR0.mat");
% x = complex_values(:,1) ;
% y = complex_values(:,2);
% 
% scatter(x,y,"filled");
% axis([-2,2,-2,2]);
% grid on;
% 
% clc;clear all;close all;
% load("ATSC64.mat")
% constellation_ATSC = cons64(:,7);
% m1=zeros(2,64);
% m1(1,:)=real(constellation_ATSC);
% m1(2,:)=imag(constellation_ATSC);
% scatter(m1(1,:),m1(2,:),"filled");
% axis([-2,2,-2,2]);
% grid on;
% 
% save square_qam64.mat m1;

clc
clear all;
num_conste_points = 64;
com_order = 0:num_conste_points-1;
% bit_seq = de2bi(0:1:num_conste_points-1,"left-msb");
% order = bi2de(bit_seq);
complex_value = load('qam64Train_SNR2complex_values.mat');
complex_value = struct2array(complex_value);
complex_value = reshape(complex_value,[],2);
complex_value_x = complex_value(:,1);
complex_value_y = complex_value(:,2);

scatter(complex_value_x,complex_value_y,"filled");
axis([-2 2 -2 2]);
grid on;
for i= 1:num_conste_points
    text (complex_value_x(i),complex_value_y(i),num2str(com_order(i)),'FontSize',14,'FontWeight','bold');
%    text (complex_value_x(i),complex_value_y(i),num2str(order(i)));
end



