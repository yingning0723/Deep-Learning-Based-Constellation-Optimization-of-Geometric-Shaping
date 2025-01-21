clc;clear all;close all;
M = 256;

load("ATSC256.mat")
constellation_ATSC = cons256(:,10);
m1=zeros(2,256);
m1(1,:)=real(constellation_ATSC);
m1(2,:)=imag(constellation_ATSC);

N = log2(M);
com_order = 0:M-1;
bit_labels = de2bi(com_order,'left-msb');
bit_labels = circshift(bit_labels,[0,-2]);
bit_labels(:,1:6) = ~bit_labels(:,1:6);
com_order_2nd = bi2de(bit_labels,'left-msb');

% figure(1)
% scatter(m1(1,:),m1(2,:),"filled");
% axis([-2,2,-2,2]);
% grid on;

% m2 = zeros(2,256);
% for i= 1:M
%     m2(1,i) = m1(1,com_order_2nd(i));
%     m2(2,i) = m1(2,com_order_2nd(i));
% %     plot(m2(i,1),m2(i,2));
% %     text (m2(i,1),m2(i,2),num2str(com_order(i)),'FontSize',12,'FontWeight','normal');
% end
 figure(2)
 scatter(m1(1,:),m1(2,:),"filled");
axis([-2,2,-2,2]);
grid on;
for i= 1:M
    text (m1(1,i),m1(2,i),num2str(com_order_2nd(i)),'FontSize',12,'FontWeight','normal');
end
%save ATSC_CoRe_256.mat m2;
%save ATSC_1st.mat m1;