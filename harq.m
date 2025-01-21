%Type I HARQ混合自动请求重传请求
clc
clear
%设置莱斯衰落信道环境参数
%coderate=1e9;%采样频率1Ghz
%K_factor=1000;%莱斯因子
%chan=Ricianchan(1/coderate,1000,K_factor,[0,0.0001/coderate,0.00001/coderate,0.000001/coderate],[0,-1000,-100,-100]);%定义莱斯衰落信道
%delay=chan.ChannelFilterDelay;%莱斯信道滤波器的时延
%chan;
max_erorr_bits = 10000;
max_frames = 100000;
max_erorr_frames = 50;


%load("qam256Train_SNR4complex_values.mat")
load("ATSC64.mat")
load("qam64_AE.mat")
constellation_AE = complex_values(:,1) + complex_values(:,2)*1i;
constellation_ATSC = cons64(:,6);

load("BG1_Z48.mat");
z = 48;
K = 22 ;
N = 46 ;
r = 1/2;

H_SV = H_SV(1:(N-K),1:N);
pcmatrix  = ldpcQuasiCyclicMatrix(z,H_SV);
cfgLDPCDec = ldpcDecoderConfig(pcmatrix,"norm-min-sum");


pun_pre = 2*z;
pun_suf = 0;
llr_pre = zeros(pun_pre,1);
llr_suf = zeros(pun_suf,1);



M=[2,2,4,4,16,16,64,64,0];% M进制调制
rate=[1/2,3/4,1/2,3/4,1/2,3/4,2/3,3/4,0];%码率为k/n的卷积码编码
r=[1.8 4.46 6.3 8.3 12.95 15.15 21.1 22.45];%门限值限定
SNR=6:0.2:6.8;
eff=(3600-96)/3600;
model=1:8;
 m=3;
% n=5;
for n=1:length(SNR)
    frame_num = 0;
    erorr_bits = 0;
    error_frames = 0;
    EbN0 = 10 ^(SNR(n) / 10);
    sigma = 1/sqrt(2* log2(M(7))* r(1) * EbN0);
%      data=[];
      data_buffer=[];
      buffer=[];
      kk=0;
      num_repeat=zeros(41,3);
%等待传输的数据 初始化
      msg = randi(2,1,K*z)-1;

      data_buffer=msg;
      buffer=[];
      window=4;
      num_nack=zeros(1,window);
      ack=[ones(1,window),zeros(1,274)];
      rxsig=[];
      rx=[];
      rx11=[];
      rx1=[];
      rx_inter=[];
      decoded=[];
      receiver=[];
      
%首先发送window个数据
 
      for i=1:window  
          i
         data_buffer(:,i)=msg(:,i);
         buffer(:,i)=data_buffer(:,i);
   
%编码调制
         %code=convo(buffer(i,:),rate(m));%卷积编码
         codeword = nrldpc_encode(H_SV,z,msg);
         %code_inter=randintrlv(code,4831);%随机交织
         %[code1,yu]=translate(code_inter,M(m));%进制转换
         %modulatedsig=modulation_model(code1,model(m)); % dpsk调制
         codeword_pun = codeword(pun_pre+1:end-pun_suf);
         tx = qammod(codeword_pun,M(7),'UnitAveragePower',true,'InputType','bit','PlotConstellation',0);
%通过莱斯衰落信道
%          modulatedsig1=[modulatedsig,zeros(1,delay)];%为时移补0
%          fadedsig=filter(chan,modulatedsig1);%调制编码信号通过莱斯衰落信道  
%通过加性高斯白噪声信道
         %rxsig(i,:)=awgn(fadedsig,SNR(n),'measured');%加性高斯白噪声
         noise(:,i) = sigma * randn(size(tx)) + sigma*randn(size(tx)) * 1i;
         rx(:,i) = tx + noise(:,i);
 %解调解码
        llr(:,i) = qamdemod(rx(:,i),M(7),"UnitAveragePower",1,"OutputType",'approxllr');
        llr_1 = [llr_pre;llr;llr_suf];
        [decode,~,~] = ldpcDecode(llr_1(1:N*z),cfgLDPCDec,50,"DecisionType","soft");%25是最大迭代次数
        decode = decode' < 0;
        receiver(:,i)=decode(:,i);
%        e1 = sum(abs(decode - msg));
%         if e1>0 
%             erorr_bits = erorr_bits + e1;
%             error_frames = error_frames + 1;
%         end
     

%          rx(i,:)=demodulation_model(rxsig(i,:),model(m));%解调
%          rx11(i,:)=rx(i,1+delay:end);%去除时移
%          rx1(i,:)=detranslate(rx11(i,:),M(m),yu);%进制转换
%          rx_inter(i,:)=randdeintrlv(rx1(i,:),4831);%解随机交织
%          decoded(i,:)=deconv(rx_inter(i,:),rate(m));%viterbi解码  
%          receiver(i,:)=decoded(i,:);
        [num(i),nErrors(i)]=biterr(receiver(:,i),buffer(:,i));%比较误比特率
        if num(i)==0
           ack(i)=1;
           ber(i)=0;
           num_nack(i)=0;
        else
           ack(i)=0;
           ber(i)=nErrors(i);
           num_nack(i)=1;
        end
      end
% length(find(ack==1))
% 继续发送数据
 
      i=window+1;
       while i<=278
          [n,i]
%     for ii=i/window:9 %第几组数据
          for jj=1:window %每组中的window个数据      
 %确定window个数据的值
            if ack(i-window)==1
              num_cuo=length(find(ack(1:i-window)==0));
              data_buffer(i,:)=data(i-num_cuo,:); 
              buffer(jj,:)=data_buffer(i,:);
            else
              num_repeat(n,m)=num_repeat(n,m)+1;%重传帧的数目
              buffer(jj,:)=data_buffer(i-window,:);
              data_buffer(i,:)=data_buffer(i-window,:);
            end
%编码调制
           code=convo(buffer(jj,:),rate(m));%卷积编码
           code_inter=randintrlv(code,4831);%随机交织
           [code1,yu]=translate(code_inter,M(m));%进制转换
           modulatedsig=modulation_model(code1,model(m)); % dpsk调制
%通过莱斯衰落信道
            modulatedsig1=[modulatedsig,zeros(1,delay)];%为时移补0
           fadedsig=filter(chan,modulatedsig1);%调制编码信号通过莱斯衰落信道
 %通过加性高斯白噪声信道
           rxsig(i,:)=awgn(fadedsig,SNR(n),'measured');%加性高斯白噪声
 %解调解码
           rx(i,:)=demodulation_model(rxsig(i,:),model(m));%解调
           rx11(i,:)=rx(i,1+delay:end);%去除时移
           rx1(i,:)=detranslate(rx(i,:),M(m),yu);%进制转换
           rx_inter(i,:)=randdeintrlv(rx1(i,:),4831);%解随机交织
            decoded(i,:)=deconv(rx_inter(i,:),rate(m));%viterbi解码  
           receiver(i,:)=decoded(i,:);
           [num(i),nErrors(i)]=biterr(receiver(i,:),data_buffer(i,:));%比较误比特率
        
 
%判断         
           if num(i)==0
             ber(i)=0;
             ack(i)=1;
             num_nack(jj)=0;
             i=i+1;
           else
             num_nack(jj)=num_nack(jj)+1;
               ber(i)=nErrors(i);
               if num_nack(jj)==3
                  ack(i)=1;
                  kk=kk+1;
                  num_nack(jj)=0;
                  i=i+1;
                else
                  ack(i)=0;
                   i=i+1;
                end
             end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
           end
       end
   ber_end(n,m)=mean(ber)
   num_succ_mean(n,m)=length(find(ack(1:278)==1))-kk;
   num_succ_rate_mean(n,m)=num_succ_mean(n,m)./(278-num_repeat(n,m));
   fer(n,m)=1-num_succ_rate_mean(n,m);
   num_time(n,m)=278./num_succ_mean(n,m);
   throughout(n,m)=num_succ_mean(n,m).*3600.*log2(M(m)).*rate(m).*eff./10^6;
% end
end  
%计算每帧平均传输时间，误帧率，系统吞吐量
for i=1:length(SNR)
%     for j=1:length(model)
       if num_time(i,3)>=3
         num_time(i,3)=3;
       else
         num_time(i,3)=num_time(i,3);
       end
%     end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)
plot(SNR,fer)
grid on
xlabel('信噪比Eb/No(dB)')
ylabel('数据帧误帧率fer')
legend('选择重传 Type I HARQ')
title('选择重传 Type I HARQ误帧率性能，QPSK调制,码率1/2')
figure(2)
plot(SNR,ber_end)
grid on
xlabel('信噪比Eb/No(dB)')
ylabel('数据帧误码率ber')
legend('FEC差错控制方式')
title('选择重传 Type I HARQ误码率性能，QPSK调制,码率1/2')   
figure(3)
plot(SNR,num_time)
grid on
xlabel('信噪比Eb/No(dB)')
ylabel('数据帧平均传输时间')
legend('选择重传 Type I HARQ')
title('选择重传 Type I HARQ平均传输时间性能，QPSK调制,码率1/2')
figure(4)
plot(SNR,throughout)
grid on
xlabel('信噪比Eb/No(dB)')
ylabel('系统吞吐量/Mbps')
legend('选择重传 Type I HARQ')
title('选择重传 Type I HARQ系统吞吐量性能，QPSK调制,码率1/2')
              
   

