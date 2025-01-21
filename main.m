clear;
% load('D:\graduation\LDPC\BG1_Z48_gpu.mat');  
% save('D:\graduation\LDPC\BG1_Z48_gpu.txt','H_SV');

M = 64;
m_bits_per_s = log2(M);


max_erorr_bits = 10000;
max_frames = 100000;
max_erorr_frames = 50;


%load("qam256Train_SNR4complex_values.mat")
load("ATSC64.mat")
load("qam64_AE.mat")
constellation_AE = complex_values(:,1) + complex_values(:,2)*1i;
constellation_ATSC = cons64(:,6);%256需要10/15码率的 为9；64qam的是1/2码率，为6

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


EbN0_db_qam = 10.4:0.2:12.2;
for i= 1:1:length(EbN0_db_qam)
frame_num = 0;
erorr_bits = 0;
error_frames = 0;
EbN0 = 10 ^(EbN0_db_qam(i) / 10);
sigma = 1/sqrt(2* log2(M)* r * EbN0);
    while(frame_num < max_frames && error_frames < max_erorr_frames)
        frame_num = frame_num + 1;  %记录帧数
        msg = randi(2,1,K*z)-1;
        codeword = nrldpc_encode(H_SV,z,msg);
        codeword_pun = codeword(pun_pre+1:end-pun_suf);
        tx = qammod(codeword_pun',M,'UnitAveragePower',true,'InputType','bit','PlotConstellation',0);
        noise = sigma * randn(size(tx)) + sigma*randn(size(tx)) * 1i;
        rx = tx + noise;
        llr = qamdemod(rx,M,"UnitAveragePower",1,"OutputType",'approxllr');
        llr = [llr_pre;llr;llr_suf];
        [decode,~,~] = ldpcDecode(llr(1:N*z),cfgLDPCDec,25,"DecisionType","soft");%25是最大迭代次数
        decode = decode' < 0;
        e1 = sum(abs(decode - msg));
        if e1>0 
            erorr_bits = erorr_bits + e1;
            error_frames = error_frames + 1;
        end
    end

    fer_qam(i) = error_frames/frame_num;
    ber_qam(i) = erorr_bits /K/z/frame_num;
    bler_qam(i) = erorr_bits/z/frame_num;
   % disp(ber_qam(i))
end


EbN0_db_AE = 9.6:0.2:11.4;
for i= 1:1:length(EbN0_db_AE)
frame_num = 0;
erorr_bits = 0;
error_frames = 0;
EbN0 = 10 ^(EbN0_db_AE(i) / 10);
sigma = 1/sqrt(2* log2(M)* r * EbN0);
    while(frame_num < max_frames && error_frames < max_erorr_frames)
        frame_num = frame_num + 1;  %记录帧数
        msg = randi(2,1,K*z)-1;
        codeword = nrldpc_encode(H_SV,z,msg);
        codeword_pun = codeword(pun_pre+1:end-pun_suf);
        tx = modulation(codeword_pun,M,0:1:M-1,constellation_AE);
        noise = sigma * randn(size(tx)) + sigma*randn(size(tx)) * 1i;
        rx = tx + noise;
        llr = df_qamdemod(rx,M,de2bi(0:1:M-1,"left-msb"),constellation_AE);
        llr = reshape(llr,1,[])';
        llr = [llr_pre;llr;llr_suf];
        [decode,~,~] = ldpcDecode(llr(1:N*z),cfgLDPCDec,25,"DecisionType","soft");
        decode = decode' < 0;
        e1 = sum(abs(decode - msg));
        if e1>0 
            erorr_bits = erorr_bits + e1;
            error_frames = error_frames + 1;
        end
    end

    fer_AE(i) = error_frames/frame_num;
    ber_AE(i) = erorr_bits /K/z/frame_num;
 %   disp(ber_AE(i))
end


EbN0_db_ATSC = 9.6:0.2:11.6;
for i= 1:1:length(EbN0_db_ATSC)
frame_num = 0;
erorr_bits = 0;
error_frames = 0;
EbN0 = 10 ^(EbN0_db_ATSC(i) / 10);
sigma = 1/sqrt(2* log2(M)* r * EbN0);
    while(frame_num < max_frames && error_frames < max_erorr_frames)
        frame_num = frame_num + 1;  %记录帧数
        msg = randi(2,1,K*z)-1;
        codeword = nrldpc_encode(H_SV,z,msg);
        codeword_pun = codeword(pun_pre+1:end-pun_suf);
        tx = modulation(codeword_pun,M,0:1:M-1,constellation_ATSC);
        noise = sigma * randn(size(tx)) + sigma*randn(size(tx)) * 1i;
        rx = tx + noise;
        llr = df_qamdemod(rx,M,de2bi(0:1:M-1,"left-msb"),constellation_ATSC);
        llr = reshape(llr,1,[])';
        llr = [llr_pre;llr;llr_suf];
        [decode,~,~] = ldpcDecode(llr(1:N*z),cfgLDPCDec,25,"DecisionType","soft");
        decode = decode' < 0;
        e1 = sum(abs(decode - msg));
        if e1>0 
            erorr_bits = erorr_bits + e1;
            error_frames = error_frames + 1;
        end
    end

    fer_ATSC(i) = error_frames/frame_num;
    ber_ATSC(i) = erorr_bits /K/z/frame_num;
%    disp(ber_ATSC(i))
end


close
figure

markersize =6 ;
linewidth = 0.9;
xlabel("Eb/N0");ylabel("BER");
semilogy(EbN0_db_qam,ber_qam,"Marker", 'o',"MarkerSize",markersize,"LineStyle","-","Color",[0.8 0.1 0.1],"LineWidth",linewidth);hold on
semilogy(EbN0_db_AE,ber_AE,"Marker", 'p',"Markersize",markersize,"LineStyle","-","Color",[0.5 0.4 0.1],"LineWidth",linewidth);hold on
semilogy(EbN0_db_ATSC,ber_ATSC,"Marker", 'square',"MarkerSize",markersize,"LineStyle","-","Color",[0.1 0.1 0.8],"LineWidth",linewidth);hold on
legend("BER-QAM256",'BER-AE256','BER-ATSC256',"Location","southwest");
%legend("BER-QAM64",'BER-AE64',"Location","southwest");
grid on
set(gca,'FontSize',14,'FontName','Times New Roman');
set(gca,'ytick',[1e-4,1e-3,1e-2,1e-1,1],'ygrid','on','gridlinestyle','-','Gridalpha',0.1)
ylim([1e-4,1])
xlabel("Eb/N0(dB)");
ylabel("BER");

hold off 
