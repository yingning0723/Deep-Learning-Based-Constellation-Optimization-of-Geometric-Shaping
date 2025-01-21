clear;
M = 64;
m_bits_per_s = log2(M);
max_erorr_bits = 10000;
max_frames = 100000;
max_erorr_frames = 50;
load("HARQ_qam64_1st.mat")
constellation_QAM = complex_values(:,1) + complex_values(:,2)*1i;
load("ATSC64.mat")
constellation_ATSC = cons64(:,7);

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


EbN0_db_qam = 7.4:0.2:9;
%EbN0_db_qam_CoRe = 6.8:0.2:8.4;

for i= 1:1:length(EbN0_db_qam)
frame_num = 0;
erorr_bits = 0;
error_frames = 0;
erorr_bits_2 = 0;
error_frames_2 = 0;

EbN0 = 10 ^(EbN0_db_qam(i) / 10);
%EbN0_CoRe = 10 ^(EbN0_db_qam_CoRe(i) / 10);

sigma = 1/sqrt(2*EbN0);
%sigma_CoRe = 1/sqrt(2*EbN0_CoRe);
tic
    while(frame_num < max_frames && error_frames < max_erorr_frames)
        frame_num = frame_num + 1;  %记录帧数
        msg = randi(2,1,K*z)-1;
        codeword = nrldpc_encode(H_SV,z,msg);
        codeword_pun = codeword(pun_pre+1:end-pun_suf);
        tx = modulation(codeword_pun,M,0:1:M-1,constellation_ATSC);
        %tx = qammod(codeword_pun',M,'UnitAveragePower',1,'InputType','bit','PlotConstellation',0);
        noise = sigma * randn(size(tx)) + sigma*randn(size(tx)) * 1i;
        rx = tx + noise;
        llr = df_qamdemod(rx,M,de2bi(0:1:M-1,"left-msb"),constellation_ATSC);
        llr = reshape(llr,1,[])';
        llr = [llr_pre;llr;llr_suf];
        [decode_1,~,~] = ldpcDecode(llr(1:N*z),cfgLDPCDec,25,"DecisionType","soft");
        decode_1 = decode_1' < 0;
        e1 = sum(abs(decode_1 - msg));
        if e1>0 
            erorr_bits = erorr_bits + e1;
            error_frames = error_frames + 1;
        end
        llr_backup = llr;
        tx = modulation(codeword_pun,M,0:1:M-1,constellation_ATSC);
        noise = sigma * randn(size(tx)) + sigma*randn(size(tx)) * 1i;
        rx = tx + noise;
        llr = df_qamdemod(rx,M,de2bi(0:1:M-1,"left-msb"),constellation_ATSC);
        llr = reshape(llr,1,[])';
        llr = [llr_pre;llr;llr_suf];
        llr_new = llr_backup + llr;
        [decode_2,~,~] = ldpcDecode(llr_new(1:N*z),cfgLDPCDec,25,"DecisionType","soft");
        decode_2 = decode_2' < 0;
        e2 = sum(abs(decode_2 - msg));
        if e2>0 
           erorr_bits_2 = erorr_bits_2 + e2;
           error_frames_2 = error_frames_2 + 1;
        end
    end
toc
    fer_qam_2(i) = error_frames_2/frame_num;
    ber_qam_2(i) = erorr_bits_2 /K/z/frame_num;

    fer_qam(i) = error_frames/frame_num;
    ber_qam(i) = erorr_bits /K/z/frame_num;
 
end
markersize =6 ;
linewidth = 0.9;
figure
xlabel("SNR");ylabel("BER");
semilogy(EbN0_db_qam,ber_qam,"Marker", 'o',"MarkerSize",markersize,"LineStyle","-","Color",[0.8 0.1 0.1],"LineWidth",linewidth);hold on
semilogy(EbN0_db_qam,ber_qam_2,"Marker", 'p',"Markersize",markersize,"LineStyle","-","Color",[0.5 0.4 0.1],"LineWidth",linewidth);hold on
%semilogy(EbN0_db_qam,ber_qam_3,"Marker", 'square',"Markersize",markersize,"LineStyle","-","Color",[0.1 0.1 0.8],"LineWidth",linewidth);hold on
hold on

% if strcmp(XPU, 'CPU')
%     legend_name2 = sprintf("BER-QAM-%s", "MATLAB Library function comm.LDPCDecoder Decoding result (CPU)");
%     % HDec = comm.LDPCDecoder(ldpcPropertyValuePairs{:});
% elseif strcmp(XPU, 'GPU')
%     legend_name2 = sprintf("BER-QAM-%s", "MATLAB Library function comm.gpu.LDPCDecoder Decoding result (GPU)");
%     % HDec = comm.gpu.LDPCDecoder(ldpcPropertyValuePairs{:});
% end
legend("QAM64(1st)","QAM64(1st) + QAM64-CoRe(2nd)","Location","southwest")
grid on
set(gca,'FontSize',14,'FontName','Times New Roman');
set(gca,'ytick',[1e-5,1e-4,1e-3,1e-2,1e-1,1],'ygrid','on','gridlinestyle','-','Gridalpha',0.1);
ylim([1e-5,1])

hold off
%%

% EbN0_db_AE = 5.6:0.2:6.4;
% for i= 1:1:length(EbN0_db_AE)
% frame_num = 0;
% erorr_bits = 0;
% error_frames = 0;
% EbN0 = 10 ^(EbN0_db_AE(i) / 10);
% sigma = 1/sqrt(2* log2(M)* r * EbN0);
%     while(frame_num < max_frames && error_frames < max_erorr_frames)
%         frame_num = frame_num + 1;  %记录帧数
%         msg = randi(2,1,K*z)-1;
%         codeword = nrldpc_encode(H_SV,z,msg);
%         codeword_pun = codeword(pun_pre+1:end-pun_suf);
%         tx = modulation(codeword_pun,M,0:1:M-1,constellation_AE);
%         noise = sigma * randn(size(tx)) + sigma*randn(size(tx)) * 1i;
%         rx = tx + noise;
%         llr = gpuArray(df_qamdemod(rx,M,de2bi(0:1:M-1,"left-msb"),constellation_AE));
%         llr = reshape(llr,1,[])';
%         llr = gather([llr_pre;llr;llr_suf]);
%         [decode,~,~] = ldpcDecode(llr(1:N*z),cfgLDPCDec,25,"DecisionType","soft");
%         %[decode,~,~] = ldpcDecode(llr(1:N*z),cfgLDPCDec,25);
%         decode = decode' < 0;
%         e1 = sum(abs(decode - msg));
%         if e1>0 
%             erorr_bits = erorr_bits + e1;
%             error_frames = error_frames + 1;
%         end
%     end
% 
%     fer_AE(i) = error_frames/frame_num;
%     ber_AE(i) = erorr_bits /K/z/frame_num;
%     disp(ber_AE(i))
% end
% 
% %%
% EbN0_db_ATSC = 5.6:0.2:6.4;
% for i= 1:1:length(EbN0_db_ATSC)
% frame_num = 0;
% erorr_bits = 0;
% error_frames = 0;
% EbN0 = 10 ^(EbN0_db_ATSC(i) / 10);
% sigma = 1/sqrt(2* log2(M)* r * EbN0);
%     while(frame_num < max_frames && error_frames < max_erorr_frames)
%         frame_num = frame_num + 1;  %记录帧数
%         msg = randi(2,1,K*z)-1;
%         codeword = nrldpc_encode(H_SV,z,msg);
%         codeword_pun = codeword(pun_pre+1:end-pun_suf);
%         tx = modulation(codeword_pun,M,0:1:M-1,constellation_ATSC);
%         noise = sigma * randn(size(tx)) + sigma*randn(size(tx)) * 1i;
%         rx = tx + noise;
%         llr = gpuArray(df_qamdemod(rx,M,de2bi(0:1:M-1,"left-msb"),constellation_ATSC));
%         llr = reshape(llr,1,[])';
%         llr = gather([llr_pre;llr;llr_suf]);
%         [decode,~,~] = ldpcDecode(llr(1:N*z),cfgLDPCDec,25,"DecisionType","soft");
%         %[decode,~,~] = ldpcDecode(llr(1:N*z),cfgLDPCDec,25);
%         decode = decode' < 0;
%         e1 = sum(abs(decode - msg));
%         if e1>0 
%             erorr_bits = erorr_bits + e1;
%             error_frames = error_frames + 1;
%         end
%     end
% 
%     fer_ATSC(i) = error_frames/frame_num;
%     ber_ATSC(i) = erorr_bits /K/z/frame_num;
%     disp(ber_ATSC(i))
% end
% 
% %%
% close
% figure
% 
% markersize =6 ;
% linewidth = 0.9;
% xlabel("Eb/N0");ylabel("BER");
% semilogy(EbN0_db_qam,ber_qam,"Marker", 'o',"MarkerSize",markersize,"LineStyle","-","Color",[0.8 0.1 0.1],"LineWidth",linewidth);hold on
% semilogy(EbN0_db_AE,ber_AE,"Marker", 'p',"Markersize",markersize,"LineStyle","-","Color",[0.5 0.4 0.1],"LineWidth",linewidth);hold on
% semilogy(EbN0_db_ATSC,ber_ATSC,"Marker", 'square',"MarkerSize",markersize,"LineStyle","-","Color",[0.1 0.1 0.8],"LineWidth",linewidth);hold on
% legend("BER-QAM",'BER-AE','BER-ATSC',"Location","southwest")
% grid on
% set(gca,'FontSize',14,'FontName','Times New Roman');
% set(gca,'ytick',[1e-4,1e-3,1e-2,1e-1,1],'ygrid','on','gridlinestyle','-','Gridalpha',0.1);
% ylim([1e-4,1])
% 
% hold off