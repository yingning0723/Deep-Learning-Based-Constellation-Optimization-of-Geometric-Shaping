% function : 产生信息矩阵并进行加噪传输，计算误比特数
function n_errbit = transmit(msg_orig,msg_tx,num_block,length_block,num_maxARQ,snr,n,k)
%                            原始信号,编码调制后信号,包数,   包长,    最大重传数,信噪比,汉明码参数
%预分配内存
msg_rx = zeros(size(msg_tx));%高斯后信号
msg_demod = zeros(size(msg_tx));%BPSK判决结果
msg_decode = zeros(num_block, length_block);%汉明解码结果
msg_decode_buffer = zeros(num_maxARQ,length_block);
 
%加噪,ARQ过程
    for row_tx = 1:num_block
        num_ARQ = 1;
        num_errbit = 0;
         %信道传输一个分组（加高斯）
        msg_rx(row_tx,:) = awgn(msg_tx(row_tx,:),snr,'measured');     
        %BPSK解调比特判决
        msg_demod(row_tx,:) = msg_rx(row_tx,:) > 0;     
         %解码
        msg_decode(row_tx,:) = decode(msg_demod(row_tx,:),n,k,'hamming');    
         %统计误比特
        num_errbit = sum((msg_decode(row_tx,:) - msg_orig(row_tx,:)) ~= 0);   
        %错误数据放入缓存
%         ffy 第i次重传的解调译码后的数据放入一个缓存中
        msg_decode_buffer(num_ARQ,:) = msg_decode(row_tx,:);                
         
        %有错重传        
        while (num_ARQ < num_maxARQ)&&(num_errbit ~= 0)     
            num_errbit = 0;  
            num_ARQ = num_ARQ + 1;
            msg_rx(row_tx,:) = awgn(msg_tx(row_tx,:),snr,'measured');                %信道传输一个分组
            msg_demod(row_tx,:) = msg_rx(row_tx,:) > 0;                              %BPSK解调比特判决
            msg_decode(row_tx,:) = decode(msg_demod(row_tx,:),n,k,'hamming');        %解码
            num_errbit = sum((msg_decode(row_tx,:) - msg_orig(row_tx,:)) ~= 0);      %统计误比特
            msg_decode_buffer(num_ARQ,:) = msg_decode(row_tx,:);                     %错误数据放入缓存
        end;
         
        %HARQ好像没有处理超过最大出错的情况？？
        %这段代码应该不需要。
        %最后一次还出错，则采用出现概率较大的符号为判决符号
        if num_ARQ == num_maxARQ  &&  num_maxARQ ~= 1   
             %对该分组的错误接收信息的各比特位出现概率进行统计
            for loop = 2:num_maxARQ                                      
                msg_decode_buffer(1,:) = msg_decode_buffer(1,:) + msg_decode_buffer(loop,:);
            end;
             %num_maxARQ为奇数时
            if mod(num_maxARQ,2) == 1                               
                for j = 1:length_block
                    msg_decode(row_tx,j) = msg_decode_buffer(1,j) > (num_maxARQ - 1)/2;
                end;
                %num_maxARQ为偶数时
            elseif mod(num_maxARQ,2) == 0                            
                for j = 1:length_block
                    if msg_decode_buffer(1,j) > num_maxARQ/2
                        msg_decode(row_tx,j) = 1;
                    elseif msg_decode_buffer(1,j) < num_maxARQ/2
                        msg_decode(row_tx,j) = 0;
                    else
                        msg_decode(row_tx,j) = randsrc(1,1,[0,1]);
                    end;
                end;
            end;
        end;
    end;
    n_errbit = sum(sum((msg_decode - msg_orig ) ~= 0));