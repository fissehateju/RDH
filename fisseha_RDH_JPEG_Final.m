% written by FIsseha Teju
% this is Experimental prove for the paper nameed "Reversible data hiding on JPEG images based on new coefficients selection strategy"

clear all

stegoJPEG = 'fishRDH.jpg';  BN = 1;

jpegOrder = [ 1 9 2 3 10 17 25 18,...   
        11 4 5 12 19 26 33 41,...
        34 27 20 13 6 7 14 21,...   
        28 35 42 49 57 50 43 36,...
        29 22 15 8 16 23 30 37,...   
        44 51 58 59 52 45 38 31,...   
        24 32 39 46 53 60 61 54,...    
        47 40 48 55 62 63 56 64];

InDir = '.\images\';  % Cover image location: makes sure it is created having JPEG images inside
JpgFileList=dir([InDir, '*.jpg']);
ListLenCover = length(JpgFileList);

result_fish_HS_PSNR  = zeros(ListLenCover,50);
result_fish_HS_filesize  = zeros(ListLenCover,50);
org_filesize = zeros(ListLenCover,50);

% for QF = 70:10:90
 for im = 1:ListLenCover
    coverJPEG = [InDir JpgFileList(im).name];

%     if QF == 70 && im == 3
%          EC = 42;
%     elseif QF == 70 && im == 4
%          EC = 25;
%     elseif QF == 70 && im == 5
%          EC = 21;
%     elseif QF == 70 && im == 6
%          EC = 29;
%     else
%          EC = 19;
%     end

    a = 1;
    EC = 8;
    for L = 8000:8000:EC*1000
        
        jpgObj = jpeg_read(coverJPEG);                                           %read the coefficients of JPEG image
        jpgCoef = jpgObj.coef_arrays{1};
        [row,col] = size(jpgCoef);
        jpgVecCoef = im2vec(jpgCoef,[8,8]);
        jpgEmbCoef = jpgVecCoef(jpegOrder(1:64),:);                                % the 1:64 frequency bands are selected for embeding
        quantize_table=jpgObj.quant_tables{1};
        quantize_table_reshape=quantize_table(jpegOrder);
        ac_coef = jpgEmbCoef(2:64,:);

        %% calculating scale factor and determine quality factor
        % Scale factor (SF) = 5000/ QF if 1 <= QF < 50
        % Scale factor (SF) = 200 - 2 * QF if 20 <= QF < 100
        % else Scale factor (SF) = 1
        
        QF50 = [16 11 10 16 24 40 51 61;
            12 12 14 19 26 58 60 55;
            14 13 16 24 40 57 69 56;
            14 17 22 29 51 87 80 62;
            18 22 37 56 68 109 103 77;
            24 35 55 64 81 104 113 92;
            49 64 78 87 103 121 120 101;
            72 92 95 98 112 100 103 99];

        for k = 10:10:90
            if k > 50
                SF = 200 - 2 * k;
            else
                SF = double(uint8(5000 / k));
            end
            
            Q = double(uint8(QF50 * (SF/100)));
            if Q == quantize_table
                QF = k;
                break;
            end
        end
        %%
        
        msgLen = L;   %%% message size L
        [row_1,col_1] = size(ac_coef);
        
%         absS = abs(ac_coef);
        
       %%  %%%%%%%%%%% block sorting    %%%%%%%%%%%%%% 
        zero_ind = ac_coef(1:end,:) == 0;
        zero_ind2 = zero_ind;
        
        for i = 1:size(zero_ind,2)
            zero_ind2(:,i) = zero_ind(:,i).* ((63:-1:1)');
        end
        
        zero_ind_sum = sum(zero_ind2);
        % tabulate(VG2_sum);
        [~, block_sorting_ind] = sort(zero_ind_sum,'descend');
        ac_coef_block_sorted = ac_coef(:,block_sorting_ind);
        
%         absS = abs(ac_coef_sorted_block);
        %%
        
        %%  %%%%%%%%%%%%%%%  index ording  %%%%%%%%%%%%%%%%%
        for i=1:63
            emb = ac_coef_block_sorted(i,:) == 1 | ac_coef_block_sorted(i,:) == -1;
            shift = ac_coef_block_sorted(i,:) > 1 | ac_coef_block_sorted(i,:) < -1;
%             D_ratio = sum(emb)/((sum(shift) + sum(emb)/2)/quantize_table_reshape(i + 1));
%             if ((sum(shift) + sum(emb)/2)/quantize_table_reshape(i + 1)) == 0
%                 D_ratio = 0;
%             end
%             ED_metric(i,:) = [ sum(emb) sum(shift) quantize_table_reshape(i + 1) D_ratio];

            ED_metric_cumsum_embedding(i,:) = cumsum(emb);
            ED_metric_cumsum_distortion(i,:) = 1 + quantize_table_reshape(i + 1)^2*(cumsum(shift)+0.5*cumsum(emb));
            ED_metric_cumsum_final(i,:) = ED_metric_cumsum_embedding(i,:)./ED_metric_cumsum_distortion(i,:);
        end
        
        mean_ED_metric_cumsum_final = mean(ED_metric_cumsum_final')';
        [sorted_mean_ED_metric_cumsum_final ind] = sortrows(mean_ED_metric_cumsum_final,-1);
        temp_variable = [sorted_mean_ED_metric_cumsum_final ED_metric_cumsum_embedding(ind,end) quantize_table_reshape(ind+1)'  ind];
        %%
              
        %%  %%%%%%%%%%%%%% section selection and embedding %%%%%%%%%%%%%%%%%%%
        %         counter=0;
        for d = 1:63
            %                 counter=counter+1;
            nonzeros_ind = find(ac_coef_block_sorted(ind(d),:) ~= 0);
            nonzeros = ac_coef_block_sorted(ind(d),nonzeros_ind);
            padding = [nonzeros zeros(1, 4096 - size(nonzeros,2))];
            nonzeros_plus_padding(ind(d),:) = padding;
        end
        
        nonzeros_padding_save = [];
        distortion_min_save=10000000000000000000000;
        
        for d = 1:63
            section = sort(ind(1:d));
            blocks_to_be_embedded_original = reshape(nonzeros_plus_padding(section,:),1,[]);
            blocks_to_be_embedded = reshape(nonzeros_plus_padding(section,:),1,[]);
            
            num_Emb = sum(sum(blocks_to_be_embedded == -1 | blocks_to_be_embedded == 1));

            %Embed blocks_to_be_embedded
            max_msLeng = log2(row*col);
            siz = max_msLeng + 63; % 63 bits will be used to record the ordered indexes used for embedding
            dc_coef_lsb = mod(jpgEmbCoef(1,1:siz),2); % looking for embedding area for side information
            
            rand('state', 100);
            msg_org = randi([0 2^BN-1],1,msgLen);
            msg = [dc_coef_lsb msg_org];
            tot_msgLen = length(msg);
            counter = 0; siz = size(blocks_to_be_embedded,2);
            for e = 1:siz
                if counter >= tot_msgLen, break, end
                if blocks_to_be_embedded(e) > 1
                    blocks_to_be_embedded(e) = blocks_to_be_embedded(e) + (2^BN-1);
                elseif blocks_to_be_embedded(e) < -1
                    blocks_to_be_embedded(e) = blocks_to_be_embedded(e) - (2^BN-1);
                elseif blocks_to_be_embedded(e) == 1
                    counter = counter+1;
                    blocks_to_be_embedded(e) = blocks_to_be_embedded(e) + msg(counter);
                elseif blocks_to_be_embedded(e) == -1
                    counter = counter+1;
                    blocks_to_be_embedded(e) = blocks_to_be_embedded(e) - msg(counter);
                end
            end
           
           distortion(d)=sum(quantize_table_reshape(ind(1:d)+1).^2*reshape(abs(blocks_to_be_embedded_original- blocks_to_be_embedded),d,[]),2);  % modified code
            %Flag for if all the payload has been embedded
            if num_Emb < tot_msgLen
                flag(d) = 0;
            else
                flag(d) = 1;
            end
            
            if flag(d)==1 && distortion_min_save > distortion(d)
                distortion_min_save=distortion(d);
                block_final_save=[];
                block_final_save=blocks_to_be_embedded;
            end
        end
        if sum(flag)==0
            %             pause
            disp('cannot embed,')
        end
        max_distortion=2*max(distortion);
        distortion(flag==0)=max_distortion;
        
        %Take the lowest distortion with the flag==1
        [~,min_index]=min(distortion);
        position = sort(ind(1:min_index));
        %%
        
        indexes = zeros(1,63);
        indexes(position) = 1;
        max_msgLen = log2(row*col);
        max_msLeng_bit = zeros(1,max_msgLen);
        message_L = de2bi(msgLen);
        max_msLeng_bit(1:length(message_L)) = message_L;
        side_info = [indexes max_msLeng_bit];
        
        dcc = jpgEmbCoef(1,:);
        dcc_mod = dcc;
        for s = 1:length(side_info)
               dcc_lsb(1,s) = mod(dcc_mod(s),2);
           if mod(dcc_mod(s),2) ~= side_info(s)
               dcc_mod(s) = dcc_mod(s) + 1;
           end 
        end
        
        %%
        
        %%% local decode
        block_final_mod = reshape(block_final_save,min_index,[]);
        ac_coef_block_sorted_mod = ac_coef_block_sorted;
        for i1 = 1:length(position)
            block = find(ac_coef_block_sorted_mod(position(i1),:) ~= 0);
            ac_coef_block_sorted_mod(position(i1),block) = block_final_mod(i1,1:length(block));
        end
        
        ac_coef_mod = zeros(size(ac_coef));
        ac_coef_mod(:,block_sorting_ind) = ac_coef_block_sorted_mod;

        stegoVecCoef = jpgVecCoef;
        stegoVecCoef(jpegOrder(2:64),:) = ac_coef_mod;
        stegoVecCoef(1,:) = dcc_mod;
        stegoCoef = vec2im(stegoVecCoef,0, 8, row/8, col/8);
        
        load(strcat('default_gray_jpeg_obj_', num2str(QF), '.mat'));
        stegoObj = default_gray_jpeg_obj;
        stegoObj.coef_arrays{1} = stegoCoef;
        
        jpeg_write(stegoObj,stegoJPEG);
        temp = dir(stegoJPEG);
        fileSize = temp.bytes;
        temp2 = dir(coverJPEG);
        ofilesize = temp2.bytes;
        diff = double(imread(coverJPEG)) - double(imread(stegoJPEG));
        x = sum(sum(diff.*diff));
        PSNR = 10*log10(512*512*255*255/x);
        
        result_fish_HS_PSNR(im,a)  = PSNR;
        result_fish_HS_filesize(im,a)  = fileSize;
        org_filesize(im,a) = ofilesize;
        
        a =a+1;
        
        %% %%%%%%%%%% decoder %%%%%%%%%%%%%
        
       fisseha_RDH_JPEG_Decoder(msg_org) 

    end

 end
save result_original_filesize org_filesize
save result_fish_HS_final result_fish_HS_PSNR result_fish_HS_filesize
