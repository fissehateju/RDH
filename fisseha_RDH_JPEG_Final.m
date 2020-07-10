% written by FIsseha Teju
% this is Experimental prove for the paper nameed "Reversible data hiding on JPEG images based on new coefficients selection strategy"

clear
coverJPEG = 'lena70.jpg';  BN = 1;  msgLen = floor(40000/BN);  stegoJPEG = 'fishRDH.jpg';  key = 100;  QF = 70; a = 1;

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
 for im=1:ListLenCover
    coverJPEG = [InDir JpgFileList(im).name];

    if QF == 70 && im == 3
         EC = 42;
    elseif QF == 70 && im == 4
         EC = 25;
    elseif QF == 70 && im == 5
         EC = 21;
    elseif QF == 70 && im == 6
         EC = 29;
    else
         EC = 19;
    end

    a = 1;
    EC = 10;
    for L = 10000:1000:EC*1000
        
        jpgObj = jpeg_read(coverJPEG);                                           %read the coefficients of JPEG image
        jpgCoef = jpgObj.coef_arrays{1};
        [row,col] = size(jpgCoef);
        jpgVecCoef = im2vec(jpgCoef,[8,8]);
        jpgEmbCoef = jpgVecCoef(jpegOrder(1:64),:);                                % the 1:64 frequency bands are selected for embeding
        quantize_table=jpgObj.quant_tables{1};
        quantize_table_reshape=quantize_table(jpegOrder);
        S = jpgEmbCoef(2:64,:);

        msgLen = L;   %%% message size L
        [row_1,col_1] = size(S);
        absS = abs(S);
        
       %%  %%%%%%%%%%% block sorting    %%%%%%%%%%%%%% 
        VG=S(1:end,:)==0;
        for i = 1:size(VG,2)
            VG2(:,i) = VG(:,i).* ((63:-1:1)');
            %             VG2(:,i) = VG(:,i);
        end
        VG2_sum=sum(VG2);
        % tabulate(VG2_sum);
        [~, VG_sort_ind]=sort(VG2_sum,'descend');
        S_sorted=S(:,VG_sort_ind);
        absS = abs(S_sorted);
        %%
        
        %%  %%%%%%%%%%%%%%%  index ording  %%%%%%%%%%%%%%%%%
        for i=1:63
            emb = S_sorted(i,:)==1|S_sorted(i,:)==-1;
            shift = S_sorted(i,:)>1|S_sorted(i,:)<-1;
            DM = sum(emb)/(sum(shift)+sum(emb)/2)/quantize_table_reshape(i+1);
            ED_metric(i,:)=[ sum(emb) sum(shift) quantize_table_reshape(i+1) DM];
            ED_metric_cumsum_embedding(i,:)=cumsum(emb);
            ED_metric_cumsum_distortion(i,:)=1+quantize_table_reshape(i+1)^2*(cumsum(shift)+0.5*cumsum(emb));
            ED_metric_cumsum_final(i,:)=ED_metric_cumsum_embedding(i,:)./ED_metric_cumsum_distortion(i,:);
        end
        
        mean_ED_metric_cumsum_final=mean(ED_metric_cumsum_final')';
        [sorted_mean_ED_metric_cumsum_final ind]=sortrows(mean_ED_metric_cumsum_final,-1);
        temp_variable=[sorted_mean_ED_metric_cumsum_final ED_metric_cumsum_embedding(ind,end) quantize_table_reshape(ind+1)'  ind];
        %%
              
        %%  %%%%%%%%%%%%%% section selection and embedding %%%%%%%%%%%%%%%%%%%
        %         counter=0;
        for d = 1:63
            %                 counter=counter+1;
            block = find(S_sorted(ind(d),:) ~= 0);
            block2 = S_sorted(ind(d),block);
            block3 = [block2 zeros(1, 4096 - size(block2,2))];
            block_final(ind(d),:) = block3;
        end
        block_final_save=[];
        distortion_min_save=10000000000000000000000;
        for d = 1:63
            sec = sort(ind(1:d));
            blocks_to_be_embedded_original=reshape(block_final(sec,:),1,[]);
            blocks_to_be_embedded=reshape(block_final(sec,:),1,[]);
            
            num_Emb = sum(sum(blocks_to_be_embedded == -1 | blocks_to_be_embedded == 1));
            
            %Embed blocks_to_be_embedded
            max_msLeng = log2(row*col);
            si = max_msLeng + 63;
            dcc_lsb = mod(jpgEmbCoef(1,1:si),2);
            
            bseed = key;
            rand('state',bseed);
            msg_org = randi([0 2^BN-1],1,msgLen);
            msg = [dcc_lsb msg_org];
            tot_msgLen = length(msg);
            counter = 0; siz = size(blocks_to_be_embedded,2);
            for e = 1:siz
                if counter >= tot_msgLen break, end
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
        for s = 1:length(side_info);
               dcc_lsb(1,s) = mod(dcc_mod(s),2);
           if mod(dcc_mod(s),2) ~= side_info(s)
               dcc_mod(s) = dcc_mod(s) + 1;
           end 
        end
        
        %%
        
        %%% local decode
        block_final_mod = reshape(block_final_save,min_index,[]);
        for i1 = 1:length(position)
            block = find(S_sorted(position(i1),:) ~= 0);
            S_sorted(position(i1),block) = block_final_mod(i1,1:length(block));
        end
        S(:,VG_sort_ind) = S_sorted;
        

        stegoVecCoef = jpgVecCoef;
        stegoVecCoef(jpegOrder(2:64),:) = S;
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
%         fish_RDH_JPEG_Decoder(msg_org)

 end
save result_original_filesize org_filesize
save result_fish_HS_final result_fish_HS_PSNR result_fish_HS_filesize

