function fish_RDH_JPEG_Decoder(msg_org)
% clear 

jpegOrder = [ 1 9 2 3 10 17 25 18   11 4 5 12 19 26 33 41   34 27 20 13 6 7 14 21   28 35 42 49 57 50 43 36,...
    29 22 15 8 16 23 30 37   44 51 58 59 52 45 38 31   24 32 39 46 53 60 61 54    47 40 48 55 62 63 56 64];
InDir = 'D:\Image_P\huang jiao shou tiaoshi JPEG\';
JpgFileList=dir([InDir, '*.jpg']);
ListLenCover = length(JpgFileList);

BN = 1; QF = 70;
coverJPEG = ('fishRDH.jpg');
jpgObj = jpeg_read(coverJPEG);                                           %read the coefficients of JPEG image
jpgCoef = jpgObj.coef_arrays{1};
[row,col] = size(jpgCoef);
jpgVecCoef = im2vec(jpgCoef,[8,8]);
jpgEmbCoef = jpgVecCoef(jpegOrder(1:64),:);                                % the 1:64 frequency bands are selected for embeding
quantize_table=jpgObj.quant_tables{1};
quantize_table_reshape=quantize_table(jpegOrder);
S_dec = jpgEmbCoef(2:64,:);

max_msLeng = log2(row*col);
side_info = max_msLeng + 63;
indexes = mod(jpgEmbCoef(1,1:side_info),2);
msgLen = bi2de(indexes(64:end));
tot_msgLen = msgLen + side_info;

a = 1;
for k=1:63
    if indexes(k) == 1
        ind_dec(a)= k;
        a = a+1;
    end
end

VG=S_dec(1:end,:)==0;
for i = 1:size(VG,2)
    VG2(:,i) = VG(:,i).* ((63:-1:1)');
    %             VG2(:,i) = VG(:,i);
end
VG2_sum=sum(VG2);
% tabulate(VG2_sum);
[~, VG_sort_ind_dec]=sort(VG2_sum,'descend');
S_sorted_dec=S_dec(:,VG_sort_ind_dec);


        for d = 1:length(ind_dec)
            block = find(S_sorted_dec(ind_dec(d),:) ~= 0);
            block2 = S_sorted_dec(ind_dec(d),block);
            block3 = [block2 zeros(1, 4096 - size(block2,2))];
            block_final(ind_dec(d),:) = block3;

        end
        block_final_dec = block_final(ind_dec,:);
        blocks_embedded_dec = reshape(block_final_dec,1,[]);

counter = 0; siz = size(blocks_embedded_dec,2);
for e = 1:siz
    if counter >= tot_msgLen break, end
    if blocks_embedded_dec(e) > 2
        blocks_embedded_dec(e) = blocks_embedded_dec(e) - 1;
    elseif blocks_embedded_dec(e) < -2
        blocks_embedded_dec(e) = blocks_embedded_dec(e) + 1;
    elseif abs(blocks_embedded_dec(e)) == 1
        counter = counter+1;
        mess(counter) = 0;
    elseif blocks_embedded_dec(e) == -2
        counter = counter+1;
        mess(counter) = 1;
        blocks_embedded_dec(e) = blocks_embedded_dec(e) + 1;
    elseif blocks_embedded_dec(e) == 2
        counter = counter+1;
        mess(counter) = 1;
        blocks_embedded_dec(e) = blocks_embedded_dec(e) - 1;
    end
end

        block_final_dec = reshape(blocks_embedded_dec,length(ind_dec),[]);
        for i1 = 1:length(ind_dec)
            block = find(S_sorted_dec(ind_dec(i1),:) ~= 0);
            S_sorted_dec(ind_dec(i1),block) = block_final_dec(i1,1:length(block));
        end
        S_dec(:,VG_sort_ind_dec) = S_sorted_dec;
        dcc_lsb_dec = mess(1:side_info);
        msg_dec = mess(side_info+1:end);
                
        dcc_dec = jpgEmbCoef(1,:);
        for s = 1:side_info;
           if mod(dcc_dec(s),2) ~= dcc_lsb_dec(s)
               dcc_dec(s) = dcc_dec(s) - 1;
           end 
        end

        stegoVecCoef = jpgVecCoef;
        stegoVecCoef(jpegOrder(2:64),:) = S_dec;
        stegoVecCoef(1,:) = dcc_dec;
        stegoCoef = vec2im(stegoVecCoef,0, 8, row/8, col/8);
        
        load(strcat('default_gray_jpeg_obj_', num2str(QF), '.mat'));
        stegoObj = default_gray_jpeg_obj;
        stegoObj.coef_arrays{1} = stegoCoef;
        
        jpeg_write(stegoObj,'fishRDH_decoded.jpg');
        check_message = sum(msg_org - msg_dec)
end
% keyboard