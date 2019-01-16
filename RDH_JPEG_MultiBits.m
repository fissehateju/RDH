% function [PSNR fileSize EC] = RDH_JPEG_MultiBits(coverJPEG,stegoJPEG,msgLen,QF,Th,key)

coverJPEG = 'lena70.jpg';  BN = 1;  msgLen = floor(40000/BN);  stegoJPEG = 'huangRDH.jpg';  key = 100;  QF = 70;


jpegOrder = [ 1 9 2 3 10 17 25 18   11 4 5 12 19 26 33 41   34 27 20 13 6 7 14 21   28 35 42 49 57 50 43 36,...
             29 22 15 8 16 23 30 37   44 51 58 59 52 45 38 31   24 32 39 46 53 60 61 54    47 40 48 55 62 63 56 64];       
InDir = 'D:\Image_P\huang jiao shou tiaoshi JPEG\images\';
JpgFileList=dir([InDir, '*.jpg']);
ListLenCover = length(JpgFileList);

result_huang_HS_PSNR = zeros(ListLenCover,50);
result_huang_HS_filesize = zeros(ListLenCover,50);

% for QF = 60:10:90
for im=1:ListLenCover
%     coverJPEG = (['compfish_lena', num2str(QF), '.jpg']); 
    coverJPEG = [InDir JpgFileList(im).name];
%     if QF == 50
%         n = 1; EC = 14;
%     elseif QF == 60
%         n = 1; EC = 16;
%     elseif QF == 70
%         n = 1; EC = 18;
%     elseif QF == 80
%         n = 2; EC = 24;
%     elseif QF == 90
%         n = 3; EC = 36;
%     end

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
    for L = 1000:1000:EC*1000
        
jpgObj = jpeg_read(coverJPEG);                                           %read the coefficients of JPEG image
jpgCoef = jpgObj.coef_arrays{1};
[row,col] = size(jpgCoef);
jpgVecCoef = im2vec(jpgCoef,[8,8]);                                                                                                  
jpgEmbCoef = jpgVecCoef(jpegOrder(1:64),:);                                % the 1:64 frequency bands are selected for embeding

S = jpgEmbCoef(2:64,:);
EC = (sum(S(:)==1)+sum(S(:)==-1));
msgLen = EC;
msgLen = L;   %%% message size L
[row_1,col_1] = size(S);
absS = abs(S); 
% maxS = max(max(absS)); 
% Ht = zeros(1+maxS,1);
% for i = 0:maxS
%     Ht(i+1) = length(find(absS(:)==i));
% end
% for i = maxS+1:-1:1
%     if Ht(i)>=msgLen
%         BN = i-1;
%         break;
%     end
% end

Hz = zeros(col_1,1);                 %save the number of zeros in each column
for k = 1:col_1
    Hz(k) = length(find(S(1:63,k)==0));
end
for k = 63:-1:0                      %the number of zeros in each block belongs to [1,63]
    posZ = find(Hz>=k);
    num = sum(sum(absS(:,posZ)==1));
    if num>=msgLen
        Tz = k;
        break;
    end
end

bseed = key;
rand('state',bseed);    
msg = randi([0 2^BN-1],1,msgLen);

bseed = key;
rand('state',bseed);    
posC = randperm(col_1);
bseed = key;
rand('state',bseed);    
posR = randperm(row_1);

msgNum = 0;
for i = 1:row_1
    if msgNum>=msgLen break; end
    for j = 1:col_1
        if Hz(posC(j))>=Tz  
            if S(posR(i),posC(j))>1
                S(posR(i),posC(j)) = S(posR(i),posC(j))+(2^BN-1);
            elseif S(posR(i),posC(j))<-1
                S(posR(i),posC(j)) = S(posR(i),posC(j))-(2^BN-1);
            elseif S(posR(i),posC(j))==1
                msgNum = msgNum+1;
                S(posR(i),posC(j))=S(posR(i),posC(j))+msg(msgNum);
                if msgNum>=msgLen break; end
            elseif S(posR(i),posC(j))==-1
                msgNum = msgNum+1;
                S(posR(i),posC(j))=S(posR(i),posC(j))-msg(msgNum);
                if msgNum>=msgLen break; end
            end        
        end
    end
end

stegoVecCoef = jpgVecCoef;
stegoVecCoef(jpegOrder(2:64),:) = S;    
stegoCoef = vec2im(stegoVecCoef,0, 8, row/8, col/8);

load(strcat('default_gray_jpeg_obj_', num2str(QF), '.mat'));
stegoObj = default_gray_jpeg_obj;
stegoObj.coef_arrays{1} = stegoCoef;

jpeg_write(stegoObj,stegoJPEG);
temp = dir(stegoJPEG);
fileSize = temp.bytes;
% temp2 = dir(coverJPEG);
% ofilesize = temp2.bytes
diff = double(imread(coverJPEG)) - double(imread(stegoJPEG));
x = sum(sum(diff.*diff));
PSNR = 10*log10(512*512*255*255/x);

result_huang_HS_PSNR(im,a)  = PSNR;
result_huang_HS_filesize(im,a)  = fileSize;
a =a+1;
    end
% end
% res_huang = (['res_huang', num2str(QF)]);
% save res_huang60 huang_PSNR huang_filesize
end
save result_huang_HS_final result_huang_HS_PSNR result_huang_HS_filesize

