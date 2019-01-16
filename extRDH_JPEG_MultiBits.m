% function extRDH_JPEG_MultiBits(coverJPEG,stegoJPEG,extStegoJPEG,msgLen,Th,Tz,key,QF)

coverJPEG = 'lena70.jpg'; stegoJPEG = 'lenaRDH.jpg'; extStegoJPEG = 'extLenaRDH.jpg'; 
BN = 2;  msgLen = floor(40016/BN);   key = 100;  Tz = 34; QF = 70;

jpegOrder = [ 1 9 2 3 10 17 25 18   11 4 5 12 19 26 33 41   34 27 20 13 6 7 14 21   28 35 42 49 57 50 43 36,...
             29 22 15 8 16 23 30 37   44 51 58 59 52 45 38 31   24 32 39 46 53 60 61 54    47 40 48 55 62 63 56 64];
          
jpgObj = jpeg_read(stegoJPEG);                                           %read the coefficients of JPEG image
jpgCoef = jpgObj.coef_arrays{1};
[row,col] = size(jpgCoef);
jpgVecCoef = im2vec(jpgCoef,[8,8]);                                                                                                  
jpgEmbCoef = jpgVecCoef(jpegOrder(1:64),:);                                % the 1:64 frequency bands are selected for embeding
% diffJpgEmbCoef = diff(jpgEmbCoef,1,2);
% diffJpgEmbCoef = jpgEmbCoef(:,1:4095)-jpgEmbCoef(:,2:4096);


S = jpgEmbCoef(2:64,:);
[row_1,col_1] = size(S);
Hz = zeros(col_1,1);                 %save the number of zeros in each column
for k = 1:col_1
    Hz(k) = length(find(S(1:63,k)==0));
end

extMsg = zeros(1,msgLen);

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
            if S(posR(i),posC(j))>2^BN
                S(posR(i),posC(j)) = S(posR(i),posC(j))-(2^BN-1);
            elseif S(posR(i),posC(j))<-2^BN
                S(posR(i),posC(j)) = S(posR(i),posC(j))+(2^BN-1);
            elseif S(posR(i),posC(j))>=1 & S(posR(i),posC(j))<=2^BN
                msgNum = msgNum+1;
                extMsg(msgNum)=S(posR(i),posC(j))-1;
                S(posR(i),posC(j))=1;    
                if msgNum>=msgLen break; end
            elseif S(posR(i),posC(j))<=-1 & S(posR(i),posC(j))>=-2^BN
                msgNum = msgNum+1;
                extMsg(msgNum)=abs(S(posR(i),posC(j))+1);   
                S(posR(i),posC(j))=-1;      
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

jpeg_write(stegoObj,extStegoJPEG);

bseed = key;
rand('state',bseed);    
msg = randint(1,msgLen,[0 2^BN-1]);
disp('extractMessageDiff')
sum(sum(abs(msg-extMsg)))

diff = imread(extStegoJPEG) - imread(coverJPEG);
x = sum(sum(diff.*diff));
psnr = 10*log10(512*512*255*255/x)