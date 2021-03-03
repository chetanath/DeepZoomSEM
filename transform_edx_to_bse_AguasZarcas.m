clearvars
close all

dzibsefile='AguasZarcasFeb2021_bse.png';
originalbsefile='~/Data/deepzoom/AguasZarcas.png';
samplenameedx=dzibsefile(1:end-8);
distorted=imread(dzibsefile);
[~,samplename,~]=fileparts(originalbsefile);
original=imread(originalbsefile);
original=imresize(original,size(distorted,1)/size(original,1));

mincontrastvalue=0.05; % =0.05; default = 0.2
minqualityvalue=0.1;

ptsOriginal  = detectFASTFeatures(original,'MinQuality',minqualityvalue,'MinContrast',mincontrastvalue);
ptsDistorted = detectFASTFeatures(distorted,'MinQuality',minqualityvalue,'MinContrast',mincontrastvalue);

[featuresOriginal,  validPtsOriginal]  = extractFeatures(original,  ptsOriginal);
[featuresDistorted, validPtsDistorted] = extractFeatures(distorted, ptsDistorted);

indexPairs = matchFeatures(featuresOriginal, featuresDistorted, 'Unique', true, 'Method', 'Approximate');

matchedOriginal  = validPtsOriginal(indexPairs(:,1));
matchedDistorted = validPtsDistorted(indexPairs(:,2));

figure;
showMatchedFeatures(original,distorted,matchedOriginal,matchedDistorted);
title('Putatively matched points (including outliers)');

[tform, inlierDistorted,inlierOriginal] = estimateGeometricTransform(...
    matchedDistorted, matchedOriginal, 'projective','MaxDistance',30, 'Confidence', 95, 'MaxNumTrials', 1e4);
% inlierDistorted = matchedDistorted(inlierIdx, :);
% inlierOriginal  = matchedOriginal(inlierIdx, :);

figure;
showMatchedFeatures(original,distorted,inlierOriginal,inlierDistorted);
title('Matching points (inliers only)');
legend('ptsOriginal','ptsDistorted');

outputView = imref2d(size(original));
recovered  = imwarp(distorted,tform,'OutputView',outputView);

figure, imshowpair(imresize(original,3),imresize(recovered,3),'falsecolor')

samplefolder='~/Documents/MATLAB/';
%%

elnames={'Fe','Ca','Mg'};
if exist([samplefolder samplenameedx '_' elnames{1} '_Map.png'],'file')~=0 && exist([samplefolder samplenameedx '_' elnames{2} '_Map.png'],'file')~=0 && exist([samplefolder samplenameedx '_' elnames{3} '_Map.png'],'file')~=0

El1=imread([samplefolder samplenameedx '_' elnames{1} '_Map.png']);
El2=imread([samplefolder samplenameedx '_' elnames{2} '_Map.png']);
El3=imread([samplefolder samplenameedx '_' elnames{3} '_Map.png']);


El1=(imwarp(El1,tform,'OutputView',outputView));
El2=(imwarp(El2,tform,'OutputView',outputView));
El3=(imwarp(El3,tform,'OutputView',outputView));

g1=1;
g2=1;
g3=1;
%RGB=cat(3,double(El1).^g1./max(max(double(El1).^g1)),double(El2).^g2./max(max(double(El2).^g2)),double(El3).^g3./max(max(double(El3).^g3)));
RGB=cat(3,El1, El2, El3);


%RGB = decorrstretch(RGB);
image(RGB)

imwrite(RGB,'RGB.png');

if exist(['dzi/' samplename '_' elnames{1} '_' elnames{2} '_' elnames{3} '_files'],'dir')==7
    
    [~,~]=system(['rm -rf dzi/' samplename '_' elnames{1} '_' elnames{2} '_' elnames{3} '_files']);
    [~,~]=system(['rm -rf dzi/' samplename '_' elnames{1} '_' elnames{2} '_' elnames{3} '.dzi']);
    
end

[~,~]=system(['vips dzsave RGB.png dzi/' samplename '_' elnames{1} '_' elnames{2} '_' elnames{3} ' --suffix .jpg[Q=100]']);
end

%%

elnames={'S','Si','O'};
if exist([samplefolder samplenameedx '_' elnames{1} '_Map.png'],'file')~=0 && exist([samplefolder samplenameedx '_' elnames{2} '_Map.png'],'file')~=0 && exist([samplefolder samplenameedx '_' elnames{3} '_Map.png'],'file')~=0
    
El1=adapthisteq(imread([samplefolder samplenameedx '_' elnames{1} '_Map.png']));
El2=adapthisteq(imread([samplefolder samplenameedx '_' elnames{2} '_Map.png']));
El3=adapthisteq(imread([samplefolder samplenameedx '_' elnames{3} '_Map.png']));


El1=imwarp(El1,tform,'OutputView',outputView);
El2=imwarp(El2,tform,'OutputView',outputView);
El3=imwarp(El3,tform,'OutputView',outputView);

g1=1;
g2=1;
g3=1;
%RGB=cat(3,double(El1).^g1./max(max(double(El1).^g1)),double(El2).^g2./max(max(double(El2).^g2)),double(El3).^g3./max(max(double(El3).^g3)));
RGB=cat(3,El1, El2, El3);


%RGB = imadjustn(RGB);
image(RGB)

imwrite(RGB,'RGB.png');

if exist(['dzi/' samplename '_' elnames{1} '_' elnames{2} '_' elnames{3} '_files'],'dir')==7
    
    [~,~]=system(['rm -rf dzi/' samplename '_' elnames{1} '_' elnames{2} '_' elnames{3} '_files']);
    [~,~]=system(['rm -rf dzi/' samplename '_' elnames{1} '_' elnames{2} '_' elnames{3} '.dzi']);
    
end

[~,~]=system(['vips dzsave RGB.png dzi/' samplename '_' elnames{1} '_' elnames{2} '_' elnames{3} ' --suffix .jpg[Q=100]']);
end
%%

elnames={'Fe','S','O'};
if exist([samplefolder samplenameedx '_' elnames{1} '_Map.png'],'file')~=0 && exist([samplefolder samplenameedx '_' elnames{2} '_Map.png'],'file')~=0 && exist([samplefolder samplenameedx '_' elnames{3} '_Map.png'],'file')~=0

El1=imread([samplefolder samplenameedx '_' elnames{1} '_Map.png']);
El2=imread([samplefolder samplenameedx '_' elnames{2} '_Map.png']);
El3=imread([samplefolder samplenameedx '_' elnames{3} '_Map.png']);


El1=(imwarp(El1,tform,'OutputView',outputView));
El2=(imwarp(El2,tform,'OutputView',outputView));
El3=(imwarp(El3,tform,'OutputView',outputView));

g1=1;
g2=1;
g3=1;
%RGB=cat(3,double(El1).^g1./max(max(double(El1).^g1)),double(El2).^g2./max(max(double(El2).^g2)),double(El3).^g3./max(max(double(El3).^g3)));
RGB=cat(3,El1, El2, El3);


%RGB = decorrstretch(RGB);
image(RGB)

imwrite(RGB,'RGB.png');


if exist(['dzi/' samplename '_' elnames{1} '_' elnames{2} '_' elnames{3} '_files'],'dir')==7
    
    [~,~]=system(['rm -rf dzi/' samplename '_' elnames{1} '_' elnames{2} '_' elnames{3} '_files']);
    [~,~]=system(['rm -rf dzi/' samplename '_' elnames{1} '_' elnames{2} '_' elnames{3} '.dzi']);
    
end

[~,~]=system(['vips dzsave RGB.png dzi/' samplename '_' elnames{1} '_' elnames{2} '_' elnames{3} ' --suffix .jpg[Q=100]']);
end
%%

elnames={'Ca','Al','Si'};
if exist([samplefolder samplenameedx '_' elnames{1} '_Map.png'],'file')~=0 && exist([samplefolder samplenameedx '_' elnames{2} '_Map.png'],'file')~=0 && exist([samplefolder samplenameedx '_' elnames{3} '_Map.png'],'file')~=0

El1=imread([samplefolder samplenameedx '_' elnames{1} '_Map.png']);
El2=imread([samplefolder samplenameedx '_' elnames{2} '_Map.png']);
El3=imread([samplefolder samplenameedx '_' elnames{3} '_Map.png']);


El1=imadjust(imwarp(El1,tform,'OutputView',outputView));
El2=imadjust(imwarp(El2,tform,'OutputView',outputView));
El3=imadjust(imwarp(El3,tform,'OutputView',outputView));

% g1=1;
% g2=1;
% g3=1;
%RGB=cat(3,double(El1).^g1./max(max(double(El1).^g1)),double(El2).^g2./max(max(double(El2).^g2)),double(El3).^g3./max(max(double(El3).^g3)));
RGB=cat(3,El1, El2, El3);


%RGB = decorrstretch(RGB);
image(RGB)

imwrite(RGB,'RGB.png');


if exist(['dzi/' samplename '_' elnames{1} '_' elnames{2} '_' elnames{3} '_files'],'dir')==7
    
    [~,~]=system(['rm -rf dzi/' samplename '_' elnames{1} '_' elnames{2} '_' elnames{3} '_files']);
    [~,~]=system(['rm -rf dzi/' samplename '_' elnames{1} '_' elnames{2} '_' elnames{3} '.dzi']);
    
end

[~,~]=system(['vips dzsave RGB.png dzi/' samplename '_' elnames{1} '_' elnames{2} '_' elnames{3} ' --suffix .jpg[Q=100]']);
end


%%
D=dir([samplefolder samplename '*Map.png']);

% For AZ: convert to lower case with "lower

for nn=1:numel(D)

    %El = extractBetween(D(nn).name,[samplefolder samplename '_'],'_Map.png');
    
    El = extractBetween(D(nn).name,[samplenameedx '_'],'_Map.png');
    
    disp(El{1})
    
    Fe=imread([samplefolder D(nn).name]);
    
    Fe_warp=adapthisteq(imwarp(Fe,tform,'OutputView',outputView));
    
    %Fe_warp=imwarp(Fe,tform,'OutputView',outputView);
    
    pxx=size(Fe_warp,1);
    pxy=size(Fe_warp,2);
    
    bits=8;
    
    pp1=gray(2^bits);
    pp2=inferno(2^bits);
    pp3=magma(2^bits);
    pp4=parula(2^bits);
    pp5=viridis(2^bits);
    pp6=plasma(2^bits);
    pp7=jet(2^bits);
    
    Fe_warp=uint8(reshape(255*pp3(Fe_warp+1,:),pxx,pxy,3));
    
    imwrite(Fe_warp,[samplename '_Map_warp.png']);
    
    if exist(['dzi/' samplename '_' El{1} '_files'],'dir')==7
        
        [~,~]=system(['rm -rf dzi/' samplename '_' El{1} '_files']);
        [~,~]=system(['rm -rf dzi/' samplename '_' El{1} '.dzi']);

    end
    
    [~,~]=system(['vips dzsave ' samplename '_Map_warp.png dzi/' samplename '_' El{1} ' --suffix .jpg[Q=100]']);
end

