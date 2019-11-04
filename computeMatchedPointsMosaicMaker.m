function [matchedcountqqkk,matchedvcountqqkk,matchedPointsqqkk,matchedPointsPrevqqkk,matchedPoints_vqqkk,matchedPointsPrev_vqqkk,brightnessvaluesqqkk,brightnessvaluesvqqkk,emission_current_ampsqqkk,matchmetricqqkk,matchmetricvqqkk] = computeMatchedPointsMosaicMaker(kk__,rows_array,cols_array,qq,imgext,imgfolder,colt_split,height,searchheight,width,searchwidth,dosavefigurematches,dosavefigurematchesv,dodisplayfeatures,minqualityvalue,mincontrastvalue,transformtype,transformtypev,mincountpoints,pointsAreColinear,searchmethod,maxratio)

kk=rows_array(kk__);

inlierOriginal=[];
inlierDistorted=[];
inlierDistortedv=[];
inlierOriginalv=[];

matchedPointsqqkk=cornerPoints;
matchedPointsPrevqqkk=cornerPoints;
matchedPoints_vqqkk=cornerPoints;
matchedPointsPrev_vqqkk=cornerPoints;

matchedvcountqqkk=0;
matchedcountqqkk=0;
matchmetricqqkk=nan;
matchmetricvqqkk=nan;
brightnessvaluesqqkk=[nan,nan];
brightnessvaluesvqqkk=[nan,nan];
emission_current_ampsqqkk=nan;

confidencelevel=90;
maxnumtrials=1e5;
maxdistance=3;

% imshow(I1); hold on; plot(points); hold off;

imname1=sprintf(['%03d_%03d.' imgext],kk-1,qq);
imname2=sprintf(['%03d_%03d.' imgext],kk,qq);

%fprintf('\n %s',imname2)

hdrname=[imgfolder imname2(1:end-4) '-' imname2(end-2:end) '.hdr'];
hdr_=readhdr(hdrname);

hdrv=hdr_.value;

emission_index=strcmp(hdr_.name,'EmissionCurrent');

emission_current_ampsqqkk=str2double(hdrv{emission_index});


searchROI_1=[1 1 searchwidth height];
searchROI_2=[width-searchwidth 1 searchwidth height];


if qq>colt_split
    qqp=qq-1;
    searchROI_1v=[1 height-searchheight width searchheight];
    searchROI_2v=[1 1 width searchheight];
else
    qqp=qq+1;
    searchROI_1v=[1 1 width searchheight];
    searchROI_2v=[1 height-searchheight width searchheight];
end

imname2v=sprintf(['%03d_%03d.' imgext],kk,qqp);

I2=imread([imgfolder imname2]);

if isfile([imgfolder imname1]) %&& kk>rows_array(1)
    
    
    I1=imread([imgfolder imname1]);
    pointsPrevious  = detectFASTFeatures(I1,'ROI',searchROI_1,'MinQuality',minqualityvalue,'MinContrast',mincontrastvalue); %,'MinContrast',0.01); %,'MinContrast',0.01); %detectSURFFeatures(I1,'ROI',[1 1 searchwidth height],'MetricThreshold',1,'NumScaleLevels',6,'NumOctaves',1);
    [featuresPrevious, pointsPrevious]  = extractFeatures(I1,  pointsPrevious, 'Method','Block');
    
    points = detectFASTFeatures(I2,'ROI',searchROI_2,'MinQuality',minqualityvalue,'MinContrast',mincontrastvalue); %,'MinContrast',0.01); %,'MetricThreshold',100,'NumScaleLevels',6);
    [features, points]  = extractFeatures(I2,  points, 'Method','Block');
    [indexPairs, matchmetric] = matchFeatures(features, featuresPrevious, 'Unique', ...
        true, 'Method', searchmethod, ...
        'MaxRatio', maxratio);
    
    if dodisplayfeatures
        fprintf('%d ',size(indexPairs,1))
    end
        
    if size(indexPairs,1)>0
        matchedPointsqqkk  = points(indexPairs(:,1));
        matchedPointsPrevqqkk = pointsPrevious(indexPairs(:,2));
        
        matchedcountqqkk=matchedPointsqqkk.Count;
        matchmetricqqkk=matchmetric;
        
    end
    
    
    
    
    
    if matchedcountqqkk>mincountpoints
        [~, inlierDistorted, inlierOriginal, ~] = ...
            estimateGeometricTransform(matchedPointsqqkk, ...
            matchedPointsPrevqqkk, transformtype, 'Confidence', confidencelevel, 'MaxNumTrials', maxnumtrials, 'MaxDistance', maxdistance); %, 'MaxDistance', round(height/100));
        if ~pointsAreColinear(inlierOriginal.Location) && ~pointsAreColinear(inlierDistorted.Location)
            try
                bbo=convhull(double(inlierOriginal.Location(:,1)),double(inlierOriginal.Location(:,2)));
                bbd=convhull(double(inlierDistorted.Location(:,1)),double(inlierDistorted.Location(:,2)));
                BW1=poly2mask(double(inlierOriginal.Location(bbo,1)),double(inlierOriginal.Location(bbo,2)),width,height);
                BW2=poly2mask(double(inlierDistorted.Location(bbd,1)),double(inlierDistorted.Location(bbd,2)),width,height);
                brightnessvaluesqqkk = [mean(double(I2(BW2))),mean(double(I1(BW1)))];
            catch
            end
        end
    end
    
    if dosavefigurematches && matchedcountqqkk>mincountpoints
        
        f=figure('Visible','off');
        showMatchedFeatures(I1,I2,inlierOriginal,inlierDistorted);
        title([imname1 ' : ' imname2],'Interpreter','None')
        saveas(f,[imname1(1:end-4) '-' imname2(1:end-4) '_match.png'])
        close(f);
        
        
    end
    
    
end


if isfile([imgfolder imname2v]) %&& qq~=cols_array(1)
    
    I1_=imread([imgfolder imname2v]);
    
    pointsv_ = detectFASTFeatures(I1_,'ROI',searchROI_1v,'MinQuality',minqualityvalue,'MinContrast',mincontrastvalue); %,'MinContrast',0.01); %,'MetricThreshold',100,'NumScaleLevels',6);
    [featuresv_, pointsv_]  = extractFeatures(I1_,  pointsv_, 'Method','Block');
    
    pointsv = detectFASTFeatures(I2,'ROI',searchROI_2v,'MinQuality',minqualityvalue,'MinContrast',mincontrastvalue); %,'MinContrast',0.01); %,'MetricThreshold',100,'NumScaleLevels',6);
    [featuresv, pointsv]  = extractFeatures(I2,  pointsv, 'Method','Block');
    
    %imshow(I1_); hold on; plot(pointsv_); hold off;
    %imshow(I1); hold on; plot(pointsv); hold off;
    
    
    [indexPairsv, matchmetricv] = matchFeatures(featuresv, featuresv_, 'Unique', ...
        true, 'Method', searchmethod, ...
        'MaxRatio', maxratio);
    
    
    if numel(indexPairsv)>0
        
        matchedPoints_vqqkk=pointsv(indexPairsv(:,1));
        matchedPointsPrev_vqqkk=pointsv_(indexPairsv(:,2));
        
        matchedvcountqqkk=matchedPointsPrev_vqqkk.Count;
        matchmetricvqqkk=matchmetricv;
        
    end
    
    
    if matchedvcountqqkk>mincountpoints
        
        [~, inlierDistortedv, inlierOriginalv, ~] = ...
            estimateGeometricTransform(matchedPoints_vqqkk, ...
            matchedPointsPrev_vqqkk, transformtypev, 'Confidence', confidencelevel, 'MaxNumTrials', maxnumtrials, 'MaxDistance', maxdistance); %, 'MaxDistance', round(height/100));
        if ~pointsAreColinear(inlierOriginalv.Location) && ~pointsAreColinear(inlierDistortedv.Location)
            try
                bbo=convhull(double(inlierOriginalv.Location(:,1)),double(inlierOriginalv.Location(:,2)));
                bbd=convhull(double(inlierDistortedv.Location(:,1)),double(inlierDistortedv.Location(:,2)));
                BW1=poly2mask(double(inlierOriginalv.Location(bbo,1)),double(inlierOriginalv.Location(bbo,2)),width,height);
                BW2=poly2mask(double(inlierDistortedv.Location(bbd,1)),double(inlierDistortedv.Location(bbd,2)),width,height);
                brightnessvaluesvqqkk = [mean(double(I2(BW2))),mean(double(I1_(BW1)))];
            catch
            end
        end
        
    end
    
    
    
    if dosavefigurematchesv && matchedvcountqqkk>mincountpoints
        
        %tft_(qq-cols_array(1))=tf_;
        f=figure('Visible','off');
        showMatchedFeatures(I1_,I2,inlierOriginalv,inlierDistortedv);
        title([imname2v ' : ' imname2],'Interpreter','None')
        %num2str(inlierOriginal.Count);
        saveas(f,[imname2v(1:end-4) '-' imname2(1:end-4) '_vmatch.png'])
        close(f);
        
        
    end
    
    
end

 %[matchedvcountqqkk,matchedcountqqkk]

if isfile([imgfolder imname2v]) && qq==cols_array(1)
    
    matchedPoints_vqqkk=cornerPoints;
    matchedPointsPrev_vqqkk=cornerPoints;
    matchedvcountqqkk=matchedPoints_vqqkk.Count;
    
end

