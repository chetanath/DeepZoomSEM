function [matchedcountqqkk,matchedvcountqqkk,matchedPointsqqkk,matchedPointsPrevqqkk,matchedPoints_vqqkk,matchedPointsPrev_vqqkk,brightnessvaluesqqkk,brightnessvaluesvqqkk,emission_current_ampsqqkk,matchmetricqqkk,matchmetricvqqkk,tformpcqqkk,peakcorrqqkk] = computeMatchedPointsMosaicMaker(kk__,rows_array,cols_array,qq,imgext,imgfolder,height,searchheight,width,searchwidth,dosavefigurematches,dosavefigurematchesv,dodisplayfeatures,mincontrastvalue,minqualityvalue,transformtype,transformtypev,mincountpoints,pointsAreColinear,searchmethod,maxratio,numoctaves,rszscale,confidencelevel,maxnumtrials,maxdistance,dopcmatching)

matchfolder='~/Data/MatlabData/mosaic_matched_points/';
tformpcqqkk=[affine2d,affine2d,affine2d,affine2d];
peakcorrqqkk=[0,0,0,0];
% if ~exist(matchfolder, 'dir')
%     mkdir(matchfolder)
% end

kk=rows_array(kk__);

% inlierOriginal=[];
% inlierDistorted=[];
% inlierDistortedv=[];
% inlierOriginalv=[];

matchedPointsqqkk=BRISKPoints;
matchedPointsPrevqqkk=BRISKPoints;
matchedPoints_vqqkk=BRISKPoints;
matchedPointsPrev_vqqkk=BRISKPoints;

matchedvcountqqkk=0;
matchedcountqqkk=0;
matchmetricqqkk=nan;
matchmetricvqqkk=nan;
brightnessvaluesqqkk=[nan,nan];
brightnessvaluesvqqkk=[nan,nan];
%emission_current_ampsqqkk=nan;

% imshow(I1); hold on; plot(points); hold off;

imname1=sprintf(['%03d_%03d.' imgext],kk-1,qq);
imname2=sprintf(['%03d_%03d.' imgext],kk,qq);

%fprintf('\n %s',imname2)

hdrname=[imgfolder imname2(1:end-4) '-' imname2(end-2:end) '.hdr'];
hdr_=readhdr(hdrname);

hdrv=hdr_.value;

emission_index=strcmp(hdr_.name,'EmissionCurrent');

emission_current_ampsqqkk=str2double(hdrv{emission_index});


searchROI_1=[1 1 searchwidth/rszscale height/rszscale];
searchROI_2=[(width-searchwidth)/rszscale 1 searchwidth/rszscale height/rszscale];

if dopcmatching
    [optimizer, metric] = imregconfig('monomodal');
    optimizer.MaximumStepLength=optimizer.MaximumStepLength/10;
end


% if qq>colt_split
qqp=qq-1;
searchROI_1v=[1 (height-searchheight)/rszscale width/rszscale searchheight/rszscale];
searchROI_2v=[1 1 width/rszscale searchheight/rszscale];
% else
%     qqp=qq+1;
%     searchROI_1v=[1 1 width searchheight];
%     searchROI_2v=[1 height-searchheight width searchheight];
% end

imname2v=sprintf(['%03d_%03d.' imgext],kk,qqp);

I2=imread([imgfolder imname2]);
if rszscale~=1
    I2=imresize(I2,1/rszscale);
end

if isfile([imgfolder imname1]) %&& kk>rows_array(1)
    
    
    I1=imread([imgfolder imname1]);
    if rszscale~=1
        I1=imresize(I1,1/rszscale);
    end
    

    if dopcmatching
        
        %[tformpcqqkk(1),peakcorrqqkk(1)] = imregcorr(imcrop(I2,searchROI_2),imcrop(I1,searchROI_1));
        %[tformpcqqkk(1),peakcorrqqkk(1)] = imregcorr(imcrop(I2,searchROI_2),...
        %    imref2d([searchROI_2(4),searchROI_2(3)+1],[searchROI_2(1),searchROI_2(1)+searchROI_2(3)],[searchROI_2(2),searchROI_2(2)+searchROI_2(4)-1]),...
        %    imcrop(I1,searchROI_1),...
        %    imref2d([searchROI_1(4),searchROI_1(3)+1],[searchROI_1(1),searchROI_1(1)+searchROI_1(3)],[searchROI_1(2),searchROI_1(2)+searchROI_1(4)-1])); 
%          tformEstimate=affine2d;
%          tformEstimate.T(3,1)=-20;
%          tformEstimate.T(3,2)=0;
%         imshowpair(imcrop(I1,searchROI_1),imcrop(I2,searchROI_2),'montage')
         [tformpcqqkk(1),peakcorrqqkk(1)] = imregcorr(imcrop(I2,searchROI_2),...
            imref2d([searchROI_2(4),searchROI_2(3)+1],[searchROI_2(1),searchROI_2(1)+searchROI_2(3)],[searchROI_2(2),searchROI_2(2)+searchROI_2(4)-1]),...
            imcrop(I1,searchROI_1),...
            imref2d([searchROI_1(4),searchROI_1(3)+1],[searchROI_1(1),searchROI_1(1)+searchROI_1(3)],[searchROI_1(2),searchROI_1(2)+searchROI_1(4)-1]),...
            transformtype);
        tformpcqqkk(1) = imregtform(imcrop(I2,searchROI_2),...
            imref2d([searchROI_2(4),searchROI_2(3)+1],[searchROI_2(1),searchROI_2(1)+searchROI_2(3)],[searchROI_2(2),searchROI_2(2)+searchROI_2(4)-1]),...
            imcrop(I1,searchROI_1),...
            imref2d([searchROI_1(4),searchROI_1(3)+1],[searchROI_1(1),searchROI_1(1)+searchROI_1(3)],[searchROI_1(2),searchROI_1(2)+searchROI_1(4)-1]),...
            'affine',optimizer, metric, 'InitialTransformation',tformpcqqkk(1));

            %,...
            %'similarity',optimizer, metric, 'InitialTransformation',tformEstimate); 
         %imwarped=imwarp(imcrop(I2,searchROI_2), tformpcqqkk(1),'OutputView',imref2d([searchROI_1(4),searchROI_1(3)+1],[searchROI_1(1),searchROI_1(1)+searchROI_1(3)],[searchROI_1(2),searchROI_1(2)+searchROI_1(4)-1])); %,'OutputView',imref2d([searchROI_2(4),searchROI_2(3)+1],[searchROI_2(1),searchROI_2(1)+searchROI_2(3)],[searchROI_2(2),searchROI_2(2)+searchROI_2(4)-1]));
         %peakcorrqqkk(1)=sum(sum(imabsdiff(imcrop(I1,searchROI_1),imwarped).^2))./numel(imwarped);         
         %imshowpair(imcrop(I1,searchROI_1),imwarped,'montage')
         %peakcorrqqkk(1)
         [msg, ~] = lastwarn;
             %if ~isempty(msg)
             %    peakcorrqqkk(1)=inf;
             %end
            %warning('off', id);
        %tformpcqqkk(1).T
        
        [tformpcqqkk(3),peakcorrqqkk(3)] = imregcorr(imcrop(I1,searchROI_1),...
            imref2d([searchROI_1(4),searchROI_1(3)+1],[searchROI_1(1),searchROI_1(1)+searchROI_1(3)],[searchROI_1(2),searchROI_1(2)+searchROI_1(4)-1]),...
            imcrop(I2,searchROI_2),...
            imref2d([searchROI_2(4),searchROI_2(3)+1],[searchROI_2(1),searchROI_2(1)+searchROI_2(3)],[searchROI_2(2),searchROI_2(2)+searchROI_2(4)-1]),transformtype); 
            [msg, ~] = lastwarn;
            %if ~isempty(msg)
            %    peakcorrqqkk(3)=inf;
            %end
            tformpcqqkk(3) = imregtform(imcrop(I1,searchROI_1),...
            imref2d([searchROI_1(4),searchROI_1(3)+1],[searchROI_1(1),searchROI_1(1)+searchROI_1(3)],[searchROI_1(2),searchROI_1(2)+searchROI_1(4)-1]),...
            imcrop(I2,searchROI_2),...
            imref2d([searchROI_2(4),searchROI_2(3)+1],[searchROI_2(1),searchROI_2(1)+searchROI_2(3)],[searchROI_2(2),searchROI_2(2)+searchROI_2(4)-1]),...
            'affine',optimizer, metric, 'InitialTransformation',tformpcqqkk(3)); 
            [msg, ~] = lastwarn;
            
            %warning('off', id)
        %tformpc.T(3,1)=tformpc.T(3,1)-( (width-searchwidth)/rszscale );
        
    end
    
    pointsPrevious  = detectBRISKFeatures(I1,'ROI',searchROI_1,'MinContrast',mincontrastvalue,'MinQuality',minqualityvalue,'NumOctaves',numoctaves);
    %pointsPrevious  = detectFASTFeatures(I1,'ROI',searchROI_1,'MinQuality',minqualityvalue,'MinContrast',mincontrastvalue); %,'MinContrast',0.01); %,'MinContrast',0.01); %detectSURFFeatures(I1,'ROI',[1 1 searchwidth height],'MetricThreshold',1,'NumScaleLevels',6,'NumOctaves',1);
    [featuresPrevious, pointsPrevious]  = extractFeatures(I1,  pointsPrevious);
    % [featuresPrevious, pointsPrevious]  = extractFeatures(I1,  pointsPrevious, 'Method','Block');
    
    points = detectBRISKFeatures(I2,'ROI',searchROI_2,'MinContrast',mincontrastvalue,'MinQuality',minqualityvalue,'NumOctaves',numoctaves);
    [features, points]  = extractFeatures(I2,  points);
    %points = detectFASTFeatures(I2,'ROI',searchROI_2,'MinQuality',minqualityvalue,'MinContrast',mincontrastvalue); %,'MinContrast',0.01); %,'MetricThreshold',100,'NumScaleLevels',6);
    %[features, points]  = extractFeatures(I2,  points, 'Method','Block');
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
        %[~, inlierDistorted, inlierOriginal, ~] = ...
        [~, matchedPointsqqkk, matchedPointsPrevqqkk, ~] = ...
            estimateGeometricTransform(matchedPointsqqkk, ...
            matchedPointsPrevqqkk, transformtype, 'Confidence', confidencelevel, 'MaxNumTrials', maxnumtrials, 'MaxDistance', maxdistance); %, 'MaxDistance', round(height/100));
        matchedcountqqkk=matchedPointsqqkk.Count;
        matchmetricqqkk=matchmetric;
        try
            if ~pointsAreColinear(matchedPointsqqkk.Location) && ~pointsAreColinear(matchedPointsPrevqqkk.Location)
                try
                    bbo=convhull(double(matchedPointsPrevqqkk.Location(:,1)),double(matchedPointsPrevqqkk.Location(:,2)));
                    bbd=convhull(double(matchedPointsqqkk.Location(:,1)),double(matchedPointsqqkk.Location(:,2)));
                    BW1=poly2mask(double(matchedPointsPrevqqkk.Location(bbo,1)),double(matchedPointsPrevqqkk.Location(bbo,2)),width/rszscale,height/rszscale);
                    BW2=poly2mask(double(matchedPointsqqkk.Location(bbd,1)),double(matchedPointsqqkk.Location(bbd,2)),width/rszscale,height/rszscale);
                    brightnessvaluesqqkk = [mean(double(I2(BW2))),mean(double(I1(BW1)))];
                catch
                end
            end
        catch
        end
        
        
        if dosavefigurematches
            
            f=figure('Visible','off');
            ax = axes;
            showMatchedFeatures(I1,I2,matchedPointsPrevqqkk,matchedPointsqqkk,'montage');
            title([imname1 ' : ' imname2],'Interpreter','None')
            saveas(f,[matchfolder imname1(1:end-4) '-' imname2(1:end-4) '_match.png'])
            close(f);
            
            
        end
        
        
    end
    
    
    
end


if isfile([imgfolder imname2v]) %&& qq~=cols_array(1)
    
    I1_=imread([imgfolder imname2v]);
    if rszscale~=1
        I1_=imresize(I1_,1/rszscale);
    end
    
    if dopcmatching
        
        %[tformpcqqkk(3),peakcorrqqkk(3)] = imregcorr(imcrop(I2,searchROI_2v),imcrop(I1_,searchROI_1v));
        %[tformpcqqkk(4),peakcorrqqkk(4)] = imregcorr(imcrop(I1_,searchROI_1v),imcrop(I2,searchROI_2v));
        
            [tformpcqqkk(2),peakcorrqqkk(2)] = imregcorr(imcrop(I2,searchROI_2v),...
            imref2d([searchROI_2v(4)+1,searchROI_2v(3)],[searchROI_2v(1),searchROI_2v(1)+searchROI_2v(3)-1],[searchROI_2v(2),searchROI_2v(2)+searchROI_2v(4)]),...
            imcrop(I1_,searchROI_1v),...
            imref2d([searchROI_1v(4)+1,searchROI_1v(3)],[searchROI_1v(1),searchROI_1v(1)+searchROI_1v(3)-1],[searchROI_1v(2),searchROI_1v(2)+searchROI_1v(4)]),transformtypev); 
            tformpcqqkk(2) = imregtform(imcrop(I2,searchROI_2v),...
            imref2d([searchROI_2v(4)+1,searchROI_2v(3)],[searchROI_2v(1),searchROI_2v(1)+searchROI_2v(3)-1],[searchROI_2v(2),searchROI_2v(2)+searchROI_2v(4)]),...
            imcrop(I1_,searchROI_1v),...
            imref2d([searchROI_1v(4)+1,searchROI_1v(3)],[searchROI_1v(1),searchROI_1v(1)+searchROI_1v(3)-1],[searchROI_1v(2),searchROI_1v(2)+searchROI_1v(4)]),...
            'affine',optimizer, metric, 'InitialTransformation',tformpcqqkk(2));
        %tformpcqqkk(1).T
            [msg, ~] = lastwarn;
            %if ~isempty(msg)
            %    peakcorrqqkk(2)=inf;
            %end
            %warning('off', id)
        
            [tformpcqqkk(4),peakcorrqqkk(4)] = imregcorr(imcrop(I1_,searchROI_1v),...
            imref2d([searchROI_1v(4)+1,searchROI_1v(3)],[searchROI_1v(1),searchROI_1v(1)+searchROI_1v(3)-1],[searchROI_1v(2),searchROI_1v(2)+searchROI_1v(4)]),...
            imcrop(I2,searchROI_2v),...
            imref2d([searchROI_2v(4)+1,searchROI_2v(3)],[searchROI_2v(1),searchROI_2v(1)+searchROI_2v(3)-1],[searchROI_2v(2),searchROI_2v(2)+searchROI_2v(4)]),transformtypev);
            [msg, ~] = lastwarn;
            %if ~isempty(msg)
            %    peakcorrqqkk(4)=inf;
            %end
            tformpcqqkk(4) = imregtform(imcrop(I1_,searchROI_1v),...
            imref2d([searchROI_1v(4)+1,searchROI_1v(3)],[searchROI_1v(1),searchROI_1v(1)+searchROI_1v(3)-1],[searchROI_1v(2),searchROI_1v(2)+searchROI_1v(4)]),...
            imcrop(I2,searchROI_2v),...
            imref2d([searchROI_2v(4)+1,searchROI_2v(3)],[searchROI_2v(1),searchROI_2v(1)+searchROI_2v(3)-1],[searchROI_2v(2),searchROI_2v(2)+searchROI_2v(4)]),...
            'affine',optimizer,metric, 'InitialTransformation',tformpcqqkk(4));

            %warning('off', id)
    end
    
    
    pointsv_  = detectBRISKFeatures(I1_,'ROI',searchROI_1v,'MinContrast',mincontrastvalue,'MinQuality',minqualityvalue,'NumOctaves',numoctaves);
    [featuresv_, pointsv_]  = extractFeatures(I1_,  pointsv_);
    %pointsv_ = detectFASTFeatures(I1_,'ROI',searchROI_1v,'MinQuality',minqualityvalue,'MinContrast',mincontrastvalue); %,'MinContrast',0.01); %,'MetricThreshold',100,'NumScaleLevels',6);
    %[featuresv_, pointsv_]  = extractFeatures(I1_,  pointsv_, 'Method','Block');
    
    pointsv = detectBRISKFeatures(I2,'ROI',searchROI_2v,'MinContrast',mincontrastvalue,'MinQuality',minqualityvalue,'NumOctaves',numoctaves);
    %pointsv = detectFASTFeatures(I2,'ROI',searchROI_2v,'MinQuality',minqualityvalue,'MinContrast',mincontrastvalue); %,'MinContrast',0.01); %,'MetricThreshold',100,'NumScaleLevels',6);
    [featuresv, pointsv]  = extractFeatures(I2,  pointsv); % 'Method','Block');
    
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
        
        %[~, inlierDistortedv, inlierOriginalv, ~] = ...
        [ttf, matchedPoints_vqqkk, matchedPointsPrev_vqqkk, ~] = ...
            estimateGeometricTransform(matchedPoints_vqqkk, ...
            matchedPointsPrev_vqqkk, transformtypev, 'Confidence', confidencelevel, 'MaxNumTrials', maxnumtrials, 'MaxDistance', maxdistance); %, 'MaxDistance', round(height/100));
        matchedvcountqqkk=matchedPointsPrev_vqqkk.Count;
        matchmetricvqqkk=matchmetricv;
        
        % Error in the the follow line when making sure points are not colinear?
        %        matchedPoints_vqqkk.Location
        %        matchedPointsPrev_vqqkk.Location
        try
            if ~pointsAreColinear(matchedPoints_vqqkk.Location) && ~pointsAreColinear(matchedPointsPrev_vqqkk.Location)
                try
                    bbo=convhull(double(matchedPointsPrev_vqqkk.Location(:,1)),double(matchedPointsPrev_vqqkk.Location(:,2)));
                    bbd=convhull(double(matchedPoints_vqqkk.Location(:,1)),double(matchedPoints_vqqkk.Location(:,2)));
                    BW1=poly2mask(double(matchedPointsPrev_vqqkk.Location(bbo,1)),double(matchedPointsPrev_vqqkk.Location(bbo,2)),width/rszscale,height/rszscale);
                    BW2=poly2mask(double(matchedPoints_vqqkk.Location(bbd,1)),double(matchedPoints_vqqkk.Location(bbd,2)),width/rszscale,height/rszscale);
                    brightnessvaluesvqqkk = [mean(double(I2(BW2))),mean(double(I1_(BW1)))];
                catch
                end
            end
        catch
        end
        
        
        if dosavefigurematchesv
            
            %tft_(qq-cols_array(1))=tf_;
            f=figure('Visible','off');
            ax = axes;
            matchedPointsPrev_vqqkk2=matchedPointsPrev_vqqkk;
            matchedPointsPrev_vqqkk2.Location(:,2)=matchedPointsPrev_vqqkk.Location(:,1);
            matchedPointsPrev_vqqkk2.Location(:,1)=size(I1_,1)-matchedPointsPrev_vqqkk.Location(:,2);
            matchedPoints_vqqkk2=matchedPoints_vqqkk;
            matchedPoints_vqqkk2.Location(:,2)=matchedPoints_vqqkk.Location(:,1);
            matchedPoints_vqqkk2.Location(:,1)=size(I1_,1)-matchedPoints_vqqkk.Location(:,2);            
            showMatchedFeatures(rot90(I1_,-1),rot90(I2,-1),matchedPointsPrev_vqqkk2,matchedPoints_vqqkk2,'montage');
            title([imname2v ' : ' imname2],'Interpreter','None')
            %num2str(inlierOriginal.Count);
            saveas(f,[matchfolder imname2v(1:end-4) '-' imname2(1:end-4) '_vmatch.png'])
            close(f);
            
            
        end
        
    end
    
    
end

%[matchedvcountqqkk,matchedcountqqkk]

if isfile([imgfolder imname2v]) && qq==cols_array(1)
    
    matchedPoints_vqqkk=BRISKPoints;
    matchedPointsPrev_vqqkk=BRISKPoints;
    matchedvcountqqkk=matchedPoints_vqqkk.Count;
    
end

if rszscale~=1
    matchedPoints_vqqkk.Location=matchedPoints_vqqkk.Location.*rszscale;
    matchedPointsPrev_vqqkk.Location=matchedPointsPrev_vqqkk.Location*rszscale;
    matchedPointsqqkk.Location=matchedPointsqqkk.Location.*rszscale;
    matchedPointsPrevqqkk.Location=matchedPointsPrevqqkk.Location.*rszscale;
end