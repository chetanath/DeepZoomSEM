clearvars;
close all;

doBSE=true;
doNN=true;
doelementmaps=true;


parentfolder='~/Data/EDAX_EDS/ProgramData/EDAX/Ryan/Cole/Mapping/Lsm/AguasZarcasFeb2021/';
subfolder='New Sample/';
multifieldfile='AguasZarcasFeb2021-20210201114032296.tmt';

endout=regexp(parentfolder,filesep,'split');
samplename=endout{end-1};

SS = edax_tmt_file([parentfolder multifieldfile]);

Id=SS.Id;
stagex=-SS.stagex;
stagey=SS.stagey;
nmaps=numel(Id);
mapnumbers=1:nmaps;



xraymapfolder=[parentfolder subfolder];

d=dir(join([xraymapfolder 'Area ' num2str(Id(2)) '/*.spd'],''));
s1 = processEDAXmaps(join([xraymapfolder 'Area ' num2str(Id(2)) '/' d.name],''));
dipr=dir(join([xraymapfolder 'Area ' num2str(Id(2)) '/map*_PhaseOverlay.ipr'],''));
iprs = readipr(join([xraymapfolder 'Area ' num2str(Id(2)) '/' dipr.name],''));
sspc = readspc(join([xraymapfolder 'Area ' num2str(Id(2)) '/' d.name(1:end-3) 'spc'],''));
microns_per_pixel_x=iprs.mppX;
microns_per_pixel_y=iprs.mppY;
imageSize=[size(s1.allmaps,1),size(s1.allmaps,2)];
fovy=microns_per_pixel_x*imageSize(1);
fovx=microns_per_pixel_y*imageSize(2);
nm_per_pixel_mean=mean([microns_per_pixel_x,microns_per_pixel_y])*1e3;
fprintf('%.3f\n',nm_per_pixel_mean)
%%
xpixel=round(1e3*stagex./microns_per_pixel_x);
ypixel=round(1e3*stagey./microns_per_pixel_y);
xrange=range(stagex);
yrange=range(stagey);
nmaps=numel(Id);



for ii=1:nmaps
    
    tforms(ii)=affine2d;
    tforms(ii).T(3,1)=xpixel(ii)-imageSize(1)/2;
    tforms(ii).T(3,2)=ypixel(ii)-imageSize(1)/2;
    
end

% Create empty pano:
for ii = 1:numel(tforms)
    [xlim(ii,:), ylim(ii,:)] = outputLimits(tforms(ii), [1 imageSize(1)], [1 imageSize(2)]);
end


% Find the minimum and maximum output limits
xMin = min(xlim(:));
xMax = max(xlim(:));

yMin = min(ylim(:));
yMax = max(ylim(:));

% Width and height of panorama.
% width  = round(xMax - xMin)+imageSize(1);
% height = round(yMax - yMin)+imageSize(2);
width  = round(xMax - xMin);
height = round(yMax - yMin);


% Initialize the "empty" panorama.
panoramabse = zeros([height width 1],'uint8');
panoramaNN = zeros([height width 3],'uint8');

%panoramaNN = zeros([height width 3]);

blender = vision.AlphaBlender('Operation', 'Binary mask', ...
    'MaskSource', 'Input port');

% Create a 2-D spatial reference object defining the size of the panorama.
% xLimits = [xMin-imageSize(1)/2 xMax+imageSize(1)/2];
% yLimits = [yMin-imageSize(2)/2 yMax+imageSize(2)/2];
xLimits = [xMin xMax];
yLimits = [yMin yMax];

panoramaView = imref2d([height width], xLimits, yLimits);
panoramaView3 = imref2d([height width 3], xLimits, yLimits);

dist_colors = distinguishable_colors(nmaps);

%parfor ii=2:nmaps+1 %7
%mapnumbers=1:nmaps;
%%
if doBSE
    for jj=1:numel(mapnumbers) %7
        
        disp(['BSE: Area ' num2str(mapnumbers(jj))])
        
        d=dir(join([xraymapfolder 'Area ' num2str(mapnumbers(jj)) '/fov*.bmp'],''));
        
        bse=imread([d(end).folder '/' d(end).name]);
        
        if size(bse,2)~=imageSize(1)
            bse=imresize(bse,imageSize(1)./size(bse,2));
        end
        
        % Transform I into the panorama.
        warpedImage = imwarp(bse, tforms(jj), 'OutputView', panoramaView);
        
        % Generate a binary mask.
        mask = imwarp(true(size(bse,1),size(bse,2)), tforms(jj), 'OutputView', panoramaView);
        
        % Overlay the warpedImage onto the panorama.
        panoramabse = step(blender, panoramabse, warpedImage, mask);
        
    end
    
    imwrite(panoramabse,[samplename '_bse.png'])
    release(blender);
    blender = vision.AlphaBlender('Operation', 'Binary mask', ...
        'MaskSource', 'Input port');
end


%%
if doNN
    for jj=1:numel(mapnumbers)
        disp(['NN: Area ' num2str(mapnumbers(jj))])
        fsz=100;
        gfac=1.5;
        NN2=permute(repmat(dist_colors(jj,:)',1,imageSize(1),imageSize(2)),[2,3,1]);
        nnx=floor(size(NN2,2)/(gfac*fsz));
        nny=floor(size(NN2,1)/(gfac*fsz));
        
        for kk=1:nnx
            for mm=1:nny
                if mean(dist_colors(jj,:))>=0.5
                    textcolor='black';
                else
                    textcolor='white';
                end
                NN2=uint8(round(255*insertText(NN2,round([kk*fsz*gfac,mm*fsz*gfac]),num2str(mapnumbers(jj)),'FontSize',fsz,'TextColor',textcolor,'AnchorPoint','Center','BoxOpacity',0)));
            end
        end
        
        warpedImageNN = imwarp(NN2, tforms(jj), 'OutputView', panoramaView3);
        maskNN = imwarp(true(imageSize(1),imageSize(2)), tforms(jj), 'OutputView', panoramaView3);
        panoramaNN = step(blender, panoramaNN, warpedImageNN, maskNN);
        
    end
    
    imwrite(panoramaNN,[samplename '_edx_locations.png']);
    release(blender);
    blender = vision.AlphaBlender('Operation', 'Binary mask', ...
        'MaskSource', 'Input port');
end
%%



if doelementmaps
    
    n_element_maps=size(s1.allmaps,3);
    Izzkk=zeros(imageSize(1),imageSize(2),numel(mapnumbers),n_element_maps);
    
    
    parfor jj=1:numel(mapnumbers) %7
        
        fprintf(['Area ' num2str(mapnumbers(jj)) '\n']);
        
        d=dir(join([xraymapfolder 'Area ' num2str(mapnumbers(jj)) '/*.spd'],''));
        
        s = processEDAXmaps(join([xraymapfolder 'Area ' num2str(mapnumbers(jj)) '/' d(end).name],''));
         
        Izzkk(:,:,jj,:)=s.allmaps./sum(s.data,3);
        
    end
    
    s=[];
    
    
    
    
    parfor kk=1:n_element_maps
        
        panorama = zeros([height width 1]);
        
        for jj=1:numel(mapnumbers) %7
            
            fprintf([s1.allelementsnames{kk} ' , Area ' num2str(mapnumbers(jj)) '\n']);
            
            d=dir(join([xraymapfolder 'Area ' num2str(mapnumbers(jj)) '/*.spd'],''));
            
            %if numel(d)==1
                              
                % Transform I into the panorama.
                warpedImage = imwarp(Izzkk(:,:,jj,kk)', tforms(jj), 'OutputView', panoramaView);
                
                % Generate a binary mask.
                mask = imwarp(true(imageSize(2),imageSize(1)), tforms(jj), 'OutputView', panoramaView);
                
                % Overlay the warpedImage onto the panorama.
                panorama = step(blender, panorama, warpedImage, mask);
            %end
            
        end
        
        %imwrite(panoramaNN,'panoramaNN.png')
        imwrite(panorama./max(panorama(:)),[samplename '_' s1.allelementsnames{kk} '_Map.png'])
        
    end
end