function OFFSET=globalbalance2(brightnessvalues,brightnessvaluesv,cols_array)
% function [XX_factor,OFFSET]=globalbalance2(brightnessvalues,brightnessvaluesv,cols_array)

threshdiff=1e-5;

polyfittype='poly22';

sz1=size(brightnessvalues,1);
sz2=size(brightnessvalues,2);

XX_offset=nan(sz1,sz2);
XX_factor=nan(sz1,sz2);

stdimg_=inf;
stdimg=realmax;

iternumber=0;

meanbrightness0=nanmean(brightnessvalues(:));

while stdimg/stdimg_<(1-threshdiff) %diffsumsq./diffsumsqp<0.95

    iternumber=iternumber+1;
    
    stdimg_=stdimg;
        
    for qq=cols_array
        
        %if qq<=colt_split
        %    qqn=qq-1;
        %else
            qqn=qq+1;
        %end
        
        firstpair=find(isfinite(brightnessvalues(qq,:,1)),1,'first');
        lastpair=find(isfinite(brightnessvalues(qq,:,1)),1,'last');
        
        
        for kk=firstpair+1:lastpair
            
            allhist1=nan(4,1);
            allhist2=nan(4,1);
            
            
            if sum(isnan(brightnessvalues(qq,kk,:)))==0
                allhist1(1)=brightnessvalues(qq,kk,1);
                allhist2(1)=brightnessvalues(qq,kk,2);
                
            end
            
            if sum(isnan(brightnessvaluesv(qq,kk,:)))==0
                allhist1(2)=brightnessvaluesv(qq,kk,1);
                allhist2(2)=brightnessvaluesv(qq,kk,2);
                
            end
            
            
            if kk+1<sz2
                if sum(isnan(brightnessvalues(qq,kk+1,:)))==0
                    allhist1(3)=brightnessvalues(qq,kk+1,2);
                    allhist2(3)=brightnessvalues(qq,kk+1,1);
                end
            end
            
            if qqn<=max(cols_array) && qqn>=min(cols_array)
                if sum(isnan(brightnessvaluesv(qqn,kk,:)))==0
                    allhist1(4)=brightnessvaluesv(qqn,kk,2);
                    allhist2(4)=brightnessvaluesv(qqn,kk,1);
                    
                end
            end
            
            isnanallhist1=isnan(allhist1);
            
            if sum(~isnanallhist1)>0

                isfinitehist=isfinite(allhist2) & isfinite(allhist1);
                XX=rms(allhist2(isfinitehist))-rms(allhist1(isfinitehist));
 
                if ~isnanallhist1(1)
                    brightnessvalues(qq,kk,1)=brightnessvalues(qq,kk,1)+XX;
                end
                
                if ~isnanallhist1(2)
                    brightnessvaluesv(qq,kk,1)=brightnessvaluesv(qq,kk,1)+XX;
                end
                
                if ~isnanallhist1(3)
                    brightnessvalues(qq,kk+1,2)=brightnessvalues(qq,kk+1,2)+XX;
                end
                
                if ~isnanallhist1(4)
                    brightnessvaluesv(qqn,kk,2)=brightnessvaluesv(qqn,kk,2)+XX;
                end
                
                if isnan(XX_offset(qq,kk))
                    XX_offset(qq,kk)=XX;
                else
                    XX_offset(qq,kk)=XX_offset(qq,kk)+XX;
                end

            end
            
        end
        
        
    end

    stdimg=nanstd(brightnessvalues(:));                              
    
    fprintf('Iter=%d, Mean Brightness=%.6e, New STD=%.6e, Old STD=%.6e, Ratio=%.3e\n',iternumber,nanmean(brightnessvalues(:)),stdimg,stdimg_,1-stdimg./stdimg_)
    
    
    % meanpixelvalue=[meanpixelvalue,mean([nanmean(nanmean(brightnessvalues(:,:,1))), nanmean(nanmean(brightnessvaluesv(:,:,1)))])];
       
end

OFFSET=inpaint_nans(XX_offset,0);
SLOPE=ones(size(OFFSET));

% SLOPE CALCULATION (NOT USED):
% isf=isfinite(brightnessvalues(:,:,1)) & isfinite(brightnessvaluesv(:,:,1));
% brightnessvalues(isf)=brightnessvalues(isf)+OFFSET(isf);
% brightnessvaluesv(isf)=brightnessvaluesv(isf)+OFFSET(isf);
% 
% [X,Y]=ndgrid(1:size(brightnessvalues,1),1:size(brightnessvalues,2));
% BV=0.5.*(brightnessvalues(:,:,1)+brightnessvaluesv(:,:,1));
% sf=fit([X(isf),Y(isf)],BV(isf),polyfittype);
% SLOPE=nan(size(OFFSET));
% SLOPE(isf)=mean(BV(isf))./sf(X(isf),Y(isf));
% SLOPE=inpaint_nans(SLOPE,0);
%SLOPE=SLOPE/mean(SLOPE(:));
%OFFSET=SLOPE.*OFFSET;






% if numel(cols_array)>1
%     
%     avgbrightness=nanmean(BV(:));
%     
%     for ii=1:size(brightnessvalues,1)
%         for jj=1:size(brightnessvalues,2)
%             shift=sf(X(ii),Y(ii));
%             if isfinite(shift)
%                 OFFSET(ii,jj)=avgbrightness-shift;
%             end
%         end
%     end
%     
% end

% OFFSET__=nanmean(BV(:))-sf(X,Y);
% OFFSET(isf)=OFFSET__(isf);


% isfXX=isfinite(OFFSET);

% brightnessvalues(isf)=brightnessvalues(isf)+OFFSET(isf);
% brightnessvaluesv(isf)=brightnessvaluesv(isf)+OFFSET(isf);
% 
% OFFSET=inpaint_nans(OFFSET+XX_offset,0);
% 
% % brightnessvalues=brightnessvalues.*meanbrightness0./nanmean(brightnessvalues(:));
% 
% %[X,Y]=ndgrid(1:size(brightnessvalues,1),1:size(brightnessvalues,2));
% BV=0.5.*(brightnessvalues(:,:,1)+brightnessvaluesv(:,:,1));
% isf=isfinite(BV(:));
% sf=fit([X(isf),Y(isf)],BV(isf),polyfittype);
% 
% XX_factor=nanmean(BV(:))./sf(X,Y);
%XX_factor(isf)=XX_factor__(isf);

% if numel(cols_array)>1
%     
%     avgbrightness=nanmean(BV(:));
%     
%     for ii=1:size(brightnessvalues,1)
%         for jj=1:size(brightnessvalues,2)
%             shift=sf(X(ii),Y(ii));
%             if isfinite(shift) && shift~=0
%                 XX_factor(ii,jj)=avgbrightness./shift;
%             end
%         end
%     end
%     
% end

%XX_factor-inpaint_nans(XX_factor,0);

%%



