function etfT = errortformfun(qq,kk,Tt2,matchedPointsPrev_v,matchedPoints_v,matchedPointsPrev,matchedPoints,colmin,colmax,rowmin,rowmax)

errortformT=nan(colmax,rowmax,2);

for qq_=(qq-1):(qq+1)
    mp_llpx=[];
    mp_llpy=[];
    llpx=[];
    llpy=[];
    %    rows_array=rowsforcol{colt==qq_}';
    if qq_==(qq-1) || qq_==(qq+1)
        rows_array=kk;
    else
        rows_array=(kk-1):(kk+1);
    end
    
    for kk_=rows_array
        
        
        %disp(['ROWS: ' num2str([qq_,kk_])]);
        
        t00=projective2d;
        if qq_>=colmin && qq_<=colmax && kk_>=rowmin && kk_<=rowmax
            if isfinite(Tt2(qq_,kk_,1,1))
                
                t00.T=squeeze(Tt2(qq_,kk_,:,:));
                
                if qq_-1>=colmin
                    if isfinite(Tt2(qq_-1,kk_,1,1))
                        t01=projective2d;
                        t01.T=squeeze(Tt2(qq_-1,kk_,:,:));
                        [llpx,llpy]        =transformPointsForward(t01,matchedPointsPrev_v{qq_,kk_}.Location(:,1),matchedPointsPrev_v{qq_,kk_}.Location(:,2));
                        [mp_llpx,mp_llpy] = transformPointsForward(t00,matchedPoints_v{qq_,kk_}.Location(:,1),matchedPoints_v{qq_,kk_}.Location(:,2));
                        errortformT(qq_,kk_,1)=sum((llpx-mp_llpx).^2+(llpy-mp_llpy).^2);
                    else
                        llpx=[];
                        llpy=[];
                        mp_llpx=[];
                        mp_llpy=[];
                    end
                end
                
                if kk_-1>=colmin
                    if isfinite(Tt2(qq_,kk_-1,1,1))
                        t02=projective2d;
                        t02.T=squeeze(Tt2(qq_,kk_-1,:,:));
                        [llpx2,llpy2]        =transformPointsForward(t02,matchedPointsPrev{qq_,kk_}.Location(:,1),matchedPointsPrev{qq_,kk_}.Location(:,2));
                        [mp_llpx2,mp_llpy2]=transformPointsForward(t00,matchedPoints{qq_,kk_}.Location(:,1),matchedPoints{qq_,kk_}.Location(:,2));
                        errortformT(qq_,kk_,2)=sum((llpx2-mp_llpx2).^2+(llpy2-mp_llpy2).^2);
                    else
                        llpx2=[];
                        llpy2=[];
                        mp_llpx2=[];
                        mp_llpy2=[];
                    end
                end
                
                
            end
        end
        
        
    end
    
end
etfT=nansum(nansum(errortformT(:,:,2)))+nansum(nansum(errortformT(:,:,1)));
