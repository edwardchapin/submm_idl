; AUTHORS
; Isaac Roseboom, igr@roe.ac.uk
;
; HISTORY
; 14OCT2011: Initial committed version

pro lstdrv_initsolveSP,sx,sy,segsrc,prf,fwhm,paxis,noisy_map,sig_map,x_pix,y_pix,seg,bl,bu,psrc,psig,pmat,amat2,sn_map,bvec,xvec,iBIC,solve=solve,sx_pix=sx_pix,sy_pix=sy_pix

snsrc=n_elements(segsrc)
paxis1=paxis[0]
paxis2=paxis[1]

   prf=psf_gaussian(fwhm=fwhm,npixel=[paxis1,paxis2])
      
      ssx=sx[segsrc]
      ssy=sy[segsrc]

      mux=max(ssx,min=mlx)
      muy=max(ssy,min=mly)
      npix=where(x_pix lt mux+20 and y_pix lt muy+20 and y_pix ge mly-20 and x_pix ge mlx-20,snpix)
      sx_pix=x_pix[npix]
      sy_pix=y_pix[npix]
      snoisy_map=noisy_map[npix]
      ssig_map=sig_map[npix]
      
    
      pmat=dblarr(snsrc,snsrc)
      bvec=dblarr(snsrc)
      amat=dblarr(snpix,snsrc)
      amat2=dblarr(snpix,snsrc)
      
    
      for s=0L,snsrc-1 do begin
         dx = -round(ssx[s])+(paxis1-1.)/2.+sx_pix
         dy = -round(ssy[s])+(paxis2-1.)/2.+sy_pix
         pindx=findgen(paxis1)-ssx[s]+round(ssx[s])
         pindy=findgen(paxis2)-ssy[s]+round(ssy[s])
         px=ssx[s]-round(ssx[s])+(paxis1-1.)/2.
         py=ssy[s]-round(ssy[s])+(paxis2-1.)/2.
         
         good = where(dx GE 0 AND dx LT paxis1 AND dy GE 0 AND dy LT paxis2, ngood, comp=bad, ncomp=nbad)
         if(ngood gt 0.5*n_elements(prf)) then begin
         make_2d,pindx,pindy,ipx2,ipy2   
         nprf=interpolate(prf,ipx2,ipy2,cubic=-0.5)
            amat[good,s]=nprf[round(dx[good]),round(dy[good])]/sqrt(ssig_map[good])
            amat2[good,s]=nprf[round(dx[good]),round(dy[good])]/(ssig_map[good])
            if(s eq 0) then begin
                arow=fltarr(ngood)+s
                acol=good
            endif else begin
                arow=[arow,fltarr(ngood)+s]
                acol=[acol,good]
            endelse
         endif
         endfor
      ; do sparse multiplication
         t1=systime(/sec)



         pmat=dblarr(snsrc,snsrc)
         sindex=sort(acol)
         ncnt=n_elements(acol)
; create pmat
         cnt=0L

         while(cnt lt ncnt) do begin
             
             wpix=acol[sindex[cnt]]
             wind=cnt
             nipix=1L
             cnt=cnt+1
             if(cnt lt ncnt) then begin
                 while(acol[sindex[cnt]] eq wpix ) do begin
                     
                     
                     
                     nipix=nipix+1L    
                     wind=[wind,cnt]
                     cnt=cnt+1L
                     if(cnt eq ncnt) then break
                 endwhile
             endif

             ipix=sindex[wind]
             irow=arow[ipix]
             make_2d,irow,irow,tPx,tPy
             make_2d,ipix,ipix,tIx,tIy
             tPx=reform(tPx,nipix^2)
             tPy=reform(tPy,nipix^2)
             tIx=reform(tIx,nipix^2)
             tIy=reform(tIy,nipix^2)
             pmat[tPx,tPy]=pmat[tPx,tPy]+amat[acol[tIx],arow[tIx]]*amat[acol[tIy],arow[tIy]]
            
         endwhile
        
;psig
         psig=dblarr(snsrc,snsrc)
         sindex=sort(acol)
         ncnt=n_elements(acol)
; create pmat
         cnt=0L

         while(cnt lt ncnt) do begin
             
             wpix=acol[sindex[cnt]]
             wind=cnt
             nipix=1L
             cnt=cnt+1
             if(cnt lt ncnt) then begin
                 while(acol[sindex[cnt]] eq wpix ) do begin
                     
                     
                     
                     nipix=nipix+1L    
                     wind=[wind,cnt]
                     cnt=cnt+1L
                     if(cnt eq ncnt) then break
                 endwhile
             endif

             ipix=sindex[wind]
             irow=arow[ipix]
             make_2d,irow,irow,tPx,tPy
             make_2d,ipix,ipix,tIx,tIy
             tPx=reform(tPx,nipix^2)
             tPy=reform(tPy,nipix^2)
             tIx=reform(tIx,nipix^2)
             tIy=reform(tIy,nipix^2)
             psig[tPx,tPy]=psig[tPx,tPy]+amat2[acol[tIx],arow[tIx]]*amat2[acol[tIy],arow[tIy]]
            
         endwhile
         


; create bvec
         bvec=dblarr(snsrc)
         for i=0L,snsrc-1 do begin
             ipix=where(arow eq i)
             if(ipix[0] ne -1) then bvec[i]=total(amat2[acol[ipix],arow[ipix]]*snoisy_map[acol[ipix]]/ssig_map[acol[ipix]])
             
         endfor
         
         
         print,'finished constructing arrays ',systime(/sec)-t1
         sn_map=snoisy_map/(ssig_map)
         xvec=dblarr(snsrc)
    

end
