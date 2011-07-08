
;+
; NAME:
; segsrc
; 
; PURPOSE:
; Take an input list of positions and SPIRE maps at 3 wavelengths and
; group those sources based on enclosing regions around them above a given S/N
; threshold. The code works by first identifying pixels which are
; above the S/N threshold, then using a flood-fill algorithm to find
; connected regions of pixels. Finally sources which are located in
; pixels below the S/N threshold are placed in the segment closest to them.
;
; CALLING SEQUENCE
; segsrc,sx,sy,pfwhm250,pfwhm350,pfwhm500,hdr250, map250,sig250 $
;                   ,hdr350, map350,sig350 $
;                   ,hdr500, map500,sig500,segsrc,bsn=bsn
; INPUT PARAMETERS
;
; sx: x pixel coordinates of input sources in the 250 micron map
; sy: y pixel coordinates of input sources in the 250 micron map
; pfwhm250: FWHM of PRF at 250 micron in pixels 
; pfwhm350: FWHM of PRF at 350 micron in pixels 
; pfwhm500: FWHM of PRF at 500 micron in pixels 
; hdr250: FITS header of 250 micron map
; map250: 250 micron map
; sig250: 250 micron noise map
; hdr350: FITS header of 350 micron map
; map350: 350 micron map
; sig350: 350 micron noise map
; hdr500: FITS header of 500 micron map
; map500: 500 micron map
; sig500: 500 micron noise map
; bsn: signal to noise threshold at which to segment map
;
; OUTPUTS
; ssrc: an integer array which contains the segment ID for each input
; source
;
; AUTHORS
; Isaac Roseboom, igr@roe.ac.uk
;
; HISTORY
; 08JUL2011: Initial committed version

PRO segsrc, sx,sy,pfwhm250,pfwhm350,pfwhm500,hdr250, map250,sig250 $
                   ,hdr350, map350,sig350 $
                   ,hdr500, map500,sig500,ssrc,bsn=bsn
                        
     
; set SN threshold below which we will segment map
if(not keyword_set(bsn)) then bsn=0.5


prf250=psf_gaussian(fwhm=pfwhm250,npix=[11,11])


; build blend map
naxis1=sxpar(hdr250,'naxis1')
naxis2=sxpar(hdr250,'naxis2')

bim=fltarr(naxis1,naxis2)
print,'build blend map'

; 250 micron SN map

bim[sx,sy]=1.
flat=fltarr(49,49)+1.
cb250=convol(bim,flat) ; only interested in regions where there are sources (need this for speed)
bad=where(cb250 gt 1)
cb250[bad]=1.
sm250=cb250*map250/sig250
; blank is anything below some arbitrary SN threshold
blank250=where(sm250 le bsn, nblank250,comp=source250)
sm250[blank250]=-100.
smap250=intarr(naxis1,naxis2)
smap250[blank250]=-1
smap250[source250]=0


;350 micron SN map




sm350=map350/sig350
; move to the same pixel scale as 250 micron 
hastrom,sm350,hdr350,smap350,chdr350,hdr250,interp=1
smap350=smap350*cb250
blank350=where(smap350 le bsn, nblank350,comp=source350)
smap350[blank350]=-1.
smap350[source350]=0

;500 micron SN map



sm500=map500/sig500
; move to the same pixel scale as 250 micron 
hastrom,sm500,hdr500,smap500,chdr500,hdr250,interp=1
smap500=smap500*cb250
blank500=where(smap500 le bsn, nblank500,comp=source500)
smap500[blank500]=-1
smap500[source500]=0



; need to make sure 'dead' (i.e. low S/N) pixels are not significant
; in any of the maps
smap=smap250*smap350*smap500
lim=where(smap ne 0)
smap[lim]=-1.

tacc=1
snum=0L

; flood fill each segment
while tacc ne 0 do begin
   ; find first unallocated pixel
   tp=where(smap eq 0,tacc)
   ; fill that pixel with current segment
   ; number
   snum=snum+1L

; try and fill surrounding pixels
   sp=[tp[0]]
   fill=where(smap[sp] eq 0.,fm)
   k=0
   stnum=snum
   while fm ne 0 do begin
      ; check for mergers
      merge=where(smap[sp] ne stnum and smap[sp] ne 0 and smap[sp] ne -1)
      if(merge[0] ne -1) then begin
      remark=where(smap eq stnum)
      stnum=smap[sp[merge[0]]]
      smap[remark]=stnum
      endif
      smap[sp[fill]]=stnum
     
      if(smap[sp[0]] ne -1)then sp=[sp[0]-1,sp]
      if(smap[sp[n_elements(sp)-1]] ne -1) then sp=[sp,sp[n_elements(sp)-1]+1]
      fill=where(smap[sp] eq 0 ,fm)      
      
    
       if fm eq 0 then begin
          
          good=where(smap[sp] eq stnum)
          sp=sp[good]+naxis1
          
          fill=where(smap[sp] ge 0,fm) 
      
       endif
      
    endwhile
  
endwhile
; find segment for each source
n_src=n_elements(sx)
segsrc=lonarr(n_src)
for i=0L,n_src-1 do begin
    segsrc[i]=smap[sx[i],sy[i]]
    rad=1
    while(segsrc[i] le 0) do begin
        xrad=sx[i]+lindgen(2*rad+1)-rad
        yrad=sy[i]+lindgen(2*rad+1)-rad
        make_2d,xrad,yrad,dx,dy
        tmap=smap[dx,dy]
       tst=where(tmap gt 0)
       if(tst[0] ne -1) then segsrc[i]=tmap[tst[0]]
       rad=rad+1
   endwhile

endfor

end
