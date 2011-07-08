;+
; NAME:
; modselLASSO
; 
; PURPOSE:
; Solve the linear problem Ax=b for x using weighted non-negative LASSO algorithm
; (Zou et al. 2006). This works via a active set algorithm, initially
; all paramaters (x) are set to zero (dormant) and then parameters are
; switched to the active set iteratively by examining the gradient of the chi^2 
;
; CALLING SEQUENCE
; modselLASSO,apmat,bvec,psrc,active,xvec,snsrc,f_src=f_src,thresh=thresh
; INPUT PARAMETERS
; apmat: matrix A, in this application this is a matrix of
; size number of sources time number of sources and encodes the
; correlations between the sources due to overlapping beams 
; 
; bvec: RHS of linear eqn Ax=b, basically the PRF convolved flux at
; the position of each of the sources 
; 
; psrc: the parameter weights, NOTE this should not be the inverse
; variance in the pixels (if you want to weight by this divide both
; amat and bvec prior to entry to this code). These weights help
; determine which parameters (sources) are made active and which are excluded by
; from the solution. In the XID code psrc is set to the phi factor
; from Roseboom et al. (2010), i.e. for source i, phi_i=(number
; density of sources brighter than f24_i)/(SPIRE beam size in this band)
; 
; snsrc: the number of parameters (sources)
;
; thresh: a maximum level for the norm of the solution, this can be
; useful if we have some idea of what the total of the solution vector
; should be (i.e. CIB constraints). If in doubt set this to a very big
; number, or leave alone
;
; quiet: set this flag to not get iteration and timing information
;
; OUTPUTS
; xvec: the solution vector
; active: the indices of which parameters were active in the final
; solution
;
; AUTHORS
; Isaac Roseboom, igr@roe.ac.uk
;
; HISTORY
; 08JUL2011: Initial committed version


pro modselLASSO,apmat,bvec,psrc,active,xvec,snsrc,thresh=thresh,quiet=quiet
ilthresh=200



; adjust pmat for weights psrc

pmat=dblarr(snsrc,snsrc)
for i=0,snsrc-1 do begin
  pmat[i,*]=apmat[i,*]/psrc[i]
endfor



spmat=sprsin(transpose(pmat))




; determine initial negative gradient
w=dblarr(snsrc)
xvec=dblarr(snsrc)
ind=indgen(snsrc)
res=bvec-matrix_multiply(pmat,xvec,/atranspose)
pres=res
w=-matrix_multiply(pmat,res)
omid=-1
dead=indgen(snsrc)
nact=0
lambda=-min(w,dind)
JB=total(abs(xvec))
if(not keyword_set(thresh)) then thresh=1e30
status=intarr(snsrc)
status[dind]=1
active=where(status eq 1,nact,comp=inactive)

xxt=matrix_multiply(pmat,pmat,/btranspose)
t1=systime(/sec)
niter=0L
id2=-1
il=0
mid=-1

while(JB lt thresh and lambda gt 1e-30 and niter lt 10.*snsrc) do begin
ttmp2=systime(/sec)
gamma=fltarr(snsrc)
tmp=active
make_2d,tmp,tmp,ina,inb
xxta=double(xxt[ina,inb])
if nact gt 1 then begin 
    ones=(fltarr(nact)+1.)
    ttmp=systime(/sec)


DEFSYSV, '!XLIB', EXISTS = usexlib
if(usexlib eq 1) then begin
; IF USING MP MACHINE AND HAVE libXID.so USE THESE LINES  
;--------------  
    ntmp=long(nact)
    blah=call_external('libXID.so','idlinvert_',$
             ntmp,xxta,value=[1,0])
    gammaA=xxta#ones

;-------------  
endif else begin
; IDL ONLY LINE
;-------------
   gammaA=la_invert(xxta)#ones
;------------
endelse



  ix=(systime(/sec)-ttmp)
endif else begin
    gammaA=1./xxta
    ix=0.
endelse
gamma[active]=gammaA


d1=-1
if(nact ne snsrc) then d1=(lambda+w[inactive])/(1.-(matrix_multiply(xxt[inactive,*],gamma)))
d2=-xvec[active]/gamma[active]

d=lambda
db=lambda
index=0
indexb=0
id1=-1

td1=where(d1 gt 1e-20 and d1 le d and ind[inactive] ne omid and status[inactive] ne -1)
if(td1[0] ne -1) then begin
 sd=sort(d1[td1])   
 d=d1[td1[sd[0]]]
 id1=inactive[td1[sd[0]]]
 index=1
 if(n_elements(td1) ge 2) then begin
     db=d1[td1[sd[1]]]
     idb=inactive[td1[sd[1]]]
     indexb=1
 endif
 
 
endif
id2=-1
td2=where(d2 gt 1e-20 and d2 lt lambda)
if(td2[0] ne -1) then begin
    sd=sort(d2[td2])   
 if(d2[td2[sd[0]]] lt d) then begin
 d=d2[td2[sd[0]]]
 id2=active[td2[sd[0]]]
 index=2
endif
if(n_elements(td2) ge 2)then begin
 if( d2[td2[sd[1]]] lt db ) then begin
     db=d2[td2[sd[1]]]
     idb=active[td2[sd[1]]]
     indexb=2
     endif
 endif
endif

; need a catch for cases where the gradient is already negative at 0
vbad=where(d2 eq 0 and gamma[active] lt 0) 
if(vbad[0] ne -1)then begin
;d=0
;id2=active[vbad]
status[active[vbad]]=0
;index=2
endif 



wtf=where(xvec lt 0)
omid=mid
if(index eq 1) then begin
    status[id1]=1
    mid=id1
endif
if(index eq 2) then begin
    status[id2]=0
    mid=id2
endif
;check for loops
 if(mid eq omid) then begin
     il=il+1
 endif else begin
     il=1
 endelse

if(il eq ilthresh) then begin
; just take best one loop and get on with it
bst=min([psrc[id1],psrc[id2]],bid)
if(bid eq 0) then begin
status[id1]=0
status[id2]=-1
endif else begin
status[id1]=-1
status[id2]=0
endelse
; do next best step

d=db

if(indexb eq 1) then status[idb]=1
if(indexb eq 2) then status[idb]=0



endif

xvec=xvec+d*gamma
active=where(status eq 1,nact,comp=inactive)
if(nact ne snsrc) then xvec[inactive]=0.
bx=fltarr(snsrc)
ttmp=systime(/sec)
tpmat=pmat[active,*]
bx=matrix_multiply(tpmat,xvec[active],/atranspose)

res=bvec-bx

w=-sprsax(spmat,res,/double)
rw=systime(/sec)-ttmp
JB=total(abs(xvec))


lambda=lambda-d

    
   
  niter=niter+1
etmp=systime(/sec)-ttmp2
if (niter mod 500 eq 0 and not keyword_set(quiet)) then print,'iter ',niter,il,ix,rw,(etmp-ix-rw),index,id1,id2,lambda,d

endwhile



xvec=xvec/psrc

if not keyword_set(quiet) then print,'finished ',(systime(/sec)-t1)

end
