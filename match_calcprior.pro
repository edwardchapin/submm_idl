;------------------------------------------------------------------------------
; match_calcprior: Calculate priors for matching based on catalogues
;
; The purpose of this routine is to take a list of sources (the
; primary catalogue) for which you would like to identify matches, and
; an external catalogue (the matching catalogue) in which you would
; like to search for those matches, and establish priors (e.g. radial
; offset distribution, typical brightnesses and colours in matching
; catalogue etc.). The method used is to search for systematic
; differences in the matching catalogue around primary catalogue
; positions compared to the rest of the catalogue.
;
; This routine is based on techniques outlined in Chapin et al. (2011)
; MNRAS 411 505.
;
; Required Inputs:
;
;   p_ra       = array of primary catalogue right ascensions in degrees
;   p_dec      = "          "     "         declinations "
;   m_ra       = array of matching catalogue right ascensions in degrees
;   m_dec      =  "         "      "         declinations "
;   m_mask     = 2d array with 1 everywhere the matching catalogue has coverage
;   m_maskhead = FITS header of mask
;   rmax       = maximum search radius for finding counterparts in arcsec
;   nr         = number of steps in search radius to consider
;
; Optional Keyword Inputs:
;
;   showplot   = If set, show plots to depict progress
;   rprior     = Override internally-estimate ideal search radius for priors
;   e1         = use this value if rprior set (see Optional Keyword Outputs)
;   sig1       = use this value if rprior set (see Optional Keyword Outputs)
;   e2         = use this value if rprior set (see Optional Keyword Outputs)
;   sig2       = use this value if rprior set (see Optional Keyword Outputs)
;   rayleigh2  = use sum of 2 rayleigh dist (see e2, sig2)
;   noradial   = If set, and rprior set, skip radial distribution measurement
;   properties = N primary catalogue entries X M properties for priors
;   propbins   = M X 2: bin_min, delta_bin for binned priors arrays
;   propdims   = M array lengths for output binned priors arrays
;   proplabels = array of string labels for each property
;   p_ind      = index to subset of primaries to analyze after mask made
;   postscript = disables plotting to x in case you're dumping to ps instead
;   verbose    = controls verbosity (0 or 1, default 1)
;
; Optional Keyword Outputs:
;
;   lr_sr      = 99% search radius for LR
;   po_bincen  = pointer array to bin centres for each property prior
;   po_fg      = pointer array to foreground (match) priors for each property
;   po_fg_err  = pointer array to errors in p_fg
;   po_bg      = pointer array to background priors for each property
;   po_bg_err  = pointer array to errors in p_bg
;   e1         = average number of matches to primaries expected
;   sig1       = sigma for Rayleigh radial distribution of matches
;   e2         = if rayleigh2 set, integral of the second component
;   sig2       = "     "       "   width of the second component
;   avmatch    = av. matches expected per primary (e1 or e1+e2 see rayleigh2)
;   q          = match prior distribution (integrates to avmatch)
;   bdist      = background number density per arcsec^2
;   success    = success flag (1=success)
;   errmsg     = string containing error message if abnormal exit
;   suggest_rmax = suggest rmax based on 95% peak of cumulative radial dist.
;   bgmask     = mask indicating where background sources were selected
;
; Dependencies:
;
;   tvscale.pro: Coyote (http://www.idlcoyote.com/programs/tvscale.pro)
;   IDL astro library: (http://idlastro.gsfc.nasa.gov/)
;
; Authors:
;   Edward Chapin, echapin@phas.ubc.ca (EC)
;
; History:
;   27APR2011: Initial version -- only radial dist, and no outputs yet (EC)
;   28APR2011: Add showplot, rprior, noradial (EC)
;              Add properties, propbins, propdims
;   29APR2011: Check all where statements generate error messages (EC)
;   02MAY2011: Combine property priors independently and return (EC)
;              Add suggest_rmax
;              Implement p_ind
;   02MAY2011: Fixed a minor bug for when only one prior is present,
;              cleaned up some of the informational i/o to allow for
;              'verbose' and 'postscript' keywords. (MZ)
;   06MAY2011: Add bgmask output (EC)
;   19MAY2011: Do better job of excess integrals vs. radius (EC)
;   08SEP2011: Fix offset bug in radial integral, add rayleigh2 (EC)
;
;------------------------------------------------------------------------------


;------------------------------------------------------------------------------
; This is a local function used for fitting rayleigh radial offset distribution
;------------------------------------------------------------------------------

function fitrayleigh, m, model_excess=model_excess
common fitrayleigh_block, all_r, all_rcen, excess, excess_err,$
   dist_lowres, dist_highres, dist_res_lookup

e1 = m[0]
sig1 = m[1]
n = n_elements(all_rcen)

chisq = 0d
all_y = dblarr(n)

; integrate model across each excess bin and calculate contribution to chi^2


; old way assumes that blocky 2d dist array is good
; approximation... it isn't
;nint = 10d
;int_bin = loggen(0,all_r[1], nint+1,/linear)
;int_dbin = int_bin[1:nint] - int_bin[0:nint-1]
;int_cen = (int_bin[1:nint] + int_bin[0:nint-1])/2d

;for i=0, n-1 do begin
;  r = all_r[i]+int_cen
;  dr = int_dbin

;  y = e1*total( (r/sig1^2d)*exp(-r^2d/(2d*sig1^2d))*dr )
;  all_y[i] = y

;  chisq = chisq + (excess[i]-y)^2d/(excess_err[i])^2d
;endfor

;all_y_old = all_y

; new way integrates across low-res pixels
da = dist_highres[1]^2d ; area of high-resolution pixel
for i=0, n-1 do begin
  ; identify pixels in high-resolution dist that match the low-resolution
  ; dist binning
  ind = where( (dist_res_lookup ge all_r[i]) and $
               (dist_res_lookup lt all_r[i+1]) )

  ; integrate the 2d Gaussian over the relevant area of dist_highres
  y = e1*total( (1/(2*!dpi*sig1^2d)) * $
                     exp( -dist_highres[ind]^2d/(2d*sig1^2d) ) * da )

  all_y[i] = y
  chisq = chisq + (excess[i]-y)^2d/(excess_err[i])^2d
endfor

model_excess = all_y

chisq = chisq/(n-2)

if (sig1 le 0) or (e1 le 0) then chisq=1000d

;print, chisq,':', m

return, chisq
end

;------------------------------------------------------------------------------
; Two Rayleigh distributions to allow a larger tail
;------------------------------------------------------------------------------

function fitrayleigh2, m, model_excess=model_excess
common fitrayleigh_block, all_r, all_rcen, excess, excess_err, $
  dist_lowres, dist_highres, dist_res_lookup

e1 = m[0]
sig1 = m[1]
e2 = m[2]
sig2 = m[3]

n = n_elements(all_rcen)

chisq = 0d
all_y = dblarr(n)

; high-resolution 2d positional distribution
pdf2d = e1*(1/(2*!dpi*sig1^2d))*exp(-dist_highres^2d/(2d*sig1^2d)) + $
        e2*(1/(2*!dpi*sig2^2d))*exp(-dist_highres^2d/(2d*sig2^2d))

; integrate across low-res pixels
da = dist_highres[1]^2d ; area of high-resolution pixel

for i=0, n-1 do begin
  ; identify pixels in high-resolution dist that match the low-resolution
  ; dist binning
  ind = where( (dist_res_lookup ge all_r[i]) and $
               (dist_res_lookup lt all_r[i+1]) )

  ; integrate the two 2d Gaussians over the relevant area of dist_highres
  y = total(pdf2d[ind]*da)

  all_y[i] = y
  chisq = chisq + (excess[i]-y)^2d/(excess_err[i])^2d
endfor

model_excess = all_y

chisq = chisq/(n-4)

if (sig1 le 0) or (e1 le 0) or (sig2 le 0) or (e2 le 0) then chisq=1000d

;if chisq le 0.92 then stop

;print, chisq,':', m

return, chisq
end

;-----------------------------------------------------------------------------
; Main routine
;-----------------------------------------------------------------------------

pro match_calcprior, p_ra, p_dec, m_ra, m_dec, m_mask, maskhead, rmax, nr, $
                     showplot=showplot, $
                     rprior=rprior, noradial=noradial, properties=properties, $
                     propbins=propbins, propdims=propdims, $
                     proplabels=proplabels, p_ind=p_ind, $
                     rayleigh2=rayleigh2, $
                     e1=e1, sig1=sig1, $
                     e2=e2, sig2=sig2, avmatch=avmatch, $
                     po_bincen=po_bincen, po_fg=po_fg, po_fg_err=po_fg_err, $
                     po_bg=po_bg, po_bg_err=po_bg_err, $
                     q=q, bdist=bdist, $
                     POSTSCRIPT=postscript,$
                     VERBOSE=verbose,SUCCESS=success,$
                     ERRMSG=errmsg, suggest_rmax=suggest_rmax, $
                     bgmask=bgmask

common fitrayleigh_block, all_r, all_rcen, excess, excess_err, $
  dist_lowres, dist_highres, dist_res_lookup

  COMPILE_OPT IDL2, STRICTARRSUBS
  success = 0b

  IF ~(N_ELEMENTS(verbose)) THEN verbose=1 ELSE verbose=verbose
  IF ~(N_ELEMENTS(postscript)) THEN postscript=0 ELSE postscript=postscript

  ; characters in greek font set
  greek_sigma = 'r'
  greek_chi = 'v'

  ; Calculate a cropped mask that only contains the region than is nonzero.
  ; We also has a border that is rmax  in size around the edges to avoid
  ; wraparound problems when we make masks.
  ind = where(m_mask)
  if ind[0] eq -1 then begin
    errmsg="Error: problem finding good regions in m_mask"
    IF verbose THEN MESSAGE,errmsg,/INFORMATIONAL
    RETURN
  endif

  width = n_elements(m_mask[*,0])
  height = n_elements(m_mask[0,*])
  xmask = ind mod width
  ymask = ind / width

  pixres = sxpar( maskhead, "CD1_1")
  if pixres eq 0 then pixres = sxpar( maskhead, "CDELT1" )
  pixres = abs( pixres )*3600d                   ; arcsec/pixel

  xmin = round(min(xmask) - rmax/pixres) > 0
  xmax = round(max(xmask) + rmax/pixres) < (width-1)
  ymin = round(min(ymask) - rmax/pixres) > 0
  ymax = round(max(ymask) + rmax/pixres) < (height-1)

  mask = m_mask[xmin:xmax,ymin:ymax]

  ; create an updated FITS header to go with cropped mask
  header = maskhead
  nx = xmax-xmin+1
  ny = ymax-ymin+1

  CRPIX1 = sxpar(header,"CRPIX1")
  CRPIX2 = sxpar(header,"CRPIX2")
  CRPIX1 = CRPIX1 - xmin
  CRPIX2 = CRPIX2 - ymin

  sxaddpar, header, 'NAXIS1', nx
  sxaddpar, header, 'NAXIS2', ny
  sxaddpar, header, 'CRPIX1', CRPIX1
  sxaddpar, header, 'CRPIX2', CRPIX2

  ; Calculate area of the entire catalogue
  pixres = sxpar( header, "CD1_1")
  if pixres eq 0 then pixres = sxpar( header, "CDELT1" )
  pixres = abs( pixres )*3600d                   ; arcsec/pixel

  ind = where(mask)
  if ind[0] eq -1 then begin
    errmsg = "Error: no good pixels in mask"
    IF verbose THEN MESSAGE,errmsg
    RETURN
  endif
  area = n_elements( ind ) * pixres^2d   ; area in arcsec^2

  ; project the catalogues on to the mask, and determine the x, y positions
  ; of sources that land within it
  cat_pix, p_ra, p_dec, p_x, p_y, header
  cat_pix, m_ra, m_dec, m_x, m_y, header

  p_onmap = where( (p_x ge 0) and (p_x lt nx) and $
                   (p_y ge 0) and (p_y lt ny) )

  if p_onmap[0] eq -1 then begin
    errmsg="Error: no primary sources with valid map locations"
    IF verbose THEN MESSAGE,errmsg
    RETURN
  endif

  m_onmap = where( (m_x ge 0) and (m_x lt nx) and $
                   (m_y ge 0) and (m_y lt ny) )

  if p_onmap[0] eq -1 then begin
    errmsg="Error: no matching sources with valid map locations"
    IF verbose THEN MESSAGE,errmsg
    RETURN
  endif

  pkeep = -1 & mkeep = -1

  if p_onmap[0] ne -1 then pkeep = where( mask[p_x[p_onmap],p_y[p_onmap]] )
  if m_onmap[0] ne -1 then mkeep = where( mask[m_x[m_onmap],m_y[m_onmap]] )

  if (pkeep[0] eq -1) or (mkeep[0] eq -1) then begin
    errmsg="Error: either no primary or matching catalogue sources land on mask"
    IF verbose THEN MESSAGE,errmsg
    RETURN
  endif

  p_x = round(p_x[p_onmap[pkeep]])
  p_y = round(p_y[p_onmap[pkeep]])

  m_x = round(m_x[m_onmap[mkeep]])
  m_y = round(m_y[m_onmap[mkeep]])
  if keyword_set(properties) then prop = properties[m_onmap[mkeep],*]

  n_p = n_elements(p_x)     ; number of primary sources on mask
  n_m = n_elements(m_x)     ; number of matching sources on mask

  ; Set up a distance array (in arcsec) with the same dimensions as the mask
  dim = max([nx,ny])
  d = shift(dist(dim),dim/2,dim/2)*pixres
  d = d[dim/2-nx/2:dim/2-nx/2+nx-1, $
         dim/2-ny/2:dim/2-ny/2+ny-1]
  d = shift(d,-nx/2,-ny/2)

  ; create masks: first cut out circles at the maximum search radius
  ; centered over the primary sources -- this will be used to
  ; establish a "background" mask to find sources from the matching
  ; catalogue that are not contaminated by the primary sources we are
  ; studying; then a "foreground" mask that is the area that is within
  ; the maximum search radius from _only one_ primary source

  ind = where( d le rmax )
  if ind[0] eq -1 then begin
    errmsg="Error: requested circle is too small"
    IF verbose THEN MESSAGE,errmsg
    RETURN
  endif

  circle = intarr(nx,ny)
  circle[ind] = 1

  bgmask = intarr(nx,ny)

  p_delta = dblarr(nx,ny)
  for i=0, n_p-1 do p_delta[p_x[i],p_y[i]] = p_delta[p_x[i],p_y[i]] + 1
  temp = round(double(fft(fft(circle,1)*fft(p_delta,1),-1))) ; convolve

  fg = where( (temp eq 1) and (mask) )
  if fg[0] eq -1 then begin
    errmsg="Error: foreground mask is empty"
    IF verbose THEN MESSAGE,errmsg
    RETURN
  endif
  fgmask = intarr(nx,ny)
  fgmask[fg] = 1

  ; now that we have fgmask we can look at a subset of the primary positions
  ; if p_ind was supplied
  if keyword_set(p_ind) then begin
    ; re-project subset of positions on to map
    cat_pix, p_ra[p_ind], p_dec[p_ind], p_x, p_y, header

    p_onmap = where( (p_x ge 0) and (p_x lt nx) and $
                     (p_y ge 0) and (p_y lt ny) )

    if p_onmap[0] eq -1 then begin
      errmsg="Error: no primary sources within p_ind with valid map locations"
      IF verbose THEN MESSAGE,errmsg
      RETURN
    endif

    if p_onmap[0] eq -1 then begin
      errmsg="Error: no matching sources within p_ind with valid map locations"
      IF verbose THEN MESSAGE,errmsg
      RETURN
    endif

    pkeep = -1

    if p_onmap[0] ne -1 then pkeep = where( mask[p_x[p_onmap],p_y[p_onmap]] )

    if pkeep[0] eq -1 then begin
      errmsg="Error: no primary catalogue source within p_ind lands on mask"
      IF verbose THEN MESSAGE,errmsg
      RETURN
    endif

    p_x = round(p_x[p_onmap[pkeep]])
    p_y = round(p_y[p_onmap[pkeep]])
    n_p = n_elements(p_x)       ; new number of primary sources
    p_delta = dblarr(nx,ny)
    for i=0, n_p-1 do p_delta[p_x[i],p_y[i]] = p_delta[p_x[i],p_y[i]] + 1
  endif

  bg = where( (temp eq 0) and (mask) )
  if bg[0] eq -1 then begin
    errmsg="Error: background mask is empty"
    IF verbose THEN MESSAGE,errmsg
    RETURN
  endif
  bgmask = intarr(nx,ny)
  bgmask[bg] = 1

  if keyword_set(showplot) then begin
     LOADCT,0
     IF ~postscript THEN BEGIN
        set_plot,'x'
        window, 0
        erase = 0
     ENDIF ELSE BEGIN
        erase = 1
     ENDELSE
    tvscale, bgmask, /keep, /noint, ERASE=erase, BACKGROUND='white'
    xyouts, 0.93, 0.5, 'Background Mask', charsize=1.5, charthick=3.0, $
            /normal, color=255,ORIENTATION=90
    xyouts, 0.04, 0.5, 'Background Mask', charsize=1.5, /normal,ORIENTATION=90

    IF ~postscript THEN window,1
    tvscale, fgmask, /keep, /noint, ERASE=erase, BACKGROUND='white'
    xyouts, 0.93, 0.5, 'Foreground Mask', charsize=1.5, charthick=3.0, $
            /normal, color=255,ORIENTATION=90
    xyouts, 0.04, 0.5, 'Foreground Mask', charsize=1.5, /normal,ORIENTATION=90
  endif

  ; calculate area and surface density + Poisson error, of sources in bg mask
  area_bg = n_elements(bg) * pixres^2d             ; arcsec^2
  ind = where( bgmask[m_x,m_y] )
  if ind[0] eq -1 then begin
    errmsg="Error: no matching sources land on background mask"
    IF verbose THEN MESSAGE,errmsg
    RETURN
  endif
  n_bg = n_elements(ind)
  rho_bg = n_bg / area_bg                       ; #/arcsec^-2
  rho_bg_err = sqrt(n_bg) / area_bg

  ; Calculate bins of search radius
  all_r = loggen( 0, rmax, nr+1, /linear)                ; bin edges
  all_dr = (all_r[1:nr] - all_r[0:nr-1])                 ; bin widths
  all_rcen = (all_r[0:nr-1] + all_r[1:nr])/2. -pixres/2  ; real radial bincen

  ; 2-D distance array, high-res distance array for integrals, and
  ; a high-res lookup table to match low-res to high-res pixels

  res_scale = 20. ; use 10 x higher resolution for integration
  dist_lowres = dist(nr*2)*all_dr[0]
  dist_highres = dist(nr*2*res_scale)*all_dr[0]/res_scale
  dist_res_lookup = rebin(dist_lowres, nr*2*res_scale, nr*2*res_scale, /sample)

  ; this shift is needed to centre things properly for the low-res pixels
  dist_res_lookup = shift(dist_res_lookup,-res_scale/2,-res_scale/2)

  ; Check for step sizes that are smaller than map resolution. The 1.01
  ; just covers cases where it's really close.
  if all_dr[0]*1.01 lt pixres then begin
    errmsg = "Radial step size" + strcompress(all_dr[0]) + $
           "(arcsec) is smaller than mask resolution " + $
           strcompress(pixres) + "(arcsec)" + $
           "increase rmax, or decrease nr."
    IF verbose THEN MESSAGE,errmsg
    RETURN
  endif

  if (keyword_set(noradial) eq 0) then begin

    ; loop over search radius, and calculate a mask that has annuli
    ; centered over every primary catalogue position, and also lands
    ; within the fgmask. Then count the number of sources (and Poisson
    ; error) within this mask, as well as the expected number based on
    ; rho_bg and the area

    counts = dblarr(nr)
    counts_err = dblarr(nr)
    expect = dblarr(nr)
    expect_err = dblarr(nr)
    excess = dblarr(nr)
    excess_err = dblarr(nr)

    all_n_full = dblarr(nr)
    all_n_clipped = dblarr(nr)

    if keyword_set(showplot) and ~postscript then window, 2

    for i=0, nr-1 do begin

      ind = where( (d ge all_r[i]) and (d lt all_r[i+1]) )
      if ind[0] eq -1 then begin
        errmsg="Error: no pixels in annulus"
        IF verbose THEN MESSAGE,errmsg
        RETURN
      endif

      annulus = intarr(nx, ny)
      annulus[ind] = 1

      amask = round(double(fft(fft(annulus,1)*fft(p_delta,1),-1))) ; convolve

      ; To calculate the average excess per primary source, we have to
      ; account for the fact that of the n_p sources, we have an
      ; effectively smaller number because we are only considering the
      ; portion of the annuli that land within both the total mask, and
      ; the fgmask (no contamination from other nearby primary sources).
      ; To do this we count the pixels in amask before and after
      ; applying fgmask, and multiply n_p by that ratio.

      ind = where(amask)
      if ind[0] eq -1 then begin
        errmsg="Error: no pixels in annulus mask"
        IF verbose THEN MESSAGE,errmsg
        RETURN
      endif
      n_full = double(n_elements(ind))
      amask = amask * fgmask
      ind = where(amask)
      if ind[0] eq -1 then begin
        errmsg="Error: no pixels in foreground-clipped annulus mask"
        IF verbose THEN MESSAGE,errmsg
      endif
      n_clipped = double(n_elements(ind))

      n_p_eff = (n_clipped / n_full) * n_p

      all_n_full[i] = n_full       ; store for calculating cumulative
      all_n_clipped[i] = n_clipped ; equivalents later

      if keyword_set(showplot) then begin
        tvscale, amask, /keep, /noint, ERASE=erase, BACKGROUND='white'
        xyouts, 0.93, 0.5, 'FG Masked Measurement Annuli', charsize=1.5, $
                charthick=3.0, /normal, color=255,ORIENTATION=90
        xyouts, 0.04, 0.5, 'FG Masked Measurement Annuli', charsize=1.5, $
                /normal,ORIENTATION=90
      endif

      ; count matching sources in annulus mask
      ind = where( amask[m_x,m_y] )
      if ind[0] ne -1 then begin
        counts[i] = n_elements(ind)
        counts_err[i] = sqrt(counts[i])
      endif

      ; how many were expected?
      ind = where(amask)
      if ind[0] eq -1 then begin
        errmsg,"Error: no pixels in annulus mask"
        IF verbose THEN MESSAGE,errmsg
        RETURN
      endif
      expect[i] = n_elements(ind)*pixres^2d * rho_bg
      expect_err[i] = (rho_bg_err/rho_bg)*expect[i]

      ; what is the average excess given effective number of sources?
      excess[i] = (counts[i]-expect[i]) / n_p_eff
      excess_err[i] = (counts_err[i]+expect_err[i])/(counts[i]-expect[i]) ; rel
      excess_err[i] = excess_err[i]*excess[i]                             ; abs

      effstr = " effective sources: " + string( n_p_eff, format='(F6.1)' )  + $
               ' / ' + strcompress(n_p,/remove_all)

      IF verbose THEN MESSAGE, 'Annulus:'+ $
         string(all_r[i],format='(I3)') + '--' + $
         string(all_r[i+1],format='(I3)') + effstr + ' excess=' + $
         string(excess[i],format='(F5.2)') + ' +/- ' + $
         string(excess_err[i],format='(F5.2)'),/INFORMATIONAL

    endfor

    ; fit rayleigh distribution to excess counts
    e1_guess = total(excess)
    ind = (where(excess eq max(excess)))[0]
    if ind[0] eq -1 then begin
      errmsg="Error: error finding guess sigma for Rayleigh distribution"
      IF verbose THEN MESSAGE,errmsg
      RETURN
    endif
    sigma_guess = all_rcen[ind]

    r_highres = loggen(0,rmax,100,/linear) ; for evaluating high-res PDF

    if keyword_set(rayleigh2) then begin
      ; --- 2-component Rayleigh distribution ---
      m_init = [e1_guess*0.5,sigma_guess,e1_guess*0.5,sigma_guess*2.]
      m_scale = 0.3*m_init

      fit = amoeba( 1d-6, function_name='fitrayleigh2', p0=m_init, $\
                     scale=m_scale )
      if fit[0] eq -1 then begin
        errmsg="Error: couldn't fit double Rayleigh distribution to excesses"
        IF verbose THEN MESSAGE,errmsg
        RETURN
      endif

      e1 = fit[0]
      sig1 = fit[1]
      e2 = fit[2]
      sig2 = fit[3]

      avmatch = e1 + e2

      ; high-res PDF
      pdf1 = e1*(r_highres/sig1^2d)*exp(-r_highres^2d/(2d*sig1^2d))
      pdf2 = e2*(r_highres/sig2^2d)*exp(-r_highres^2d/(2d*sig2^2d))

      pdf = pdf1 + pdf2

      ; low-res
      chisq = fitrayleigh2(fit,model_excess=model_excess)

    endif else begin
    ; --- single-component rayleigh distribution ---
      m_init = [e1_guess,sigma_guess]
      m_scale = 0.3*m_init

      fit = amoeba( 1d-6, function_name='fitrayleigh', p0=m_init, $\
                    scale=m_scale )
      if fit[0] eq -1 then begin
        errmsg="Error: couldn't fit Rayleigh distribution to excesses"
        IF verbose THEN MESSAGE,errmsg
        RETURN
      endif

      e1 = fit[0]
      sig1 = fit[1]

      avmatch = e1

      ; high-res sampled version of the smooth distribution
      pdf = e1*(r_highres/sig1^2d)*exp(-r_highres^2d/(2d*sig1^2d))

      ; ...and also the low-resolution binned version from which chi^2
      ; was calculated
      chisq = fitrayleigh(fit,model_excess=model_excess)
    endelse

    ; note that all_r refers to low-resolution pixel edges. For
    ; example, if the first range is from 0 to 1, it really contains
    ; sources that land within a box from -0.5 to +0.5 centered over
    ; the primary position. Similarly, the next "annulus" defined as 1
    ; to 2 is really radii 0.5 to 1.5. We defined r_cen earlier so
    ; that it lies at the centres of these shifted bins, i.e. 0, 1,
    ; 2...

    if keyword_set(showplot) then begin
       IF ~postscript THEN window,3
      plot, all_rcen, excess, psym=5, $
            xtitle='Search Radius (arcsec)', $
            ytitle='Excess / bin', xstyle=1, xrange=[0,rmax], $
            charsize=1.5
      errplot, all_rcen, excess-excess_err, excess+excess_err
      oplot,[0,rmax],[0,0]

      oplot, r_highres, pdf*all_dr[0] ; account for bin size
      if keyword_set(rayleigh2) then begin
        oplot, r_highres, pdf1*all_dr[0], linestyle=1
        oplot, r_highres, pdf2*all_dr[0], linestyle=1
      endif

      oplot_hist, all_r-pixres/2., model_excess


      modelstr = '!4'+greek_sigma+'!3!D1!n='+string(sig1,format='(F4.1)') + $
                 ' E!d1!n='+string(e1,format='(F3.1)')
      if keyword_set(rayleigh2) then $
        modelstr = modelstr +  ' !4'+greek_sigma+'!3!D2!n='+ $
                   string(sig2,format='(F4.1)') + $
                   ' E!d2!n='+string(e2,format='(F3.1)')

      modelstr = modelstr + $
                 ' !4'+greek_chi+'!3!u2!N='+string(chisq,format='(F4.1)')

      xyouts, 0.4, 0.85, modelstr, charsize=1.5, charthick=!p.thick, /normal
    endif

    ; now calculate the cumulative distribution
    cumexcess = dblarr(nr)      ; binned cumulative dist
    cumexcess_err = dblarr(nr)  ; binned cumulative error
    cumexcess_model = dblarr(nr); forward-modeled binned cumulative dist
    cumr = all_r[1:nr]-pixres/2.; radii corresponding to cumexcess

    for i=0, nr-1 do begin
      ; calculate n_full and n_clipped for circles centered over primary
      ; sources by adding up values for sub-annuli, and then n_p_eff

      n_full = total(all_n_full[0:i])
      n_clipped = total(all_n_clipped[0:i])

      n_p_eff = (n_clipped / n_full) * n_p

      ; now calculate cumulative excess and its error
      counts_cum = total(counts[0:i])
      expect_cum = total(expect[0:i])

      counts_cum_err = sqrt(counts_cum)
      expect_cum_err = (rho_bg_err/rho_bg)*expect_cum

      cumexcess[i] = (counts_cum - expect_cum) / n_p_eff
      cumexcess_err[i] = (counts_cum_err + expect_cum_err) / $
                         (counts_cum - expect_cum )    ; relative error
      cumexcess_err[i] = cumexcess_err[i]*cumexcess[i] ; absolute error

      ; the chunky forward-model cumulative excess
      cumexcess_model[i] = total(model_excess[0:i])
    endfor

    ; identify a maximum search radius from the cumulative counts that
    ; can be used in future calls to match_calcprior -- the radius
    ; at which we first get 95% of the peak, and then go 50% larger

    maxcum = max(cumexcess)
    ind = where( cumexcess ge 0.95*maxcum )
    if ind[0] ne -1 then begin
      suggest_rmax = cumr[ind[0]]
    endif

    ; increase by 50% and round to multiple of step size in r
    suggest_rmax = suggest_rmax*1.5
    suggest_rmax = round(suggest_rmax / all_dr[0])*all_dr[0]

    ; smooth fitted model
    cum_highres = pdf*0
    cumr_highres = r_highres - (r_highres[1]/2.)
    for i=0, n_elements(cum_highres)-1 do $
      cum_highres[i] = total( pdf[0:i] )

    ; account for bin size
    cum_highres = cum_highres * (r_highres[1]-r_highres[0])

    if keyword_set(showplot) then begin
      IF ~postscript THEN window, 4
      plot, [0], [0], xtitle='Search Radius (arcsec)', $
            ytitle='Cumulative Excess (<R)', charsize=1.5, $
            xrange=[0,rmax], xstyle=1, yrange=[0,max(cumexcess)]

      oplot, cumr, cumexcess, psym=5
      errplot, cumr, cumexcess-cumexcess_err, cumexcess+cumexcess_err

      oplot, cumr_highres, cum_highres
      oplot_hist, all_r, cumexcess_model

    endif

    ; identify search radius in cumulative counts that maximizes SNR of
    ; the excess. This will define the search radius used to calculate
    ; other priors (assuming that they are uncorrelated with distance)

    cum_snr = cumexcess / cumexcess_err
    maxind = (where( cum_snr eq max(cum_snr) ))[0]
    if maxind[0] eq -1 then begin
      errmsg="Error: couldn't find max SNR in cumulative radial distribution"
      IF verbose THEN MESSAGE,errmsg
      RETURN
    endif
    IF verbose THEN MESSAGE, "Best SNR at search radius of " + $
                             string(all_rcen[maxind],format='(F5.1)') + $
                             " arcsec",$
                             /INFORMATIONAL
  endif else begin
    IF verbose THEN MESSAGE,$
       "Skipping measurement of radial offset distribution...",$
       /INFORMATIONAL
  endelse

  ; Work out avmatch
  if keyword_set(rayleigh2) then avmatch = e1 + e2 $
  else avmatch = e1

  ; use externally supplied rprior
  if keyword_set( rprior ) then begin
    maxind = where( all_rcen le rprior )
    maxind = (maxind[n_elements(maxind)-1])[0]
    if maxind[0] eq -1 then begin
      errmsg="Error: couldn't find radial bin matching rprior"
      IF verbose THEN MESSAGE,errmsg
      RETURN
    endif

    IF verbose THEN MESSAGE,$
       "Will measure priors based on supplied properties using " + $
       "search radius of "+ $
       string(all_rcen[maxind],format='(F5.1)') + " arcsec",$
       /INFORMATIONAL
  endif

  ; If properties supplied, work out other priors at best SNR radius
  if keyword_set(properties) and keyword_set(propbins) and $
    keyword_set(propdims) then begin

    ; How many properties were supplied?
    nprop = n_elements(propdims)

    ; Work out the region that is within the requested search radius of
    ; the primary catalogue positions, but also on the foreground mask

    ind = where( d lt all_r[maxind+1] )
    if ind[0] eq -1 then begin
      errmsg="Error: region for prior calculation is empty (too small?)"
      IF verbose THEN MESSAGE,errmsg
      RETURN
    endif
    circle = intarr(nx,ny)
    circle[ind] = 1

    pmask = round(double(fft(fft(circle,1)*fft(p_delta,1),-1))) ; convolve
    ind = where(pmask and fgmask)
    if ind[0] eq -1 then begin
      errmsg="Error: prior mask is empty"
      IF verbose THEN MESSAGE,errmsg
      RETURN
    endif
    pmask = pmask*0
    pmask[ind] = 1
    area_p = n_elements(where(pmask)) * pixres^2d ; arcsec^2

    if keyword_set(showplot) then begin
      IF ~postscript THEN window, 5
      tvscale, pmask, /keep, /noint, ERASE=erase, BACKGROUND='white'
      xyouts, 0.93, 0.5, 'Prior Mask', charsize=1.5, charthick=3.0, $
              /normal, color=255,ORIENTATION=90
      xyouts, 0.04, 0.5, 'Prior Mask', charsize=1.5, /normal,ORIENTATION=90
    endif

    ; set up arrays for storing the 1d marginalized prior probabilities
    p_bincen = ptrarr(propdims)
    p_fg = ptrarr(propdims)
    p_fg_err = ptrarr(propdims)
    p_bg = ptrarr(propdims)
    p_bg_err = ptrarr(propdims)

    ; loop over prior
    for i=0, nprop-1 do begin
       if keyword_set(proplabels) then $
          label = proplabels[i] else $
             label = 'Property'+strcompress(i+1)

      MESSAGE, "--- Calculating priors for "+label,/INFORMATIONAL
      ; Set up arrays for foreground and background counts as a function
      ; of property bins
      prop_bgcounts = dblarr( propdims[i] )
      prop_bgcounts_err = dblarr( propdims[i] )
      prop_fgcounts = dblarr( propdims[i] )
      prop_fgcounts_err = dblarr( propdims[i] )



      prop_bins = loggen( propbins[i,0], $
                          propbins[i,0]+propbins[i,1]*(propdims[i]+1), $
                          propdims[i]+1, /linear )
      prop_bincen = (prop_bins[0:propdims[i]-1] + $
                     prop_bins[1:propdims[i]])/2.
      dp = propbins[i,1]

      ; loop over bins and count sources with properties in range
      ; for the background and prior masks
      for j=0, propdims[i]-1 do begin
        ind = where( bgmask[m_x,m_y] and $
                     (prop[*,i] ge prop_bins[j]) and $
                     (prop[*,i] lt prop_bins[j+1]) )
        if ind[0] ne -1 then begin
          prop_bgcounts[j] = n_elements(ind)
          prop_bgcounts_err[j] = sqrt(prop_bgcounts[j])
        endif

        ind = where( pmask[m_x,m_y] and $
                     (prop[*,i] ge prop_bins[j]) and $
                     (prop[*,i] lt prop_bins[j+1]) )
        if ind[0] ne -1 then begin
          prop_fgcounts[j] = n_elements(ind)
          prop_fgcounts_err[j] = sqrt(prop_fgcounts[j])
        endif
      endfor

      ; Now work out the prior as the excess in each bin
      expect_p = (prop_bgcounts/area_bg) * area_p
      expect_p_err = (prop_bgcounts_err/prop_bgcounts) * expect_p

      prior_fg = prop_fgcounts - expect_p
      prior_fg_err = sqrt( prop_fgcounts_err^2d + expect_p_err^2d )

      ; normalize the foreground and background distributions to calculate
      ; probability densities for plotting
      scale = total(prior_fg*dp)
      probability_fg = prior_fg / scale
      probability_fg_err = prior_fg_err / scale

      scale = total(prop_bgcounts*dp)
      probability_bg = prop_bgcounts / scale
      probability_bg_err = prop_bgcounts_err / scale

      ; store prior results for this property -- ensuring that we
      ; have positive values by using smoothed values to help fill-in.

      bg_smooth = smooth(probability_bg>0,5)
      ind = where( probability_bg lt max(probability_bg)*0.01 )
      if ind[0] ne -1 then probability_bg[ind] = bg_smooth[ind]
      ind = where( probability_bg gt 0, complement=fix )
      if fix[0] ne -1 then probability_bg[fix] = min(probability_bg[ind])
      probability_bg = probability_bg / total(probability_bg*dp)

      fg_smooth = smooth(probability_fg>0,5)
      ind = where( probability_fg lt max(probability_fg)*0.01 )
      if ind[0] ne -1 then probability_fg[ind] = fg_smooth[ind]
      ind = where( probability_fg gt 0, complement=fix )
      if fix[0] ne -1 then probability_fg[fix] = min(probability_fg[ind])
      probability_fg = probability_fg / total(probability_fg*dp)

      p_bincen[i] = ptr_new( prop_bincen )
      p_fg[i] = ptr_new( probability_fg>0 )
      p_fg_err[i] = ptr_new( probability_fg_err )
      p_bg[i] = ptr_new( probability_bg>0 )
      p_bg_err[i] = ptr_new( probability_bg_err )

      ; plot the prior probability densities
      if keyword_set(showplot) then begin
         IF ~postscript THEN window,6+i

        plot, prop_bincen, probability_bg, linestyle=2, charsize=1.5, $
              ytitle='Probability Density', xtitle=label, $
              yrange=[0,max([probability_bg,probability_fg])]
        errplot, prop_bincen, probability_bg-probability_bg_err, $
                 probability_bg+probability_bg_err

        oplot, prop_bincen, probability_fg
        errplot, prop_bincen, probability_fg-probability_fg_err, $
                 probability_fg+probability_fg_err

        plots, 0.8+[0,0.03], 0.9*[1.,1.], /normal
        xyouts,0.84, 0.9, 'Matches', /normal

        plots, 0.8+[0,0.03], 0.84*[1.,1.],linestyle=2, /normal
        xyouts,0.84, 0.84, 'Background', /normal
      endif

    endfor

    ; Now we have calculated all of the individual priors. Assume that
    ; they are uncorrelated and generate the joint foreground and
    ; background priors in the format expected by cat_lr:
    ;
    ;   q is the distribution of matches, and it is normalized to integrate
    ;   to avmatch
    ;
    ;   bdist is the background number density of sources per arcsec^2

    q = dblarr(propdims)
    bdist = dblarr(propdims)

    ntot = n_elements(q)     ; total number of elements in q and bdist
    coord = intarr(nprop)    ; coordinates along all dimensions
    sublen = lonarr(nprop)   ; products of dims less than i

    sublen[0] = 1
    for i=1, nprop-1 do sublen[i] = sublen[i-1]*propdims[i-1]

    for i=0l, ntot-1 do begin
      ; turn this absolute index into coordinates along each of the axes
      remainder = i
      for j=0, nprop-1 do begin
        coord[nprop-j-1] = remainder / sublen[nprop-j-1]
        remainder = remainder mod sublen[nprop-j-1]
      endfor

      ; evaluate q and bdist at these coordinates
      q[i] = 1.
      for j=0, nprop-1 do q[i] = q[i] * (*p_fg[j])[coord[j]]

      bdist[i] = 1.
      for j=0, nprop-1 do bdist[i] = bdist[i] * (*p_bg[j])[coord[j]]
    endfor

    ; Finally, ensure that neither bdist nor q are zero anywhere, and
    ; normalize q by avmatch, and bdist by rho_bg
    ind = where(q)
    minq = min(q[ind])
    ind = where( q eq 0)
    q = q * avmatch

    ind = where(bdist)
    minbdist = min(bdist[ind])
    ind = where( bdist eq 0)
    bdist = bdist * rho_bg

    ; free pointers if they are not being returned to the caller
    if arg_present( po_bincen ) then po_bincen = p_bincen $
    else ptr_free, p_bincen

    if arg_present( po_fg ) then po_fg = p_fg $
    else ptr_free, p_fg

    if arg_present( po_fg_err ) then po_fg_err = p_fg_err $
    else ptr_free, p_fg_err

    if arg_present( po_bg ) then po_bg = p_bg $
    else ptr_free, p_bg

    if arg_present( po_bg_err) then po_bg_err = p_bg_err $
    else ptr_free, p_bg_err
  endif

end
