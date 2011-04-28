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
;   noradial   = If set, and rprior set, skip radial distribution measurement
;   properties = N primary catalogue entries X M properties for priors
;   propbins   = M X 2: bin_min, delta_bin for binned priors arrays
;   propdims   = M array lengths for output binned priors arrays

;
; Optional Keyword Outputs:
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
;   28APR2011: Add showplot, rprior, noradial
;              Add properties, propbins, propdims (EC)
;
;------------------------------------------------------------------------------


;------------------------------------------------------------------------------
; This is a local function used for fitting rayleigh radial offset distribution
;------------------------------------------------------------------------------

function fitrayleigh, m, model_excess=model_excess
common fitrayleigh_block, all_r, all_rcen, excess, excess_err

scale = m[0]
sigma = m[1]
n = n_elements(all_rcen)

nint = 10d
int_bin = loggen(0,all_r[1], nint+1,/linear)
int_dbin = int_bin[1:nint] - int_bin[0:nint-1]
int_cen = (int_bin[1:nint] + int_bin[0:nint-1])/2d

chisq = 0d
all_y = dblarr(n)

; integrate model across each excess bin and calculate contribution to chi^2
for i=0, n-1 do begin
  r = all_r[i]+int_cen
  dr = int_dbin

  y = scale*total( (r/sigma^2d)*exp(-r^2d/(2d*sigma^2d))*dr )
  all_y[i] = y

  chisq = chisq + (excess[i]-y)^2d/(excess_err[i])^2d
endfor

model_excess = all_y

chisq = chisq/(n-2)

if (sigma le 0) or (scale le 0) then chisq=1000d

;print, chisq,':', m

return, chisq
end

;-----------------------------------------------------------------------------
; Main routine
;-----------------------------------------------------------------------------

pro match_calcprior, p_ra, p_dec, m_ra, m_dec, m_mask, maskhead, rmax, nr, $
                     showplot=showplot, $
                     rprior=rprior, noradial=noradial, $
                     properties=properties, propbins=propbins, propdims=propdims

common fitrayleigh_block, all_r, all_rcen, excess, excess_err

  ; characters in greek font set
  greek_sigma = 'r'
  greek_chi = 'v'

  ; Calculate a cropped mask that only contains the region than is nonzero.
  ; We also has a border that is rmax  in size around the edges to avoid
  ; wraparound problems when we make masks.
  ind = where(m_mask)
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
  area = n_elements( where(mask) ) * pixres^2d   ; area in arcsec^2

  ; project the catalogues on to the mask, and determine the x, y positions
  ; of sources that land within it
  cat_pix, p_ra, p_dec, p_x, p_y, header
  cat_pix, m_ra, m_dec, m_x, m_y, header

  p_onmap = where( (p_x ge 0) and (p_x lt nx) and $
                   (p_y ge 0) and (p_y lt ny) )

  m_onmap = where( (m_x ge 0) and (m_x lt nx) and $
                   (m_y ge 0) and (m_y lt ny) )

  pkeep = -1 & mkeep = -1

  if p_onmap[0] ne -1 then pkeep = where( mask[p_x[p_onmap],p_y[p_onmap]] )
  if m_onmap[0] ne -1 then mkeep = where( mask[m_x[m_onmap],m_y[m_onmap]] )

  if (pkeep[0] eq -1) or (mkeep[0] eq -1) then begin
    print, "Error: either no primary or matching catalogue sources land on mask"
    stop
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
  circle = intarr(nx,ny)
  circle[ind] = 1

  bgmask = intarr(nx,ny)

  p_delta = dblarr(nx,ny)
  for i=0, n_p-1 do p_delta[p_x[i],p_y[i]] = p_delta[p_x[i],p_y[i]] + 1
  temp = round(double(fft(fft(circle,1)*fft(p_delta,1),-1))) ; convolve

  fg = where( (temp eq 1) and (mask) )
  fgmask = intarr(nx,ny)
  fgmask[fg] = 1

  bg = where( (temp eq 0) and (mask) )
  bgmask = intarr(nx,ny)
  bgmask[bg] = 1

  if keyword_set(showplot) then begin
    set_plot,'x'
    window, 0
    tvscale, bgmask, /keep, /noint
    xyouts, 0.02, 0.95, 'Background Mask', charsize=1.5, charthick=3.0, $
            /normal, color=255
    xyouts, 0.02, 0.95, 'Background Mask', charsize=1.5, /normal


    window, 1
    tvscale, fgmask, /keep, /noint
    xyouts, 0.02, 0.95, 'Foreground Mask', charsize=1.5, charthick=3.0, $
            /normal, color=255
    xyouts, 0.02, 0.95, 'Foreground Mask', charsize=1.5, /normal
  endif

  ; calculate area and surface density + Poisson error, of sources in bg mask
  area_bg = n_elements(where(bgmask)) * pixres^2d                  ; arcsec^2
  n_bg = n_elements(where(bgmask[m_x,m_y]))
  rho_bg = n_bg / area_bg            ; #/arcsec^-2
  rho_bg_err = sqrt(n_bg) / area_bg

  ; Calculate bins of search radius
  all_r = loggen( 0, rmax, nr+1, /linear)       ; bin edges
  all_dr = (all_r[1:nr] + all_r[0:nr-1])        ; bin widths
  all_rcen = (all_r[0:nr-1] + all_r[1:nr])/2.   ; bin centres

  if (keyword_set(noradial) eq 0) and (keyword_set(rprior) eq 0) then begin

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

    if keyword_set(showplot) then window, 2

    for i=0, nr-1 do begin

      ind = where( (d ge all_r[i]) and (d lt all_r[i+1]) )
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

      n_full = double(n_elements(where(amask)))
      amask = amask * fgmask
      n_clipped = double(n_elements(where(amask)))

      n_p_eff = (n_clipped / n_full) * n_p

      all_n_full[i] = n_full       ; store for calculating cumulative
      all_n_clipped[i] = n_clipped ; equivalents later

      if keyword_set(showplot) then begin
        tvscale, amask, /keep, /noint
        xyouts, 0.02, 0.95, 'FG Masked Measurement Annuli', charsize=1.5, $
                charthick=3.0, $
                /normal, color=255
        xyouts, 0.02, 0.95, 'FG Masked Measurement Annuli', charsize=1.5, $
                /normal
      endif

      ; count matching sources in annulus mask
      ind = where( amask[m_x,m_y] )
      if ind[0] ne -1 then begin
        counts[i] = n_elements(ind)
        counts_err[i] = sqrt(counts[i])
      endif

      ; how many were expected?
      expect[i] = n_elements(where(amask))*pixres^2d * rho_bg
      expect_err[i] = (rho_bg_err/rho_bg)*expect[i]

      ; what is the average excess given effective number of sources?
      excess[i] = (counts[i]-expect[i]) / n_p_eff
      excess_err[i] = (counts_err[i]+expect_err[i])/(counts[i]-expect[i]) ; rel
      excess_err[i] = excess_err[i]*excess[i]                             ; abs

      effstr = " effective sources: " + string( n_p_eff, format='(F5.1)' )  + $
               ' / ' + strcompress(n_p,/remove_all)

      print, 'annulus:'+ $
             string(all_r[i],format='(I3)') + '--' + $
             string(all_r[i+1],format='(I3)') + effstr + ' excess=' + $
             string(excess[i],format='(F5.2)') + ' +/- ' + $
             string(excess_err[i],format='(F5.2)')

    endfor

    ; fit rayleigh distribution to excess counts
    scale_guess = total(excess)
    ind = (where(excess eq max(excess)))[0]
    sigma_guess = all_rcen[ind]
    m_init = [scale_guess,sigma_guess]
    m_scale = 0.3*m_init
    fit = amoeba( 1d-6, function_name='fitrayleigh', p0=m_init, scale=m_scale )
    a = fit[0]
    sigma = fit[1]

    ; high-res sampled version of the smooth distribution
    r_highres = loggen(0,rmax,100,/linear)
    pdf = (r_highres/sigma^2d)*exp(-r_highres^2d/(2d*sigma^2d))
    scale=a*(all_dr[0])         ; since total #'s given in each bin

    ; ...and also the low-resolution binned version from which chi^2 was
    ; calculated
    chisq = fitrayleigh(fit,model_excess=model_excess)

    if keyword_set(showplot) then begin
      window,3
      plot, all_rcen, excess, psym=5, xtitle='Search Radius (arcsec)', $
            ytitle='Excess / bin', xstyle=1, xrange=[0,rmax], $
            charsize=1.5
      errplot, all_rcen, excess-excess_err, excess+excess_err
      oplot,[0,rmax],[0,0]
      oplot, r_highres, pdf*scale
      oplot, all_rcen, model_excess, psym=10

      modelstr = '!4'+greek_sigma+'!3!Dr!n='+string(sigma,format='(F4.1)')
      modelstr = modelstr+' E='+string(a,format='(F3.1)')+ $
                 ' !4'+greek_chi+'!3!u2!N='+string(chisq,format='(F3.1)')

      xyouts, 0.6, 0.85, modelstr, charsize=1.5, charthick=!p.thick, /normal
    endif

    ; now calculate the cumulative distribution
    cumexcess = dblarr(nr)      ; cumulative dist
    cumexcess_err = dblarr(nr)  ; cumulative error
    cumr = all_r[1:nr]          ; radii corresponding to cumexcess

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
    endfor

    ; smooth fitted model
    cum_highres = pdf*0
    cumr_highres = r_highres - (r_highres[1]/2.)
    for i=0, n_elements(cum_highres)-1 do $
      cum_highres[i] = total( a*pdf[0:i] )

    cum_highres = cum_highres * (r_highres[1]-r_highres[0])

    if keyword_set(showplot) then begin
      window, 4
      plot, cumr, cumexcess, psym=5, xtitle='Search Radius (arcsec)', $
            ytitle='Cumulative Excess (<R)', charsize=1.5, $
            xrange=[0,rmax], xstyle=1

      errplot, cumr, cumexcess-cumexcess_err, cumexcess+cumexcess_err
      oplot, cumr_highres, cum_highres
    endif

    ; identify search radius in cumulative counts that maximizes SNR of
    ; the excess. This will define the search radius used to calculate
    ; other priors (assuming that they are uncorrelated with distance)

    cum_snr = cumexcess / cumexcess_err
    maxind = where( cum_snr eq max(cum_snr) )

    print, "Best SNR at search radius of " + $
           string(all_rcen[maxind],format='(F5.1)') + " arcsec"
  endif else begin
    print, "Skipping measurement of radial offset distribution..."

    maxind = where( all_rcen le rprior )
    maxind = (maxind[n_elements(maxind)-1])[0]

    print, "Will measure priors based on supplied properties using " + $
           "search radius of ", $
           string(all_rcen[maxind],format='(F5.1)') + " arcsec"
  endelse

  ; If properties supplied, work out other priors at best SNR radius
  if keyword_set(properties) and keyword_set(propbins) and $
    keyword_set(propdims) then begin

    ; How many properties were supplied?
    nprop = n_elements(propdims)

    ; Work out the region that is within the requested search radius of
    ; the primary catalogue positions, but also on the foreground mask

    ind = where( d lt all_r[maxind+1] )
    circle = intarr(nx,ny)
    circle[ind] = 1

    pmask = round(double(fft(fft(circle,1)*fft(p_delta,1),-1))) ; convolve
    ind = where(pmask and fgmask)
    pmask = pmask*0
    pmask[ind] = 1
    area_p = n_elements(where(pmask)) * pixres^2d ; arcsec^2

    if keyword_set(showplot) then begin
      window, 5
      tvscale, pmask, /keep, /noint
      xyouts, 0.02, 0.95, 'Prior Mask', charsize=1.5, charthick=3.0, $
              /normal, color=255
      xyouts, 0.02, 0.95, 'Prior Mask', charsize=1.5, /normal
    endif

    for i=0, nprop-1 do begin
      print, "--- Calculating prior for property ", strcompress(i+1), $
             " / ", strcompress(nprop), " ---"
      ; Set up arrays for foreground and background counts as a function
      ; of property bins
      prop_bgcounts = dblarr( propdims[i] )
      prop_bgcounts_err = dblarr( propdims[i] )
      prop_fgcounts = dblarr( propdims[i] )
      prop_fgcounts_err = dblarr( propdims[i] )

      prop_bins = loggen( propbins[i,0], $
                          propbins[i,0]+propbins[i,1]*(propdims[i]+1), $
                          propdims[i]+1, /linear )
      prop_bincen = (prop_bins[0:propdims[i]-1] + prop_bins[1:propdims[i]])/2.
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

      ; normalize the foreground and background distributions
      scale = total(prior_fg*dp)
      prior_fg = prior_fg / scale
      prior_fg_err = prior_fg_err / scale

      scale = total(prop_bgcounts*dp)
      prior_bg = prop_bgcounts / scale
      prior_bg_err = prop_bgcounts_err / scale

      window,6+i
      plot, prop_bincen, prior_bg, linestyle=2, charsize=1.5, $
            ytitle='Probability Density'
      errplot, prop_bincen, prior_bg-prior_bg_err, prior_bg+prior_bg_err

      oplot, prop_bincen, prior_fg
      errplot, prop_bincen, prior_fg-prior_fg_err, prior_fg+prior_fg_err
    endfor

  endif

end
