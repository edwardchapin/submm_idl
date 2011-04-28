;------------------------------------------------------------------------------
; cat_pval: Find catalogue matches and calculate P-values
;
; This routine does the Downes et al. (1986) MNRAS 218 31 logarithmic
; correction to give false identification probabilities based on the
; surface density of the objects in the matching catalogue as a
; function of brightness, and given a search radius.
;
; Required Inputs:
;
;   ra       = array of input source right ascensions in hours
;   dec      = "          "     "    declinations     in degrees
;   catra    = array of matching catalogue RA's in hours
;   catdec   =  "         "     "          DEC's in degrees
;   catflux  = brightness of each source in matching catalogue
;   catarea  = area in deg^2 covered by the matching catalogue
;   insr     = input search radius in arcsec (scalar or vector)
;   maxmatch = maximum number of matches per source
;
; Required Outputs:
;
;   matches  = N x maxmatch x 3 array giving index into matching
;              cat/distance/P. The indices are -1 in each slot where no
;              match found. Error message generated if maxmatch too small.
;   P_c      = critical Poisson chance probability
;   N_p      = number density deg^-2 at "plate limit" (entire catalogue)
;
; Optional Keyword Inputs:
;
;   If the following are specified, variable background source rates
;   will be calculated as a function of position, overriding the
;   catarea variable above.
;
;   catmask  = 2d map with 1's indicating locations to be considered in cat
;   cathead  = FITS header corresponding to catmask
;   catsr    = radius in cat (arcsec) in which to measure background counts
;;
; Optional Keyword Outputs:
;
;   smcounts  = smooth counts image if catmask etc. were specified
;
; Authors:
;   Edward Chapin, echapin@phas.ubc.ca (EC)
;
; History:
;   01APR2009: Initial version (EC)
;   26MAY2009: Added catmask/cathead/catsr (EC)
;
;------------------------------------------------------------------------------

pro cat_pval, ra, dec, catra, catdec, catflux, catarea, insr, maxmatch, $
              matches, P_c, N_p, catmask=catmask, cathead=cathead, $
              catsr=catsr, smcounts=smcounts

  ; dimensions of arrays
  n = n_elements(ra)
  ncat = n_elements(catra)

  ; make array of search radii from scalar if necessary
  if n_elements(insr) eq 1 then sr = dblarr(n)+insr $
  else sr = insr

  ; area enclosed in search rad. deg^2
  sa  = !DPI*(double(sr)/3600d)^2d

  ; if requested, set up variables for handling a variable surface density
  if keyword_set(catmask) then begin
    print, "Analyzing masked region..."

    ; backup the catalogue before we select a subset
    old_catarea = catarea
    old_catra = catra
    old_catdec = catdec
    old_catflux = catflux

    ; Determine which entries from the catalogue land on the mask
    nx = sxpar(cathead,"NAXIS1")
    ny = sxpar(cathead,"NAXIS2")
    CDELT1 = sxpar(cathead,"CD1_1")
    CDELT2 = sxpar(cathead,"CD2_2")
    if CDELT1 eq 0 then begin
      CDELT1 = sxpar(cathead,"CDELT1")
      CDELT2 = sxpar(cathead,"CDELT2")
    endif
    CRPIX1 = sxpar(cathead,"CRPIX1")
    CRPIX2 = sxpar(cathead,"CRPIX2")
    CRVAL1 = sxpar(cathead,"CRVAL1")
    CRVAL2 = sxpar(cathead,"CRVAL2")
    pixres = abs( cdelt1 )*3600d ; arcsec/pixel

    radec2tan, CRVAL1, CRVAL2, catra*15d, catdec, ra_off, dec_off
    x = round(ra_off/CDELT1 + CRPIX1 - 1)
    y = round(dec_off/CDELT2 + CRPIX2 - 1)

    good = where( (x ge 0) and (x lt nx) and $
                  (y gt 0) and (y lt ny ) )

    onmask = where( catmask[x[good],y[good]] eq 1 )
    if onmask[0] ne -1 then begin
      catra = catra[good[onmask]]
      catdec = catdec[good[onmask]]
      catflux = catflux[good[onmask]]
      x = x[good[onmask]]
      y = y[good[onmask]]
    endif else begin
      print, "Error: no matching catalogue entries land on mask"
      stop
    endelse

    ; while we have the projected catalogue positions handy, create a
    ; map of the counts

    counts = dblarr(nx,ny)
    for i=0, n_elements(x)-1 do begin
      counts[x[i],y[i]] = counts[x[i],y[i]] + 1
    endfor

    ; Next, work out the effective area of the masked region centered
    ; over every position in the map within catsr. Accomplish this
    ; by convolving the mask with a disc of this radius to sum pixels.
    ; First we embed in a larger map with a border at least the size of
    ; the search radius to avoid wrap-around problems in the convolution.

    sz = max( [nx, ny] )            ; longest dimension of map
    border = round(catsr/pixres)+5  ; border to put around the edge
    dim = sz + border*2
    dim = dim + (dim mod 2)         ; ensure that we have an even number

    bigmask = fltarr( dim, dim )    ; larger array in which to embed mask
    bigmask[dim/2-nx/2:dim/2-nx/2+nx-1,dim/2-ny/2:dim/2-ny/2+ny-1] = catmask

    d = dist(dim)                   ; disc for doing the convolution
    disc = fltarr(dim, dim)
    disc[where( d le catsr/pixres )] = 1

    bigarea = float(fft( fft(disc,1)*fft(bigmask,1),-1))
    smarea = bigarea[dim/2-nx/2:dim/2-nx/2+nx-1,dim/2-ny/2:dim/2-ny/2+ny-1]
    smarea = smarea*(pixres/3600d)^2d ; convert to deg^2 from pixels

    ; Work out the effective area for each catalogue entry by
    ; projecting source positions on to the map

    radec2tan, CRVAL1, CRVAL2, ra*15d, dec, ra_off, dec_off
    x = round(ra_off/CDELT1 + CRPIX1 - 1)
    y = round(dec_off/CDELT2 + CRPIX2 - 1)

    ind = where(catmask[x,y] eq 0)
    if ind[0] ne -1 then begin
      print, "Warning: at least one source does not land on mask"
      stop
    endif

    catarea = smarea[x,y]

    ; work out the surface density of the catalogue vs. position
    ; by convolving with disc and dividing by effective area. This is the
    ; number density at the "plate limit"

    bigcounts = fltarr(dim,dim)
    bigcounts[dim/2-nx/2:dim/2-nx/2+nx-1,dim/2-ny/2:dim/2-ny/2+ny-1] = counts
    smcounts = float(fft( fft(disc,1)*fft(bigcounts,1),-1))
    smcounts = smcounts[dim/2-nx/2:dim/2-nx/2+nx-1,dim/2-ny/2:dim/2-ny/2+ny-1]
    smcounts = catmask*smcounts/smarea

    N_p = smcounts[x,y]         ; number density deg^-2 at "plate limit"
    P_c = 1-exp(-N_p*sa)        ; critical poisson chance prob.
  endif else begin
    N_p = ncat / catarea        ; number density deg^-2 at "plate limit"
    P_c = 1-exp(-N_p*sa)        ; critical poisson chance prob.
  endelse

  ; loop over sources and search for counterparts

  matches = dblarr( n, maxmatch, 3 )-1 ; index, distance, P

  for i=0l, n-1 do begin
    gcirc, 1, ra[i], dec[i], catra, catdec, d

    ; identify sources within the search radius
    ind = where( d le sr[i] )
    if ind[0] ne -1 then begin
      nmatch = n_elements(ind)

      if  nmatch gt maxmatch then begin
        print, "cat_pval: Error, found", nmatch," matches, > ", maxmatch
        stop
      endif

      ; Calculate P-values for each counterpart
      for j=0, nmatch-1 do begin

        if keyword_set(catmask) then begin ; Variable area/surface density
          matches[i,j,0] = good[onmask[ind[j]]] ; index to matching catalogue
          matches[i,j,1] = d[ind[j]] ; distance to counterpart in arcsec

          ; surface density of sources brighter than match
          N_star = n_elements( where( (catflux ge catflux[ind[j]]) and $
                                      (d le catsr) ) )
          N_star = N_star / catarea[i]
        endif else begin                   ; Constant surface density
          matches[i,j,0] = ind[j]          ; index to matching catalogue
          matches[i,j,1] = d[ind[j]]       ; distance to counterpart in arcsec

          ; surface density of sources brighter than match
          N_star = n_elements(where(catflux ge catflux[ind[j]]))
          N_star = N_star / catarea
        endelse

        ; Uncorrected false identification P-value
        P_star = 1 - exp(-!DPI*N_star*(d[ind[j]]/3600d)^2d)

        ; Perform the logarithmic correction as needed
        if P_star ge P_c[i] then E = P_c[i] $
        else E = P_star*(1 + alog(P_c[i]/P_star))
        P = 1d - exp(-E)
        matches[i,j,2] = P

      endfor

    endif

  endfor

  ; restore catalogue
  if keyword_set(catmask) then begin
    catarea = old_catarea
    catra = old_catra
    catdec = old_catdec
    catflux = old_catflux
  endif

end
