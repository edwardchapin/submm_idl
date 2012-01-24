;------------------------------------------------------------------------------
; findpeaks.pro: identifies peaks in a map
;
; This routine identifies a list of local maxima in a 2-d map.  You can
; optionally choose a threshold, and a minimum distance to reject many
; sources that are close to eachother. If a PSF is supplied it is also
; possible to calculation locations with sub-pixel precision.
;
; Required Inputs:
;
;  map     = 2d array
;
; Optional Keyword Inputs:
;
;  thresh  = threshold in same units as map
;  mindist = minimum distance between sources
;  header  = FITS header for astrometry
;  cntrd   = set to 0,0-centered PSF to calculate x_cen, y_cen, ra_cen, dec_cen
;
; Optional Keyword Outputs:
;
;   x      = x pixel coordinates of peaks
;   y      = y "
;   z      = map values at peaks
;   ra     = right ascension in degrees if header specified
;   dec    = declination      "    "
;   cenx   = calculate centroid x-pixel value (sub-pixel resolution)
;   ceny   = "                  y-pixel
;   cenra  = "                  RA if header specified
;   cendec = "                  Dec     "
;
; Dependencies:
;
;   IDL astro library: (http://idlastro.gsfc.nasa.gov/)
;
; Authors:
;
;   Edward Chapin, echapin@phas.ubc.ca (EC)
;
; History:
;
;   06MAY2011: Initial committed version (EC)
;   22AUG2011: Added cntrd/cenx/ceny/cenra/cendec (EC)
;   24JAN2012: Don't return list if nothing above threshold (EC)
;
;------------------------------------------------------------------------------

pro findpeaks, map, x, y, z, mindist=mindist, thresh=thresh, header=header, $
               ra=ra, dec=dec, cntrd=cntrd, cenx=cenx, ceny=ceny, $
               cenra=cenra, cendec=cendec

  xpix = n_elements(map[*,0])
  ypix = n_elements(map[0,*])

  ; x- and y- offsets to all of the pixels in the perimeter of a 3x3 box
  ; centered over the pixel we're examining
  per_xoff = [-1, 0, 1, 1, 1, 0,-1,-1]
  per_yoff = [ 1, 1, 1, 0,-1,-1,-1, 0]

  nsource = 0l

  ; loop over all pixels except for the edges

  print, "Finding local maxima..."

  srcflag = bytarr(xpix,ypix)

  for i=1l, xpix-2 do begin
    for j=1l, ypix-2 do begin

      ; identify indices for the perimeter in the map centered over this
      ; pixel
      per_xind = i+per_xoff
      per_yind = j+per_yoff
      perimeter = map[per_xind, per_yind]

      ; the pixel is a source if it is the local maximum
      if( map[i,j] gt max(perimeter)) then begin
        srcflag[i,j] = 1b
      endif
    endfor
  endfor

  srcind = where(srcflag eq 1b, nsource)
  if srcind[0] ne -1 then begin
    x = srcind mod xpix
    y = srcind / xpix
  endif

  z = map[x,y]
  havesrc = 1

  ; threshold the sources
  if keyword_set(thresh) then begin
    ind = where(z ge thresh)
    if ind[0] ne -1 then begin
      x = x[ind]
      y = y[ind]
      z = z[ind]
    endif else begin
      delvarx, x, y, z
      havesrc = 0
    endelse
  endif

  ; multiple source rejection

  if keyword_set(mindist) and havesrc then begin
    print, "Rejecting multiple nearby sources..."

    goodflag = intarr(n_elements(x))+1

    ; loop over all sources that we found
    for i=0, n_elements(x)-1 do if(goodflag[i]) then begin

      ; calculate dist^2 to all other sources from this source
      d = (x[i]-x)^2d + (y[i]-y)^2d

      ; flag all sources that are closer than the minimum distance, that
      ; haven't been rejected already
      ind = where( (d le mindist^2d) and (goodflag) )

      if ind[0] ne -1 then begin
        ; if we found some other sources nearby, reject all but the brightest
        badind = where( z[ind] lt max(z[ind]) )

        ;if n_elements(badind) gt 1 then badind = badin[1:n_elements(badins)-1]
        if badind[0] ne -1 then goodflag[ind[badind]] = 0
      endif
    endif

    ; update the list to contain only the good ones
    ind = where(goodflag)
    if ind[0] ne -1 then begin
      x = x[ind]
      y = y[ind]
      z = z[ind]
    endif
  endif

  ; calculate RA and Dec for each source

  if keyword_set(header) and havesrc then begin
    cat_pix, ra, dec, x, y, header, /pix2cat
  endif

  ; refine positions by fitting a Gaussian, fit to the supplied
  ; PSF, to the vicinity of each identified peak

  if keyword_set(cntrd) and havesrc then begin

    print, "Calculating sub-pixel peak locations..."

    ; shift PSF to the centre of the map and fit a Gaussian to it. Just
    ; use the central 5x5 pixels to get a rough idea

    psf = cntrd

    nx_psf = n_elements(psf[*,0])
    ny_psf = n_elements(psf[0,*])

    psf = shift(psf, nx_psf/2, ny_psf/2)

    sz = (5 < nx_psf) < ny_psf

    psf_centre = psf[nx_psf/2-sz/2:nx_psf/2-sz/2+sz-1, $
                     ny_psf/2-sz/2:ny_psf/2-sz/2+sz-1]

    fit = gauss2dfit( psf_centre, coeff )

    a = coeff[2]
    b = coeff[3]
    fwhm = mean( (2.35/sqrt(2))*[a,b] )

    ; then run over all of the peaks and fit central 3x3
    n = n_elements(x)
    cenx = x*1d
    ceny = y*1d

    for i=0l, n-1 do begin
      ; place into image that is 2*fwhm x 2*fwhm
      sz_img = 2d*round(fwhm)

      ; only look at sources far enough from the edges
      if( (x[i] ge sz_img/2) and $
          (x[i] lt xpix-sz_img/2) and $
          (y[i] ge sz_img/2) and $
          (y[i] lt ypix-sz_img/2) ) then begin

        img = map[x[i]-sz_img/2:x[i]-sz_img/2+sz_img-1, $
                  y[i]-sz_img/2:y[i]-sz_img/2+sz_img-1 ]

        cntrd, img, sz_img/2, sz_img/2, xc, yc, fwhm, /silent

        if( (xc ne -1) and (yc ne -1) and $
            (xc ge (sz_img/2-1.0)) and (xc lt (sz_img/2+1.0)) and $
            (yc ge (sz_img/2-1.0)) and (yc lt (sz_img/2+1.0)) ) then begin

          cenx[i] = xc+x[i]-sz_img/2
          ceny[i] = yc+y[i]-sz_img/2
        endif
      endif
    endfor

    if keyword_set(header) then begin
      cat_pix, cenra, cendec, cenx, ceny, header, /pix2cat
    endif

  endif

end
