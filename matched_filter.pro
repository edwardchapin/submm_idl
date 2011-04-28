;------------------------------------------------------------------------------
; matched_filter.pro
;
; Calculate the "matched filter" for identifying point-sources in a
; map with instrumental (white) and point-source confusion noise
; components. For a complete description of this algorithm, see
; Appendix A in Chapin et al. 2010, MNRAS, submitted,
; arXiv:1003.2647.
;
; If a signal and noise map are supplied this routine can be used to
; do the noise-weighted cross-correlation with the matched filter. If
; there is a user-supplied PSF it must be be square, and centered at
; (0,0). If a scalar is supplied it indicates the FWHM (in pixels) of
; a Gaussian PSF to be generated internally. The signal and noise maps
; can have arbitrary dimensions but they should be the same.
;
; Required Inputs:
;   psf      = 2d square array containing PSF (centre pixel at 0,0)
;                                 OR
;              a scalar indicating FWHM of Gaussian PSF in
;              pixels. Will default to 128x128 pixels if no map provided.
;   white    = typical white noise RMS level in a pixel
;                                 OR
;              if not set estimate value from noise map (see below)
;   conf     = confusion RMS noise level in a pixel (same units as white)
;
; Optional Keyword Inputs:
;   map      = signal map if you want to smooth it with matched filter
;   noise    = noise map (optional) if you want to smooth with matched filter,
;              and/or which to calculate white automatically
;   pad      = if set, pad borders with this number of pixels
;
; Optional Keyword Outputs:
;   filt     = the matched filter (real space)
;   effective= the effective point source profile (psf correlated with filt)
;   sm_map   = smoothed map (if map supplied)
;   sm_noise = smoothed noise map (if map and noise supplued)
;   psf_out  = output PSF if generated internally
;   white_out= output white noise if calculated automatically
;
; Authors:
;   Edward Chapin, echapin@phas.ubc.ca (EC)
;
; History:
;   26APR2010: Initial Version (EC)
;   27AUG2010: Guess white noise level, padding (EC)
;------------------------------------------------------------------------------

pro matched_filter, psf, white, conf, $
                    filt=filt, effective=effective, $
                    map=map, noise=noise, $
                    sm_map=sm_map, sm_noise=sm_noise, $
                    psf_out=psf_out, hits=hits, pad=pad

  if n_elements(psf) eq 1 then begin
    ; generate a peak-normalized PSF if needed

    if( keyword_set(map) ) then begin
      n = max( [n_elements(map[*,0]), n_elements(map[0,*])] )
    endif else begin
      n = 128
    endelse

    d = dist(n)
    fwhm = psf
    psf = exp( -d^2d/(2d*(fwhm/2.35d)^2d) )
    psf = psf/max(psf)
  endif else begin
    ; Get dimensions from the PSF
    n = n_elements(psf[*,0])
  endelse

  ; work out the white noise level if needed
  if keyword_set( white ) eq 0 then begin
    ; Create a 1/sigma^2 map and calculate the mean for a central region
    weightmap = 1/noise^2d
    nx = n_elements(weightmap[*,0])
    ny = n_elements(weightmap[0,*])

    sz = fix(nx*0.1) ; look at the central 20% x 20% of the map

    white = 1/sqrt( mean(weightmap[nx/2-sz:nx/2+sz,ny/2-sz:ny/2+sz]) )
  endif

  ; The white noise level in Fourier space is related to the real-space
  ; noise in a pixel, and the size of the map

  nsamp = double(n)^2d
  white = sqrt(nsamp*white^2d) ; white noise level -- amplitude in Fourier Space
  p_w = white^2d               ; power spectrum level

  ; Confusion noise power spectrum is estimated to have the same shape
  ; as the PSF, but normalized to the specified confusion noise

  psf_fft = fft(psf,1)
  p_beam = abs(psf_fft)^2d
  scale_confusion = conf/stdev(psf)
  p_c = scale_confusion^2d*p_beam

  ; the matched filter is the psf divided by the noise power spectrum in
  ; Fourier space.

  p_noise = p_w + p_c                    ; total noise power spectrum
  filt_fft = psf_fft/p_noise             ; Fourier space matched filter
  filt = double( fft(filt_fft, -1) )     ; real space matched filter

  ; next we work out the normalization empirically: we want to be able
  ; to filter a map with the PSF and get the same peak amplitude. Since
  ; we are cross-correlating the map with the PSF, remember to take
  ; the transpose (in real space) of one of the arguments before taking
  ; the product in Fourier space.

  map_filt = double( fft( fft(transpose(psf),1) * filt_fft, -1 ) ) / $
             total(filt^2d)

  scale_filt = max(map_filt)

  filt_fft = filt_fft * scale_filt
  filt = filt * scale_filt
  effective = map_filt * scale_filt

  ; do map smoothing here if requested

  if keyword_set(map) then begin

    ; pad signal and noise maps
    if keyword_set(pad) then begin
      nx = n_elements(map[*,0])
      ny = n_elements(map[0,*])

      newmap = dblarr(nx+2*pad,ny+2*pad)
      newmap[pad:pad+nx-1,pad:pad+ny-1] = map
      map = newmap

      if keyword_set(noise) then begin
        newnoise = dblarr(nx+2*pad,ny+2*pad) + max(noise)
        newnoise[pad:pad+nx-1,pad:pad+ny-1] = noise
        noise = newnoise
      endif
    endif

    ; create a PSF that matches the map dimensions
    nx = n_elements(map[*,0])
    ny = n_elements(map[0,*])

    flt = dblarr(ny,nx)            ; swap because we will take the transpose
    filt = shift(filt,n/2,n/2)

    xsz = min( [nx,n] )
    ysz = min( [ny,n] )

    flt[ny/2-ysz/2:ny/2-ysz/2+ysz-1, nx/2-xsz/2:nx/2-xsz/2+xsz-1] = $
      filt[n/2-ysz/2:n/2-ysz/2+ysz-1, n/2-xsz/2:n/2-xsz/2+xsz-1]

    filt = shift(filt,-n/2,-n/2)   ; original
    flt = shift(flt,-ny/2,-nx/2)   ; matched to dimensions of map

    if keyword_set(noise) then begin
      ; noise-weighted cross-correlation
      weightmap = 1/noise^2.

      smoothweightflux = double(fft( fft(map*weightmap,1) * $
                                     fft(transpose(flt),1), -1 ))

      smoothweightmap = double(fft( fft(weightmap,1) * $
                                    fft(transpose(flt)^2d,1), -1 ))

      sm_map = smoothweightflux/smoothweightmap
      sm_noise = sqrt(1./smoothweightmap)

    endif else begin
      ; simple cross-correlation
      sm_map = double(fft( fft(map,1)*fft(transpose(flt),1), -1 )) / $
               total(flt^2d)
    endelse

    ; undo padding
    if keyword_set(pad) then begin
      nx = nx - 2*pad
      ny = ny - 2*pad

      map = map[pad:pad+nx-1,pad:pad+ny-1]
      sm_map = sm_map[pad:pad+nx-1,pad:pad+ny-1]

      if keyword_set(noise) then begin
        noise = noise[pad:pad+nx-1,pad:pad+ny-1]
        sm_noise = sm_noise[pad:pad+nx-1,pad:pad+ny-1]
      endif
    endif
  endif

  ; other return values
  psf_out = psf
  white_out = white
end
