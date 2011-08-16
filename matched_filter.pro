;------------------------------------------------------------------------------
; matched_filter.pro
;
; Calculate the "matched filter" for identifying point-sources in a
; map with instrumental (white) and point-source confusion noise
; components. For a complete description of this algorithm, see
; Appendix A in Chapin et al. (2011) MNRAS 411 505.
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
;   16AUG2011: Handle arbitrary dims, fix broken behaviour asymmetric PSFs (EC)
;------------------------------------------------------------------------------

pro matched_filter, psf, white, conf, $
                    filt=filt, effective=effective, $
                    map=map, noise=noise, $
                    sm_map=sm_map, sm_noise=sm_noise, $
                    psf_out=psf_out, hits=hits, pad=pad, $
                    white_out=white_out

  ; obtain dimensions from map if supplied, otherwise use a small
  ; hard-wired (square) map size for the filter
  if keyword_set(psf) then begin
    nx = n_elements(map[*,0])
    ny = n_elements(map[0,*])
  endif else begin
    nx = 128
    ny = nx
  endelse

  if n_elements(psf) eq 1 then begin
    ; generate a peak-normalized PSF

    maxdim = max( [nx,ny] )     ; remember the longest dimension

    d = shift(dist(maxdim),maxdim/2,maxdim/2)
    d = d[maxdim/2-nx/2:maxdim/2-nx/2+nx-1, $
          maxdim/2-ny/2:maxdim/2-ny/2+ny-1]
    d = shift(d,-nx/2,-ny/2) ; this now has possibly non-square aspect ratio

    fwhm = psf
    thepsf = exp( -d^2d/(2d*(fwhm/2.35d)^2d) )
    thepsf = thepsf/max(thepsf)
  endif else begin
    ; otherwise user supplied the psf

    thepsf = psf
  endelse

  ; work out the white noise level if needed
  if keyword_set( white ) eq 0 then begin
    ; Create a 1/sigma^2 map and calculate the mean for a central region
    weightmap = 1/noise^2d

    sz = fix(nx*0.1) ; look at the central 20% x 20% of the map

    white = 1/sqrt( mean(weightmap[nx/2-sz:nx/2+sz,ny/2-sz:ny/2+sz]) )
  endif

  ; The white noise level in Fourier space is related to the real-space
  ; noise in a pixel, and the size of the map

  nsamp = double(nx)*double(ny)
  white = sqrt(nsamp*white^2d) ; white noise level -- amplitude in Fourier Space
  p_w = white^2d               ; power spectrum level

  ; Confusion noise power spectrum is estimated to have the same shape
  ; as the PSF, but normalized to the specified confusion noise

  psf_fft = fft(thepsf,1)
  p_beam = abs(psf_fft)^2d
  scale_confusion = conf/stdev(thepsf)
  p_c = scale_confusion^2d*p_beam

  ; the matched filter is the psf divided by the noise power spectrum in
  ; Fourier space.

  p_noise = p_w + p_c                    ; total noise power spectrum
  filt_fft = psf_fft/p_noise             ; Fourier space matched filter
  filt = double( fft(filt_fft, -1) )     ; real space matched filter

  ; next we work out the normalization empirically: we want to be able
  ; to filter a map with the PSF and get the same peak
  ; amplitude. Since we are cross-correlating rather that convolving
  ; the map with the PSF, we take the complex conjugate of the PSF
  ; before taking the product in Fourier space.

  map_filt = double( fft( conj(fft(thepsf,1)) * filt_fft, -1 ) ) / $
             total(filt^2d)

  scale_filt = max(map_filt)

  filt_fft = filt_fft * scale_filt
  filt = filt * scale_filt
  effective = map_filt * scale_filt

  ; do map smoothing here if requested

  if keyword_set(map) then begin

    ; pad signal and noise maps if required
    if keyword_set(pad) then begin
      themap = dblarr(nx+2*pad,ny+2*pad)
      themap[pad:pad+nx-1,pad:pad+ny-1] = map

      if keyword_set(noise) then begin
        thenoise = dblarr(nx+2*pad,ny+2*pad) + max(noise)
        thenoise[pad:pad+nx-1,pad:pad+ny-1] = noise
      endif
    endif else begin
      themap = map
      if keyword_set(noise) then thenoise = noise
    endelse

    ; create a PSF that matches the map dimensions
    mapx = n_elements(themap[*,0])
    mapy = n_elements(themap[0,*])

    flt = dblarr(mapx,mapy)

    flt[mapx/2-nx/2:mapx/2-nx/2+nx-1, mapy/2-ny/2:mapy/2-ny/2+ny-1] = $
      shift(filt, nx/2, ny/2)

    flt = shift(flt,-mapx/2,-mapy/2)

    if keyword_set(noise) then begin
      ; noise-weighted cross-correlation
      weightmap = 1/thenoise^2.

      smoothweightflux = double(fft( fft(themap*weightmap,1) * $
                                     conj(fft(flt,1)), -1 ))

      smoothweightmap = double(fft( fft(weightmap,1) * $
                                    conj(fft(flt^2d,1)), -1 ))

      sm_map = smoothweightflux/smoothweightmap
      sm_noise = sqrt(1./smoothweightmap)

    endif else begin
      ; simple cross-correlation
      sm_map = double(fft( fft(themap,1)*conj(fft(flt,1)), -1 )) / $
               total(flt^2d)
    endelse

    ; undo padding
    if keyword_set(pad) then begin
      mapx = mapx - 2*pad
      mapy = mapy - 2*pad

      sm_map = sm_map[pad:pad+mapx-1,pad:pad+mapy-1]

      if keyword_set(noise) then begin
        sm_noise = sm_noise[pad:pad+mapx-1,pad:pad+mapy-1]
      endif
    endif
  endif

  ; other return values
  psf_out = thepsf
  white_out = white
end
