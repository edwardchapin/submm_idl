;------------------------------------------------------------------------------
;+
;NAME
; matched_filter
;PURPOSE
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
;INPUTS
;   psf      = 2d square array containing PSF (centre pixel at 0,0)
;                                 OR
;              a scalar indicating FWHM of Gaussian PSF in
;              pixels. Will default to 128x128 pixels if no map provided.
;   conf     = confusion RMS noise level in a pixel (same units as map)
;
;OPTIONAL INPUTS
;   map      = signal map if you want to smooth it with matched filter
;   white_in = typical white noise RMS level in a pixel
;   noise    = noise map (optional) if you want to smooth with matched filter,
;              and/or which to calculate white automatically
;   pad      = if set, pad borders with this number of pixels
;
;OPTIONAL OUTPUTS
;   filt     = the matched filter (real space)
;   effective= the effective point source profile (psf correlated with filt)
;   sm_map   = smoothed map (if map supplied)
;   sm_noise = smoothed noise map (if map and noise supplued)
;   psf_out  = psf actually used
;   white_out= white noise actually used
;
; Authors:
;   Edward Chapin, echapin@phas.ubc.ca (EC)
;   Alex Conley, alexander.conley@colorado.edu (AC)
;   Gaelen Marsde, gmarsden@phas.ubc.ca (GM)
;
; History:
;   26APR2010: Initial Version (EC)
;   27AUG2010: Guess white noise level, padding (EC)
;   11JUL2011: have it not modify input arguments (AC)
;   16AUG2011: Handle arbitrary dims, fix broken behaviour asymmetric PSFs (EC)
;   07SEP2011: Merge SMAP patched (GM)
;-
;------------------------------------------------------------------------------

PRO matched_filter, psf, conf, WHITE_IN=white, $
                    FILT=filt, EFFECTIVE=effective, $
                    MAP=map, NOISE=noise, $
                    SM_MAP=sm_map, SM_NOISE=sm_noise, $
                    PSF_OUT=psf_out, HITS=hits, PAD=pad,$
                    WHITE_OUT=white_out
  COMPILE_OPT IDL2

  IF N_ELEMENTS(psf) EQ 0 THEN MESSAGE,"No user provided psf information"

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
    ;; generate a peak-normalized PSF if needed

    d = DIST(nx,ny)

    fwhm = psf
    psf_out = exp( -d^2d/(2d*(fwhm/2.35d)^2d) )
    psf_out = psf_out/max(psf_out)
  endif else begin
    ; otherwise user supplied the psf

    psf_out = psf
  endelse

  ;; work out the white noise level if needed
  IF N_ELEMENTS(white) EQ 0 THEN BEGIN
     ;; Create a 1/sigma^2 map and calculate the mean for a central
     ;; region
     IF N_ELEMENTS(noise) EQ 0 THEN MESSAGE,"Need noise to compute white"
     weightmap = 1/noise^2d
     nx = n_elements(weightmap[*,0])
     ny = n_elements(weightmap[0,*])

     sz = fix(nx*0.1) ;; look at the central 20% x 20% of the map

     white_out = 1/sqrt( mean(weightmap[nx/2-sz:nx/2+sz,ny/2-sz:ny/2+sz]) )
  ENDIF ELSE BEGIN
     white_out = white
     IF ~ FINITE(white_out) THEN MESSAGE,"Invalid (non-finite) white value"
     IF white_out LE 0.0 THEN MESSAGE,"Invalid (non-positive) white value"
  ENDELSE

  ;; The white noise level in Fourier space is related to the real-space
  ;; noise in a pixel, and the size of the map

  nsamp = double(nx)*double(ny)
  p_w = nsamp * white_out^2d               ; power spectrum level

  ;; Confusion noise power spectrum is estimated to have the same shape
  ;; as the PSF, but normalized to the specified confusion noise

  psf_fft = fft(psf_out,1)
  p_beam = abs(psf_fft)^2d
  scale_confusion = conf/stdev(psf_out)
  p_c = scale_confusion^2d*p_beam

  ;; the matched filter is the psf divided by the noise power spectrum in
  ;; Fourier space.

  p_noise = p_w + p_c                    ; total noise power spectrum
  filt_fft = psf_fft/p_noise             ; Fourier space matched filter
  filt = double( fft(filt_fft, -1) )     ; real space matched filter

  ;; next we work out the normalization empirically: we want to be able
  ;; to filter a map with the PSF and get the same peak
  ;; amplitude. Since we are cross-correlating rather that convolving
  ;; the map with the PSF, we take the complex conjugate of the PSF
  ;; before taking the product in Fourier space.

  map_filt = double( fft( conj(fft(psf_out,1)) * filt_fft, -1 ) ) / $
             total(filt^2d)

  scale_filt = max(map_filt)

  filt_fft = filt_fft * scale_filt
  filt = filt * scale_filt
  effective = map_filt / scale_filt

  ;; do map smoothing here if requested
  IF N_ELEMENTS(map) NE 0 THEN BEGIN

     ;; pad signal and noise maps
     IF N_ELEMENTS(pad) NE 0 THEN BEGIN

        working_map = dblarr(nx+2*pad,ny+2*pad)
        working_map[pad:pad+nx-1,pad:pad+ny-1] = map

        IF N_ELEMENTS(noise) NE 0 THEN BEGIN
           working_noise = dblarr(nx+2*pad,ny+2*pad) + max(noise)
           working_noise[pad:pad+nx-1,pad:pad+ny-1] = noise
        ENDIF
     ENDIF ELSE BEGIN
        working_map = map
        IF N_ELEMENTS(noise) NE 0 THEN working_noise = noise
     ENDELSE

     ; create a PSF that matches the map dimensions
     mapx = n_elements(working_map[*,0])
     mapy = n_elements(working_map[0,*])

     flt = dblarr(mapx,mapy)

     flt[mapx/2-nx/2:mapx/2-nx/2+nx-1, mapy/2-ny/2:mapy/2-ny/2+ny-1] = $
        shift(filt, nx/2, ny/2)

     flt = shift(flt,-mapx/2,-mapy/2)

     IF N_ELEMENTS(noise) NE 0 THEN BEGIN
        ;;Apply smoothing to noise as well
        ;; noise-weighted cross-correlation
        weightmap = 1/working_noise^2.

        smoothweightflux = double(fft( fft(working_map*weightmap,1) * $
                                       conj(fft(flt,1)), -1 ))

        smoothweightmap = double(fft( fft(weightmap,1) * $
                                      conj(fft(flt^2d,1)), -1 ))

        sm_map = smoothweightflux/smoothweightmap
        sm_noise = sqrt(1./smoothweightmap)

     endif else begin
        ; simple cross-correlation
        sm_map = double(fft( fft(working_map,1)*conj(fft(flt,1)), -1 )) / $
                 total(flt^2d)
     endelse

     ;; undo padding
     IF N_ELEMENTS(pad) NE 0 THEN BEGIN
        mapx = mapx - 2*pad
        mapy = mapy - 2*pad

        IF N_ELEMENTS(sm_map) NE 0 THEN $
           sm_map = sm_map[pad:pad+mapx-1,pad:pad+mapy-1]

        IF N_ELEMENTS(sm_noise) NE 0 THEN $
           sm_noise = sm_noise[pad:pad+mapx-1,pad:pad+mapy-1]

     ENDIF
  ENDIF

END
