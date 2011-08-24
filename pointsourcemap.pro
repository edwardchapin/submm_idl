;------------------------------------------------------------------------------
; pointsourcemap.pro: produce a map of point sources
;
; This routine creates a map of point sources given a PSF, locations,
; and scale factors (brightnesses). Sub-pixel positioning is used.
;
; Required Inputs:
;
;  nx       = x-size of the map
;  ny       = y-size "   "
;  psf      = 2d map of PSF (arbitrary dimensions, but centered at 0,0)
;  xpos     = x pixel locations of sources (floating point)
;  ypos     = y  "     "
;  scale    = PSF scale factors for each source
;
; Optional Outputs:
;
;  deltamap = 2d map of delta functions prior to smoothing
;
; Return Value:
;
;  Map of sources smoothed by the supplied PSF
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
;   23AUG2011: Initial committed version (EC)
;
;------------------------------------------------------------------------------

function pointsourcemap, nx, ny, psf, xpos, ypos, scale, deltamap=deltamap

  ; which points are good
  good = where( (xpos ge -0.5) and (xpos lt (nx-0.5)) and $
                (ypos ge -0.5) and (ypos lt (ny-0.5)) )

  if good[0] ne -1 then begin
    ; Use cloud-in-cell sampling to place delta functions down -- this
    ; scheme will spread flux among adjacent pixels if xpos/ypos are
    ; not integer valued. Note that the 0.5 offsets are needed since CIC
    ; calls 0 the lower edge of the first pixel, whereas we interpret
    ; 0 as the centre of the first pixel.

    deltamap = cic( scale, xpos+0.5, nx, ypos+0.5, ny, /isolated )

    ; place the PSF in a map of the same dimensions and then convolve
    xpsf = n_elements(psf[*,0])
    ypsf = n_elements(psf[0,*])

    width = min( [xpsf, nx] )
    height = min( [ypsf, ny] )

    bigpsf = dblarr(nx, ny)
    bigpsf[ nx/2-width/2:nx/2-width/2+width-1, $
            ny/2-height/2:ny/2-height/2+height-1] = $
      (shift( psf, xpsf/2, ypsf/2 ))[ xpsf/2-width/2:xpsf/2-width/2+width-1, $
                                      ypsf/2-height/2:ypsf/2-height/2+height-1]
    bigpsf = shift( bigpsf, -nx/2, -ny/2 )

    map = fft( fft(deltamap,1) * fft(bigpsf,1), -1 )
    map = double( map )
  endif else begin
    print, "Warning: no sources land on the map"
    map = dblarr( nx, ny )
  endelse

  return, map
end
