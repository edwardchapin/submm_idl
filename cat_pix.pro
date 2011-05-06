;------------------------------------------------------------------------------
; cat_pix: Convert between spherical and pixel coordinates
;
; The purpose of this routine is to convert between spherical
; coordinates and pixel coordinates for a map given its FITS header
; (or the inverse transformation).
;
; Required Inputs/Outputs:
;
;   lon     = input array of longitudes in degrees
;   lat     = input array of latitudes in degrees
;   x       = output x pixel coordinates
;   y       = output y pixel coordinates
;   header  = FITS header for transformation
;
; Optional Keyword Inputs:
;
;   pix2cat = if set convert x,y --> lon,lat instead
;
; Optional Keyword Outputs:
;
;   pixres  = pixel resolution (absolute value of CDELT1 -- assumed square)
;   cdelt1  = CDELT1 from header (or CD1_1 if CDELT1 doesn't exist)
;   cdelt2  = CDELT2 from header (or CD2_2 if CDELT2 doesn't exist)
;   crpix1  = CRPIX1 from header
;   crpix2  = CRPIX2 from header
;   crval1  = CRVAL1 from header
;   crval2  = CRVAL2 from header
;
; Dependencies:
;
;   IDL astro library: (http://idlastro.gsfc.nasa.gov/)
;
; Authors:
;   Edward Chapin, echapin@phas.ubc.ca (EC)
;
; History:
;   27APR2011: Initial committed version (EC)
;   28APR2011: Switch to IDL astro libs for coordinate conversion (EC)
;
;------------------------------------------------------------------------------

pro cat_pix, lon, $           ; Lon. in degrees
             lat, $           ; Lat. in degrees
             x, $             ; x pixel coordinate
             y, $             ; y pixel coordinate
             header, $        ; FITS header for transformation
             pix2cat=pix2cat,$; if set, convert x,y -> lon,lat, otherwise does
                              ; the inverse by default
             pixres=pixres, $ ; return projection information if requested
             cdelt1=cdelt1, $
             cdelt2=cdelt2, $
             crpix1=crpix1, $
             crpix2=crpix2, $
             crval1=crval1, $
             crval2=crval2

  ; Get projection from the header
  CDELT1 = sxpar(header,"CD1_1")
  CDELT2 = sxpar(header,"CD2_2")
  if CDELT1 eq 0 then begin
    CDELT1 = sxpar(header,"CDELT1")
    CDELT2 = sxpar(header,"CDELT2")
  endif
  CRPIX1 = sxpar(header,"CRPIX1")
  CRPIX2 = sxpar(header,"CRPIX2")
  CRVAL1 = sxpar(header,"CRVAL1")
  CRVAL2 = sxpar(header,"CRVAL2")

  if keyword_set(pix2cat) then begin
    ; convert x,y --> lon, lat
    xyad, header, x, y, lon, lat
  endif else begin
    ; convert lon, lat --> x,y
    adxy, header, lon, lat, x, y
  endelse

  pixres = abs(CDELT1)

end
