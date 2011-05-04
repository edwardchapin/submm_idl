;------------------------------------------------------------------------------
; cat_lr: Find catalogue matches using likelihood ratios
;
; This routine is based on Sutherland & Saunders (1992) MNRAS 259 413.
; 'bdist' is the background differential counts as a function of M
; properties. 'q' is the prior distribution for matches as a function
; of M properties.  In both cases, these arrays are estimates of the
; densities, so that to obtain absolute numbers of sources, the sum of
; these times the M-dimensional volume element must be used. Also note
; that such an integral over q should give the total expected match
; fraction (i.e. in the range 0 -- 1, although it might be higher if
; their are multiple identifications expected on average).
;
; Required Inputs:
;
;   p_ra       = array of primary source right ascensions in degrees
;   p_dec      = "          "     "       declinations     "    "
;   m_ra       = array of matching catalogue RA's in degrees
;   m_dec      =  "         "     "          Dec's "   "
;   sr         = maximum search radius (arcsec)
;   maxmatch   = maximum number of matches per source
;   sigma      = sigma for positional uncertainties (arcsec)
;   properties = flux density, colour etc. (N cat entries x M properties)
;   q          = prior match dist of M properties. (M axes of lengths dims)
;   bdist      = prior background dist of M properties (M axes of lengths dims)
;   dims       = length of each  of the M dimensions
;   bins       = M x 2: bin_min, dbin for each dimension
;
; Optional Keyword Outputs:
;
;   matches    = N x maxmatch x 4 array giving index into matching
;                cat/distance/LR/reliability. The indices are -1 in each
;                index slot where no match found. Error message generated if
;                maxmatch too small.
;
; Authors:
;   Edward Chapin, echapin@phas.ubc.ca (EC)
;
; History:
;   29MAY2009: Initial version (EC)
;   02MAY2011: added a comment (MZ)
;   04MAY2011: Switch to using degrees for RA (EC)
;
;------------------------------------------------------------------------------

pro cat_lr, p_ra, p_dec, m_ra, m_dec, sr, maxmatch, sigma, $
            properties, q, bdist, dims, bins, matches=matches

  ; dimensions of arrays
  n = n_elements(p_ra)          ; number of objects to find counterparts for
  ncat = n_elements(m_ra)       ; elements in matching catalogue
  ndims = n_elements(dims)      ; number of priors

  ; convert right_ascensions to hours
  p_ra_h = p_ra/15d
  m_ra_h = m_ra/15d

  ; integrate q to get the total match fraction
  vol = 1.  ; total volume of an element in q
  for i=0, ndims-1 do vol = vol*bins[i,1]
  matchfrac = total(q)*vol

  ; loop over primary sources and search for counterparts
  matches = dblarr( n, maxmatch, 4 )-1 ; index, distance, LR, R

  for i=0l, n-1 do begin
    gcirc, 1, p_ra_h[i], p_dec[i], m_ra_h, m_dec, d

    ; identify sources within the search radius
    ind = where( d le sr )
    if ind[0] ne -1 then begin
      nmatch = n_elements(ind)

      if  nmatch gt maxmatch then begin
        print, "cat_lr: Error, found", nmatch," matches, > ", maxmatch
        stop
      endif

      ; Calculate likelihood ratios for each object
      for j=0, nmatch-1 do begin

        vals = properties[ind[j],*]     ; properties of this catalogue entry

        ; figure out which bin of q and n the current source lands in
        which = intarr(ndims) ; multi-dimensional index
        for k=0, ndims-1 do begin
          which[k] = fix( (vals[k] - bins[k,0])/bins[k,1] )
          if which[k] lt 0 then which[k] = 0
          if which[k] ge dims[k] then which[k] = dims[k]-1
        endfor

        whichi = 0l ; calculate the equivalent 1-d index into the arrays
        dimlen = 1l
        for k=0, ndims-1 do begin
          whichi = whichi + which[k] * dimlen
          dimlen = dimlen * dims[k]
        endfor

        lr = q[whichi]*exp( -d[ind[j]]^2d/(2d*sigma^2d) ) / $
             (2d*!DPI*sigma^2d*bdist[whichi])

        matches[i,j,0] = ind[j]     ; index of source
        matches[i,j,1] = d[ind[j]]  ; distance to source
        matches[i,j,2] = lr         ; likelihood ratio
      endfor

      ; now that we have likelihood ratios, calculate reliabilities
      tot_lr = total( matches[i,0:nmatch-1,2] )
      for j=0, nmatch-1 do begin
        matches[i, j, 3] = matches[i,j,2] / (tot_lr + (1-matchfrac))
      endfor
    endif
  endfor


end
