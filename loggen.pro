;------------------------------------------------------------------------------
;
; Create an array of n logarithmically spaced points between two
; limits
;
; Inputs:
;
; minval  = lower bound
; maxval  = upper bound
; npoints = # of points
;
; Optional:
;
; double  = set keyword if we want a double array instead of floats
; linear  = force the array to be linearly spaced
;
; Output: The array
;
;------------------------------------------------------------------------------

function loggen, minval, maxval, npoints, double=double, linear=linear

  if keyword_set(double) then points = dindgen(npoints)/(npoints-1) $
  else points = findgen(npoints)/(npoints-1)

  if keyword_set(linear) then $
    return, (maxval - minval)*points + minval $
  else $
    return, 10^( (alog10(maxval/minval))*points + alog10(minval) )

end
