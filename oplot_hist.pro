; use oplot line segments to do histogram with different line styles

pro oplot_hist, bins, data, _EXTRA=extra
  
  n = n_elements(data)

  for i=0l, n-1 do begin
    ; flat part of bin
    oplot, bins[i:i+1], data[i]*[1.0,1.0], _EXTRA=extra

    ; verticle step
    if i ne n-1 then oplot, bins[i+1]*[1.0,1.0], [data[i],data[i+1]], $
                            _EXTRA=extra
  endfor

end
