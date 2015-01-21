compile_opt idl2, strictarrsubs

base_execstring = './wp 420.0 ../tests/data/gals_Mr19.ff f BINFILE PIMAX NTHREADS > xx'
binfile = '../tests/bins'
timings_file = 'timings_openmp_Mr19.txt'
min_nthreads = 1
max_nthreads = 20
ntries = 5
pimax = 40.0

readcol, binfile, rmin, rmax, format = 'D,D', /silent
if (findfile(timings_file))[0] eq '' then begin
   openw, lun, timings_file, /get_lun
   printf, lun, "# rmax = ", max(rmax), " pimax = ", string(pimax, format = '(F0.2)')
   printf, lun, "#######################################"
   printf, lun, "#  Iteration     Nthreads     Time     "
   printf, lun, "#######################################"
   for ithread = max_nthreads, min_nthreads, -1 do begin
      execstring = str_replace(base_execstring, 'NTHREADS', strn(ithread))
      execstring = str_replace(execstring, 'BINFILE', binfile)
      execstring = str_replace(execstring, 'PIMAX', strn(pimax))
      for itry = 0, ntries-1 do begin
         t0 = systime(/seconds)
         spawn, execstring, dummy, dummy1
         t1 = systime(/seconds)
         printf, lun, itry+1, ithread, t1-t0, format = '(I10," ", I10," ",G12.4)'
         flush, lun
      endfor
   endfor
   free_lun, lun
endif

;;; Now read in the data
readcol, timings_file, iteration, nthreads, time, format = 'L,D,D', comment = '#', /silent

totnthreads = max_nthreads-min_nthreads+1
mean_times  = dblarr(totnthreads)
sigma_times = dblarr(totnthreads)
fixed_x = lindgen(totnthreads) + min_nthreads
ind = where(nthreads eq 1, cnt)
if cnt eq 0 then stop
base_time = median(time[ind])
speedup = dblarr(totnthreads)
speedup_error = dblarr(totnthreads)

for ithread = min_nthreads, max_nthreads do begin
   ind = where(nthreads eq ithread, cnt)
   if cnt le 1 then stop

   mean_times[ithread-min_nthreads] = median(time[ind])
   sigma_times[ithread-min_nthreads] = mad(time[ind])

   xx = base_time/time[ind]
   speedup[ithread-min_nthreads] = median(xx)
   speedup_error[ithread-min_nthreads] = mad(xx)

endfor

position = [0.2, 0.2, 0.9, 0.9]
size = 700
window, xsize = size, ysize = size
xrange = [0, 10]
yrange = xrange
cgplot, [0], /nodata, $
        xtitle = 'Number of Threads', ytitle = 'Median Speedup', $
        xrange = xrange, yrange = yrange, $
        position = position

oploterror, fixed_x, speedup, speedup_error, psym = -16, symsize = 2, errcolor = 'dodgerblue', color = 'dodgerblue'
cgplots, fixed_x, fixed_x, line = 2, thick = 4, color = 'red', noclip = 0
end
