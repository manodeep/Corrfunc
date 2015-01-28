compile_opt idl2, strictarrsubs

codes = ['../bin/DD', '../bin/DDrppi', '../bin/wp']
codestring = ['DD', 'DDrppi', 'wp']
linestyle = [0, 0, 0]
symbols = [1, 2, 4]
colors = ['red', 'dodgerblue', 'green']
legendstring = [tex2idl("$\xi(r)$    "), tex2idl("$\xi(r_p,\pi)$"), tex2idl("$w_p(r_p)$ ")]
generate_eps = 0

partfile = '../tests/data/gals_Mr19.ff'
base_execstrings = ['PARTFILE f  PARTFILE f BINFILE NTHREADS  > xx', $
                    'PARTFILE f  PARTFILE f BINFILE PIMAX NTHREADS > xx', $
                    '420.0 PARTFILE f BINFILE PIMAX NTHREADS > xx']
binfile = '../tests/bins'
timings_file = 'timings_Mr19_openmp.txt'
pimax = 40.0
ntries = 5
min_nthreads = 1
max_nthreads = 16
NumSteps = max_nthreads-min_nthreads+1
ncodes = n_elements(codes)

readcol, binfile, rmin, rmax, format = 'D,D', /silent
if (findfile(timings_file))[0] eq '' then begin
   openw, lun, timings_file, /get_lun
   printf, lun, "# rmax = ", max(rmax), " pimax = ", pimax
   printf, lun, "#################################################"
   printf, lun, "#  Iteration     Nthreads     Time       Code    "
   printf, lun, "#################################################"
   for nthreads = max_nthreads, min_nthreads, -1 do begin
      for icode = 0, ncodes-1 do begin
         execstring = codes[icode] + " " + base_execstrings[icode]
         execstring = str_replace(execstring, 'BINFILE', binfile)
         execstring = str_replace(execstring, 'PIMAX', strn(pimax))
         execstring = str_replace(execstring, 'PARTFILE', partfile)
         execstring = str_replace(execstring, 'NTHREADS', strn(nthreads))
         for itry = 0, ntries-1 do begin
            t0 = systime(/seconds)
            spawn, execstring, dummy, dummy1
            t1 = systime(/seconds)
            printf, lun, itry+1, nthreads, t1-t0, codestring[icode], format = '(I10," ", I10," ",G12.4,A10)'
            flush, lun
         endfor
      endfor
   endfor
   free_lun, lun
endif

readcol, timings_file, iteration, numthreads, cpu_time, code_string, format = 'L,L,D,A', /silent
timings = dblarr(ncodes, NumSteps)
scatter = dblarr(ncodes, NumSteps)

speedup       = dblarr(ncodes, NumSteps)
speedup_error = dblarr(ncodes, NumSteps)

xdata = reverse(dindgen(NumSteps) + min_nthreads)
for icode = 0, ncodes-1 do begin

   base_ind = where(numthreads eq 1 and code_string eq codestring[icode], cnt)
   if cnt ne ntries then stop
   base_time = median(cpu_time[base_ind])

   for ipart = 0, NumSteps-1 do begin
      ind = lindgen(ntries) + icode*ntries + ipart*ncodes*ntries
      timings[icode, ipart] = median(cpu_time[ind])
      scatter[icode, ipart] = mad(cpu_time[ind])

      xx =  base_time[0]/cpu_time[ind]
      speedup[icode, ipart] =  median(xx)
      speedup_error[icode, ipart] =  mad(xx)

   endfor
endfor

xrange = minmax(xdata)
xrange[0] -= 1 > 0
xrange[1] += 1

yrange =  xrange

xticklen =  0.04
yticklen =  0.04
xtitle =  'Nthreads'
ytitle =  'Median Speedup'
position =  [0.2,  0.2,  0.9,  0.9]
symsize = 4  

if generate_eps eq 0 then begin
   thick =  3
   xthick =  4
   ythick =  4
   size =  800
   window,  0,  xsize =  size,  ysize =  size
endif else begin
   set_plot,  'ps'
   size =  12
   thick =  6
   xthick =  6
   ythick =  6
   device,   decomposed =   1
   !p.font =   1
   psfname =  str_replace(timings_file,  '.txt',  '.eps')
   device,   filename =   psfname,   /encap,   /color, /cmyk, bits_per_pixel =   8,   xsize =   size,   ysize =   size,   /in,   /helvetica,  /tt_font,   /bold,   $
             xoffset =   0,   yoffset =   0
   
endelse

plot,  [0],  xrange =  xrange,  yrange =  yrange,  /nodata,  $
             xthick =  xthick,  ythick =  ythick,  xticklen =  xticklen,  yticklen =  yticklen,  $
             xtitle =  xtitle,  ytitle =  ytitle,  position =  position,  thick =  thick

for icode =  0,  ncodes-1 do begin
   this_timings = reform(timings[icode,  *])
   this_scatter = reform(scatter[icode,  *]) > 0.0001
   this_speedup = reform(speedup[icode, *])
   this_speedup_error = reform(speedup_error[icode, *]) > 0.0001
   oploterror,  xdata,  this_speedup, this_speedup_error,  color =  colors[icode],  line =  linestyle[icode],  thick =  thick,  $
                errcolor =  colors[icode],  errthick =  thick, psym = symbols[icode], symsize = symsize
   cgplots, xdata, this_speedup, line = linestyle[icode], thick = thick, color = colors[icode], noclip = 0
endfor
cgplots, xdata, xdata, line = 2, thick = thick, color = 'black'
al_legend, legendstring, psym = symbols, color = colors, /top, /left, symsize = symsize, $
           textcolor = colors

plot,  [0],  xrange =  xrange,  yrange =  yrange,  /nodata,  $
             xthick =  xthick,  ythick =  ythick,  xticklen =  xticklen,  yticklen =  yticklen,  $
             xtitle =  xtitle,  ytitle =  ytitle,  position =  position,  thick =  thick, /noerase

if !d.name eq 'PS' then begin
   device,  /close
   set_plot,  'x'
   @idl_reset_graphics
endif

;;; print a latex table style data
texfname =  str_replace(timings_file,  '.txt',  '.tex')
openw, lun, texfname, /get_lun
for i = NumSteps-1, 0, -1 do begin
   printf, lun, xdata[i], format = '(I10," ",$)'
   for icode = 0, ncodes-1 do begin
      printf, lun, 100.0*(speedup[icode, i]/xdata[i]), format = '(" & ",F12.2," ",$)'
   endfor
   printf, lun, '\\'
endfor
free_lun, lun

end
