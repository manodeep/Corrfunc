function myfunct,  X,  P
  return,  p[0] + p[1]*X^p[2]
end


compile_opt idl2, strictarrsubs

codes = ['../bin/DD', '../bin/DDrppi', '../bin/wp']
codestring = ['DD', 'DDrppi', 'wp']
linestyle = [0, 0, 0]
symbols = [1, 2, 4]
colors = ['red', 'dodgerblue', 'green']
legendstring = [tex2idl("$\xi(r)$    "), tex2idl("$\xi(r_p,\pi)$"), tex2idl("$w_p(r_p)$ ")]
generate_eps = 1

partfile = '../tests/data/gals_Mr19.ff'
base_execstrings = ['PARTFILE f  PARTFILE f BINFILE NTHREADS  > xx', $
                    'PARTFILE f  PARTFILE f BINFILE PIMAX NTHREADS > xx', $
                    '420.0 PARTFILE f BINFILE PIMAX NTHREADS > xx']
binfile = './bins1 '
timings_file = 'timings_Mr19_rmax.txt'
pimax = 40.0
ntries = 5
min_rmax = 5.0
max_rmax = 50.0
NumSteps = 10

rbinsize = (max_rmax-min_rmax)/(NumSteps-1.0)
target_rmax = dindgen(NumSteps)*rbinsize + min_rmax
nrpbins = 14
nthreads = 4

ncodes = n_elements(codes)

if (findfile(timings_file))[0] eq '' then begin
   openw, lun, timings_file, /get_lun
   printf, lun, "# pimax = ", pimax, " nthreads = ", nthreads
   printf, lun, "#############################################"
   printf, lun, "#  Iteration     rmax     Time       Code    "
   printf, lun, "#############################################"
   for irmax = 0, NumSteps-1 do begin
      this_rmax = target_rmax[irmax]
      execstring = "../bin/logbins 0.1 " + strn(this_rmax) + " " + strn(nrpbins) + " > " + binfile
      spawn, execstring
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
            printf, lun, itry+1, this_rmax, t1-t0, codestring[icode], format = '(I10," ", G10.4," ",G10.4," ",A10)'
            flush, lun
         endfor
      endfor
   endfor
   free_lun, lun
endif

readcol, timings_file, iteration, cpu_time, code_string, format = 'L,X,D,A', /silent
timings = dblarr(ncodes, NumSteps)
scatter = dblarr(ncodes, NumSteps)

xdata = target_rmax
for icode = 0, ncodes-1 do begin
   for ipart = 0, NumSteps-1 do begin
      ind = lindgen(ntries) + icode*ntries + ipart*ncodes*ntries
      timings[icode, ipart] = median(cpu_time[ind])
      scatter[icode, ipart] = mad(cpu_time[ind])
   endfor
endfor

xrange = minmax(xdata)
xrange[0] *= 0.8
xrange[1] *= 2.0

yrange =  minmax([timings-scatter,  timings+scatter]) > 0.01
yrange[0] *= 0.8
yrange[1] *= 1.2


xticklen =  0.04
yticklen =  0.04
xtitle =  'R!Dmax!N'
ytitle =  'runtime [seconds]'
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

parinfo =   replicate({value:0.D,   fixed:0,   limited:[0,  0],   $
                       limits:[0.D,  0]},  3)
parinfo[*].limited[*] =  1
parinfo[0].limits =  [0.0001,  1.0]
parinfo[1].limits =  [1d-10,  100.0]
parinfo[2].limits =  [0.0001,  4.0]
parinfo[*].value =  [0.2d,  0.5,  2.0]


plot,  [0],  xrange =  xrange,  yrange =  yrange,  /nodata,  $
             xthick =  xthick,  ythick =  ythick,  xticklen =  xticklen,  yticklen =  yticklen,  $
             xtitle =  xtitle,  ytitle =  ytitle,  position =  position,  thick =  thick, /xlog, /ylog

reliable_N = 0 > 0
for icode =  0,  ncodes-1 do begin
   this_timings =  reform(timings[icode,  *])
   this_scatter =  reform(scatter[icode,  *]) > 0.0001
   oploterror,  xdata,  this_timings,  this_scatter,  color =  colors[icode],  line =  linestyle[icode],  thick =  thick,  $
                errcolor =  colors[icode],  errthick =  thick, psym = symbols[icode], symsize = symsize
   cgplots, xdata, this_timings, line = linestyle[icode], thick = thick, color = colors[icode], noclip = 0
   
   result =  mpfitfun('MYFUNCT',  xdata[reliable_N:*],  this_timings[reliable_N:*],  this_scatter[reliable_N:*],  parinfo =  parinfo,  yfit =  yfit, /quiet)
   cgplots,  xdata[reliable_N:*],  yfit,  line =  2,  thick =  thick,  color =  colors[icode], noclip = 0
   textstring =  tex2idl("$\propto$") + ' r!U' + string(result[2],  format =  '(F3.1)') + '!N'
   legendstring[icode] +=  + " ( " + textstring  + " ) "
endfor

al_legend, legendstring, psym = symbols, color = colors, /top, /left, symsize = symsize, $
           textcolor = colors

plot,  [0],  xrange =  xrange,  yrange =  yrange,  /nodata,  $
             xthick =  xthick,  ythick =  ythick,  xticklen =  xticklen,  yticklen =  yticklen,  $
             xtitle =  xtitle,  ytitle =  ytitle,  position =  position,  thick =  thick, /noerase, /xlog, /ylog

if !d.name eq 'PS' then begin
   device,  /close
   set_plot,  'x'
   @idl_reset_graphics
endif


end
