function myfunct,  X,  P
  return,  p[0] + p[1]*X^p[2]
end


function read_fortran,  lun,  nelements,  size,  var
  compile_opt idl2,  strictarrsubs
  skip1 =  0L
  skip2 =  0L

  if n_elements(var) ne nelements then stop

  readu,  lun,  skip1
  readu,  lun,  var
  readu,  lun,  skip2

  if skip1 ne skip2 or skip1 ne size*nelements then stop

  return,  var

end

pro write_fortran,  lun,  nelements,  size,  var
  compile_opt idl2,  strictarrsubs
  bytes = long(size*nelements)
  
  writeu,  lun,  bytes
  writeu,  lun,  var
  writeu,  lun,  bytes

end


compile_opt idl2, strictarrsubs

codes = ['../../bin/DD', '../../bin/DDrppi', '../../bin/wp', '../../bin/xi']
codestring = ['DD', 'DDrppi', 'wp', 'xi']
linestyle = [0, 0, 0, 0]
symbols = [1, 2, 4, 6]
colors = ['red', 'dodgerblue', 'green', 'cyan']
legendstring = [tex2idl("$DD(r)    $"), tex2idl("$DD(r_p,\pi)$"), tex2idl("$w_p(r_p)    $"), tex2idl("$\xi(r)        $")]
generate_eps = 1

partfile = '../tests/data/gals_Mr19.ff'
base_execstrings = ['PARTFILE f  PARTFILE f BINFILE NTHREADS  > xx', $
                    'PARTFILE f  PARTFILE f BINFILE PIMAX NTHREADS > xx', $
                    '420.0 PARTFILE f BINFILE PIMAX NTHREADS > xx', $
                    '420.0 PARTFILE f BINFILE NTHREADS > xx']

binfile = '../tests/bins'
timings_file = 'timings_Mr19_numpart.txt'
pimax = 40.0
ntries = 5
nthreads = 1
NumPartSteps = 10
MinNumPartFactor = 1d-2


ncodes = n_elements(codes)
openr, partlun, partfile, /get_lun
idat = lonarr(5)
fdat = fltarr(9)
znow = 0.0
idat = read_fortran(partlun, 5, 4, idat)
fdat = read_fortran(partlun, 9, 4, fdat)
znow = read_fortran(partlun, 1, 4, znow)

Ngal = idat[1]
x = fltarr(Ngal)
y = fltarr(Ngal)
z = fltarr(Ngal)

x = read_fortran(partlun, Ngal, 4, x)
y = read_fortran(partlun, Ngal, 4, y)
z = read_fortran(partlun, Ngal, 4, z)

close, partlun

logbinsize = (alog10(Ngal)-alog10(Ngal*MinNumPartFactor))/(NumPartSteps-1)
logTargetNumParts = dindgen(NumPartSteps)*logbinsize + alog10(Ngal*MinNumPartFactor)
TargetNumParts =  long(10.0d^logTargetNumParts) < Ngal

readcol, binfile, rmin, rmax, format = 'D,D', /silent, comment = '#'
if (findfile(timings_file))[0] eq '' then begin
   openw, lun, timings_file, /get_lun
   printf, lun, "# rmax = ", max(rmax), " pimax = ", pimax, " nthreads = ", nthreads
   printf, lun, "#################################################"
   printf, lun, "#  Iteration     NumPart     Time       Code    "
   printf, lun, "#################################################"
   for ipart = 0, NumPartSteps-1 do begin

      this_Ngal = TargetNumParts[ipart] 
      if this_Ngal lt Ngal then begin
         tmp_partfile = 'tmp_particles.ff'
         idat[1] = this_Ngal       
         random_indices = cgrandomindices(Ngal, this_Ngal)
         new_x = x[random_indices]
         new_y = y[random_indices]
         new_z = z[random_indices]
         openw, partlun, tmp_partfile
         write_fortran, partlun, 5, 4, idat
         write_fortran, partlun, 9, 4, fdat
         write_fortran, partlun, 1, 4, znow
         
         write_fortran, partlun, this_Ngal, 4, new_x
         write_fortran, partlun, this_Ngal, 4, new_y
         write_fortran, partlun, this_Ngal, 4, new_z
         
         close, partlun
         print, this_Ngal, tmp_partfile, format = '("Wrote ", I0," elements to file ",A0)'
      endif else begin
         execstring = "rm -f " + tmp_partfile
         spawn, execstring
         tmp_partfile = partfile
      endelse

      for icode = 0, ncodes-1 do begin
         execstring = codes[icode] + " " + base_execstrings[icode]
         execstring = str_replace(execstring, 'BINFILE', binfile)
         execstring = str_replace(execstring, 'PIMAX', strn(pimax))
         execstring = str_replace(execstring, 'PARTFILE', tmp_partfile)
         execstring = str_replace(execstring, 'NTHREADS', strn(nthreads))
         for itry = 0, ntries-1 do begin
            t0 = systime(/seconds)
            spawn, execstring, dummy, dummy1
            t1 = systime(/seconds)
            printf, lun, itry+1, this_Ngal, t1-t0, codestring[icode], format = '(I10," ", I10," ",G12.4,A10)'
            flush, lun
         endfor
      endfor
   endfor
   free_lun, lun
endif
free_lun, partlun

readcol, timings_file, iteration, numpart, cpu_time, codestring, format = 'L,L,D,A', comment = '#'
timings = dblarr(ncodes, NumPartSteps)
scatter = dblarr(ncodes, NumPartSteps)
xdata = TargetNumParts

for icode = 0, ncodes-1 do begin
   for ipart = 0, NumPartSteps-1 do begin
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
yrange[1] *= 2.0

xticklen =  0.04
yticklen =  0.04
xtitle =  '# of galaxies'
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
parinfo[*].limited[*] =  0
parinfo[0].limits =  [0.0001,  1.0]
parinfo[1].limits =  [1d-10,  100.0]
parinfo[2].limits =  [0.0001,  4.0]
parinfo[*].value =  [0.2d,  0.5,  2.0]

plot,  [0],  xrange =  xrange,  yrange =  yrange,  /nodata,  $
             xthick =  xthick,  ythick =  ythick,  xticklen =  xticklen,  yticklen =  yticklen,  $
             xtitle =  xtitle,  ytitle =  ytitle,  position =  position,  thick =  thick,  /ylog, /xlog

reliable_N = 0 > 0
for icode =  0,  ncodes-1 do begin
   this_timings =  reform(timings[icode,  *])
   this_scatter =  reform(scatter[icode,  *]) > 0.0001
   oploterror,  xdata,  this_timings,  this_scatter,  color =  colors[icode],  line =  linestyle[icode],  thick =  thick,  $
                errcolor =  colors[icode],  errthick =  thick, psym = symbols[icode], symsize = symsize
   cgplots, xdata, this_timings, line = linestyle[icode], thick = thick, color = colors[icode], noclip = 0
   
   result =  mpfitfun('MYFUNCT',  xdata[reliable_N:*],  this_timings[reliable_N:*],  this_scatter[reliable_N:*],  parinfo =  parinfo,  yfit =  yfit, /quiet)
   cgplots,  xdata[reliable_N:*],  yfit,  line =  2,  thick =  thick,  color =  colors[icode], noclip = 0
   textstring =  tex2idl("$\propto$") + ' N!U' + string(result[2],  format =  '(F3.1)') + '!N'
   legendstring[icode] +=  + " ( " + textstring  + " ) "
endfor

al_legend, legendstring, psym = symbols, color = colors, /top, /left, symsize = symsize, $
           textcolor = colors

plot,  [0],  xrange =  xrange,  yrange =  yrange,  /nodata,  $
             xthick =  xthick,  ythick =  ythick,  xticklen =  xticklen,  yticklen =  yticklen,  $
             xtitle =  xtitle,  ytitle =  ytitle,  position =  position,  thick =  thick,  /ylog, /xlog, /noerase

if !d.name eq 'PS' then begin
   device,  /close
   set_plot,  'x'
   @idl_reset_graphics
endif


end
