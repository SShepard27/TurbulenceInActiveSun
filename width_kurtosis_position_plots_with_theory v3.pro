;
; subroutines for theoretical curves
; 
function turb_integrand, x
  common params,r_pass,alpha_pass
  return, 2*x/(1.+2.*r_pass*x^(6.-alpha_pass))
end
;
function line_profile, x
  common params,r_prof,alpha_prof
  result = exp(-qromo('turb_integrand',0,x))
  return, result
end
;
function w2_profile, x
  result2=x^2. * line_profile(x)
  return,result2
end
;
function w4_profile, x
  result4=x^4. * line_profile(x)
  return,result4
end
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;                                              ;
hard = 1 ; 0 for screen plot, 1 for hard copy  ;
;                                              ;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;                                                                ;
pos_step = 3 ; interval between slit positions used for analysis ;
;                                                                ;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; name of folder where data FITS file is located, and name of the file
;
folder_name='c:\users\agr75681\Documents\IDL\IRIS\'
file_name='data_for_gordon_SiIV1394_oct252014.fits'
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;                                                                  ;
; properties of the spectral line under observation (e.g., Si IV)  ;
;                                                                  ;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
t_formation=7.e4 ;formation temperature for Si IV (K)
ion_mass=28.*1.67e-24 ; mass of Si ion in grams
rest_wavelength=1394. ; rest wavelength in Angstroms
v_thermal = sqrt(1.38e-16*t_formation/ion_mass)
;
inst_width = 0.01214 ; instrument width in Angstroms
;
; read in data
;
wavelength=readfits(folder_name+file_name,header,exten_no=0)
intensity=readfits(folder_name+file_name,header,exten_no=1)
time=readfits(folder_name+file_name,header,exten_no=2)
;
; determine size of arrays
;
n_lambda=n_elements(wavelength)
n_time=n_elements(time)
n_pos=30 ; number of positions along slit
intensity_pixel=fltarr(n_time,n_pos) ; lightcurves for each pixel
intensity_overall=fltarr(n_time) ; for overall lightcurve
;
; define arrays for normalized moments
;
delta_lambda=fltarr(n_lambda,n_time,n_pos)
raw_moment=fltarr(n_time,n_pos,5)
norm_moment=fltarr(n_time,n_pos,5)
sdev=fltarr(n_time,n_pos)
v_excess=fltarr(n_time,n_pos)
red_kurt=fltarr(n_time,n_pos)
;
; determine lightcurves for each pixel and for the overall flare
;
for i_time=0,n_time-1 do begin
  for i_pos=0,n_pos-1 do begin
    intensity_pixel(i_time,i_pos)=total(intensity(*,i_time,i_pos))
  endfor
  intensity_overall(i_time)=total(intensity_pixel(i_time,*))
endfor
;
; arrays for Gaussian fit coefficients and fitted line profile
;
coeff_fit=fltarr(n_time,n_pos,4)
int_gauss=fltarr(n_lambda,n_time,n_pos)
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;                            ;
; main (slit position) loop  ;
;                            ;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
for i_pos=0,n_pos-1,pos_step do begin
  ;
  print,'Slit Position = ',i_pos ; just to keep track of what's going on
  ;
  ; get some colors
  ;
  white=FSC_COLOR("White", !D.Table_Size-1)
  black=FSC_COLOR("Black", !D.Table_Size-2)
  blue=FSC_COLOR("Blue", !D.Table_Size-3)
  red=FSC_COLOR("Red", !D.Table_Size-4)
  green=FSC_COLOR("Green", !D.Table_Size-5)
  charcoal=FSC_COLOR("Charcoal", !D.Table_Size-6)
  yellow=FSC_COLOR("Yellow", !D.Table_Size-7)
  cyan=FSC_COLOR("Cyan", !D.Table_Size-8)
  navy=FSC_COLOR("Navy", !D.Table_Size-9)
  gray=FSC_COLOR("Gray", !D.Table_Size-10)
  axiscolor=white ; axis color for screen prints
  backcolor=black ; background color for screen prints
  ;
  if (hard eq 0) then window,i_pos/pos_step else begin
    !p.charsize=1.5
    !p.charthick=3
    !x.thick=3
    !y.thick=3
    !z.thick=3
    !p.thick=3
    axiscolor=black ; axis color for hard copy
    backcolor=white ; background color for hard copy
    set_plot,'ps'
    device,/color,filename=folder_name+'width_and_kurtosis_position_'+string(format='(i2)',i_pos)+'.eps'
  endelse
  ;
  ; do Gaussian fits to each line profile
  ; 
  for i_time=0,n_time-1 do begin
    gfit=gaussfit(wavelength,intensity(*,i_time,i_pos),coeff,nterms=4)
    for k = 0,3 do coeff_fit(i_time,i_pos,k)=coeff(k)
     for i_lambda=0,n_lambda-1 do int_gauss(i_lambda,i_time,i_pos)=coeff(3) + coeff(0) * exp(-(wavelength(i_lambda)-coeff(1))^2./(2.*coeff(2)^2.))
  endfor
  ;
  ; range of wavelengths for moment calculations (eliminates noisy data at extreme wavelengths)
  ;
  i_lambda_low=5
  i_lambda_high=30
  ;
  ; calculate wavelength moments (0 through 4)
  ;
  for i_time=0,n_time-1 do begin
    for i_lambda=i_lambda_low,i_lambda_high do begin
      delta_lambda(i_lambda,i_time,i_pos)=wavelength(i_lambda)-coeff_fit(i_time,i_pos,1) ; deviation from centroid position
    endfor
    for n_mom=0,4 do begin
      for i_lambda=i_lambda_low,i_lambda_high do begin
        raw_moment(i_time,i_pos,n_mom)=raw_moment(i_time,i_pos,n_mom)+(delta_lambda(i_lambda,i_time,i_pos))^n_mom * intensity(i_lambda,i_time,i_pos) 
      endfor
    endfor
    ;
    ; normalized moments
    ;
    for n_mom=0,4 do begin
      norm_moment(i_time,i_pos,n_mom)=raw_moment(i_time,i_pos,n_mom)/raw_moment(i_time,i_pos,0)
    endfor
    ;
    ; excess width (in units of v_thermal)
    ;
    sdev(i_time,i_pos)=sqrt(norm_moment(i_time,i_pos,2)) ; sdev in Angstrom units
    ;
    ; excess velocity width (total - instrument - thermal, in quadrature), in units of v_thermal
    ;   
    v_excess(i_time,i_pos)=sqrt(((sdev(i_time,i_pos)/rest_wavelength)*3.e10/v_thermal)^2. - $
      ((inst_width/rest_wavelength)*3.e10/v_thermal)^2. - 1.)
    ;
    ; reduced kurtosis (dimensionless)
    ;
    red_kurt(i_time,i_pos)=norm_moment(i_time,i_pos,4)/sdev(i_time,i_pos)^4. - 3.
  endfor
  ;
  ; calculate theoretical curves for (v_excess,kurtosis) for different values of (R,alpha)
  ;
  ; define arrays to be used
  ;
  r_arr=findgen(21)/2.
  alpha_arr=findgen(13)/2.
  common params,r,alpha
  f=fltarr(n_elements(r_arr),n_elements(alpha_arr),201)
  second_moment=fltarr(n_elements(r_arr),n_elements(alpha_arr))
  fourth_moment=fltarr(n_elements(r_arr),n_elements(alpha_arr))
  sdev_theory=fltarr(n_elements(r_arr),n_elements(alpha_arr))
  sdev_excess_theory=fltarr(n_elements(r_arr),n_elements(alpha_arr))
  kurt_theory=fltarr(n_elements(r_arr),n_elements(alpha_arr))
  w=findgen(201)/10.-10.
  ;
  ; calculate moments using subroutines "line_profile" (and associated sub-subroutine "turb_integrand"),
  ;  "w2_profile", and "w4_profile"
  ;
  for ialpha=0,n_elements(alpha_arr)-1 do begin
    alpha=alpha_arr(ialpha)
    for ir=0,n_elements(r_arr)-1 do begin
      r=r_arr(ir)
      for i=0,200 do begin
        f(ir,ialpha,i)=line_profile(abs(w(i)))
      endfor
    endfor
  endfor
  ;
  for ialpha=0,n_elements(alpha_arr)-1 do begin
    alpha=alpha_arr(ialpha)
    for ir=0,n_elements(r_arr)-1 do begin
      r=r_arr(ir)
      second_moment(ir,ialpha)= qromo('w2_profile',0.,10.)/qromo('line_profile',0.,10.)
      fourth_moment(ir,ialpha) = qromo('w4_profile',0.,10.)/qromo('line_profile',0.,10.)
      sdev_theory(ir,ialpha)=sqrt(2.*second_moment(ir,ialpha))
      sdev_excess_theory(ir,ialpha)=sqrt(sdev_theory(ir,ialpha)^2.-1.)
      kurt_theory(ir,ialpha)=4.*fourth_moment(ir,ialpha)/sdev_theory(ir,ialpha)^4.-3. ; subtract 3 so that Gaussian has zero kurtosis
    endfor
  endfor
  ;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;                                            ;
  i_time_min=200                               ;
  i_time_max=300                               ;
  ;                                            ;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;
  ; plot scatterplot of (excess width, kurtosis)
  ;  
  plot,v_excess(*,i_pos),red_kurt(*,i_pos),xtitle='Excess Width (in units of v!dth!n)',$
    xrange=[0,10],yrange=[-2,6],psym=0,/nodata,background=backcolor,color=axiscolor,$
    position=[0.2,0.15,0.8,0.95]
  xyouts,-0.8,0.3,orientation=90,'Reduced Kurtosis'
  xyouts,5,5,'Slit Position = '+string(format='(i2)',i_pos),color=axiscolor
  xyouts,12.7,1.8,'t(s)',orientation=90
  xyouts,-3,2,'R',orientation=90
  xx=findgen(i_time_max-i_time_min+1)
  loadct,74
  cc=bytscl(xx,top=!D.n_colors-1)
  for i_time=i_time_min,i_time_max do plots,v_excess(i_time,i_pos),red_kurt(i_time,i_pos),$
    psym=6,color=cc(i_time-i_time_min)
  ;
  ; Add colorbars for time (right) and value of R (left)
  ;
  colorbar,ncolors=256,/vertical,/right,position=[0.85,0.15,0.88,0.95],range=[i_time_min,i_time_max],$
    divisions=5
  loadct,2
  colorbar,ncolors=256,/vertical,/invertcolors,position=[0.05,0.15,0.08,0.95],$
    divisions=5,ticknames=[10,8,6,4,2,0],range=[max(r_arr),0]
  ;
  ; overplot theoretical curves
  ;
  xxx=findgen(21)
  ccc=bytscl(xxx,top=!D.n_colors-1)
  for ir=0,20,1 do oplot,sdev_excess_theory(ir,*),kurt_theory(ir,*),color=ccc(ir)
  ;
  xyouts,8,-1,'!7a!3=0'
  oplot,[1,3],[-0.4,-0.4]
  oplot,[1,1],[-0.4,-0.2]
  oplot,[3,3],[-0.4,-0.2]
  xyouts,1.2,-0.9,'!7a!3 = 8'
  ;
  if (hard eq 1) then begin
    !p.charsize=1
    !p.charthick=1
    !p.thick=1
    !x.thick=1
    !y.thick=1
    !z.thick=1
    axiscolor=white
    backcolor=black
    device,/close_file
    set_plot,'win'
  endif ; restore plot parameters to defaults
  ;
endfor
;
end