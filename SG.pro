;对光谱曲线进行sg滤波操作，先转成tiff再操作
PRO test, Nleft, Nright, Order, Degree
  ; REFERENCE: HELP about SAVGOL
  ; Specify the following parameters!!!
  Nleft=5   ; lt nb/2 ?
  Nright=5  ; lt nb/2 ?
  Order=0
  Degree=3

  ; OPEN FILE & read data
  spectrum_dir='C:\Users\Jiancan\Desktop\'
  spectrum_files=file_search(spectrum_dir,'*.csv',count=numfiles)
 
  fdata = fltarr(62)
  for i=0,numfiles-1 do begin
    fdata=read_csv(spectrum_files[i])
    print,fdata
   ;Savitzky-Golay with XX, Nth degree polynomial: 
   savgolFilter = SAVGOL(Nleft, Nright, Order, Degree)
   SGdata = CONVOL(fdata, savgolFilter, /EDGE_TRUNCATE)
   outname='C:\Users\Jiancan\Desktop\'+FILE_BASENAME(spectrum_files[i],'.csv')+'_sg.csv'
   write_csv,outname,SGdata
   print,outname
  endfor
END
