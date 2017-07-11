;对光谱曲线进行sg滤波操作
PRO SGFilter, Nleft, Nright, Order, Degree 
  ; REFERENCE: HELP about SAVGOL
  ; Specify the following parameters!!! 
  Nleft=5   ; lt nb/2 ?
  Nright=5  ; lt nb/2 ?
  Order=0
  Degree=3 
  
  ; OPEN FILE & read data
  spectrum_dir='E:\....'   
  spectrum_files=file_search(image_dir,'*.csv',count=numfiles)
  ENVI_OPEN_FILE,spectrum_files, r_fid=fid 
  ENVI_FILE_QUERY, fid, dims=dims, ns=ns, nl=nl, nb=nb
  for i=0,numfiles-1 do begin
    data=read_csv(spectrum_files[i])
  endfor
  
  ; Savitzky-Golay with XX, Nth degree polynomial: 
  savgolFilter = SAVGOL(Nleft, Nright, Order, Degree)
  SGdata = CONVOL(transpose(data), savgolFilter, /EDGE_TRUNCATE)
  SGdata = transpose(SGdata)
  
  outname='E:\...'+FILE_BASENAME(spectrum_files[i],'.csv')+'_sg.csv'  

END
