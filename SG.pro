PRO sgtxt
  
  ;原始数据文件路径，一列反射率数据，保存为txt文件
  file = 'C:\Users\xxxx\Desktop\22\55555.txt'
  
  ;获取数据行数line
  lines = 0L
  OPENR, lun, file, /get_lun
  WHILE ~EOF(lun) DO BEGIN
    str = ''
    READF, lun, str
    lines++
  ENDWHILE
  FREE_LUN, lun
  
  ;输出行数
  PRINT, '原始数据lines:   ', lines
  
  ;获取原始数组Data
  Data = FLTARR(lines)
  OPENR, lun, file, /get_lun  
  READF, lun, Data
  FREE_LUN, lun
  
  ;输出原始数据
  PRINT, Data
 
  ;SG滤波器参数设置 Savitzky-Golay with XX, Nth degree polynomial: 
  Nleft=5   ; lt nb/2 ?
  Nright=5  ; lt nb/2 ?
  Order=0
  Degree=3
  savgolFilter = SAVGOL(Nleft, Nright, Order, Degree)
  
  ;对原始数组进行滤波卷积
  SGdata = CONVOL(data, savgolFilter, /EDGE_TRUNCATE)
  
  ;输出结果
  help, SGdata
  print,'经过SG滤波后的结果:  ', SGdata
  outname='C:\Users\xxxx\Desktop\22\'+FILE_BASENAME(file,'.txt')+'_sg.txt'
  write_csv,outname,SGdata
  print,outname
  
END
