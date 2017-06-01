;tif格式转成ASCII码格式无非是想得到矩阵，转成csv文件即可，csv文件可用txt读取。原理是读写TIFF格式，文件写出时用csv格式

pro readtif
  image_dir='F:\DATA\AMSR2_monthly\2016LW_China\06v\';ordinary tif directory 
  image_files=file_search(image_dir,'*.tif',count=numfiles);files name filter condition 
  for i=0,numfiles-1 do begin
    data=read_tiff(image_files[i])
    
    calculate=data/100.0;像元值除以100，DN值转成温度K值
    ;get the file name insdead of file path,removing the '.tif' suffix
    outname='F:\DATA\AMSR2_monthly\2016LW_China_ascii\06v\'+FILE_BASENAME(image_files[i],'.tif')+'.txt'
    
    write_csv,outname,data
    print,outname
  endfor
end
