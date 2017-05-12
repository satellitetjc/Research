PRO reverse_xiaotu
  COMPILE_OPT idl2
  ENVI,/restore_base_save_files
  envi_batch_init,LOG_FILE='batch.log'
  image_dir='E:\tan jiancan\tm_data\SAMPLE_20160415_LC8\xiaotu\field\';ordinary tif directory 
  image_files=file_search(image_dir,'*.TIF',count=numfiles);files name filter condition 
  
  for i=0,numfiles-1 do begin
    image_file=image_files[i]
    
    ;调用envi二次开发模式读取文件，书P431
    envi_open_file,image_file,r_fid=fid
    envi_file_query,fid,ns=ns,nl=nl,nb=nb,dims=dims,data_type=data_type,sname=sname
    mapinfo=ENVI_GET_MAP_INFO(fid=fid)
    dataarr=fltarr(ns,nl,nb)
    
    for j=0,nb-1 do begin
      data= envi_get_data(fid=fid,dims=dims,pos=j)
      fliphorzing=reverse(data,1);水平翻转
      dataarr[*,*,j]=fliphorzing
    endfor
    
    out_name='E:\tan jiancan\tm_data\SAMPLE_20160415_LC8\REVERSE\fliphorzing\field\'+string(i+31)+'.TIF'
    ;文件名加31
    ;二进制方式输出,书P435
    openw,lun,out_name,/get_lun
    writeu,lun,dataarr
    free_lun,lun
    ;写出文件的头文件信息
    envi_setup_head,fname=out_name,$
    ns=ns,nl=nl,nb=nb,data_type=data_type,$
    interleave =0,$
    offset =0,map_info=mapinfo,/write,/open 
  endfor
   
 end
