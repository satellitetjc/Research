pro TIF2CSV_BANDS
  COMPILE_OPT idl2
  ENVI,/restore_base_save_files
  envi_batch_init,LOG_FILE='batch.log'
  
  tif_dir='E:\tan jiancan\tm_data\SAMPLE_20160415_LC8\xiaotu\water\';ordinary tif directory
  tif_files=file_search(tif_dir,'*.tif',count=numfiles);files name filter condition
  
  for i=0,numfiles-1 do begin
    tif_file=tif_files[i]
  
    envi_open_file,tif_file,r_fid=fid
    ;envi_select,r_fid=fid,dims=dims,pos=pos
    envi_file_query,fid,ns=ns,nl=nl,nb=nb,dims=dims,$
    xstart=xstart,ystart=ystart,$
    data_type=data_type,sname=sname
    mapinfo=ENVI_GET_MAP_INFO(fid=fid)    
      
    data0=ENVI_GET_DATA(fid=fid,dims=dims,pos=0)
    data1=ENVI_GET_DATA(fid=fid,dims=dims,pos=1)
    data2=ENVI_GET_DATA(fid=fid,dims=dims,pos=2)
    data=[[data0],[data1],[data2]]
        
    out_name='E:\tan jiancan\tm_data\SAMPLE_20160415_LC8\IDL_XIAOTU_SAMPLE_ASCII\water\'+FILE_BASENAME(tif_files[i],'.tif')+'.txt'
   
    ;二进制方式输出,书P435
    openw,lun,out_name,/get_lun
    writeu,lun,data
    free_lun,lun
    ;写出文件的头文件信息
    envi_setup_head,fname=out_name,$
    ns=ns,nl=nl,nb=nb,data_type=data_type,$
    interleave =0,$
    offset =0,xstart=xstart,ystart=ystart,$
    map_info=mapinfo,/write,/open
     
    write_csv,out_name,data
    ;help,data
    print,out_name
  endfor
end
