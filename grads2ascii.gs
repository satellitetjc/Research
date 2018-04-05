*==========================================================
*-------*--------*--------*-------*--------*------
* 函数说明:
* grads2ascii 函数用于把GrADS读进来的数据输出到ASCII文件，
* 如.txt, .csv等.
*--------------------*-------------------*-----------
* 使用方法：
* ①使用前将 grads2ascii.gs  1.9版本GrADS放到安装目录的lib目录下，2.0版本放到*  OpenGrADS\Contents\Resources\Scripts目录下；
* ②使用此函数前需用GrADS打开需要转换的数据文件；
* ③用法：
* grads2ascii expr filename [format numl numb [u]]
* expr 要输出变量的变量名
* filename 变量输出的ASCII文件路径及文件名
* format 输出数据的格式,与C语言类似,默认是%g
* numl 每个记录(行)输出的数据个数，默认是8
* numb 每两个数据间插入的空格数,默认是1
* u 输出'Undef'代替数据中的缺省值
*--------------------*--------------------*------------
* 编程思路：
* grads2ascii 主要使用到了'set gxout print' 将产生的ASCII变量
* 输出到缓存,并用临时变量temp存储起来,接着使用write函数把临时
* 变量写到ASCII文件.其中使用了
* 'set prnopts format numl numb <u>'
* 控制变量输出的格式,用法可参照
* http://www.iges.org/grads/gadoc/gradcomdsetprnopts.html
*---------------------*--------------------*------------
*
* HUANG Yongjie(SYSU)
* Nov 25 2011 
*兰溪的程序已经写得很清楚怎么用的了，来个测试gs吧：
*'reinit '
*'sdfopen H:\test\slp.2010.nc'
*'grads2ascii slp H:\test\slp.txt %8.0f 144 1'
*'reinit'
*;
*-------*--------*---------*----------*----------*------
*==========================================================


function grads2ascii(args)

*=================
* 截取变量
*=================
expr = subwrd(args,1)
filename = subwrd(args,2)
format = subwrd(args,3)
numl = subwrd(args,4)
numb = subwrd(args,5)
u = subwrd(args,6)


*=================================
* 当参数为空时屏幕输出函数用法
*=================================
if(filename='')
say ''
say ' grads2ascii'
say ' 将GrADS读进来的数据输出到ASCII文件，用法:'
say ''
say ' grads2ascii expr filename [format numl numb [u]]'
say ' expr 要输出变量的变量名.'
say ' filename 变量输出的ASCII文件路径及文件名.'
say ' format 输出数据的格式,与C语言类似,默认是%g.'
say ' numl 每个记录(行)输出的数据个数，默认是8.'
say ' numb 每两个数据间插入的空格数,默认是1.'
say ' u 输出'Undef'代替数据中的缺省值.'
say ''
say ' grads2ascii version "1.0" of Nov 25 2011 01:07:00 say ' Copyright (c) HUANG Yongjie (SYSU) ^_^'
say ''
return 
endif 

*===================================
* 使用临时变量temp存储要输出的数据
*===================================
'set gxout print'

*--若数据格式参数不为空则设置输出数据格式
if(format != '')
'set prnopts 'format' 'numl' 'numb' 'u''
endif 

*---存储输出数据---
'd ' expr
if (rc !=0); return -1; endif 
temp = result 


*==========================================================
* 使用 sublin 获取每一行，并使用 write 写入ASCII文件中
*==========================================================

*--------------------------------------------------------
* 第一行为：Printing Grid -- XXX Values -- Undef = XXX
* 其中XXX表示数值，一般不把这行写入文件中
*--------------------------------------------------------
i = 1; line = sublin(temp,i)
num = subwrd(line,4)
undef = subwrd(line,9)

*--从第二行开始把数据写入ASCII文件中--
i = 2; line = sublin(temp,i)
while(line !='')
rc = write(filename,line)
if (subwrd(rc,1) !=0); return -2; endif 
i = i + 1
line = sublin(temp,i)
endwhile 

*-----关闭文件------
res = close(filename)
if(res != 0); return -3; endif 

say ' 'num' Values have been written to 'filename'!'
say ' Thanks to 兰溪! 兰溪哥哥亲亲！'

return 
