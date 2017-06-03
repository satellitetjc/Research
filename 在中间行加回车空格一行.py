#会出现一个小问题：找不到路径，在spyder按ctrl+shift+alt+p在按run修改路径
import sys, string, os

dir = "F:\Training"
files = os.listdir(dir)
for f in files:
        if os.path.splitext(f)[1] == '.txt':
                
                Output_Workspace = "F:\Training\mmp"
                basename = os.path.splitext(f)[0]
                Output_file = Output_Workspace + os.sep + basename + ".txt"

                if os.path.exists(Output_file) == False:
                         fin=open(f)
                         a=fin.readlines()
                         fout=open(Output_file,'w')
                         #should write 'w'
                         #for i in range(0,len(a)-1,3):
                         #b='\n'.join(a[3:])
                         b=''.join(a[0:32])+'\n'+''.join(a[32:64])+'\n'+''.join(a[64:96])
                         #在中间行加回车空格一行
                         fout.write(b)
                         fin.close()
                         fout.close()
