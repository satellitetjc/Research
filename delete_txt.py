import sys, string, os

dir = "F:\Training"
files = os.listdir(dir)
for f in files:
        if os.path.splitext(f)[1] == '.txt':
                
                Output_Workspace = "F:\Training\delete_wet_land"
                basename = os.path.splitext(f)[0]
                Output_file = Output_Workspace + os.sep + basename + ".txt"

                if os.path.exists(Output_file) == False:
                         fin=open(f)
                         a=fin.readlines()
                         fout=open(Output_file,'w')
                         #should write 'w'
                         b=''.join(a[6:])
                         fout.write(b)
                         fin.close()
                         fout.close()
                

