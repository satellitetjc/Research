fid0 = fopen('10BT.txt','r');
fid1 = fopen('11BT.txt','r');
fid2 = fopen('12BT.txt','r');
fid3 = fopen('13BT.txt','r');
fid4 = fopen('14BT.txt','r');
T10= fscanf(fid0,'%g ',[1 inf]);
T11= fscanf(fid1,'%g ',[1 inf]);
T12 = fscanf(fid2,'%g ',[1 inf]);
T13 = fscanf(fid3,'%g ',[1 inf]);
T14 = fscanf(fid4,'%g ',[1 inf]);
y = [T10;T11;T12; T13;T14];
fid = fopen('ASTER1B5.txt', 'wt');
fprintf(fid, '%6.3f %6.3f %6.3f %6.3f %6.3f\n', y);
fclose(fid);



fid1 = fopen('LST.txt','r');
fid2 = fopen('E11.txt','r');
fid3 = fopen('E12.txt','r');
fid4 = fopen('E13.txt','r');
fid5 = fopen('E14.txt','r');
l = fscanf(fid1,'%g ',[1 inf]);
e11 = fscanf(fid2,'%g ',[1 inf]);
e12 = fscanf(fid3,'%g ',[1 inf]);
e13 = fscanf(fid4,'%g ',[1 inf]);
e14 = fscanf(fid5,'%g ',[1 inf]);
y = [l;e11;e12;e13;e14];
fid = fopen('ASTER_product.txt', 'wt');
fprintf(fid, '%6.3f %6.3f %6.3f %6.3f %6.3f\n', y);
fclose(fid);


////////////////////////////////////////////////////
fid1 = fopen('modisretrieval.txt','r');
y= fscanf(fid1,'%g %g %g %g',[4 inf]);
LST=rand(50,50);
E29=rand(50,50);
E31=rand(50,50);
E32=rand(50,50);
[LST;E29; E31;E32]=y;
fid = fopen('LST.txt', 'wt');
fprintf(fid, '%6.3f\n', y[1]);
fclose(fid);


fid1 = fopen('result_10000_LST.txt','r');
LST = fscanf(fid1,'%g ',[100 100]);
lst=LST';
fid = fopen('LST_10000.txt', 'wt');
fprintf(fid, '%6.3f \n', lst);
fclose(fid);

////////////////////////////////////////////////////////
fid1 = fopen('T11.txt','r');
fid2 = fopen('T12.txt','r');
fid3 = fopen('T13.txt','r');
fid4 = fopen('T14.txt','r');
T11= fscanf(fid1,'%g ',[1 inf]);
T12 = fscanf(fid2,'%g ',[1 inf]);
T13 = fscanf(fid3,'%g ',[1 inf]);
T14 = fscanf(fid4,'%g ',[1 inf]);
fid5 = fopen('LST.txt','r');
fid6 = fopen('E11.txt','r');
fid7 = fopen('E12.txt','r');
fid8 = fopen('E13.txt','r');
fid9 = fopen('E14.txt','r');
lst = fscanf(fid1,'%g ',[1 inf]);
e11 = fscanf(fid2,'%g ',[1 inf]);
e12 = fscanf(fid3,'%g ',[1 inf]);
e13 = fscanf(fid4,'%g ',[1 inf]);
e14 = fscanf(fid5,'%g ',[1 inf]);
y = [T11;T12; T13;T14;lst;e11;e12;e13;e14];
fid = fopen('ASTERTD.txt', 'wt');
fprintf(fid, '%6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f\n', y);
fclose(fid);

将数组写到单独的文件


fid = fopen('lst11.txt','r');
y= fscanf(fid,'%g %g %g %g %g',[5 inf]);
fid1 = fopen('lstresult.txt', 'wt');
fid2 = fopen('e11.txt', 'wt');
fid3 = fopen('e12.txt', 'wt');
fid4 = fopen('e13.txt', 'wt');
fid5 = fopen('e14.txt', 'wt');
fprintf(fid1, '%6.3f \n', y(1,:));
fprintf(fid2, '%6.3f \n', y(2,:));
fprintf(fid3, '%6.3f \n', y(3,:));
fprintf(fid4, '%6.3f \n', y(4,:));
fprintf(fid5, '%6.3f \n', y(5,:));
fclose(fid);
fclose(fid1);
fclose(fid2);
fclose(fid3);
fclose(fid4);
fclose(fid5);


转换成矩阵
fid1 = fopen('lstresult.txt','r');
LST = fscanf(fid1,'%g ',[654 699]);
lst=LST';
dlmwrite('LST_457146.txt',lst, '\t' );
fclose(fid1);

fid1 = fopen('E11.txt','r');
E11 = fscanf(fid1,'%g ',[654 699]);
e11=E11';
dlmwrite('E11_457146.txt',e11, '\t' );
fclose(fid1);


fid1 = fopen('E12.txt','r');
E12 = fscanf(fid1,'%g ',[654 699]);
e12=E12';
dlmwrite('E12_457146.txt',e12, '\t' );
fclose(fid1);


fid1 = fopen('E13.txt','r');
E13 = fscanf(fid1,'%g ',[654 699]);
e13=E13';
dlmwrite('E13_457146.txt',e13, '\t' );
fclose(fid1);

fid1 = fopen('E14.txt','r');
E14 = fscanf(fid1,'%g ',[654 699]);
e14=E14';
dlmwrite('E14_457146.txt',e14, '\t' );
fclose(fid1); 
