
sufilep=~/Dropbox/toto_y_yo/snap_sp.su
sufilev=~/Dropbox/toto_y_yo/snap_svz.su

sfsegyread <$sufilep tfile=xxx.rsf endian=0 su=y>wavefield.rsf && \
echo n1=251 d1=4 o1=0 n2=601 d2=4 o2=0 n3=476 d3=0.004 o3=0 >> wavefield.rsf
bin=`sfin wavefield.rsf |grep in|tail -1|sed -e "s/\"/ /g" |awk '{print $2}'`
mv $bin wave.dat&& echo in=./wave.dat >> wavefield.rsf 



sfsegyread <$sufilev tfile=xxx.rsf endian=0 su=y>wavefieldv.rsf && \
echo n1=251 d1=4 o1=0 n2=601 d2=4 o2=0 n3=476 d3=0.004 o3=0 >> wavefieldv.rsf
bin=`sfin wavefieldv.rsf |grep in|tail -1|sed -e "s/\"/ /g" |awk '{print $2}'`
mv $bin wavev.dat&& echo in=./wavev.dat >> wavefieldv.rsf 


# compare wavefields:
sfcat wavefield.rsf wavefieldv.rsf axis=5 |sfwindow n3=1 f3=100 |sfgrey gainpanel=e|xtpen
