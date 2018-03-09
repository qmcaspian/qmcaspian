xm2 = ∑i xm*xm;
xp2 = ∑i xp*xp;

ym2 = ∑i ym*ym;
yp2 = ∑i yp*yp;

zm2 = ∑i zm*zm;
zp2 = ∑i zp*zp;

xmym = ∑i xm*ym;
xmyp = ∑i xm*yp;

xmzm = ∑i xm*zm;
xmzp = ∑i xm*zp;

xpym = ∑i xp*ym;
xpyp = ∑i xp*yp;

xpzm = ∑i xp*zm;
xpzp = ∑i xp*zp;

ymzm = ∑i ym*zm;
ymzp = ∑i ym*zp;

ypzm = ∑i yp*zm;
ypzp = ∑i yp*zp;

where
   xm = si[0]-ti[0], xp = si[0]+ti[0],
   ym = si[1]-ti[1], yp = si[1]+ti[1],
   zm = si[2]-ti[2], zp = si[2]+ti[2];

M = [ [xm2+ym2+zm2,   ypzm-ymzp,   xmzp-xpzm,   xpym-xmyp],
      [  ypzm-ymzp, xm2+yp2+zp2,   xmym-xpyp,   xmzm-xpzp],
      [  xmzp-xpzm,   xmym-xpyp, xp2+ym2+zp2,   ymzm-ypzp],
      [  xpym-xmyp,   xmzm-xpzp,   ymzm-ypzp, xp2+yp2+zm2] ];

val_vec = Jacobi(M);  // Eigen-values and -vectors
val = val_vec[0], vec = val_vec[1];
small = 0;
for(i = 1; i < val.length; i ++)
  if(val[i] < val[small]) small=i;
return rotationQ( vec[small] );