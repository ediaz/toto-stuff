clear all;
clc;
format short e;
epsv = [10e-4 10e-6 10e-8];
diary problem3.txt
for i=1:3
  epsilon = epsv(i);
  A = matrixp3A(epsilon);

  [hatQ,hatR] = cgs(A);
  QTQ = hatQ'*hatQ;

  errAcgs = norm(A-hatQ*hatR);
  errIcgs = norm(eye(size(QTQ))-QTQ);


  [hatQ,hatR] = mgs(A);
  QTQ = hatQ'*hatQ;

  errAmgs = norm(A-hatQ*hatR);
  errImgs = norm(eye(size(QTQ))-QTQ);
  fprintf('epsilon = %e\n',epsilon);
  fprintf('   ||A - QR||:\n   cgs = %e   mgs = %e \n',errAcgs,errAmgs);
  fprintf('   ||I - Q^T*Q||:\n   cgs = %e   mgs = %e \n\n\n',errIcgs,errImgs);
end

diary off;
