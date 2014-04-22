clc; clear;
format long e;
H_1 = [2 2 3;
       2 4 5;
       0 1 2];

H_2 = [1 2 3;
       1 0 1;
       0 -2 -2];

H=H_1;
for i=1:50
  [Q,R] = house(H);
  H = R*Q;
end
H1 = H

H=H_2;
for i=1:50
  [Q,R] = house(H);
  H = R*Q;
end
H2 = H
