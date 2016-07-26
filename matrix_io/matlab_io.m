clear all
close all
clf

K = mmread('K.mtx');
C = mmread('C.mtx');
M = mmread('M.mtx');

f  = 4;
om = 2*pi*f;

mat = K + 1i*om*C - om^2*M;
spy(mat); % Have fun...
