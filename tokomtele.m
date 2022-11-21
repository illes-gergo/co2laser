clc;clear;close all;
cd ./21-Nov-2022_09.33.28/;

M = readmatrix("PumpSpec.txt");
fasz = diff(M,1,1);
plot(fasz(:,1));