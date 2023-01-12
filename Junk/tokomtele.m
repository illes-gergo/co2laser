clc;clear;close all;

%M = readmatrix("21-Nov-2022_12.07.29/PumpSpec-4000_um.txt");
%M2 = readmatrix("21-Nov-2022_14.46.15/PumpSpec-4000_um.txt");

M = readmatrix("21-Nov-2022_12.07.29/efficSH.txt");
M2 = readmatrix("21-Nov-2022_14.46.15/efficSH.txt");


plot(M(:,1),M(:,2),M2(:,1),M2(:,2));