function [year, incidence, immigration, T]=loadTBCanada()

%year = [2008:2020];
%incidence = [14.6 14.4 14.1 14.7 14.6 15 14.3 15 15.5 15 14.8 15.9 14.3]; %Actual TB rate from 2008 - 2020

year = [2010:2020];
incidence = [14.1 14.7 14.6 15 14.3 15 15.5 15 14.8 15.9 14.3]; %Actual TB rate from 2010 - 2020
immigration = [259110 260036 263101 267924 240763 323192 272707 303325 313601 284157 226309];

T = [1054 1108 1112 1153 1110 1175 1231 1319 1315 1427 1303]; %TB surveillance 
% https://open.canada.ca/data/en/dataset/4dbb9bff-022d-4aab-a11d-0a2e1b0afaad
end 