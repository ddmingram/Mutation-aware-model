%%% Script that plots data for the 'loss of switching capacity' of toggle
%%% switches, as analysed in "doi.org/10.1101/2023.04.08.536106".
%%% Code author: Duncan Ingram, 2023.

set(0,'DefaultLineLineWidth',3);
set(0,'defaultAxesFontSize',16);

%%% Time values used by all dependent variables
t_values = [0,1,2,3,3.5,4,4.5,5,5.25,5.5,...
           5.75,6,6.25,6.5,6.75,7,7.25,7.5,7.75,8,...
           8.25,8.5,8.75,9,9.25,9.5,9.75,10,10.25,10.5,...
           10.75,11,11.25,11.5,12,12.25,12.5,13,13.1,14,...
           15,16,17,18,19,19.5,20,20.25,20.5,20.75,...
           21,21.25,21.5,21.75,22,22.25,22.35,32];

%%% Values for when there are no mutations
TS_nomut = ones(1,length(t_values))*544;

%%% Values for an s2 framework with alpha_E=1e5
TS_s2_1e5 = [737,756,779,806,821,837,856,876,887,898,...
             910,924,937,951,967,984,1002,1020,1040,1063,...
             1087,1112,1141,1172,1206,1244,1286,1333,1387,1447,...
             1517,1602,1699,1820,2161,2436,2809,4851,5850,NaN,...
             NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN...
             NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN];

%%% Values for an s3 framework with alpha_I=5e3
TS_s3_1e5_5e3 = [843,900,965,1051,1104,1168,1243,1341,1403,1472,...
                 1549,1640,1766,1910,2101,2411,2851,4146,3300,2642,...
                 2128,1729,1418,1173,982,832,716,625,554,498,...
                 455,422,395,375,348,339,332,322,321,316,...
                 318,325,337,356,388,414,451,477,509,552,...
                 608,692,816,1037,1504,3416,7709,NaN];

             
%% Plots

close
set(gcf, 'Position',  [500, 200, 450, 325])

plot(t_values, TS_nomut, 'k-')
xlim([0 32])
xlabel('Time of inhibitor/h')
ylim([0 5000])
ylabel('Inhibitor for 50% switch')


%%

close
set(gcf, 'Position',  [500, 200, 450, 325])

plot(t_values, TS_nomut, 'k--', 'LineWidth', 1.5)
hold on
plot(t_values, TS_s2_1e5, 'k-')
xlim([0 32])
xlabel('Time of inhibitor/h')
ylim([0 5000])
ylabel('Inhibitor for 50% switch')

%%

close
set(gcf, 'Position',  [500, 200, 450, 325])

plot(t_values, TS_nomut, 'k--', 'LineWidth', 1.5)
hold on
plot(t_values, TS_s3_1e5_5e3, 'k-')
xlim([0 32])
xlabel('Time of inhibitor/h')
ylim([0 5000])
ylabel('Inhibitor for 50% switch')