function db = CCDTdatabase
% function db = CCDTdatabase
%   Database of subjects performing the CCDT in the EMU.
%
%   DR 04/2019

% subject, session
i = 1;
db{i,1} = 'HUP069'; db{i,2} = '04Oct17'; i=i+1;
db{i,1} = 'HUP133'; db{i,2} = '01Feb17'; i=i+1; % 500, 1000, 1500ms DTs
db{i,1} = 'HUP133'; db{i,2} = '02Feb17'; i=i+1; % 500, 1000, 1500ms DTs
db{i,1} = 'HUP136'; db{i,2} = '07Mar17'; i=i+1; % 500, 1000, 1500ms DTs
db{i,1} = 'HUP139'; db{i,2} = '29Apr17'; i=i+1; % 20 null channels
db{i,1} = 'HUP140'; db{i,2} = '18May17'; i=i+1;
db{i,1} = 'HUP142'; db{i,2} = '30Jun17'; i=i+1; % grid (ch2-49)
db{i,1} = 'HUP143'; db{i,2} = '19Jul17'; i=i+1;
db{i,1} = 'HUP145'; db{i,2} = '15Aug17'; i=i+1; 
db{i,1} = 'HUP146'; db{i,2} = '23Aug17'; i=i+1;
db{i,1} = 'HUP150'; db{i,2} = '02Oct17'; i=i+1;
db{i,1} = 'HUP152'; db{i,2} = '10Nov17'; i=i+1;
db{i,1} = 'HUP153'; db{i,2} = '14Nov17'; i=i+1;
db{i,1} = 'HUP154'; db{i,2} = '30Nov17'; i=i+1;
db{i,1} = 'HUP157'; db{i,2} = '20Dec17'; i=i+1;
db{i,1} = 'HUP160'; db{i,2} = '26Jan18'; i=i+1; % eeg/eog channels; 4 null channels
db{i,1} = 'HUP165'; db{i,2} = '20Mar18'; i=i+1;
db{i,1} = 'HUP168'; db{i,2} = '23Apr18'; i=i+1;
db{i,1} = 'HUP171'; db{i,2} = '18Jun18'; i=i+1;
db{i,1} = 'HUP178'; db{i,2} = '17Sep18'; i=i+1; % two sessions same day; events for both
db{i,1} = 'HUP179'; db{i,2} = '29Oct18'; i=i+1; % eeg/eog
db{i,1} = 'HUP179'; db{i,2} = '01Nov18'; i=i+1; % eeg/eog
db{i,1} = 'HUP181'; db{i,2} = '03Dec18'; i=i+1; % two sessions same day; events for both; eeg/eog channels (also RE11 not processed since not consecutive with other RE channels)
db{i,1} = 'HUP182'; db{i,2} = '17Dec18'; i=i+1; % two sessions same day; events for both; eeg/eog channels
db{i,1} = 'HUP187'; db{i,2} = '13Mar19'; i=i+1; % two sessions same day; events for both; eeg/eog channels
db{i,1} = 'HUP191'; db{i,2} = '17Apr19'; i=i+1; % eeg/eog channels
db{i,1} = 'HUP191'; db{i,2} = '18Apr19'; i=i+1; % eeg/eog channels
