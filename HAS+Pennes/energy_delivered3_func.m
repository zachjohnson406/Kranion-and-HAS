function success = energy_delivered3_func(pathTreat,pathOut)
% % extract power and energy over time from Insightec treatment exports

% this creates a shell script that extracts the data from the log files and 
% writes .txt files for each sonication that are opened in Matlab

% clear all;
% close all;

% Select "SonicationSummary.xml" from treatment export 
%[filename directoryname] = uigetfile('*.*','Pick the SonicationSummary.xml file from Treatment Export');
computer = 0;
success = 1;        %for now it always works
cd(pathTreat)
filename = 'SonicationSummary.xml';
directoryname = [pathTreat,'/'];

% need to check if there is more than 1 Csa_Brain file
CPC_Csa_Brain_files = dir(strcat(directoryname,'CPCFiles/Csa_Brain*'));
% find Csa_Brain file that is the largest - so far this worked
filesize = vertcat(CPC_Csa_Brain_files.bytes); % first make a vector of field in struct for max command to work
index = find([CPC_Csa_Brain_files.bytes] == max(filesize));
%index = find([CPC_Csa_Brain_files.bytes] ==
%max(CPC_Csa_Brain_files.bytes)); this only worked for two files, not
%three; don't know why
Csa_Brain_filename = CPC_Csa_Brain_files(index).name;

%addpath(directoryname,strcat(directoryname,'CPCFiles/'))
CSAfilename = fullfile(directoryname,'CPCFiles',Csa_Brain_filename); %CPC/CSA_Brain_xxx.txt log file
SSfilename = fullfile(directoryname,filename); % SonicationSummary.xml
XDiniFilename = fullfile(directoryname,'CPCFiles/Site/Ini','Xd_7207.ini'); % ini file containing efficiency information for different frequencies

% extract data from SonicationSummary.xml file
xmlDoc = xmlread(SSfilename);
rootNode = xmlDoc.getDocumentElement;
entries = rootNode.getChildNodes;
%entries.getLength % not sure why there are 5 entires, but #3 contains the data
theData = entries.item(3)
% theData = entries.item(3).getChildNodes; % just both in one line
% xmlwrite(theData) %display all the data in the xml file

% Get the first "Entry"'s children
theData = entries.item(3).getChildNodes;
% Iterate over the nodes to find the "sontime"
% once there are no more siblings, "node" will be empty
node = theData.getFirstChild;
ii=1;
while ~isempty(node)
    if strcmpi(node.getNodeName,'z:row') %z:row is the nodename for meaningful data, skip rest
        sonTimeString(ii) = node.getAttribute('sontime')
        sonicationIndex(ii) = str2num(node.getAttribute('sonicationindex'));
        %sonTime(ii) = str2num(node.getAttribute('sontime'));
        measuredEnergy(ii) = str2num(node.getAttribute('measuredenergy'));
        measuredPower(ii) = str2num(node.getAttribute('measuredpower'));
        actualFrequency(ii) = str2num(node.getAttribute('frequency')); 
        ii = ii+1;
    end
    node = node.getNextSibling;
end

% now find the sonications times and power information in the CSA log file
sonTimeCharacterArray = char(sonTimeString);
%sonTimeCharacterArrayPlus1 = num2str(str2num(sonTimeCharacterArray(:,9:12))+1);
sonTimeCharacterArrayPlus1 = sonTimeCharacterArray(:,9:12);
tt = str2num(sonTimeCharacterArray(:,9:12))+1;
for jj=1:size(tt,1)
    if(tt(jj) >=1000)
        sonTimeCharacterArrayPlus1(jj,:) = num2str(tt(jj));
    else
        sonTimeCharacterArrayPlus1(jj,:) = ['0',num2str(tt(jj))];
    end
end
sonDate = sonTimeCharacterArray(1,1:8);
sonHrMinSec = strcat(sonTimeCharacterArray(:,9:10),':',sonTimeCharacterArray(:,11:12),':',sonTimeCharacterArray(:,13:14))
sonHrMin = strcat(sonTimeCharacterArray(:,9:10),':',sonTimeCharacterArray(:,11:12));
sonHrMinPlus1 = strcat(sonTimeCharacterArrayPlus1(:,1:2),':',sonTimeCharacterArrayPlus1(:,3:4));
sonMin = strcat(sonTimeCharacterArray(:,11:12));
sonMinPlus1 = strcat(sonTimeCharacterArrayPlus1(:,3:4));
% % this is not working yet
% x1 = 'grep "Hw Server Status: <SONICATION>" ';
% x2 = ' > poweroutput.txt';
% command = strcat(x1, CSAfilename,x2)
% [status,cmdout] = unix(command)
cd(pathOut);
% create unix script to run on commandline
fileID = fopen('run_log_extract.sh','w');
% first grap every poweroutput from CSA log and dump into one file
formatSpec0 = '#!/bin/sh\n\n\n';
fprintf(fileID,formatSpec0);
%formatSpec0 = 'grep -e  "FREQ" -m 10 %s > efficiency.txt \n';
%fprintf(fileID,formatSpec0, XDiniFilename);
formatSpec0 = 'grep -e  " Power <" %s > power_output.txt \n';
fprintf(fileID,formatSpec0, CSAfilename);
% now remove all < and > in poweroutput.txt
% this works on the mac, but not on the computeboxes
if(computer ~=1)
    formatSpec0 = 'sed "/%s<//gw power_outputx.txt" power_output.txt >tmp\n';
    fprintf(fileID,formatSpec0,'\</s/\');
    formatSpec01 = 'sed "/%s>//gw power_output.txt" power_outputx.txt >tmp\n\n';
    fprintf(fileID,formatSpec01,'\>/s/\');
else
    % this works on the computeboxes
    formatSpec0 =  'sed "/ </s/ </ /gw power_outputx.txt" power_output.txt >tmp\n';
    fprintf(fileID,formatSpec0);
    formatSpec01 = 'sed "/> /s/> / /gw power_output.txt" power_outputx.txt >tmp\n\n';
    fprintf(fileID,formatSpec01);
end

for ii = 1:max(sonicationIndex) 
    formatSpec1 = 'grep -e "%s:" -e "%s:" power_output.txt > power_son%i.txt\n';
    fprintf(fileID,formatSpec1,sonHrMin(ii,:),sonHrMinPlus1(ii,:),ii);
    formatSpec2 = 'sed "s/:/ /gw power_son%ix.txt" power_son%i.txt >tmp\n';
    fprintf(fileID,formatSpec2,ii,ii);
    formatSpec3 = 'awk ''{ if($2 == "%s" || $2 == "%s") {print}}'' power_son%ix.txt > power_son%i.txt\n\n';
    fprintf(fileID,formatSpec3,sonMin(ii,:),sonMinPlus1(ii,:),ii,ii);
    ii=ii+1;
end
fprintf(fileID,'rm power_son*x.txt\n');
fprintf(fileID,'rm power_outputx.txt\n');
fclose(fileID);
%%
% Now run everything in run_log_extract.sh 
unix('chmod +x run_log_extract.sh');
unix('./run_log_extract.sh');
%%this code was added by Zach Johnson to deal with mismatched clocks
%%between .xml and .csa files. Reads information directly from .csa file
%%read into power_output.txt
filename = 'power_output.txt';  
fileID = fopen(filename, 'r');
data = fread(fileID, '*char')';  
fclose(fileID);



lines = strsplit(data, '\n');
offset = lines{9}(1:5);
hrOffset = str2double(sonHrMin(1, 1:2)) - str2double(offset(1:2));
minOffset = str2double(sonHrMin(1, 4:5)) - str2double(offset(4:5));

for i = 1:max(sonicationIndex)
    newHr = str2double(sonHrMin(i, 1:2)) - hrOffset;
    newMin = str2double(sonHrMin(i, 4:5)) - minOffset;

    % Adjust minutes and hours if minutes exceed 60 or go below 0
    if newMin < 0
        newMin = newMin + 60;
        newHr = newHr - 1;
    elseif newMin >= 60
        newMin = newMin - 60;
        newHr = newHr + 1;
    end
    
    % Adjust hours to stay within 12-hour format
    if newHr < 1
        newHr = newHr + 12; % Wrap around if hours go below 1
    elseif newHr > 12
        newHr = newHr - 12; % Wrap around if hours go beyond 12
    end
    newMinPlus1 = newMin + 1;
    % Ensure proper formatting for hours and minutes
    if newHr < 10
        newHr = ['0', num2str(newHr)];
    else
        newHr = num2str(newHr);
    end
    
    if newMin < 10
        newMin = ['0', num2str(newMin)];
    else
        newMin = num2str(newMin);
    end


    if newMinPlus1 < 10
        newMinPlus1 = ['0', num2str(newMinPlus1)];
    else
        newMinPlus1 = num2str(newMinPlus1);
    end

    newTime = [newHr, ':', newMin];
    sonHrMin(i, :) = newTime; 

    newTimePlus1 = [newHr, ':', newMinPlus1];
    sonHrMinPlus1(i, :) = newTimePlus1; 
end

for ii = 1:max(sonicationIndex)
     sonMin(ii,:) = sonHrMin(ii,4:end);
     sonMinPlus1(ii,:) = sonHrMinPlus1(ii,4:end);
end

fileID = fopen('run_log_extract.sh','w');
for ii = 1:max(sonicationIndex) 
    formatSpec1 = 'grep -e "%s:" -e "%s:" power_output.txt > power_son%i.txt\n';
    fprintf(fileID,formatSpec1,sonHrMin(ii,:),sonHrMinPlus1(ii,:),ii);
    formatSpec2 = 'sed "s/:/ /gw power_son%ix.txt" power_son%i.txt >tmp\n';
    fprintf(fileID,formatSpec2,ii,ii);
    formatSpec3 = 'awk ''{ if($2 == "%s" || $2 == "%s") {print}}'' power_son%ix.txt > power_son%i.txt\n\n';
    fprintf(fileID,formatSpec3,sonMin(ii,:),sonMinPlus1(ii,:),ii,ii);
    ii=ii+1;
end
fprintf(fileID,'rm power_son*x.txt\n');
fprintf(fileID,'rm power_outputx.txt\n');
fclose(fileID);



%%
% these are the frequencies and efficiencies from /CPC/Site/Ini/Xd_7207.ini
% file from treatment export; they should not change until we have a
% different transducer so I'm not reading it new every time
freq = [335000 630000 640000 650000 660000 670000 680000 690000 700000];
efficiency = [0.6 0.59 0.73 0.76 0.78 0.83 0.81 0.85 0.77];

% now there is one file per sonication: e.g. power_son1.txt
cd(pathOut)
ii = 1;
for ii = 1:max(sonicationIndex) 
    % read sonication text file stripped log file
    filename_power = sprintf('power_son%i.txt',ii);
     % add check if file is empty - this happens when sonication was halted
    s = dir(filename_power);
    if s.bytes == 0
        disp('Sonication halted or stopped before power delivered.');
        disp('No .mat file created for this sonication');
    else
        powerTable = readtable(filename_power);
        timevector = table2array(powerTable(:,1:4));
        timeVectorSec = timevector(:,2)*60 + timevector(:,3) + timevector(:,4)/1000;
        timeDuration = timeVectorSec - timeVectorSec(1);
        timeDiffVector = diff(timeVectorSec); % taking difference shortens length of vector
        timeDiffVector(length(timeDiffVector)+1) = 0; % add zero element to vector to make it old length; power is zero then anyway
        power_elec = table2array(powerTable(:,8)); % this is electric power
        % correct with efficiency to get acoustic power
        %efficiency = 0.78; %@FREQ = 660000 MHz % may want to extract efficiciency from ini-file for other frequencies
        power_acou = power_elec * efficiency(find(freq==actualFrequency(ii)));
        energyFromFile = table2array(powerTable(:,10)).*efficiency(find(freq==actualFrequency(ii)));
        %energyDelta = power_acou(1:length(power)-1).*timeDiffVector;
        energyCalculated_acou = cumsum(power_acou.*timeDiffVector);
        energyCalculated_elec = cumsum(power_elec.*timeDiffVector);

       % save variables
        filenameSave = sprintf('son%i_var',ii);
            disp([filenameSave,' timevector(1,:)=',num2str(timevector(1,:))])  % check if the time stamps make sense
        cd(pathOut)    
        save(filenameSave,'timeDuration','power_acou','energyCalculated_acou','energyFromFile');
        %power_and_time{ii} = 
%         figure(20+ii)
%         yyaxis left;
%         plot(timeDuration,power_acou,'c');
%         ylabel('Power');
%         yyaxis right;
%         plot(timeDuration,energyCalculated_acou,'b'); hold on
%         %plot(timeDuration,energyCalculated_elec,'r');hold on
%         plot(timeDuration,energyFromFile,'g'); 
%         ylabel('Energy');
%         hold off
%         xlabel('time (sec)');
    end 
    ii=ii+1;
end

