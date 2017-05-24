%% EEG preprocessing

clear all; restoredefaultpath; addpath('C:\work\spm12_6470\'); spm('defaults', 'eeg');  
anapath='C:\work\numcum_analysis\';

addpath(genpath('C:\Users\Berni\Dropbox\numcum\scripts'));
anadir=[anapath 'proc10b_home\']; mkdir(anadir); cd(anadir);

subjects={'01';'02';'03';'04';'05';'06';'07';'08';'09';'10';'11';'12';'13';'14';'15';'16';'17';'18';'19';'20';'21';'22';'23';'24'}; %N=24

tic;
% [check below if and where badchannels are excluded]
badchans=  {'TP7','F5','',''    % sub01
            '','','',''         % sub02
            '','','',''         % sub03
            '','','',''         % sub04
            'Oz','','',''       % sub05
            'F4','','',''       % sub06      
            '','','',''         % sub07         
            '','','',''         % sub08
            'F4','','',''       % sub09
            '','','',''         % sub10
            'P8','P4','',''     % sub11
            '','','',''         % sub12 
            '','','',''         % sub13 
            '','','',''         % sub14 
            'P7','','',''       % sub15 
            'P4','','',''       % sub16 
            '','','',''         % sub17  
            'F1','PO3','O1',''  % sub18
            '','','',''         % sub19
            'F4','FC2','',''    % sub20
            'TP8','','',''      % sub21 
            '','','',''         % sub22 
            'O1','','',''       % sub23 
            'P1','','',''};     % sub24 

correyes=subjects;  
blinkthr=repmat(5,length(subjects),1); 
numcomps=ones(length(subjects),1); 
 
for sub=1:length(subjects)
     
%% convert 

    S = [];
    S.dataset = [anapath '\eeg\nce' subjects{sub} '.bdf'];
    S.outfile = [anadir 'nc' subjects{sub} '.mat'];
    D = spm_eeg_convert(S);
     
%% Downsample

    S = [];
    S.D = [anadir D.fname];
    S.fsample_new=256; 
    S.prefix='d';
    D = spm_eeg_downsample(S)
    delete([S.D(1:end-3),'mat']);
    delete([S.D(1:end-3),'dat']);
   
%% Reference/Montage
    
    %D=spm_eeg_load([anadir 'dnc' subjects{sub} '.mat']); % temp load D
    S = []; 
    S.D = [anadir D.fname];
    tra=zeros(66,D.nchannels);
    badind=D.indchannel(badchans(str2num(subjects{sub}),:))
    ngood=64-length(badind);
    tra(1:64,1:64)=eye(64)-1/ngood; %average reference of 64 EEG channels, ecluding visually identified badchannels
    tra(:,badind)=0; % exclude badchannels from referencing 
    for i=1:length(badind)
        tra(badind(i),badind(i))=1; % but reference badchannels to avg of others 
    end
    tra(65,[D.indchannel('Fp2') 67])=[1 -1]; %turn eye channels 67 and Fp2 into bipolar (pseudo) VEOG channel (65);
    tra(66,65:66)=[1 -1]; %turn eye channels 65&66 into bipolar HEOG channel (66);
   % figure; imagesc(tra); title(['montage matrix sub' subjects{sub}]);
   % xlabel('old channels'); ylabel('new channels'); colorbar;  % for checking
    S.montage.tra=tra;
    S.montage.labelnew=[D.chanlabels(1:64) 'VEOG' 'HEOG'];
    S.montage.labelorg=D.chanlabels;
    S.keepothers=0; 
    [D, montage] = spm_eeg_montage(S);
    D=D.badchannels(badind, 1);
    save(D);
   
%% high-pass filter
 
    S = [];
    S.D = [anadir D.fname];
    S.band = 'high';
    S.freq = 0.5;% 1 0.5;  % 0.1;
    % S.order = 3;
    D = spm_eeg_filter(S);
    delete([S.D(1:end-3),'mat']); % remove unused files
    delete([S.D(1:end-3),'dat']);

%% correct ocular artefacts
    %D=spm_eeg_load([anadir 'fMdnc' subjects{sub} '.mat']);
    if ismember(subjects{sub}, correyes)   
                
        matlabbatch={};
        matlabbatch{1}.spm.meeg.source.headmodel.D = {[anadir D.fname];};
        matlabbatch{1}.spm.meeg.source.headmodel.val = 1;
        matlabbatch{1}.spm.meeg.source.headmodel.comment = '';
        matlabbatch{1}.spm.meeg.source.headmodel.meshing.meshes.template = 1;
        matlabbatch{1}.spm.meeg.source.headmodel.meshing.meshres = 2;
        if ismember(subjects{sub}, dontelpos)  
            matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregdefault = 1;
        else
            matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(1).fidname = 'fidt9';
            matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(1).specification.select='lpa';
            matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(2).fidname = 'fidnz';
            matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(2).specification.select='nas';
            matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(3).fidname = 'fidt10';
            matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(3).specification.select='rpa';
            matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregspecify.useheadshape=1;
        end
        matlabbatch{1}.spm.meeg.source.headmodel.forward.eeg = 'EEG BEM';
        matlabbatch{1}.spm.meeg.source.headmodel.forward.meg = 'Single Shell';
        spm_jobman('serial', matlabbatch);
        parkD=spm_eeg_load([anadir D.fname]); 
        
%% this for blinks 
      
        S = [];
        S.D = [anadir D.fname];
        S.mode = 'mark';
        S.badchanthresh = 1;
        S.methods.channels = {'VEOG'};
        S.methods.fun = 'eyeblink';
        S.methods.settings.threshold = blinkthr(str2num(subjects{sub}));
        S.methods.settings.excwin = 0;
        S.append = 0;
        S.prefix = 'a';
        D = spm_eeg_artefact(S);

        S = [];
        S.D = [anadir D.fname];
        S.timewin = [-500 500];
        S.trialdef.conditionlabel = 'Eyeblink';
        S.trialdef.eventtype = 'artefact_eyeblink';
        S.trialdef.eventvalue = 'VEOG';
        S.trialdef.trlshift = 0;
        S.bc = 1;
        S.prefix = 'eyeblink';
        S.eventpadding = 0;
        D = spm_eeg_epochs(S);
        delete([S.D(1:end-3),'mat']); % remove unused files
        delete([S.D(1:end-3),'dat']);

        % define I
        S = [];
        S.D = [anadir D.fname];;
        S.mode = 'svd';
        S.ncomp = 1;
        S.timewin = [-Inf Inf];
        D = spm_eeg_spatial_confounds(S);

        % define II
        matlabbatch={}
        matlabbatch{1}.spm.meeg.preproc.sconfounds.D = {[anadir parkD.fname]};
        matlabbatch{1}.spm.meeg.preproc.sconfounds.mode{1}.spmeeg.conffile = {[anadir D.fname]};
        spm_jobman('serial', matlabbatch);
        
%% correct
        matlabbatch={}
        matlabbatch{1}.spm.meeg.preproc.correct.D = {[anadir parkD.fname]};
        matlabbatch{1}.spm.meeg.preproc.correct.mode = 'berg';
        matlabbatch{1}.spm.meeg.preproc.correct.prefix = 'T';
        spm_jobman('serial', matlabbatch);
        
      %  delete([S.D(1:end-3),'mat']); % remove eyeblink files
      %  delete([S.D(1:end-3),'dat']);
%         delete([anadir parkD.fname]); % remove unused files
%         delete([parkD.fnamedat]);
    end

 %% low-pass filter
  
    D=spm_eeg_load([anadir 'T' parkD.fname]); %temp load
 
    S = [];
    S.D = [anadir D.fname];
    S.band = 'low';
    S.freq = 45;  
    D = spm_eeg_filter(S);
    delete([S.D(1:end-3),'mat']); % remove unused files
    delete([S.D(1:end-3),'dat']);       
    
end

runtime=toc/3600;    



 