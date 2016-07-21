function demo
% Demo of the M2DP descriptor described in the following paper:
%
% Li He, Xiaolong Wang and Hong Zhang, M2DP: A Novel 3D Point Cloud 
% Descriptor and Its Application in Loop Closure Detection, IROS 2016.
%
% Li He, Dept. of Computing Science, University of Alberta
% lhe2@ualberta.ca

%% 0. Please run mex CountVote2D.c in Matlab
% mex CountVote2D.c

%% 1. Read a point cloud
% 000000.bin is from KITTI dataset of sequence 00.
% http://www.cvlibs.net/datasets/kitti/eval_odometry.php
[data, numData] = ReadBinData('000000.bin');

% display the input data
figure(10);grid on;hold on
pcshow(data);title('Input point cloud');

%% 2. M2DP
% desM2DP: descriptor of data; A: 2D signatures of data.
tstart = tic;
[desM2DP, A] = M2DP(data);
telapsed = toc(tstart);
disp(['Processing time of M2DP on ' num2str(numData) ' points: ' num2str(telapsed) ' seconds']);



function [data, numData] = ReadBinData(nameDataFile)

fid = fopen(nameDataFile,'r');
if fid==-1
    disp('Cannot find input file\n');
    data=[];
    numData=[];
    return;
end

% get the size of file
fseek(fid,0,'eof');
fsize = ftell(fid);
% number points = total bytes / 4 per float32 / 4 float32 per row
numData = floor(fsize/16);

% now, read data
fseek(fid,0,'bof');
data = fread(fid,fsize,'float32');
fclose(fid);

% reshape data
data = reshape(data,[4, numData]);
data = data';
% ignoring the last feature, laser intensity
data = data(:,1:3);