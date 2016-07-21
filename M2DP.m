function [desM2DP,A] = M2DP(data)
% Multiview 2D projection (M2DP) descriptor. 
%
% Input:
%       data        n*3     Point cloud. Each row is [x y z]
% Output:
%       desM2DP     192*1   M2DP descriptor of the input cloud data
%       A           64*128  Signature matrix
%
% Introduction:
% M2DP is a global descriptor of input point cloud. Details of M2DP can be 
% found in the following paper:
%
% Li He, Xiaolong Wang and Hong Zhang, M2DP: A Novel 3D Point Cloud 
% Descriptor and Its Application in Loop Closure Detection, IROS 2016.
%
% Li He, Dept. of Computing Science, University of Alberta
% lhe2@ualberta.ca


%% 1. Initialization
% key parameter
% number of bins in theta, the 't' in paper
numT = 16;
% number of bins in rho, the 'l' in paper
numR = 8;
% number of azimuth angles, the 'p' in paper
numP = 4;
% number of elevation angles, the 'q' in paper
numQ = 16;

% rotation invariant
data = PCARotationInvariant(data);

% Azimuthe list
azimuthList = linspace(-pi/2,pi/2,numP);

% Elevation list
elevationList = linspace(0,pi/2,numQ);

% get the farthest point distance
rho2 = sum(data.^2,2);
maxRho = sqrt(max(rho2));

% main function, get the signature matrix A
A = GetSignatureMatrix(azimuthList, elevationList, data, numT, numR, maxRho);

% run SVD on A and use [u1,v1] as the final output
[u,s,v] = svd(A);
desM2DP = [u(:,1);v(:,1)];


function A = GetSignatureMatrix(azimuthList, elevationList, data, numT, numR, maxRho)

% signature matrix A
A = zeros(length(azimuthList)*length(elevationList),numT*numR);
% plane index
n = 1;

% theta list
thetaList = linspace(-pi,pi,numT+1);
% rho list
rhoList = linspace(0,sqrt(maxRho),numR+1);
rhoList = rhoList.^2;
rhoList(end) = rhoList(end)+.001; % make sure all points in bins

% loop on azimuth
for p=1:length(azimuthList)
    azm = azimuthList(p); % pick one azimuth
    
    % loop on evevation
    for q=1:length(elevationList)
        elv = elevationList(q); % pick one elevation
        
        % normal vector vecN of the selected 2D plane
        [x,y,z] = sph2cart(azm,elv,1);
        vecN = [x y z];
        
        % distance of vector [1,0,0] to the surface with normal vector vecN
        h = [1 0 0]*vecN';
        
        % a new vector, c = h*vecN, so that vector [1,0,0]-c is the
        % projection of x-axis onto the plane with normal vector vecN
        c = h*vecN;
        
        % x-axis - c, the projection
        px = [1 0 0]-c;
        
        % given the normal vector vecN and the projected x-axis px, the
        % y- axis is cross(vecN,px)
        py = cross(vecN,px);
        
        % projection of data onto space span{px,py}
        pdata = [data*px' data*py'];
        
        % represent data in polar coordinates
        [theta,rho] = cart2pol(pdata(:,1),pdata(:,2));
        
        % main function, count points in bins, use eigher a) C codes or b)
        % Matlab built-in function
        % a) C codes
        bin = CountPoint(theta, thetaList, rho, rhoList);
%         % b) Matlab codes
%         bin = histcounts2(theta',rho',thetaList,rhoList)/size(data,1);
%         bin = bin(:);
        
        % record the sigature of the n-th plane
        A(n,:) = bin';
        n = n+1;
    end
end


function data = PCARotationInvariant(data)
% 3D rotate input data so that x-axis and y-axis are the 1st and 2nd PCs of
% data, respectively

n = size(data,1);
% shift data so that origin is the mean of data
% mean of data
md = mean(data);
data = data-repmat(md,n,1);

pc = pca(data);

X = pc(:,1)'*data';
Y = pc(:,2)'*data';
Z = pc(:,3)'*data';
data = [X' Y' Z'];


function bin = CountPoint(theta, thetaList, rho, rhoList)
% Count points in bins
% use C code
bin = CountVote2D(theta, rho,thetaList', rhoList');

