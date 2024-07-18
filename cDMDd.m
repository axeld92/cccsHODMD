function [Xreconst,growthRates,frequencies,amplitudes,Un] = cDMDd(X,d,Time,e1,e2);
% % Inputs: - V: Data matrix (snapshots matrix)
%           - d: Parameter d of the DMD-d method
%           - Time: vector of the time values. It should start with 0.
%           - r: number of modes retained in each of the two dimension
%           reductions.
% % Outputs:- Vreconst: DMD-d reconstruction of the original data
%           - the following three are self-explanatory
%           - Un: Matrix of the DMD modes
%
%%% Application of compression to the DMDd algorithm according to the
%%% description of compressed DMD by Kutz, Brunton, Brunton and Proctor.

% Written by Axel Dullak

% Compute the resolution of the time vector
dt = Time(2)-Time(1);
[m,n] = size(X);
p = ceil(m/100); %Compress aproximately 100fold
rng default
C = randn(p,m); %Compression matrix
Y = C*X; %Compressed data
% C and Psi (Sparse basis) are incoherent 


% Step one: SVD the V matrix
[U,S,V] = svd(Y,'econ');
Svec = diag(S);
NormSvec = norm(Svec);
rat1 = 10;
ii = 1;
while rat1>e1
     rat1 = norm(Svec(ii:end))/NormSvec;
     ii = ii + 1;
end

r1 = ii;
fprintf('cutoff for the first dimension reduction is %d\n', r1);

Ur = U(:,1:r1);
Sr = S(1:r1,1:r1);
Vr = V(:,1:r1);

% Reduced snapshot matrix

Xhat = Sr*Vr'; %This is how they calculate it

[~,K] = size(Xhat); % K is the number of snapshots

% Calculate the enlarged reduced snapshots matrix

Xstar = zeros(r1*d,K-d+1);
 for k = 1:K-d+1
     Xstar(:,k) =  reshape(Xhat(:,k:k+d-1),[],1);
 end


% SVD the matrix Vstar ---> Vhatstar

[U1star,Sstar,U2star] = svd(Xstar,'econ');
Sstarvec = diag(Sstar);
NormSstarvec = norm(Sstarvec);
ii = 1;
r2 = 10;
while r2 > e2
    r2 = norm(Sstarvec(ii:end))/NormSstarvec;
    ii = ii+1;
end
r2 = ii;
fprintf('cutoff for the second dimension reduction is %d\n', r2);
 
% Compute the reduced matrices

U1star = U1star(:,1:r2);
Sstar = Sstar(1:r2,1:r2);
U2star = U2star(:,1:r2);

%Comupute the reduced-enlarged-reduced matrix
Xstarhat = Sstar*U2star'; %hatT1 in the original code

Y = Xstarhat;
% Produce the matrices from time 1 to end-1 and from time 2 to end

Y1 = Y(:,1:end-1); 
Y2 = Y(:,2:end);

% SVD to compute the reduced modified approximate Koopman operator

[U3star,S1star,U4star] = svd(Y1);

% Compute the approx reduced K-operator

Rstarhat = Y2*U4star*pinv(S1star)*U3star';

% Get the eigenvalues and eigenvectors

[Qstarhat,M] = eig(Rstarhat);

mus = diag(M);
logmu = log(mus)/dt;
growthRates = real(logmu);
frequencies = imag(logmu);

% Transform Qstarhat to Qstar (reduced-enlarged-reduced to enlarged-reduced)

Qstar = U1star*Qstarhat;

% Get Qhat from the first Nhat (r) components of Qstar

Qhat = Qstar(1:r1,:);

% unhat are the rescaled versions of qnhats, the columns of Qhat. I put them
% in a matrix Unhat to differentiate from all the other matrices named
% Usomething

Unhat = bsxfun(@rdivide, Qhat, vecnorm(Qhat));

Un = Ur*Unhat; % These are the DMD modes!
Un = pinv(C)*Un; % Uncompress the Un
%Now compute the amplitudes with the first snapshot

x1 = X(:,1);

amplitudes = Un\x1;

% Reconstruct the signal

time_dynamics = zeros(r2,K);

for iter = 1:K
    time_dynamics(:,iter) = amplitudes.*exp(logmu*Time(iter));
end

Xreconst = Un*time_dynamics;

%Xreconst = real(Xreconst); %Could be necessary to take only the real part
end
