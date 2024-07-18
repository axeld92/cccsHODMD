function [Vreconst,growthRates,frequencies,amplitudes,Un] = DMDd(V,d,Time,e1,e2)
% % Inputs: - V: Data matrix (snapshots matrix)
%           - d: Parameter d of the DMD-d method
%           - Time: vector of the time values. It should start with 0.
%           - r: number of modes retained in each of the two dimension
%           reductions.
% % Outputs:- Vreconst: DMD-d reconstruction of the original data
%           - the following three are self-explanatory
%           - Un: Matrix of the DMD modes
%
% % This function is the implementation of a modification of the DMDd 
% algotrithm, as presented by Vega and LeClainche. The difference is that 
% the amplitudes of the modes are computed by right-
% multiplying the pinv of the DMD modes and first snapshot. Therfore this
% code assumes that Time(1) = 0. This method of computing the amplitudes is
% the same as the one presented for the standard DMD algorithm in the book
% by Kutz, Brunton, Brunton and Proctor. It results on a much faster
% computation of the DMD-d than the one presented on the book by Vega and
% LeClainche.

% Written by Axel Dullak

% Compute the resolution of the time vector
dt = Time(2)-Time(1);

% Step one: SVD the V matrix
[U,S,T] = svd(V,'econ');
Svec = diag(S);
NormSvec = norm(Svec); % Calculate the norm-2 of the vector of singular values
rat1 = 10; % Initialize the rat1 (ratio 1) variable, value must be > e1
ii = 1; % Initialize counter

while rat1>e1
     rat1 = norm(Svec(ii:end))/NormSvec; %calculate the cutoff point 
                                % where the norm of the 
     ii = ii + 1;
end

r1 = ii;
fprintf('cutoff for the first dimension reduction is %d\n', r1);




Ur = U(:,1:r1);
Sr = S(1:r1,1:r1);
Tr = T(:,1:r1);

% Reduced snapshot matrix
% Vhat = Ur'*V; %hatT in original code
Vhat = Sr*Tr'; %This is how they calculate it

[~,K] = size(Vhat); % K is the number of snapshots

% Calculate the enlarged reduced snapshots matrix

Vstar = zeros(r1*d,K-d+1);
 for k = 1:K-d+1
     Vstar(:,k) =  reshape(Vhat(:,k:k+d-1),[],1);
 end

% Vtilde is the tildeT in the original code
% SVD the matrix Vstar ---> Vhatstar

[U1star,Sstar,U2star] = svd(Vstar,'econ');
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
 
 
U1star = U1star(:,1:r2);
Sstar = Sstar(1:r2,1:r2);
U2star = U2star(:,1:r2);

%Comupute the reduced-enlarged-reduced matrix
Vstarhat = Sstar*U2star'; %hatT1 in the original code

% Produce the matrices from time 1 to end-1 and from time 2 to end
V1 = Vstarhat(:,1:end-1); 
V2 = Vstarhat(:,2:end);

% SVD to compute the reduced modified approximate Koopman operator

[U3star,S1star,U4star] = svd(V1);

% Compute the approx reduced K-operator

Rstarhat = V2*U4star*pinv(S1star)*U3star';

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

%Now compute the amplitudes with the first snapshot

v1 = V(:,1);

amplitudes = Un\v1;

% Reconstruct the signal

time_dynamics = zeros(r2,K);

for iter = 1:K
    time_dynamics(:,iter) = amplitudes.*exp(logmu*Time(iter));
end

Vreconst = Un*time_dynamics;

%Vreconst = real(Vreconst);
end
