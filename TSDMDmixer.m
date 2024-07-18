readmats; % Read Matrices.

%% Perform DMDd
dt = 2/30;
Time = (0:m-1)*dt; 
e = 1E-8;e1=e;
d = 4; 

[XrecDMD,GrowthRateX,FrequencyX,AmplitudeX,PhiXdmd] = DMDd(X,d,Time,e,e1);
%%
rPhiX = real(PhiXdmd);

% Delete redundant modes that come from taking only the real part (many
% modes come in conjugate pairs).

[~, index] = uniquetol(rPhiX',1E-5,'ByRows',true);
Psi_r = rPhiX(:, index);

%Construct the taylored base 
PsiPsir = Psi_r*Psi_r';
%% QR decomposition to get the optimal sensors in the pivots

[Q,R,pivot] = qr(PsiPsir','vector');
%%

sensnum = 360; % select a number of sensors
% construct the measurement matrix for the number of sensors
C = zeros(size(PsiPsir(:,1:sensnum)'));
for i = 1:sensnum
    C(i,pivot(i)) = 1;
end

% Find the indices of the nonzero elements of the measurement matrix
indxs = [];
indys = [];
for i = 1:sensnum
    C_temp = reshape(C(i,:)',mm,nn);
    [indxs(i),indys(i)] = find(C_temp);
end

% Plot the position of the sensors over the mean field.
xs = x(indxs);
ys = y(indys);

hfig = figure;
hold on
contourf(x,y,reshape(mean(X,2),mm,nn)','LineStyle','none')
colorbar
plot(xs,ys,'.k','MarkerSize',10,'LineWidth',2);
axis equal
hold on
plot(xwall,ywall,"ws",'MarkerFaceColor','w') %This is only for the mixer, so
% that the baffles are included in the plot
axis equal
title('Location of sensors')
picturewidth = 20; % set this parameter and keep it forever
set(findall(hfig,'-property','FontSize'),'FontSize',13)
set(findall(hfig,'-property','Box'),'Box','off') % optional
set(findall(hfig,'-property','LineWidth'),'LineWidth',1.5) % optional
set(hfig,'Units','Inches');

pos = get(hfig,'Position');
set(hfig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% Save the figure as a PNG file
    pngFileName = [name '_sensorposition.png'];  % Using 'snap' twice in the name for demonstration
    print(pngFileName, '-dpng', '-r600'); % '-r300' sets the resolution to 300 DPI
  
    %close(gcf);
%% tailored sensing of the validation data
Theta = C*Psi_r;
Ynot = C*Xnot;
Y = C*X;
for i = 1:10
    aa(:,i) = Theta\Ynot(:,i);
end
Xnotrec = Psi_r*aa;
% for i = 1:m

%% Plot of original and tailored-sensed snapshot

snapshot = 491;

hfig = figure;
subplot(1,2,1)
contourf(x,y,reshape(Xnot(:,1),mm,nn)','LineStyle','none')
title('Original snapshot ',num2str(snapshot));
axis equal
colorbar
hold on
plot(xwall,ywall,"ws",'MarkerFaceColor','w')
axis equal
caxis manual
%caxis(climits)
subplot(1,2,2)
contourf(x,y,reshape(Xnotrec(:,1),mm,nn)','LineStyle','none')
axis equal
titleString = ['Reconstruction of the snapshot ', num2str(snapshot), ' using ', num2str(sensnum), ' sensors'];
title(titleString)
colorbar
hold on
plot(xwall,ywall,"ws",'MarkerFaceColor','w')
axis equal
caxis manual
%caxis(climits)


picturewidth = 20; % set this parameter and keep it forever
set(findall(hfig,'-property','FontSize'),'FontSize',13)
set(findall(hfig,'-property','Box'),'Box','off') % optional
set(findall(hfig,'-property','LineWidth'),'LineWidth',1.5) % optional
%set(hfig,'Units','Inches');
%pos = get(hfig,'Position');
%set(hfig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

%% Create plot of reconstruction error vs number of sensors.
sensnums = 10:5:750;
clear recError
for jj = 1:length(sensnums)

sensnum = sensnums(jj);
C = zeros(size(PsiPsir(:,1:sensnum)'));
for i = 1:sensnum
    C(i,pivot(i)) = 1;
end

indxs = [];
indys = [];
for i = 1:sensnum
    C_temp = reshape(C(i,:)',mm,nn);
    [indxs(i),indys(i)] = find(C_temp);
end

Theta = C*Psi_r;
Ynot = C*Xnot;
Y = C*X;
for i = 1:10
    aa(:,i) = Theta\Ynot(:,i);
end
Xnotrec = Psi_r*aa;

Xrec = Psi_r*aa;
snapshot = 491;

recError(jj) = 100*norm(Xnot(:,1)-Xnotrec(:,1),"fro")/norm(Xnot(:,1),"fro");
end

hfig = figure;
plot(sensnums,recError)
ylabel('Percent error')
xlabel('Number of sensors')

picturewidth = 20; % set this parameter and keep it forever
set(findall(hfig,'-property','FontSize'),'FontSize',13)
set(findall(hfig,'-property','Box'),'Box','off') % optional
set(findall(hfig,'-property','LineWidth'),'LineWidth',1.5) % optional
set(hfig,'Units','Inches');
pos = get(hfig,'Position');
set(hfig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

% Save the figure as a PNG file
    pngFileName = [name '_percenterror.png'];  % Using 'snap' twice in the name for demonstration
    print(pngFileName, '-dpng', '-r600'); % '-r300' sets the resolution to 300 DPI
