function varargout = Field_Tomogram_Reconstruction(varargin)
% FIELD_TOMOGRAM_RECONSTRUCTION MATLAB code for Field_Tomogram_Reconstruction.fig
%      FIELD_TOMOGRAM_RECONSTRUCTION, by itself, creates a new FIELD_TOMOGRAM_RECONSTRUCTION or raises the existing
%      singleton*.
%
%      H = FIELD_TOMOGRAM_RECONSTRUCTION returns the handle to a new FIELD_TOMOGRAM_RECONSTRUCTION or the handle to
%      the existing singleton*.
%
%      FIELD_TOMOGRAM_RECONSTRUCTION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FIELD_TOMOGRAM_RECONSTRUCTION.M with the given input arguments.
%
%      FIELD_TOMOGRAM_RECONSTRUCTION('Property','Value',...) creates a new FIELD_TOMOGRAM_RECONSTRUCTION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Field_Tomogram_Reconstruction_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Field_Tomogram_Reconstruction_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Field_Tomogram_Reconstruction

% Last Modified by GUIDE v2.5 08-Aug-2020 21:59:03

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @Field_Tomogram_Reconstruction_OpeningFcn, ...
    'gui_OutputFcn',  @Field_Tomogram_Reconstruction_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before Field_Tomogram_Reconstruction is made visible.
function Field_Tomogram_Reconstruction_OpeningFcn(hObject, eventdata, handles, varargin)
global radVal;global clickToken;global signVal;global projectionToken;
handles.output = hObject;
guidata(hObject, handles);
radVal=3;
% axes(handles.axes1);
% a = imread('coins.png');
% imagesc(a),axis image
set(handles.radiusVal,'string',num2str(radVal));
set(handles.drawBtn,'value',1);

clickToken=0;
signVal=1;
projectionToken=1;

handles.output = hObject;

set(handles.analysisOption2D,'userdata',1);
set(handles.analysisOption3D,'userdata',1);
set(handles.sliceSelection,'userdata',1);
set(handles.sampleSelect2D,'userdata',1);
set(handles.sampleSelect2D,'string',1);
set(handles.sampleSelect3D,'userdata',1);
set(handles.sampleSelect3D,'string',1);
set(handles.bgNum,'string',num2str(1));
set(handles.n_m,'string',num2str(1.337));
set(handles.overallPhase,'string',num2str(1));
set(handles.deltaPhase,'string',num2str(0.05));

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Field_Tomogram_Reconstruction wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Field_Tomogram_Reconstruction_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

function folderPath_Callback(hObject, eventdata, handles)
global folderPath
if get(hObject,'Value')==1
    folderPath = uigetdir('D:','Select the folder where data are saved.');
    if folderPath~=0
        set(handles.pathIndicator,'string',folderPath);
        cd(folderPath);
        samplePath=dir('sample*.mat');
        set(handles.sampleNum,'string',num2str(length(samplePath)));
        set(handles.sampleNum,'userdata',samplePath);
        bgPath=dir('bg*');
        bgNumSel=str2double(get(handles.bgNum,'string'));
        load(bgPath(bgNumSel).name);
        
        set(handles.wavelength,'String',num2str(lambda));
        set(handles.pixelSize,'String',num2str(res));
        set(handles.magnification,'String',num2str(4.8/res));
        set(handles.na_value,'String',num2str(NA));
        set(handles.bgPathInd,'string',strcat(folderPath,'\',bgPath(bgNumSel).name));
        savePath=dir('field_retrieval');
        if length(savePath)<1
            mkdir(folderPath,'field_retrieval')
        end
        set(handles.savePathInd,'string',strcat(folderPath,'\','field_retrieval'));
    end
end

function sampleSelect3D_Callback(hObject, eventdata, handles)
set(handles.sampleSelect3D,'userdata',str2double(get(hObject,'String')));

function sampleSelect3D_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function sampleSelect2D_Callback(hObject, eventdata, handles)
set(handles.sampleSelect2D,'userdata',str2double(get(hObject,'String')));

function sampleSelect2D_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function analysisOption2D_SelectionChangeFcn(hObject, eventdata, handles)
switch hObject
    case handles.singleSample2D
        set(handles.analysisOption2D,'userdata',1);
    case handles.analyzeAll2D
        set(handles.analysisOption2D,'userdata',1:length(get(handles.sampleNum,'userdata')));
end

function analysisOption3D_SelectionChangeFcn(hObject, eventdata, handles)
switch hObject
    case handles.singleSample3D
        set(handles.analysisOption3D,'userdata',1);
    case handles.analyzeAll3D
        set(handles.analysisOption3D,'userdata',1:length(get(handles.sampleNum,'userdata')));
end

function analyze2D_Callback(hObject, eventdata, handles)
global imField;global Mask;global posX;global posY;
set(handles.analyze2DMsg,'string','Field Retrieval Start');
drawnow

bgpath=get(handles.bgPathInd,'string');
folderPath=get(handles.pathIndicator,'string');
spath=get(handles.savePathInd,'string');
n_m=str2double(get(handles.n_m,'String'));

cd(folderPath);
load(bgpath);
img=double(squeeze(tomogMap(:,:,round(size(tomogMap,3)/3)-1)));
%     img=double(imread(bglist(1).name,49));
img=squeeze(img);
ii=length(img);

ZP=ii;
r=round(ZP*res*NA/lambda)+20;
yr = r;

f_dx=[];
f_dy=[];
for bgNum=1:size(tomogMap,3)
    img=double(squeeze(tomogMap(:,:,bgNum)));
    %     subplot(121),  imagesc(((img))),axis image,colormap('jet')
    [xSize ySize]=size(img);
    Fimg = fftshift(fft2(img))/(xSize*ySize); %FFT
    
    [f_x,f_y]=find(Fimg==max(max(Fimg(:, round(ii*0.05):round(ii*0.45) ))));
    f_dx=[f_dx;f_x];
    f_dy=[f_dy;f_y];
end

mi=mean(f_dx);
mj=mean(f_dy);
mi=round(mi-ii/2-1); mj=round(mj-ii/2-1);

c1mask = ~(mk_ellipse(r-20,yr-20,ZP,ZP));
c3mask = circshift(c1mask,[mi mj]);

Fbg=zeros(round(2*r),round(2*r),size(tomogMap,3),'single');
for bgNum=1:size(tomogMap,3)
    img=double(squeeze(tomogMap(:,:,bgNum)));
    [xSize ySize]=size(img);
    Fimg = fftshift(fft2(img))/(xSize*ySize); %FFT
    Fimg = Fimg.*c3mask;
    Fimg = circshift(Fimg,-[round(mi) round(mj)]);
    Fimg = Fimg(ii/2-r+1:ii/2+r,ii/2-r+1:ii/2+r);
    sizeFimg=length(Fimg);
    Fimg = ifft2(ifftshift(Fimg))*(sizeFimg^2);
    Fbg(:,:,bgNum)=Fimg;
end
%%
sampleNum=length(get(handles.analysisOption2D,'userdata'));
sampleList=get(handles.sampleNum,'userdata');

if sampleNum==1
    sampleNumList=get(handles.sampleSelect2D,'userdata');
else
    sampleNumList=1:sampleNum;
end

for sampleNumVal=sampleNumList(1):sampleNumList(end)
    cd(folderPath);
    load(sampleList(sampleNumVal).name);
    retPhase=zeros(round(2*r)-4,round(2*r)-4,size(tomogMap,3),'single');
    retAmplitude=zeros(round(2*r)-4,round(2*r)-4,size(tomogMap,3),'single');
    for iter=1:size(tomogMap,3)
        img=(double(squeeze(tomogMap(:,:,iter))));
        Fimg = fftshift(fft2(img))/(ii^2); %FFT
        Fimg = Fimg.*c3mask;
        Fimg = circshift(Fimg,-[round(mi) round(mj)]);
        Fimg = Fimg(ii/2-r+1:ii/2+r,ii/2-r+1:ii/2+r);
        sizeFimg=length(Fimg);
        Fimg = (ifft2(ifftshift(Fimg))*(sizeFimg^2));
        
        Fimg=Fimg./squeeze(Fbg(:,:,iter));
        Fimg=Fimg(3:end-2,3:end-2);
        retAmplitude(:,:,iter)=abs(Fimg);
        
        p=((unwrap2(double(angle(Fimg)))));
        [p,coeffVal,~]=phaseCompensation(p,1);
        p=p-mode(round(p(:)*100)/100);
        %%
        retPhase(:,:,iter)=p;
        
        if sampleNum==1
            set(handles.analyze2DMsg,'string',strcat('Field Retrieval ',sprintf('%04.2f',100*(iter/size(tomogMap,3))),'% Done'));
        else
            set(handles.analyze2DMsg,'string',strcat('Field Retrieval ',sprintf('%04.2f',100*((sampleNumVal-1)*size(tomogMap,3)+iter)/(size(tomogMap,3)*length(sampleNumList))),'% Done'));
        end
        drawnow
    end
    [xx yy frame]=size(retPhase);
    
    fig=figure;
    for kk=1:6
        subplot(3,3,kk),imagesc(squeeze(retPhase(:,:,round(size(tomogMap,3)/6*kk))),[-2 2]),axis image,axis off
    end
    subplot(3,3,[7 9]),imagesc(squeeze(retPhase(:,round(yy/2),:)),[-2 2])
    colormap('jet')
    fileName=strcat('Field_',sampleList(sampleNumVal).name);
    
    cd(spath)
    save(strcat(fileName(1:end-4),'.mat'),'retAmplitude','retPhase','xSize','f_dx','f_dy','NA','lambda','res','ZP','n_m');
    print(fig,'-dpng',strcat(fileName(1:end-4),'.png'))
    
    close(fig)
    set(handles.frameNumShow,'sliderstep',[1/frame 10/frame]);
    axes(handles.plot2D);
    imagesc(squeeze(retPhase(:,:,round(size(retPhase,3)/3)-1)),[-2 2])
    imField=squeeze(retPhase(:,:,round(size(retPhase,3)/3)-1));
end

Mask=zeros(size(imField));
[posX,posY]=meshgrid((1:size(Mask,2)),(1:size(Mask,1)));
set(handles.analyze2DMsg,'string',strcat('Field Retrieval Done'));

function analyze3D_Callback(hObject, eventdata, handles)
global Reconimg;global res3;global res4;global excludeFrame;global Mask;global imField;
global signPhase;
set(handles.analyze3DMsg,'string','Tomogram Reconstruction Start');
drawnow

spath=get(handles.savePathInd,'string');
sampleList=get(handles.sampleNum,'userdata');

sampleNum=length(get(handles.analysisOption3D,'userdata'));
if sampleNum==1
    sampleNumList=get(handles.sampleSelect3D,'userdata');
else
    sampleNumList=1:sampleNum;
end

for sampleNumVal=sampleNumList(1):sampleNumList(end)
    cd(spath)
    load(strcat('Field_',sampleList(sampleNumVal).name));
    [xx, yy, frame]=size(retPhase);
    crop_size=xx;
    f_dx2=f_dx-mean(f_dx(:)); % subtract maxpoint
    f_dy2=f_dy-mean(f_dy(:));
    original_size=xSize;
    
    n_s=n_m+0.035;
    ZP=round(1.2*xx/2)*2;
    crop_factor=crop_size/original_size;
    res2=res/crop_factor;
    padd_factor=ZP/crop_size;
    kres=1/(res*ZP)*crop_factor;
    f_dx2=f_dx2*padd_factor;f_dy2=f_dy2*padd_factor;
    k0_x=kres*f_dx2; % for actual value, multiply resolution
    k0_y=kres*f_dy2;
    k0=1/lambda;
    k0_z=real(sqrt((n_m*k0)^2-(k0_x).^2-(k0_y).^2)); % magnitude of absolute value is k0
    
    set(handles.n_m,'userdata',n_m);
    %%
    if isempty(excludeFrame)
        excludeFrame=[];
    end
    for kkk=1:frame
        p2=squeeze(retPhase(:,:,kkk));
        if sum(isnan(p2(:)))
            excludeFrame=[excludeFrame,kkk];
        end
    end
    %%
    excludeFrame=unique(excludeFrame);
    
    % Range=[54 539];
    ZP2=round(512*1.0/2)*2;
    ZP3=round(256*1.0/2)*2;
    res3=res2*ZP/ZP2;
    res4=res2*ZP/ZP3;
    frameList=(1:frame);
    frameList(excludeFrame)=[];
    IterNum=100;
    TomoParam=struct('n_m',n_m,'ZP',ZP,'ZP2',ZP2,'ZP3',ZP3,'xx',xx,'yy',yy,...
        'f_dx2',f_dx2,'f_dy2',f_dy2,'NA',NA,'lambda',lambda,...
        'k0',k0,'k0_x',k0_x,'k0_y',k0_y,'k0_z',k0_z,'kres',kres,'frameList',frameList,...
        'res2',res2,'res3',	res3,'res4',res4,'signPhase',signPhase);
    %% Tomogram Reconstruction
    [Reconimg, ORytov]=ODTReconstruction(retAmplitude,retPhase,TomoParam);
    Reconimg=real(ODTIteration(gpuArray((Reconimg)),ORytov,TomoParam,IterNum));
    ZP4=round(size(Reconimg,1)/(2*padd_factor));
    ZP5=round(size(Reconimg,3)/(4*padd_factor));
    
    for zz=1:size(Reconimg,3)
        Reconimg(:,:,zz)=conv2(Reconimg(:,:,zz),fspecial('disk',0.7),'same');
    end
    
    Reconimg=real(Reconimg(end/2-ZP4+1:end/2+ZP4,end/2-ZP4+1:end/2+ZP4,end/2-ZP5+1:end/2+ZP5));
    ZP4=round(size(Reconimg,1));
    ZP5=round(size(Reconimg,3));
    
    fig=figure;
    %%
    subplot(221),imagesc(((1:ZP4)-ZP4/2)*res3,((1:ZP5)-ZP5/2)*res4,((squeeze(Reconimg(end/2,:,:)))'),[n_m-0.005 n_s]),axis image,colorbar
    subplot(222),imagesc(((1:ZP4)-ZP4/2)*res3,((1:ZP5)-ZP5/2)*res4,((squeeze(Reconimg(:,end/2,:)))'),[n_m-0.005 n_s]),axis image,colorbar
    subplot(223),imagesc(((1:ZP4)-ZP4/2)*res3,((1:ZP4)-ZP4/2)*res3,((squeeze(Reconimg(:,:,end/2)))),[n_m-0.005 n_s]),axis image,colorbar
    subplot(224),imagesc(max((Reconimg),[],3),[n_m-0.005 n_s]),axis image, axis off
    colormap('jet')
    
    axes(handles.plot3D);
    imagesc((squeeze(Reconimg(:,:,end/2))),[n_m-0.005 n_s]),axis image,colorbar,axis off
    %     set(handles.plot3D,'userdata',Reconimg);
    
    fileName=strcat('Tomogram_',sampleList(sampleNumVal).name);
    save(strcat(fileName(1:end-4),'.mat'),'n_m','Reconimg','res3','res4','lambda');
    print(fig,'-dpng',strcat(fileName(1:end-4),'.png'))
    close(fig)
    
    [x,y,z]=size(Reconimg);
    set(handles.xPosition,'userdata',x);
    set(handles.yPosition,'userdata',y);
    set(handles.zPosition,'userdata',z);
    stringMsg=get(handles.positionMsg,'string');
    positionVal=[round(x/2) round(y/2) round(z/2)];
    set(handles.positionMsg,'userdata',positionVal);
    stringMsg{1}=strcat('x : ',num2str(positionVal(1)));
    stringMsg{2}=strcat('y : ',num2str(positionVal(2)));
    stringMsg{3}=strcat('z : ',num2str(positionVal(3)));
    set(handles.positionMsg,'string',stringMsg);
    
    set(handles.analyze3DMsg,'string',strcat('Tomogram Reconstruction ',sprintf('%04.2f',100*sampleNumVal/length(sampleNumList)),'% Done'));
    drawnow
end

function xBackward_Callback(hObject, eventdata, handles)
moveBackward(handles,1);

function yBackward_Callback(hObject, eventdata, handles)
moveBackward(handles,2);

function zBackward_Callback(hObject, eventdata, handles)
moveBackward(handles,3);

function xForward_Callback(hObject, eventdata, handles)
xSize=get(handles.xPosition,'userdata');
moveForward(handles,1,xSize);

function yForward_Callback(hObject, eventdata, handles)
ySize=get(handles.yPosition,'userdata');
moveForward(handles,2,ySize);

function zForward_Callback(hObject, eventdata, handles)
zSize=get(handles.zPosition,'userdata');
moveForward(handles,3,zSize);

function moveBackward(handles,coord)
coordList=['x';'y';'z'];
positionVal=get(handles.positionMsg,'userdata');
stringMsg=get(handles.positionMsg,'string');
if positionVal(coord)-1>0
    positionVal(coord)=positionVal(coord)-1;
    set(handles.positionMsg,'userdata',positionVal);
    stringMsg{coord}=strcat(coordList(coord),' : ',num2str(positionVal(coord)));
    set(handles.positionMsg,'string',stringMsg);
end
plot3Dslice(handles,get(handles.sliceSelection,'userdata'));

function moveForward(handles,coord,coordSize)
coordList=['x';'y';'z'];
positionVal=get(handles.positionMsg,'userdata');
stringMsg=get(handles.positionMsg,'string');
if positionVal(coord)<coordSize
    positionVal(coord)=positionVal(coord)+1;
    set(handles.positionMsg,'userdata',positionVal);
    stringMsg{coord}=strcat(coordList(coord),' : ',num2str(positionVal(coord)));
    set(handles.positionMsg,'string',stringMsg);
end
plot3Dslice(handles,get(handles.sliceSelection,'userdata'));

function frameNumShow_Callback(hObject, eventdata, handles)
retPhase=get(handles.showPhaseMap,'userdata');
[xx yy frame]=size(retPhase);
sliderVal=get(hObject,'value');
sliderVal=ceil(frame*sliderVal);
if sliderVal==0
    sliderVal=1;
end
axes(handles.plot2D);
imagesc(squeeze(retPhase(:,:,sliderVal)),[-3 3]),title(strcat(num2str(sliderVal),' frame')), axis off, axis image

function frameNumShow_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function sliceSelection_SelectionChangeFcn(hObject, eventdata, handles)
switch hObject
    case handles.xySlice
        set(handles.sliceSelection,'userdata',1);
        plot3Dslice(handles,1);
    case handles.xzSlice
        set(handles.sliceSelection,'userdata',2);
        plot3Dslice(handles,2);
    case handles.yzSlice
        set(handles.sliceSelection,'userdata',3);
        plot3Dslice(handles,3);
end

function plot3Dslice(handles,slice)
global Reconimg;global res3;global res4;
ZP4=round(size(Reconimg,1));
ZP5=round(size(Reconimg,3));

n_m=get(handles.n_m,'userdata');
% Reconimg=get(handles.plot3D,'userdata');
axes(handles.plot3D);
positionVal=get(handles.positionMsg,'userdata');
if slice==1
    imagesc(((1:ZP4)-ZP4/2)*res3,((1:ZP4)-ZP4/2)*res3,conv2((squeeze(Reconimg(:,:,positionVal(3)))),fspecial('disk',1),'same'),[n_m-0.01 n_m+0.06]),colorbar,colormap('jet'), axis off,axis image,title('xy slice')
elseif slice==2
    imagesc(((1:ZP4)-ZP4/2)*res3,((1:ZP5)-ZP5/2)*res4,conv2(squeeze(Reconimg(:,positionVal(2),:))',fspecial('disk',1),'same'),[n_m-0.01 n_m+0.06]),colorbar,colormap('jet'), axis off,axis image,title('xz slice')
else
    imagesc(((1:ZP4)-ZP4/2)*res3,((1:ZP5)-ZP5/2)*res4,conv2(squeeze(Reconimg(positionVal(1),:,:))',fspecial('disk',1),'same'),[n_m-0.01 n_m+0.06]),colorbar,colormap('jet'), axis off,axis image,title('yz slice')
end

function n_m_Callback(hObject, eventdata, handles)
n_m=get(handles.n_m,'string');
set(handles.n_m,'userdata',str2num(n_m));

function n_m_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function bgNum_Callback(hObject, eventdata, handles)

function bgNum_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function exportTiff_Callback(hObject, eventdata, handles)
spath=get(handles.savePathInd,'string');
cd(spath)
tomoList=dir('Tomogram*.mat');

set(handles.analyze3DMsg,'string',strcat('Tiff Export Start'));
for tomoNum=1:length(tomoList)
    load(tomoList(tomoNum).name)
    fileName=tomoList(tomoNum).name;
    fileName=fileName(1:end-10);
    Reconimg=uint16(real(Reconimg)*10000);
    outputFileName = strcat(fileName,'.tif');
    for K=1:size(Reconimg,3)
        imwrite(fliplr(Reconimg(:, :, K)), outputFileName, 'WriteMode', 'append');
    end
    set(handles.analyze3DMsg,'string',strcat('Tiff Export',sprintf('%04.2f',100*tomoNum/length(tomoList)),'% Done'));
    drawnow
end
set(handles.analyze3DMsg,'string',strcat('Tiff Export Done'));


function overallPhase_Callback(hObject, eventdata, handles)
global delPhase;global meanPhase;global phaseLength;
figure(10),
plot(meanPhase,'or'),hold on
plot(delPhase,'og'),hold on
ylim([-2 2])
line([0 phaseLength+1],[1 1]*str2double(get(handles.overallPhase,'string')),'Color','red','LineStyle','--')
line([0 phaseLength+1],[1 1]*str2double(get(handles.deltaPhase,'string')),'Color','green','LineStyle','--')
hold off
inspectField_Callback(hObject, eventdata, handles)

function overallPhase_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function deltaPhase_Callback(hObject, eventdata, handles)
global delPhase;global meanPhase;global phaseLength;
figure(10),
plot(meanPhase,'or'),hold on
plot(delPhase,'og'),hold on
ylim([-2 2])
line([0 phaseLength+1],[1 1]*str2double(get(handles.overallPhase,'string')),'Color','red','LineStyle','--')
line([0 phaseLength+1],[1 1]*str2double(get(handles.deltaPhase,'string')),'Color','green','LineStyle','--')
hold off
inspectField_Callback(hObject, eventdata, handles)

function deltaPhase_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function inspectField_Callback(hObject, eventdata, handles)
global excludeFrame;global delPhase;global meanPhase;global phaseLength;
excludeFrame=[];
set(handles.analyze3DMsg,'string','Field Inspection Start');
drawnow

spath=get(handles.savePathInd,'string');
sampleList=get(handles.sampleNum,'userdata');

sampleNum=length(get(handles.analysisOption3D,'userdata'));
if sampleNum==1
    sampleNumList=get(handles.sampleSelect3D,'userdata');
else
    sampleNumList=1:sampleNum;
end

for sampleNumVal=sampleNumList(1):sampleNumList(end)
    cd(spath)
    load(strcat('Field_',sampleList(sampleNumVal).name));
    meanPhase=mean(squeeze(mean(abs(retPhase),1)));
    
    figure(10),
    plot(meanPhase,'or'),hold on
    excludeFrame=[excludeFrame, find(abs(meanPhase)>str2double(get(handles.overallPhase,'string')))];
    delPhase=meanPhase-circshift(meanPhase,1);
    figure(10),
    plot(delPhase,'og'),hold on
    ylim([-2 2])
    excludeFrame=[excludeFrame, find(abs(delPhase)>str2double(get(handles.deltaPhase,'string')))];
    line([0 size(retPhase,3)+1],[1 1]*str2double(get(handles.overallPhase,'string')),'Color','red','LineStyle','--')
    line([0 size(retPhase,3)+1],[1 1]*str2double(get(handles.deltaPhase,'string')),'Color','green','LineStyle','--')
end
hold off
figure(10),
phaseLength=size(retPhase,3);


function mouseMove (hObject, eventdata, handles)
global Mask;global posX;global posY;global radVal;
global clickToken;global signVal;
global imField;global signPhase;

C = get (gca, 'CurrentPoint');
MaskSelected=((posX-C(1,1)).^2+(posY-C(1,2)).^2)<(radVal.^2);
h=imagesc(gca,imField,[-2 2]);

% get(handles.drawBtn,'value')
% h=imagesc(a);
colormap('jet')
axis image;

if clickToken==0
    set(h,'AlphaData',1-(Mask+signVal*MaskSelected*0.5)+0.25);
else
    Mask=(Mask+signVal*MaskSelected)>0;
    set(h,'AlphaData',(1-Mask+0.25));
    tempPhaseMask=imField.*Mask;
    signPhase=mean(tempPhaseMask(Mask>0));
    signPhase=(signPhase>0);
end

function enddraw (hObject, eventdata, handles)
set(gcf,'WindowButtonMotionFcn',[])
set(gcf,'WindowButtonUpFcn',[])
global clickToken;
set (gcf, 'WindowButtonMotionFcn', @mouseMove);
clickToken=0;

function radiusVal_Callback(hObject, eventdata, handles)
global radVal;
radVal=str2double(get(handles.radiusVal,'string'));
%
function radiusVal_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function figure1_WindowButtonDownFcn(hObject, eventdata, handles)
global clickToken;global signVal;global imField;
% C = get (gcf, 'CurrentPoint')
C1 = get(handles.plot2D,'CurrentPoint');
C1=C1(1,1:2);

if (C1(1)>0)&&(C1(1)<size(imField,1))&&(C1(2)>0)&&(C1(2)<size(imField,2))
    switch hObject
        case handles.eraseBtn
            signVal=-1;
        case handles.drawBtn
            signVal=1;
    end
    if clickToken==0
        clickToken=1;
        set(gcf,'WindowButtonMotionFcn',@mouseMove)
        set(gcf,'WindowButtonUpFcn',@enddraw)
    end
end

function figure1_WindowButtonUpFcn(hObject, eventdata, handles)
% set(gcf,'WindowButtonUpFcn',@enddraw)

function maskSelection_SelectionChangedFcn(hObject, eventdata, handles)
global signVal;
switch hObject
    case handles.eraseBtn
        signVal=-1;
    case handles.drawBtn
        signVal=1;
end

function figure1_WindowScrollWheelFcn(hObject, eventdata, handles)
global radVal;
radVal=round(max(0.1,radVal+(eventdata.VerticalScrollCount)*(eventdata.VerticalScrollAmount)/5)*10)/10;
set(handles.radiusVal,'string',num2str(radVal));

function stopMasking_Callback(hObject, eventdata, handles)
set(gcf,'WindowButtonMotionFcn',[])