function varargout=FRETscript(mode,str)
map=importdata('map.mat');
camcal=importdata('camcal.mat');
filelist=dir([str,'*.tif']);

parameters.sigma=[0.9,1.1];
parameters.pixelsize=117;
parameters.uncert_threshold=[20,20];
parameters.LLR_threshold=[Inf,Inf];
parameters.photon_threshold=[100,100];
parameters.ming_sigma=[0.5,0.5];
parameters.max_sigma=[1.5,1.7];
parameters.ecc_threshold=[0.5,0.5];
parameters.maxproj=false;
parameters.nchannels=2;
parameters.refchannel=2;
parameters.colocradius=2;
parameters.nnradius=10;
parameters.st=[0,30];
parameters.nocuda=false;
parameters.baseline=400;
parameters.counts_per_photon=2.17;

n=10000000;
auxdata.spots=struct('frame',{zeros(n,1,'uint32'),zeros(n,1,'uint32'),[]},'coords',{zeros(n,2,'single'),zeros(n,2,'single'),[]},'photons',{zeros(n,1,'single'),zeros(n,1,'single'),zeros(n,1,'single')},'filter_idx',{false(n,1),false(n,1),false(n,1)},'idx',{0,0,[]});
auxdata.colocmatrix=struct('photons',zeros(n,3,'single'),'coords',zeros(n,4,'single'),'filter_idx',false(n,1),'idx',0,'signal_photons',0,'FRET_photons',0);
auxdata.FRETidx=0;
for k=1:numel(filelist)
    fprintf(1,'%u out of %u\n',k,numel(filelist));
    filename=filelist(k).name;
    info=imfinfo(filename);
    nframes=length(info);
    firstimg=imread(filename,'info',info,'index',1)';
    switch mode
        case 'EStrace'
            if exist([filename,'.temp.spots.mat'],'file')
                photons(k)=importdata([filename,'.temp.spots.mat']);
                continue
            end
            acq={zeros([size(firstimg),nframes/2],'single'),zeros([size(firstimg),nframes/2],'single')};
            for kk=1:nframes
                if mod(kk,2)
                    acq{1}(:,:,(kk+1)/2)=single(imread(filename,'info',info,'index',kk))';
                else
                    acq{2}(:,:,kk/2)=single(imread(filename,'info',info,'index',kk))';
                end
            end
            %calculate maximum projection
            rtFRET(cat(3,max(acq{1},[],3),max(acq{2},[],3)),1);
            
            K=genker(max(parameters.sigma));
            %acq{1}=bsxfun(@rdivide,bsxfun(@minus,acq{1},camcal.offset),camcal.gain);
            acq{1}=(acq{1}-400)/2.17;
            acq{2}=(acq{2}-400)/2.17;
            %acq{2}=bsxfun(@rdivide,bsxfun(@minus,acq{2},camcal.offset),camcal.gain);
            acq{3}=acq{1}(map.subcoords{2}(1):map.subcoords{2}(1)+map.subcoords{2}(3)-1,map.subcoords{2}(2):map.subcoords{2}(2)+map.subcoords{2}(4)-1,:);
            acq{2}=acq{2}(map.subcoords{2}(1):map.subcoords{2}(1)+map.subcoords{2}(3)-1,map.subcoords{2}(2):map.subcoords{2}(2)+map.subcoords{2}(4)-1,:);
            acq{1}=acq{1}(map.subcoords{1}(1):map.subcoords{1}(1)+map.subcoords{1}(3)-1,map.subcoords{1}(2):map.subcoords{1}(2)+map.subcoords{1}(4)-1,:);
            for kk=1:size(acq{1},3)
                acq{1}(:,:,kk)=acq{3}(:,:,kk)+interp2q(acq{1}(:,:,kk),map.xmeshrl,map.ymeshrl);
            end
            for kk=1:3
                bgacq{kk}=convn(acq{kk},K.bgmatrix,'same')/K.bgnz;
                acq{kk}=convn(acq{kk},K.signalmatrix,'same')-bgacq{kk}*K.signalnz;
            end
            photons(k).coords=auxdata.colocmatrix.coords(1:auxdata.colocmatrix.idx,:);
            coords=auxdata.colocmatrix.coords(1:auxdata.colocmatrix.idx,[3,4]);
            xcoords_AA=repmat(round(coords(:,1)'),[size(acq{2},3),1]);
            ycoords_AA=repmat(round(coords(:,2)'),[size(acq{2},3),1]);
            coords=auxdata.colocmatrix.coords(1:auxdata.colocmatrix.idx,[1,2]);
            xcoords_DD_DA=repmat(round(coords(:,1)'),[size(acq{2},3),1]);
            ycoords_DD_DA=repmat(round(coords(:,2)'),[size(acq{2},3),1]);
            frame=1:size(acq{1},3);
            frame=repmat(frame',[1,size(coords,1)]);
            photons(k).DD_DA=acq{1}(sub2ind(size(acq{1}),double(xcoords_DD_DA),double(ycoords_DD_DA),double(frame))); %DD+DA photons
            photons(k).AA=acq{2}(sub2ind(size(acq{2}),double(xcoords_AA),double(ycoords_AA),double(frame))); %AA photons
            photons(k).DA=acq{3}(sub2ind(size(acq{3}),double(xcoords_AA),double(ycoords_AA),double(frame))); %DA photons
            temp=photons(k);
            save([filename,'.temp.spots.mat'],'temp')
            auxdata.spots(1).idx=0;
            auxdata.spots(2).idx=0;
            auxdata.colocmatrix.idx=0;
        case 'ES'
            acq=zeros([size(firstimg),2],'single');
            for kk=1:nframes/2
                acq(:,:,1)=single(imread(filename,'info',info,'index',2*kk-1))';
                acq(:,:,2)=single(imread(filename,'info',info,'index',2*kk))';
                auxdata.FRETidx=auxdata.FRETidx+1;
                rtFRET(acq,auxdata.FRETidx);
            end
        case 'E'
            for kk=1:nframes
                acq=single(imread(filename,'info',info,'index',kk))';
                auxdata.FRETidx=auxdata.FRETidx+1;
                rtFRET(acq,auxdata.FRETidx);
            end
    end
end

switch mode
    case 'EStrace'
        varargout{1}=photons;
        reviewFRETtrace(str,[photons.DD_DA],[photons.DA],[photons.AA],auxdata.colocmatrix.signal_photons,auxdata.colocmatrix.FRET_photons);
        return
end

for kl=1:2
    auxdata.spots(kl).photons=auxdata.spots(kl).photons(1:auxdata.spots(kl).idx);
    auxdata.spots(kl).coords=auxdata.spots(kl).coords(1:auxdata.spots(kl).idx,:);
    auxdata.spots(kl).frame=auxdata.spots(kl).frame(1:auxdata.spots(kl).idx);
    auxdata.spots(kl).filter_idx=auxdata.spots(kl).filter_idx(1:auxdata.spots(kl).idx);
end
switch mode
    case 'ES'
        auxdata.spots(3).photons=auxdata.spots(3).photons(1:auxdata.spots(2).idx);
        auxdata.spots(3).filter_idx=auxdata.spots(3).filter_idx(1:auxdata.spots(2).idx);
    case 'E'
        auxdata.spots(3).photons=auxdata.spots(3).photons(1:auxdata.spots(1).idx);
        auxdata.spots(3).filter_idx=auxdata.spots(3).filter_idx(1:auxdata.spots(1).idx);
end
auxdata.colocmatrix.photons=auxdata.colocmatrix.photons(1:auxdata.colocmatrix.idx,:);
auxdata.colocmatrix.coords=auxdata.colocmatrix.coords(1:auxdata.colocmatrix.idx,:);
auxdata.colocmatrix.filter_idx=auxdata.colocmatrix.filter_idx(1:auxdata.colocmatrix.idx,:);
auxdata.colocmatrix.mode=mode;
reviewFRET(str,pwd,auxdata.spots,auxdata.colocmatrix,parameters);


    function rtFRET(acq,acqidx)
        parameters.colocradius=1.5;
        parameters.nnradius=10;
        
        spots=FRETspotsearch1(acq,acqidx,map,camcal,parameters);
        %filter new data with current threshold
        
        if size(acq,3)==1
            parameters.nchannels=1;
            spots(1).filter_idx=true(size(spots(1).photons));%spots(1).photons>auxdata.colocmatrix.signal_photons;
            spots=spotfun1('nnfilter',spots,parameters);
            spots(3).filter_idx=true(size(spots(3).photons));%spots(3).photons>auxdata.colocmatrix.FRET_photons;
            idx1=auxdata.spots(1).idx;
            n1=numel(spots(1).frame);
            auxdata.spots(1).frame(idx1+1:idx1+n1)=spots(1).frame;
            auxdata.spots(1).coords(idx1+1:idx1+n1,:)=spots(1).coords;
            auxdata.spots(1).photons(idx1+1:idx1+n1)=spots(1).photons;
            auxdata.spots(3).photons(idx1+1:idx1+n1)=spots(3).photons;
            auxdata.spots(1).filter_idx(idx1+1:idx1+n1)=spots(1).filter_idx;
            auxdata.spots(3).filter_idx(idx1+1:idx1+n1)=spots(3).filter_idx;
            auxdata.spots(1).idx=idx1+n1;
            return
        end
        
        spots(1).filter_idx=true(size(spots(1).photons));%spots(1).photons>auxdata.colocmatrix.signal_photons;
        spots(2).filter_idx=true(size(spots(2).photons));%spots(2).photons>auxdata.colocmatrix.signal_photons;
        spots=spotfun1('nnfilter',spots,parameters);
        colocmatrix=colocfun('spots',spots,parameters);
        colocmatrix=colocmatrix{1,2};
        colocmatrix.photons(:,3)=spots(3).photons(colocmatrix.spotidx(:,2));
        colocmatrix.filter_idx=true(size(colocmatrix.photons(:,3)));%colocmatrix.photons(:,3)>auxdata.colocmatrix.FRET_photons;
        
        %merge information
        idx1=auxdata.spots(1).idx;
        idx2=auxdata.spots(2).idx;
        n1=numel(spots(1).frame);
        n2=numel(spots(2).frame);
        auxdata.spots(1).frame(idx1+1:idx1+n1)=spots(1).frame;
        auxdata.spots(2).frame(idx2+1:idx2+n2)=spots(2).frame;
        auxdata.spots(1).coords(idx1+1:idx1+n1,:)=spots(1).coords;
        auxdata.spots(2).coords(idx2+1:idx2+n2,:)=spots(2).coords;
        auxdata.spots(1).photons(idx1+1:idx1+n1)=spots(1).photons;
        auxdata.spots(2).photons(idx2+1:idx2+n2)=spots(2).photons;
        auxdata.spots(3).photons(idx2+1:idx2+n2)=spots(3).photons;
        auxdata.spots(1).filter_idx(idx1+1:idx1+n1)=spots(1).filter_idx;
        auxdata.spots(2).filter_idx(idx2+1:idx2+n2)=spots(2).filter_idx;
        auxdata.spots(1).idx=idx1+n1;
        auxdata.spots(2).idx=idx2+n2;
        idx=auxdata.colocmatrix.idx;
        n=size(colocmatrix.coords,1);
        auxdata.colocmatrix.photons(idx+1:idx+n,:)=colocmatrix.photons;
        auxdata.colocmatrix.coords(idx+1:idx+n,:)=colocmatrix.coords;
        auxdata.colocmatrix.photons(idx+1:idx+n,:)=colocmatrix.photons;
        auxdata.colocmatrix.filter_idx(idx+1:idx+n,:)=colocmatrix.filter_idx;
        auxdata.colocmatrix.idx=idx+n;
    end

end

function photons=getphotons(acq,camcal,parameters)
K=genker(max(parameters.sigma));
min_photons_over_bg=4*K.signalnz;
spacing=(size(K.bgmatrix,1)-1)/2;
alexacq=bsxfun(@rdivide,bsxfun(@minus,single(acq),camcal.offset),camcal.gain);
%isolate image subsets
acq{3}=alexacq(map.subcoords{2}(1):map.subcoords{2}(1)+map.subcoords{2}(3)-1,map.subcoords{2}(2):map.subcoords{2}(2)+map.subcoords{2}(4)-1,1:2:end);
if size(alexacq,3)>1
    acq{2}=alexacq(map.subcoords{2}(1):map.subcoords{2}(1)+map.subcoords{2}(3)-1,map.subcoords{2}(2):map.subcoords{2}(2)+map.subcoords{2}(4)-1,2:2:end);
end
acq{1}=alexacq(map.subcoords{1}(1):map.subcoords{1}(1)+map.subcoords{1}(3)-1,map.subcoords{1}(2):map.subcoords{1}(2)+map.subcoords{1}(4)-1,1:2:end);
%generate summed FRET image
s=size(acq{1});
acq{1}=acq{3}+interp2q(acq{1},map.xmeshrl,map.ymeshrl);
%convolution and peak finding
for k=1:3
    if isempty(acq{k})
        continue
    end
    bgacq=conv2(acq{k},K.bgmatrix,'same')/K.bgnz;
    varacq=conv2(acq{k}.^2,K.bgmatrix,'same')/K.bgnz-bgacq.^2;
    acq{k}=conv2(acq{k},K.signalmatrix,'same')-bgacq*K.signalnz;
    if k<3
        if k==1
            acq{k}(~isfinite(acq{k}))=-Inf;
        end
        
        BW=(acq{k}>=max(min_photons_over_bg,4*K.signalnz*sqrt(1/K.signalnz+1/K.bgnz)*sqrt(varacq)));
        st=regionprops(BW,acq{k},{'WeightedCentroid','MaxIntensity'});
        ds=double(struct2dataset(st));
        [~,I]=sort(ds);
        ds=ds(I,:);
        sizeds=size(ds);
        


        spots(k).coords=ds(1:2,:);
        spots(k).frame=frameidx*ones(sizeds(1),1);
    end
    
        spots(k).photons=ds(:,end);
end
end

function acq=interp2q(acq,x,y)
s=size(acq);
x_floor=floor(x);
y_floor=floor(y);
x=x-x_floor;
y=y-y_floor;
oor=x_floor<1 | y_floor<1 | x_floor+1>s(1) | y_floor+1>s(2) | ~isfinite(x_floor+y_floor);
x_floor(oor)=1;
y_floor(oor)=1;
%convert to linear index
idx=x_floor+(y_floor-1)*s(1);
acq=acq(idx).*(1-x).*(1-y)+acq(idx+1).*x.*(1-y)+acq(idx+s(1)).*(1-x).*y+acq(idx+s(1)+1).*x.*y;
acq(oor)=NaN;
end