function spots=FRETspotsearch1centroid(alexacq,frameidx,map,camcal,parameters)
K=genker(max(parameters.sigma));
min_photons_over_bg=4*K.signalnz;
spacing=(size(K.bgmatrix,1)-1)/2+10; %additional spacing of 10 px due to rotation between channels
alexacq=bsxfun(@rdivide,bsxfun(@minus,single(alexacq),camcal.offset),camcal.gain);
%isolate image subsets
acq{3}=alexacq(map.coords{2}(1):map.coords{2}(1)+map.coords{2}(3)-1,map.coords{2}(2):map.coords{2}(2)+map.coords{2}(4)-1,1:2:end);
if size(alexacq,3)>1
    acq{2}=alexacq(map.coords{2}(1):map.coords{2}(1)+map.coords{2}(3)-1,map.coords{2}(2):map.coords{2}(2)+map.coords{2}(4)-1,2:2:end);
end
acq{1}=alexacq(map.coords{1}(1):map.coords{1}(1)+map.coords{1}(3)-1,map.coords{1}(2):map.coords{1}(2)+map.coords{1}(4)-1,1:2:end);
%generate summed FRET image
s=size(acq{1});
acq{1}=acq{3}+interp2(acq{1}',map.xmeshrl,map.ymeshrl,'linear',0);
%convolution and peak finding
for k=1:3
    if isempty(acq{k})
        continue
    end
    bgacq=conv2(acq{k},K.bgmatrix,'same')/K.bgnz;
    varacq=conv2(acq{k}.^2,K.bgmatrix,'same')/K.bgnz-bgacq.^2;
    acq{k}=conv2(acq{k},K.signalmatrix,'same')-bgacq*K.signalnz;
    if k<3
        
        
        BW=(acq{k}>=max(min_photons_over_bg,4*K.signalnz*sqrt(1/K.signalnz+1/K.bgnz)*sqrt(varacq)));
        st=regionprops(BW,acq{k},{'WeightedCentroid','MaxIntensity'});
        ds=double(struct2dataset(st));
        [~,I]=sort(ds);
        ds=ds(I,:);
        sizeds=size(ds);
        
     
%         idxmatrix=(acq{k}==imdilate(acq{k},K.signalmatrix))&(acq{k}>=max(min_photons_over_bg,4*K.signalnz*sqrt(1/K.signalnz+1/K.bgnz)*sqrt(varacq)));
%         idx=find(idxmatrix);
%         x=rem(idx-1,s(1))+1;
%         y=(idx-x)/s(1)+1;
%         idx=x>spacing & x<=s(1)-spacing & y>spacing & y<=s(2)-spacing;
%         spots(k).coords=[x(idx),y(idx)];
%         spots(k).frame=frameidx*ones(sum(idx),1);
        spots(k).coords=[ds(1:2,:)];
        spots(k).frame=frameidx*ones(sizeds(1),1);
    end
    
        spots(k).photons=ds(:,end);
%     spots(k).photons=acq{k}(idxmatrix);
%     spots(k).photons=spots(k).photons(idx);
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