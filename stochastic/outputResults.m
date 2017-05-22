function [f] = outputResults(unitCell,extrudedUnitCell,result,timesteps,opt)

% options for saving simulation data
date = opt.date;
time = opt.time;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PLOT INFO FOR IMPLEMENTATION OF NEW GEOMETRIES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
col=hsv(32);
colt=col([12,26,31,20],:);
colt(1:2,:)=colt(1:2,:)+(1-colt(1:2,:))*0.3;
colt(3:4,:)=colt(3:4,:)+(1-colt(3:4,:))*0.1;
colt=[colt; col(2,:)+(1-col(2,:))*0.85];

colt(4,:)=[0,153,255]/255;
colt(5,:)=[225,225,225]/255;

%col=hsv(32);
%col=col+(1-col)*0;
%colt=col([32,31,8,13],:)
%colt=[colt; col(2,:)+(1-col(2,:))*0.85];
    
viewCoor=[sind(opt.AZ) -cosd(opt.AZ) sind(opt.EL)];
opt.tranPol=0.5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PREPARE PLOTTING UNDEFORMED CONFIGURATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Check if output folder is required, and create it if it doesn't exist
nameFolder=[pwd,'/Results/',date];
if or(strcmp(opt.saveFig,'on'),strcmp(opt.saveMovie,'on'))
    if exist(nameFolder)==0
        mkdir(nameFolder)
    end
end

nref=size(unitCell.l,1);
if nref==0
    extrudedUnitCell.ref=[];
end

%prepare lattice vectors for plotting periodic structures and lattice
for ne=1:length(unitCell.Polyhedron)
    unitCell.Polyhedron(ne).latVec=detLatVec(unitCell.lhat,opt);
    unitCell.PolyhedronNew(ne)=unitCell.Polyhedron(ne);
end
extrudedUnitCell.latVec=detLatVec(unitCell.l,opt); 
%Prepare polyhdra and solid faces for efficient plotting
for ne=1:length(unitCell.Polyhedron)
    plotunitCell.Init(ne)=prepEffPlot(unitCell.Polyhedron(ne),viewCoor);
    plotunitCell.InitNew(ne)=prepEffPlot(unitCell.PolyhedronNew(ne),viewCoor);
end
plotextrudedUnitCell=prepEffPlot(extrudedUnitCell,viewCoor);

%Plot solid face with 100% transparency
f=figure('Position', [0 0 800 800]); hold on
for nc=1:size(extrudedUnitCell.latVec,1)
    for i=3:10
        c=(plotextrudedUnitCell.polFace(i).normal*viewCoor')>0;
        hs{nc,i}=patch('Faces',plotextrudedUnitCell.polFace(i).nod,'Vertices',plotextrudedUnitCell.lat(nc).coor,'facecolor','flat','facevertexCData',c*colt(4,:)+abs(1-c)*colt(5,:),'facealpha',1.0,'edgealpha',1.0);
    end
end

%Plot polyhedra with 100% transparency
for ne=1:length(unitCell.Polyhedron)
    for nc=1:size(unitCell.Polyhedron(1).latVec,1)
        for i=3:10
            [c, plotunitCell.Init(ne).polFace(i).indexe, b1]=intersect(plotunitCell.Init(ne).polFace(i).index,unitCell.Polyhedron(ne).extrude);
            [c, plotunitCell.Init(ne).polFace(i).indexs, b1]=intersect(plotunitCell.Init(ne).polFace(i).index,unitCell.Polyhedron(ne).solidify);
            hie{ne,nc,i}=patch('Faces',plotunitCell.Init(ne).polFace(i).nod(plotunitCell.Init(ne).polFace(i).indexe,:),'Vertices',plotunitCell.Init(ne).lat(nc).coor,'facecolor','flat','facevertexCData',interp1([1,max(2,length(unitCell.Polyhedron))],colt([1,2],:),(ne)),'facealpha',0.0,'edgealpha',0.0);
            his{ne,nc,i}=patch('Faces',plotunitCell.Init(ne).polFace(i).nod(plotunitCell.Init(ne).polFace(i).indexs,:),'Vertices',plotunitCell.Init(ne).lat(nc).coor,'facecolor','flat','facevertexCData',interp1([1,max(2,length(unitCell.Polyhedron))],colt([1,2],:),(ne)),'facealpha',0.0,'edgealpha',0.0);
        end
    end
end
%Set axis
axis tight
xlim=1.1*get(gca,'xlim');
ylim=1.1*get(gca,'ylim');
zlim=1.1*get(gca,'zlim');
set(gca,'xlim',xlim,'ylim',ylim,'zlim',zlim);
hl2=plotOpt(opt);

opt.xlim=xlim;
opt.ylim=ylim;
opt.zlim=zlim;
if strcmp(opt.plot,'result')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %PREPARE PLOTTING OF DEFORMED CONFIGURATION
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if result.numMode>0
        %prepare reference vectors for plotting periodic structures  
        refVec=detrefVec(unitCell,result,opt,extrudedUnitCell.ref);

        %Determine node coordinates and orientation of faces for each step in
        %the mode
        maxAxis=[0 0 0];
        minAxis=[0 0 0];
        Imax=opt.interval+1;
        Imin=-1;
        for nMode=1:result.numMode
            for nc=1:size(extrudedUnitCell.latVec,1)
                plotextrudedUnitCell.mode(nMode).lat(nc).coor=extrudedUnitCell.node+result.deform(nMode).V+ones(size(extrudedUnitCell.node,1),1)*(extrudedUnitCell.latVec(nc,:)-refVec(nMode).val(nc,:));
                ma=max(plotextrudedUnitCell.mode(nMode).lat(nc).coor);
                mi=min(plotextrudedUnitCell.mode(nMode).lat(nc).coor);
                maxAxis(maxAxis<ma)=ma(maxAxis<ma);
                minAxis(minAxis>mi)=mi(minAxis>mi);
            end
            for i=3:10
                plotextrudedUnitCell.mode(nMode).polFace(i).normal(1,:)=[0 0 0];
                for j=1:size(plotextrudedUnitCell.polFace(i).nod,1)
                    n1a=plotextrudedUnitCell.mode(nMode).lat(1).coor(plotextrudedUnitCell.polFace(i).nod(j,2),1:3)-plotextrudedUnitCell.mode(nMode).lat(1).coor(plotextrudedUnitCell.polFace(i).nod(j,1),1:3);
                    n2a=plotextrudedUnitCell.mode(nMode).lat(1).coor(plotextrudedUnitCell.polFace(i).nod(j,3),1:3)-plotextrudedUnitCell.mode(nMode).lat(1).coor(plotextrudedUnitCell.polFace(i).nod(j,1),1:3);
                    n3a=cross(n1a,n2a);
                    plotextrudedUnitCell.mode(nMode).polFace(i).normal(j,:)=n3a;
                end
            end
        end
        %Update axis
        xlim=[minAxis(1),maxAxis(1)];
        ylim=[minAxis(2),maxAxis(2)];
        zlim=[minAxis(3),maxAxis(3)];
        set(gca,'xlim',xlim,'ylim',ylim,'zlim',zlim);
    end
end

if strcmp(opt.plot,'info')
    %First make solid face with 100% transparency
    for nc=1:size(extrudedUnitCell.latVec,1)
        for i=3:10
            set(hs{nc,i},'facealpha',0,'edgealpha',0);
        end
    end
    
    %PLOT POLYHEDRA TEMPLATE
    [f1,hs,hie,his] = copyFigure(unitCell,extrudedUnitCell,opt,hs,hie,his);
    for nc=1:size(unitCell.Polyhedron(1).latVec,1)
        for ne=1:length(unitCell.Polyhedron)
            for i=3:10                  
                set(hie{ne,nc,i},'facealpha',opt.tranPol,'edgealpha',opt.tranPol);
                set(his{ne,nc,i},'facealpha',opt.tranPol,'edgealpha',opt.tranPol);
            end               
        end
    end
    close(f)
    %PLOT ALL INFORMATION
    f=figure('Position', [0 0 800 800]);hold on;
    for ne=1:length(unitCell.Polyhedron)
        for i=3:10
            patch('Faces',plotunitCell.Init(ne).polFace(i).nod(plotunitCell.Init(ne).polFace(i).indexe,:),'Vertices',plotunitCell.Init(ne).lat(1).coor,'facecolor','flat','facevertexCData',interp1([1,max(2,length(unitCell.Polyhedron))],colt([1,2],:),(ne)),'facealpha',opt.tranPol,'edgealpha',opt.tranPol);
            patch('Faces',plotunitCell.Init(ne).polFace(i).nod(plotunitCell.Init(ne).polFace(i).indexs,:),'Vertices',plotunitCell.Init(ne).lat(1).coor,'facecolor','flat','facevertexCData',interp1([1,max(2,length(unitCell.Polyhedron))],colt([1,2],:),(ne)),'facealpha',opt.tranPol,'edgealpha',opt.tranPol);
        end
    end
     
    for ne=1:length(unitCell.Polyhedron)
        coorCenter=sum(unitCell.Polyhedron(ne).node)/size(unitCell.Polyhedron(ne).node,1);
        text(coorCenter(1),coorCenter(2),coorCenter(3),num2str(ne),'fontsize',30,'color','g')
        for i=1:size(unitCell.Polyhedron(ne).node,1)
            coor=unitCell.Polyhedron(ne).node(i,:);
            coorText=[coor(1)*0.85+coorCenter(1)*0.15,coor(2)*0.85+coorCenter(2)*0.15,coor(3)*0.85+coorCenter(3)*0.15];
            line([coor(1),coorText(1)],[coor(2),coorText(2)],[coor(3),coorText(3)],'color','k','linestyle',':')
            plot3(coor(1),coor(2),coor(3),'*k')
            text(coorText(1),coorText(2),coorText(3),num2str(i),'fontsize',20)
        end
        for i=1:size(unitCell.Polyhedron(ne).face,1)
            coor=sum(unitCell.Polyhedron(ne).node(unitCell.Polyhedron(ne).face{i},:))/length(unitCell.Polyhedron(ne).face{i});
            coorText=[coor(1)*0.85+coorCenter(1)*0.15,coor(2)*0.85+coorCenter(2)*0.15,coor(3)*0.85+coorCenter(3)*0.15];
            line([coor(1),coorText(1)],[coor(2),coorText(2)],[coor(3),coorText(3)],'color','k','linestyle',':')
            plot3(coor(1),coor(2),coor(3),'*b')
            text(coorText(1),coorText(2),coorText(3),num2str(i),'fontsize',20,'color','b')
        end
    end 
    hl2=plotOpt(opt); 
    axis tight
    %EXTRUDED STATE
    [f4,hs,hie,his] = copyFigure(unitCell,extrudedUnitCell,opt,hs,hie,his);   
    for nc=1:size(extrudedUnitCell.latVec,1)
        for i=3:10
            set(hs{nc,i},'facealpha',1,'edgealpha',1)         
        end
    end
    for ne=1:length(unitCell.Polyhedron)
        for nc=1:size(unitCell.Polyhedron(1).latVec,1)
            for i=3:10
                set(hie{ne,nc,i},'facealpha',0,'edgealpha',0);
                set(his{ne,nc,i},'facealpha',0,'edgealpha',0);
            end
        end
    end
    set(gca,'xlim',xlim,'ylim',ylim,'zlim',zlim);
    %EXTRUDED STATE 'INFO
    [f5,hs,hie,his] = copyFigure(unitCell,extrudedUnitCell,opt,hs,hie,his); hold on 
    for nc=1:size(extrudedUnitCell.latVec,1)
        for i=3:10
            set(hs{nc,i},'facealpha',0.75,'edgealpha',.75)         
        end
    end
    coorCenter=sum(extrudedUnitCell.node)/size(extrudedUnitCell.node,1);    
    for i=1:size(extrudedUnitCell.node,1)
            coor=extrudedUnitCell.node(i,:);
            coorText=[coor(1)*0.85+coorCenter(1)*0.15,coor(2)*0.85+coorCenter(2)*0.15,coor(3)*0.85+coorCenter(3)*0.15];
            line([coor(1),coorText(1)],[coor(2),coorText(2)],[coor(3),coorText(3)],'color','k','linestyle',':')
            plot3(coor(1),coor(2),coor(3),'*k')
            text(coorText(1),coorText(2),coorText(3),num2str(i),'fontsize',20)
        end
    set(gca,'xlim',xlim,'ylim',ylim,'zlim',zlim);    
    %EXTRUDED STATE 'INFO' ANGLE NUMBERS
    [f6,hs,hie,his] = copyFigure(unitCell,extrudedUnitCell,opt,hs,hie,his); hold on 
    for nc=1:size(extrudedUnitCell.latVec,1)
        for i=3:10
            set(hs{nc,i},'facealpha',0.75,'edgealpha',.75)         
        end
    end
    coorCenter=sum(extrudedUnitCell.node)/size(extrudedUnitCell.node,1);    
    for i=1:size(extrudedUnitCell.nodeHingeEx,1)
            coor1=extrudedUnitCell.node(extrudedUnitCell.nodeHingeEx(i,1),:);
            coor2=extrudedUnitCell.node(extrudedUnitCell.nodeHingeEx(i,2),:);
            coor=coor1/2+coor2/2;
            coorText=[coor(1)*0.85+coorCenter(1)*0.15,coor(2)*0.85+coorCenter(2)*0.15,coor(3)*0.85+coorCenter(3)*0.15];
            line([coor(1),coorText(1)],[coor(2),coorText(2)],[coor(3),coorText(3)],'color','k','linestyle',':')
            plot3(coor(1),coor(2),coor(3),'*k')
            text(coorText(1),coorText(2),coorText(3),num2str(i),'fontsize',20)
        end
    set(gca,'xlim',xlim,'ylim',ylim,'zlim',zlim);
end

if strcmp(opt.plot,'result')
        %First make solid face with 100% transparency
        for nc=1:size(extrudedUnitCell.latVec,1)
            for i=3:10
                set(hs{nc,i},'facealpha',0,'edgealpha',0);
            end
        end
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %PLOT SPACE-FILLING ASSEMBLY POLYHEDRA
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for nc=1:size(unitCell.Polyhedron(1).latVec,1)
            for ne=1:length(unitCell.Polyhedron)
                for i=3:10
                    set(hie{ne,nc,i},'facealpha',opt.tranPol,'edgealpha',opt.tranPol);
                    set(his{ne,nc,i},'facealpha',opt.tranPol,'edgealpha',opt.tranPol);
                end
                set(gca,'xlim',xlim,'ylim',ylim,'zlim',zlim); 
            end
        end
        axis equal
        printHigRes(f,opt,'polyhedra',nameFolder)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %PLOT MOVING OF POLYHEDRA
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %New coordinates
        for ne=1:length(unitCell.Polyhedron)
            for nc=1:size(unitCell.Polyhedron(1).latVec,1)
                plotunitCell.Init(ne).lat(nc).coorNew(:,1)=plotunitCell.InitNew(ne).lat(nc).coor(:,1)+extrudedUnitCell.latVec(nc,1)-unitCell.PolyhedronNew(ne).latVec(nc,1);
                plotunitCell.Init(ne).lat(nc).coorNew(:,2)=plotunitCell.InitNew(ne).lat(nc).coor(:,2)+extrudedUnitCell.latVec(nc,2)-unitCell.PolyhedronNew(ne).latVec(nc,2);
                plotunitCell.Init(ne).lat(nc).coorNew(:,3)=plotunitCell.InitNew(ne).lat(nc).coor(:,3)+extrudedUnitCell.latVec(nc,3)-unitCell.PolyhedronNew(ne).latVec(nc,3);
            end
        end
        %Update position
        for ne=1:length(unitCell.Polyhedron)
            for nc=1:size(unitCell.Polyhedron(1).latVec,1)
                for i=3:10
                    set(hie{ne,nc,i},'vertices',plotunitCell.Init(ne).lat(nc).coorNew);
                    set(his{ne,nc,i},'vertices',plotunitCell.Init(ne).lat(nc).coorNew);
                end 
            end
        end
        set(gca,'xlim',xlim,'ylim',ylim,'zlim',zlim);
        %printHigRes(f,opt,'Polyhedra_Packing_Expanded',nameFolder)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %PLOT SELECTED FACES TO EXTRUDE, SOLIDIFY AND REMOVE
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for nc=1:size(unitCell.Polyhedron(1).latVec,1)
            for ne=1:length(unitCell.Polyhedron)
                for i=3:10
                    if length(plotunitCell.Init(1).polFace(i).indexs)~=0
                        c=(plotunitCell.Init(ne).polFace(i).normal(plotunitCell.Init(ne).polFace(i).indexs,:)*viewCoor')>0;
                        set(his{ne,nc,i},'facealpha',1,'edgealpha',1,'facevertexCData',(c*colt(4,:)+abs(1-c)*colt(5,:)));
                    end
                end
            end
        end
        set(gca,'xlim',xlim,'ylim',ylim,'zlim',zlim);  
        printHigRes(f,opt,'facetype',nameFolder)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %PLOT EXTRUSION
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for nc=1:size(extrudedUnitCell.latVec,1)
            for i=3:10
                set(hs{nc,i},'facealpha',1,'edgealpha',1)      
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %REMOVE POLYHEDRA
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for ne=1:length(unitCell.Polyhedron)
            for nc=1:size(unitCell.Polyhedron(1).latVec,1)
                for i=3:10
                    set(hie{ne,nc,i},'facealpha',0,'edgealpha',0);
                    set(his{ne,nc,i},'facealpha',0,'edgealpha',0);
                end
            end
        end
        set(gca,'xlim',xlim,'ylim',ylim,'zlim',zlim);
        printHigRes(f,opt,'undeformed',nameFolder)     

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %PLOT MODES INDIVIDUALLY
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ant = [];
        for nMode=1:result.numMode
            ant = addText(ant, nMode, timesteps);
            for nc=1:size(extrudedUnitCell.latVec,1)
                for i=3:10     
                    viewCoor=get(gca,'view');
                    viewCoor=[sind(viewCoor(1)) -cosd(viewCoor(1)) sind(viewCoor(2))];
                    c=(plotextrudedUnitCell.mode(nMode).polFace(i).normal*viewCoor'>0);
                    
                    set(hs{nc,i},'Vertices',plotextrudedUnitCell.mode(nMode).lat(nc).coor,'facecolor','flat','facevertexCData',c*colt(4,:)+abs(1-c)*colt(5,:),'facealpha',1.0);
                end
            end

            printGif(opt,nMode,f,nameFolder, time);
            %printHigRes(f,opt,['deformed_'],nameFolder)  
            
            drawnow limitrate % faster updates
        end
        drawnow
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%USED FUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function latVec=detLatVec(perVec,opt)
    latVec=[];
    nref=size(perVec,1);
    switch nref
        case 0
            latVec=[0 0 0];                 
        case 1
            for i=1:opt.plotPer
                latVec=[latVec; (i-1)*perVec(1,:)];
            end
        case 2
            for j=1:opt.plotPer
                for i=1:opt.plotPer
                    latVec=[latVec; (i-1)*perVec(1,:)+(j-1)*perVec(2,:)];
                end
            end
        case 3
            for k=1:opt.plotPer
                for j=1:opt.plotPer
                    for i=1:opt.plotPer
                        latVec=[latVec; (i-1)*perVec(1,:)+(j-1)*perVec(2,:)+(k-1)*perVec(3,:)];
                    end
                end
            end     
    end  
    
function refVec=detrefVec(unitCell,result,opt,ref)
    for nMode=1:length(result.deform)
        refVec(nMode).val=[];
        switch size(unitCell.l,1);
            case 0
                refVec(nMode).val=[0 0 0];                    
            case 1
                for i=1:opt.plotPer
                    refVec(nMode).val=[refVec(nMode).val; (i-1)*result.deform(nMode).V(ref(1),:)];
                end
            case 2
                for j=1:opt.plotPer
                    for i=1:opt.plotPer
                        refVec(nMode).val=[refVec(nMode).val; (i-1)*result.deform(nMode).V(ref(1),:)+(j-1)*result.deform(nMode).V(ref(2),:)];
                    end
                end
            case 3
                for k=1:opt.plotPer
                    for j=1:opt.plotPer
                        for i=1:opt.plotPer
                            refVec(nMode).val=[refVec(nMode).val; (i-1)*result.deform(nMode).V(ref(1),:)+(j-1)*result.deform(nMode).V(ref(2),:)+(k-1)*result.deform(nMode).V(ref(3),:)];
                        end
                    end
                end
        end            
    end
    
function plotg=prepEffPlot(som,viewCoor)
    for i=3:10
        plotg.polFace(i).nod=[];
        plotg.polFace(i).index=[];
    end
    for i=1:length(som.face)
        l=length(som.face{i});
        plotg.polFace(l).nod=[plotg.polFace(l).nod;som.face{i}];
        plotg.polFace(l).index=[plotg.polFace(l).index; i];
    end

    %COORDINATES AND FACE NORMAL POLYHEDRON
    for nc=1:size(som.latVec,1)
        plotg.lat(nc).coor=som.node+ones(size(som.node,1),1)*(som.latVec(nc,:));
        if isfield(som,'nodeNew')
            plotg.lat(nc).coorNew=som.nodeNew+ones(size(som.nodeNew,1),1)*(som.latVec(nc,:));
        end
    end
    for i=3:10
        plotg.polFace(i).normal(1,:)=[0,0,0];
        for j=1:size(plotg.polFace(i).nod,1)
            n1a=plotg.lat(nc).coor(plotg.polFace(i).nod(j,2),1:3)-plotg.lat(nc).coor(plotg.polFace(i).nod(j,1),1:3);
            n2a=plotg.lat(nc).coor(plotg.polFace(i).nod(j,3),1:3)-plotg.lat(nc).coor(plotg.polFace(i).nod(j,1),1:3);
            n3a=cross(n1a,n2a);
            plotg.polFace(i).normal(j,:)=n3a;
        end
    end    
    
function hl2=plotOpt(opt)
    lighting flat; 
    axis equal; 
    axis off; 
    view(opt.AZ,opt.EL);
    hl2=camlight(-15,30); 
    set(gcf,'color','w');
    hl2.Style='infinite';
    material dull;
    set(gca,'cameraviewanglemode','manual');
    %set(gca,'outerposition',[0 0 1.0 1.0])
    
function printGif(opt,nMode,f,nameFolder,time)
    pause(1/opt.frames)
    if strcmp(opt.saveMovie,'on')
        name=[nameFolder,'/',strcat(opt.template, '_', time)];
        switch opt.plot
            case 'modes'
                name=[name,'_Modes_'];
        end
        switch opt.periodic
            case 'on'
                name=[name,'_pcb_',num2str(opt.plotPer),'.gif'];
            case 'off'
                name=[name,'.gif'];
        end
        if opt.safeMovieAntiAlias==0
            frame = getframe(f.Number);
        else
            myaa(opt.safeMovieAntiAlias) % what's myaa?
            pause(0.01) % changed, was 0.03
            frame = getframe(f.Number+1);
        end


        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        if nMode==1;
            imwrite(imind,cm,name,'gif', 'Loopcount',inf,'delaytime',0);
        else
            imwrite(imind,cm,name,'gif','WriteMode','append','delaytime',0);
        end
        if opt.safeMovieAntiAlias~=0
            close(f.Number+1)
        end
    end
    
function printHigRes(f,opt,nam,nameFolder)
    pause(1/opt.frames)
    switch opt.saveFig
        case 'on'
            name=[nameFolder,'/',opt.template,'_',num2str(opt.plotPer),'_',nam];
            figpos=getpixelposition(f); %dont need to change anything here
            resolution=get(0,'ScreenPixelsPerInch'); %dont need to change anything here
            set(f,'paperunits','inches','papersize',figpos(3:4)/resolution,...
            'paperposition',[0 0 figpos(3:4)/resolution]); %dont need to change anything here
            print(f,name,'-dpng',['-r',num2str(opt.figDPI)],'-opengl') %save file
    end

function [f,hs,hie,his] = copyFigure(unitCell,extrudedUnitCell,opt,hs,hie,his)
    f=figure('Position', [0 0 800 800]);
    ax = axes;
    for nc=1:size(extrudedUnitCell.latVec,1)
        for i=3:10
            hs{nc,i}=copyobj(hs{nc,i},ax);
        end
    end
    for nc=1:size(unitCell.Polyhedron(1).latVec,1)
        for ne=1:length(unitCell.Polyhedron)
            for i=3:10
                hie{ne,nc,i}=copyobj(hie{ne,nc,i},ax);
                his{ne,nc,i}=copyobj(his{ne,nc,i},ax);
            end
        end
    end
    set(gca,'xlim',opt.xlim,'ylim',opt.ylim,'zlim',opt.zlim);
    hl2=plotOpt(opt);

    
function ant=addText(ant, nMode, timesteps)
    timestep_txt = strcat('iteration=', num2str(timesteps(nMode)));
    if isempty(ant)
        ant = annotation('textbox', [0.1, 0, 0.8, 0.1],...
            'String', timestep_txt,...
            'LineStyle', 'none',...
            'HorizontalAlignment', 'center');
    else
        set(ant, 'string', timestep_txt)
    end