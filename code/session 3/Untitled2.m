% SI SNR
TX = [35.5845	36.59946312	37.32580338	38.27087806];
TY=[12.84361998	12.5831808	11.65860655	11.83958886];
MWF_TS = [16.10784653	13.9506287	12.85623061	12.56780414];
MVDR = [13.15712074	10.24961724	9.761718842	9.746179091];
MWF = [14.54639209	12.15157407	11.50110777	11.24994712];
rev = [ 0.3 0.61 1 1.5];
figure;
scatter(rev,TX,'bd');hold on
scatter(rev,TY,'ro');hold on
scatter(rev,MWF_TS,'gx');hold on
scatter(rev,MVDR,'ys');hold on
scatter(rev,MWF,'mh');hold on
plot(rev,TX,'b',rev,TY,'r',rev,MWF_TS,'g',rev,MVDR,'y',rev,MWF,'m');
grid on
ylabel('\Delta SI-SNR');xlabel('Reverberation time T_{60}')
legend('TX','TY','MWF-TS','MVDR','MWF');

%% SD

TX = [-2.5365	-2.088595544	-5.707866112	-6.816831946];
TY=[-1.525045737	-1.882834393	-3.509968279	-4.0017852887];
MWF_TS = [-1.809007558	-1.829196714	-3.866680214	-4.295345695];
MVDR = [-0.998869251	-1.352074972	-2.36346076	-2.479337352];
MWF = [-1.471092028	-1.649877328	-3.331805139	-3.561520455];
rev = [ 0.3 0.61 1 1.5];
figure;
scatter(rev,TX,'bd');hold on
scatter(rev,TY,'ro');hold on
scatter(rev,MWF_TS,'gx');hold on
scatter(rev,MVDR,'ys');hold on
scatter(rev,MWF,'mh');hold on
plot(rev,TX,'b',rev,TY,'r',rev,MWF_TS,'g',rev,MVDR,'y',rev,MWF,'m');
grid on
ylabel('\Delta SD');xlabel('Reverberation time T_{60}')
legend('TX','TY','MWF-TS','MVDR','MWF');
%% PESQ
TX = [0.7331 0.670112133 0.547781229 0.4180305];
TY=[-0.014835	0.001072407	-0.020002484	-0.073816776];
MWF_TS = [0.016083121	0.021450043	0.005393505	-0.0521207];
MVDR = [0.059222102	0.036840916	0.050082326	0.007158637];
MWF = [0.18	0.081570029	0.079578519	0.014341474];
rev = [ 0.3 0.61 1 1.5];
figure;
scatter(rev,TX,'bd');hold on
scatter(rev,TY,'ro');hold on
scatter(rev,MWF_TS,'gx');hold on
scatter(rev,MVDR,'ys');hold on
scatter(rev,MWF,'mh');hold on
plot(rev,TX,'b',rev,TY,'r',rev,MWF_TS,'g',rev,MVDR,'y',rev,MWF,'m');
grid on
ylabel('\Delta PESQ');xlabel('Reverberation time T_{60}')
legend('TX','TY','MWF-TS','MVDR','MWF');
%% STOI
TX = [0.3269	0.327075891	0.330213299	0.341039141];
TY=[0.204019538	0.119662622	0.067437486	0.062716232];
MWF_TS = [0.256726437	0.160233238	0.124401907	0.102325903];
MVDR = [0.273428319	0.161486353	0.128489504	0.111021029];
MWF = [0.271222625	0.164325844	0.125776705	0.130200012];
rev = [ 0.3 0.61 1 1.5];
figure;
scatter(rev,TX,'bd');hold on
scatter(rev,TY,'ro');hold on
scatter(rev,MWF_TS,'gx');hold on
scatter(rev,MVDR,'ys');hold on
scatter(rev,MWF,'mh');hold on
plot(rev,TX,'b',rev,TY,'r',rev,MWF_TS,'g',rev,MVDR,'y',rev,MWF,'m');
grid on
ylabel('\Delta STOI');xlabel('Reverberation time T_{60}')
legend('TX','TY','MWF-TS','MVDR','MWF');


