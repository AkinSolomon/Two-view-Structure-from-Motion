%PLOTMATCH: Takes in images and produces the matched points
%across the images using their closest sift descriptors
%the plot contains the frames in each image, with a line drawn between the
%matches.

function plotmatch(Ia,Ib,f1,f2,matchindx)
    hold on
    montage([Ia Ib])
    hold on
    h1 = vl_plotframe(f1(:,matchindx));
    h2=  vl_plotframe(f1(:,matchindx));
    set(h1,'color','k','linewidth',3);
    set(h2,'color','y','linewidth',2);
    title('SIFT matched features')
    
    f2shift=f2;
    f2shift(1,:)=f2(1,:)+ size(Ia,2);
    h1 = vl_plotframe(f2shift(:,matchindx));
    h2 =  vl_plotframe(f2shift(:,matchindx));
    set(h1,'color','k','linewidth',3);
    set(h2,'color','y','linewidth',2);
    title('SIFT matched features')
    
    for i=1:length(matchindx)
        x=[f1(1,i) f2shift(1,i)];
        y=[f1(2,i) f2shift(2,i)];
        line(x,y)
    end
    hold off
end