function RedoTicks(h)
Chil=get(h,'Children');
for i = 1:length(Chil)
    if strcmp(Chil(i).Type,'axes')
        set(Chil(i),'TickDir','out','TickLength',[0.03 0.03],'box','off') 
    end
    if strcmp(Chil(i).Type,'legend')
        set(Chil(i),'box','off')
    end
end
end