function setAspectRatio()

  % axis equal;
  h = get(gca,'DataAspectRatio'); 
  if h(3)==1
    set(gca,'DataAspectRatio',[1 1 1/max(h(1:2))]);
  else
    set(gca,'DataAspectRatio',[1 1 h(3)]);
  end       

return
end	