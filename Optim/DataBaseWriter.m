function DataBaseWriter(DataBaseName,z,Aop,Iop,ATHz,ETHz,ASH,ISH)
try
h5create(DataBaseName,"/"+num2str(z*1e6)+"/Aop",size(Aop));
catch
warning("Data field already exists!")    
end
try
h5create(DataBaseName,"/"+num2str(z*1e6)+"/Eop",size(Iop));
catch
warning("Data field already exists!")    
end
try
h5create(DataBaseName,"/"+num2str(z*1e6)+"/ATHz",size(ATHz));
catch
warning("Data field already exists!")    
end
try
h5create(DataBaseName,"/"+num2str(z*1e6)+"/ETHz",size(ETHz));
catch
warning("Data field already exists!")    
end
try
h5create(DataBaseName,"/"+num2str(z*1e6)+"/ASH",size(ASH));
catch
warning("Data field already exists!")    
end
try
h5create(DataBaseName,"/"+num2str(z*1e6)+"/ESH",size(ISH));
catch
warning("Data field already exists!")    
end
h5write(DataBaseName,"/"+num2str(z*1e6)+"/Aop",Aop);
h5write(DataBaseName,"/"+num2str(z*1e6)+"/Eop",Iop);
h5write(DataBaseName,"/"+num2str(z*1e6)+"/ATHz",ATHz);
h5write(DataBaseName,"/"+num2str(z*1e6)+"/ETHz",ETHz);
h5write(DataBaseName,"/"+num2str(z*1e6)+"/ASH",ASH);
h5write(DataBaseName,"/"+num2str(z*1e6)+"/ESH",ISH);
end