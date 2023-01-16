function DataBaseEnder(DataBaseName,z,t,nu,effic,efficSH)
try
h5create(DataBaseName,"/z",size(z));
catch
warning("Data field already exists!")    
end
try
h5create(DataBaseName,"/effic",size(effic));
catch
warning("Data field already exists!")    
end
try
h5create(DataBaseName,"/efficSH",size(efficSH));
catch
warning("Data field already exists!")    
end
try
h5create(DataBaseName,"/t",size(t));
catch
warning("Data field already exists!")    
end
try
h5create(DataBaseName,"/nu",size(nu));
catch
warning("Data field already exists!")    
end
h5write(DataBaseName,"/z",z);
h5write(DataBaseName,"/effic",effic);
h5write(DataBaseName,"/efficSH",efficSH);
h5write(DataBaseName,"/t",t);
h5write(DataBaseName,"/nu",nu);
end
