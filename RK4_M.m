function varargout = RK4_M(f,step,t0,y0,t_final)
T=t0:step:t_final;
Y = zeros(length(T),size(y0,2),size(y0,3));
Y(1,:,:) = y0;
for ii = 2:length(T)
k1 = f(T(ii-1),Y(ii-1,:,:));
k2 = f(T(ii-1)+step/2,Y(ii-1,:,:)+k1*step/2);
k3 = f(T(ii-1)+step/2,Y(ii-1,:,:)+k2*step/2);
k4 = f(T(ii-1)+step,Y(ii-1,:,:)+k3*step);
Y(ii,:,:) =Y(ii-1,:,:) + 1/6*(k1+2*k2+2*k3+k4)*step;
disp(ii/length(T)*100);
end
if nargout == 1
    varargout{1} = Y;
elseif nargout == 2
    varargout{1} = T;
    varargout{2} = Y;
else
    error("Nem megfelelő számú kimenet! (1 vagy 2)");
end

end