clear

load approx_dirac_cauchy.dat  
load approx_dirac_caulor.dat  
load approx_dirac_gaussi.dat  
load approx_dirac_sincfc.dat  
load approx_dirac_triang.dat

ntimes=length(approx_dirac_cauchy);
dt=approx_dirac_cauchy(2,1)-approx_dirac_cauchy(1,1);
time=approx_dirac_cauchy(:,1);
nomega=(ntimes-1)/2;
timewidth=time(ntimes)-time(1);
for iomega=1:nomega+1
  omega(iomega) = 2*pi*(iomega-1)/timewidth ;
  frequency(iomega) = (iomega-1)/timewidth ; %! iomega cycles per sample
  if (iomega>1)
   period(iomega) = timewidth/real(iomega-1) ;
  end
%  period(1) has no meaning (corresponds to the mean) 
end

spec_cauchy=fft(approx_dirac_cauchy(:,2));
spec_caulor=fft(approx_dirac_caulor(:,2));
spec_gaussi=fft(approx_dirac_gaussi(:,2));
spec_sincfc=fft(approx_dirac_sincfc(:,2));
spec_triang=fft(approx_dirac_triang(:,2));

figure(1);clf
plot(time,approx_dirac_cauchy(:,2),'r')
hold on
plot(time,approx_dirac_caulor(:,2),'b')
plot(time,approx_dirac_gaussi(:,2),'k--')
plot(time,approx_dirac_sincfc(:,2),'m')
plot(time,approx_dirac_triang(:,2),'r-.')
legend('cauchy','cauchy lorentz','gauss','sinc','triangle')
xlabel('period [s]')
ylabel('power spectrum')

figure(2);clf
loglog(1./frequency(1:nomega+1),abs(spec_cauchy(1:nomega+1)).^2,'r')
hold on
loglog(1./frequency(1:nomega+1),abs(spec_caulor(1:nomega+1)).^2,'b')
loglog(1./frequency(1:nomega+1),abs(spec_gaussi(1:nomega+1)).^2,'k--')
loglog(1./frequency(1:nomega+1),abs(spec_sincfc(1:nomega+1)).^2,'m')
loglog(1./frequency(1:nomega+1),abs(spec_triang(1:nomega+1)).^2,'r-.')
legend('cauchy','cauchy lorentz','gauss','sinc','triangle')
xlabel('period [s]')
ylabel('power spectrum')

