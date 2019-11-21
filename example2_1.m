clear; clc; close all

x1 = [0,0]; x2 = [0,10]; x3 = [10,0]; x4 = [10,10];
X = [x1;x2;x3;x4]'; x = [2,3]'; L = size(X,2); 
d = (sqrt(sum((x*ones(1,L)-X).^2,1))).';

m = 0; iterN = 30;
for dB = -10:5:60
    sigma2 = (d.^2)/(10.^(dB/10));
    for i = 1:1000
        r = d + randn(L,1).*sqrt(sigma2);    
        for k = 1:iterN
            [H, H1] = hessian_nls(X,x,r,sigma2);
            [g,g1] = grad_nls(X,x,r,sigma2);
            x = x - inv(H)*g;
            x_nr(i,:) = x;           
        end
        
        for k = 1:iterN
            [H, H1] = hessian_nls(X,x,r,sigma2);
            [g,g1] = grad_nls(X,x,r,sigma2);
            x = x - inv(H1)*g1;
            x_gn(i,:) = x;           
        end
    end
    res_1 = (x_nr(:,1) - 2).^2  + (x_nr(:,2) - 3).^2;
    res_2 = (x_gn(:,1) - 2).^2  + (x_gn(:,2) - 3).^2;
    m = m + 1;
    result_1(m,:) = mean(res_1);
    result_2(m,:) = mean(res_2);
end

dBrange = -10:5:60;
plot(dBrange, 10*log10(result_1), 'k+', dBrange, 10*log10(result_2), 'ko'); 
xlabel('SNR (dB)'); ylabel('Mean Square Position Error (dB)');
legTitle = legend('ML', 'NLS'); title(legTitle, 'Approaches');saveas(gcf, 'MSPE.png')