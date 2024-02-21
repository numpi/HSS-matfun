addpath('util')
rng(1);
hssoption('threshold', 1e-8);
hssoption('block-size', 256);
do_print = 1;
debug = 0;

sizes = 2.^[9:17]; 
l = length(sizes); 
deltas = 2e-2 * 2.^[0 : -1 : 1 - l];

phi = 3;
fun = @(M) inv_sqrt(M);
errors = zeros(l, 2); % first column error of sqrtm, second column error of dac
times = zeros(l, 2);  % first column time of sqrtm, second column time of dac
ranks = zeros(l, 2);  % first column hssrank of A, second column hssrank of f(A)
times_dense = zeros(l, 1);
band = zeros(l, 1);   % bandwidth of A
fprintf('size & band & time \\CKR & time \\CKM & time Rat & time Dense & err \\CKR & err \\CKM & err Rat   \\\\ \n \\hline \n');
for j = 1:l
	n = sizes(j);
	delta = deltas(j);
	s = sort(rand(n, 1));

	% Construct the precision matrix A
	A = sparse(n, n);
    
	for b = 1:n
		x = abs(s(1+b:end) - s(1:end-b)) < delta;
		if  sum(x) == 0
			break
		else
			A = A + spdiags([x, zeros(n - b, 1)], -b, n, n);
			b = b + 1;
		end
	end
	A = A + A';
	D = phi * A * ones(n, 1) + ones(n, 1);
	A = spdiags(D, 0, n, n) - phi * A;
	nrmfA = 1.0; %inv_sqrt(eigs(A, 1, 'smallestabs','MaxIterations', 1e5))
	
	band(j) = bandwidth(A);
	hssA = hss(A); 
	hssA = (hssA + hssA')/2;

	ranks(j, 1) = hssrank(hssA);
    opts.treshold=256;
    opts.issymm=1;
    H = full_to_hss2(full(A),opts);
    b==eigs(A,1,'largestabs', 'MaxIterations', size(A,1));
    a=eigs(A,1,'smallestabs', 'MaxIterations', size(A,1));
    eps=1e-8;
    Phi=@(x)(((2.*x-(b+a))/(b-a))-sqrt(((2.*x-(b+a))/(b-a)).^2-1));
    m=ceil(log(8*log2(n/opts.treshold)/(sqrt(a)*eps*abs(Phi(0))))*log(16*b/a)/pi^2);
    poles = poles_Markov_functions(a,b,-Inf,0,m);
    opts.deflationTol=hssoption('threshold');
	tic
   % profile on
	sA = hss2_funm_symm_telescop(H,@(x)1./sqrt(real(x)),poles,opts);
	%profsave
    %profile off
	times(j, 1) = toc;
    
    poles = [0 inf];
	tic
    % profile on
	sA2 = hss_fun_dac_band_hermitian(A, fun, poles, debug, nrmfA, 0);
	%profsave
    %profile off
	times(j, 2) = toc;
	ranks(j, 2) = hssrank(sA2);
	
    
	if n <= 8192 
        [FH3, timeRat] = rational_evaluation(hss(A),@(x) 1./sqrt(x),eig(A));
		fA = full(A); fA = (fA + fA')/2;
		tic
		fA = inv_sqrt(fA); 
		times_dense(j) = toc;
		errors(j, 1) = norm(fA - full(sA), 'fro') / norm(fA, 'fro');
		errors(j, 2) = norm(fA - full(sA2), 'fro') / norm(fA, 'fro');
        errRat = norm(fA - full(FH3), 'fro') / norm(fA, 'fro');
		if do_print
			fprintf('%d & %d & %.2f & %.2f & %.2f & %.2f & %1.2e & %1.2e & %1.2e \\\\ \n', n, band(j), times(j, 1), times(j, 2),timeRat, times_dense(j), errors(j, 1), errors(j, 2),errRat);
		end
	else
		if do_print
			fprintf('%d & %d & %.2f & %.2f & &  &  \\\\ \n', n, band(j), times(j, 1), times(j, 2));
		end
	end
% pause
end


%----------------- Auxiliary function --------------
function F = inv_sqrt(A)
    A = (A + A')/2;
	[V, D] = eig(A, 'vector');
	D = 1./sqrt(D);
	F = V * diag(D) * V';
end
%------------------------------------------------

function poles = poles_Markov_functions(a,b,alpha,beta,m)
%eigenvalues in [a,b] integration domain (alpha,beta)
%[1] Beckermann Reichel

phi = @(x) (((2.*x-(b+a))/(b-a))-sqrt(((2.*x-(b+a))/(b-a)).^2-1));
psi =  @(x) (b-a)/4*(x+1./x) + (b+a)/2; % inverse of phi
if abs(alpha)==Inf
    kappa = -1/phi(beta);
else
    kappa = (phi(alpha)-phi(beta))/(1-phi(alpha)*phi(beta));
end
% k = ((1-sqrt(1-kappa^2))/kappa)^2;
k = ( kappa / ( 1 + sqrt(1-kappa^2) ) )^2;			% version with less cancellation
T1 = @(x)(1+phi(beta)*x)/(x+phi(beta)); % corresponds to T_1^{-1} in [1]
T2 = @(x) (sqrt(k).*x-1)./(x-sqrt(k));    % corresponds to T_1^{-1} in [1]
T = @(x) T1(T2(x));                       % corresponds to T^{-1} in [1]
K = ellipke(k^2);		
% K = ellipke(k);		
poles = zeros(1,m);
for i = 1:m
    poles(i) = sqrt(k) * ellipj(K*(m+1-2*i)/m,k^2); 		
    % poles(i) = sqrt(k) * ellipj(K*(m+1-2*i)/m, k); 		
	% poles(i) corresponds to \hat w_i of [1]
    poles(i) = T(poles(i));
    poles(i) = psi(poles(i));

end
end

