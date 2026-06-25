classdef KuramotoRing
% Shared helper routines for the Kuramoto examples.

methods(Static)

function A = nearestNeighbourRingAdjacency(n)

    A = zeros(n,n);

    for i = 1:n
        A(i,mod(i,n)+1)   = 1;
        A(i,mod(i-2,n)+1) = 1;
    end
end

function Phi = lockingResidualPolynomial(sv,cv,omega,A,K)

    n = length(sv);
    Omega = cell(n,1);

    for i = 1:n

        coupling = 0*sv(1);

        for j = 1:n
            if A(i,j) ~= 0
                coupling = coupling + A(i,j)*(sv(j)*cv(i) - cv(j)*sv(i));
            end
        end

        Omega{i} = omega(i) + (K/n)*coupling;
    end

    Phi = 0*sv(1);

    for i = 1:n
        ip = mod(i,n)+1;
        Phi = Phi - (Omega{ip} - Omega{i})^2;
    end

    Phi = Phi/n;
end

function f = polynomialVectorField(sv,cv,omega,A,K)

    n = length(sv);
    expr = cell(2*n,1);

    for i = 1:n

        coupling = 0*sv(1);

        for j = 1:n
            if A(i,j) ~= 0
                coupling = coupling + A(i,j)*(sv(j)*cv(i) - cv(j)*sv(i));
            end
        end

        thetaDot = omega(i) + (K/n)*coupling;

        expr{i}   = cv(i)*thetaDot;
        expr{n+i} = -sv(i)*thetaDot;
    end

    f = vertcat(expr{:});
end

function [V,coeffV] = phaseDifferenceAuxiliary(sv,cv,degV,useReflectionSymmetry)

    n = length(sv);

    if degV <= 0
        V = 0*sv(1);
        coeffV = [];
        return
    end

    degQ = floor(degV/2);

    if degQ <= 0
        V = 0*sv(1);
        coeffV = [];
        return
    end

    V = 0*sv(1);
    coeffV = [];

    if ~useReflectionSymmetry

        for i = 1:n

            ip = mod(i,n)+1;

            a = cv(i)*cv(ip) + sv(i)*sv(ip);
            b = sv(ip)*cv(i) - cv(ip)*sv(i);

            mon = KuramotoRing.monomialsTwoVars(a,b,degQ);
            mon = KuramotoRing.removeConstantMonomial(mon);

            ci = sdpvar(length(mon),1);

            V = V + ci.'*mon;
            coeffV = [coeffV; ci]; 
        end

        V = V/n;
        return
    end

    used = false(n,1);

    for i = 1:n

        if used(i)
            continue
        end

        r = KuramotoRing.reflectedEdgeIndex(i,n);

        ip = mod(i,n)+1;

        a = cv(i)*cv(ip) + sv(i)*sv(ip);
        b = sv(ip)*cv(i) - cv(ip)*sv(i);

        mon = KuramotoRing.monomialsTwoVars(a,b,degQ);
        mon = KuramotoRing.removeConstantMonomial(mon);

        ci = sdpvar(length(mon),1);

        V = V + ci.'*mon;
        coeffV = [coeffV; ci]; 

        used(i) = true;

        if r ~= i

            rp = mod(r,n)+1;

            ar = cv(r)*cv(rp) + sv(r)*sv(rp);
            br = sv(rp)*cv(r) - cv(rp)*sv(r);

            monr = KuramotoRing.monomialsTwoVars(ar,br,degQ);
            monr = KuramotoRing.removeConstantMonomial(monr);

            V = V + ci.'*monr;
            used(r) = true;
        end
    end

    V = V/n;
end

function mon = monomialsTwoVars(a,b,degQ)

    mon = [];

    for totalDeg = 0:degQ
        for powA = totalDeg:-1:0
            powB = totalDeg - powA;
            mon = [mon; a^powA*b^powB]; 
        end
    end
end

function mon = removeConstantMonomial(mon)

    if length(mon) > 1
        mon = mon(2:end);
    else
        mon = [];
    end
end

function r = reflectedEdgeIndex(i,n)

    r = mod(n - i - 1,n) + 1;
end

function baseCliques = buildEdgeDifferenceBaseCliques(n)

    baseCliques = cell(n,1);

    for i = 1:n
        baseCliques{i} = KuramotoRing.cyclicNodes(i-1,4,n);
    end
end

function chordalCliques = chordalExtensionCliques(n,baseCliques)

    G = false(n,n);

    for a = 1:numel(baseCliques)
        C = baseCliques{a};

        for ii = 1:numel(C)
            for jj = ii+1:numel(C)
                G(C(ii),C(jj)) = true;
                G(C(jj),C(ii)) = true;
            end
        end
    end

    G(1:n+1:end) = false;

    remaining = true(1,n);
    elimCliques = {};

    for step = 1:n 

        remNodes = find(remaining);

        if numel(remNodes) == 1
            v = remNodes(1);
        else
            degs = zeros(size(remNodes));

            for k = 1:numel(remNodes)
                vtmp = remNodes(k);
                degs(k) = nnz(G(vtmp,remNodes));
            end

            [~,idx] = min(degs);
            v = remNodes(idx);
        end

        neigh = remNodes(G(v,remNodes));
        C = sort([v neigh]);

        elimCliques{end+1} = C; 

        for ii = 1:numel(neigh)
            for jj = ii+1:numel(neigh)
                G(neigh(ii),neigh(jj)) = true;
                G(neigh(jj),neigh(ii)) = true;
            end
        end

        remaining(v) = false;
    end

    keep = true(1,numel(elimCliques));

    for a = 1:numel(elimCliques)

        Ca = elimCliques{a};

        for b = 1:numel(elimCliques)

            if a == b
                continue
            end

            Cb = elimCliques{b};

            if numel(Ca) <= numel(Cb) && all(ismember(Ca,Cb))
                if numel(Ca) < numel(Cb) || a > b
                    keep(a) = false;
                    break
                end
            end
        end
    end

    chordalCliques = elimCliques(keep);
end

function [F,coeffAll] = chordalTorusSOSConstraint(D,h,sv,cv,d,chordalCliques)

    n = length(sv);
    x = [sv; cv];

    F = [];
    coeffAll = [];

    degRho = max(0,d-2);
    eqPart = 0*x(1);

    for i = 1:n

        cliqueIdx = KuramotoRing.findSmallestCliqueContainingNode(chordalCliques,i);
        C = chordalCliques{cliqueIdx};
        z = [sv(C); cv(C)];

        [rhoi,ci] = polynomial(z,degRho);

        eqPart = eqPart + rhoi*h(i);
        coeffAll = [coeffAll; ci]; 
    end

    sosPart = 0*x(1);

    for a = 1:numel(chordalCliques)

        C = chordalCliques{a};
        z = [sv(C); cv(C)];

        [sigma,cs] = polynomial(z,d);

        F = [F, sos(sigma)]; 

        sosPart = sosPart + sigma;
        coeffAll = [coeffAll; cs]; 
    end

    resid = D - eqPart - sosPart;
    coeffResid = coefficients(resid,x);

    F = [F, coeffResid == 0];
end

function idx = findSmallestCliqueContainingNode(cliques,node)

    candidates = [];

    for a = 1:numel(cliques)
        if ismember(node,cliques{a})
            candidates(end+1) = a; 
        end
    end

    if isempty(candidates)
        error('No chordal clique contains node %d.',node);
    end

    sizes = cellfun(@numel,cliques(candidates));
    [~,j] = min(sizes);

    idx = candidates(j);
end

function nodes = cyclicNodes(startIdx,L,n)

    nodes = zeros(1,L);

    for k = 1:L
        nodes(k) = mod(startIdx+k-2,n)+1;
    end

    [~,ia] = unique(nodes,'stable');
    nodes = nodes(sort(ia));
end

function [Kc,thetaBest,info] = firstLockedCouplingByBranchSearch(n,omega)

    S = cumsum(omega(:));

    Kflow = n*(max(S) - min(S))/2;
    Klo = max(0,Kflow*(1 - 1e-10));
    Khi = max(1,Kflow*1.2 + 1);

    opts.numAlphaGrid = 5000;
    opts.resTol = 1e-8;

    thetaBest = [];
    info = struct();

    found = false;

    for tries = 1:40 

        [found,thetaTmp,alphaTmp,branchTmp,targetTmp,resTmp] = ...
            KuramotoRing.feasibleLockedAtK(Khi,n,S,opts);

        if found
            break
        end

        Khi = 1.5*Khi + 1;
    end

    if ~found
        Kc = nan;
        info.found = false;
        return
    end

    bestTheta = thetaTmp;
    bestAlpha = alphaTmp;
    bestBranch = branchTmp;
    bestTarget = targetTmp;
    bestRes = resTmp;

    for it = 1:60 

        Kmid = 0.5*(Klo + Khi);

        [found,thetaTmp,alphaTmp,branchTmp,targetTmp,resTmp] = ...
            KuramotoRing.feasibleLockedAtK(Kmid,n,S,opts);

        if found
            Khi = Kmid;
            bestTheta = thetaTmp;
            bestAlpha = alphaTmp;
            bestBranch = branchTmp;
            bestTarget = targetTmp;
            bestRes = resTmp;
        else
            Klo = Kmid;
        end
    end

    Kc = Khi;
    thetaBest = bestTheta;

    info.found = true;
    info.alpha = bestAlpha;
    info.branch = bestBranch;
    info.target = bestTarget;
    info.closureResidual = bestRes;
    info.Kflow = Kflow;
end

function [found,thetaBest,alphaBest,branchBest,targetBest,resBest] = ...
    feasibleLockedAtK(K,n,S,opts)

    found = false;
    thetaBest = [];
    alphaBest = nan;
    branchBest = [];
    targetBest = nan;
    resBest = inf;

    alphaMin = max(S - K/n);
    alphaMax = min(S + K/n);

    if alphaMin > alphaMax
        return
    end

    alphaGrid = linspace(alphaMin,alphaMax,opts.numAlphaGrid);
    numBranches = 2^n;

    for bint = 0:numBranches-1

        branch = bitget(bint,1:n);
        gVals = nan(size(alphaGrid));

        for ia = 1:numel(alphaGrid)
            gVals(ia) = KuramotoRing.closureSum(alphaGrid(ia),K,n,S,branch);
        end

        gMin = min(gVals);
        gMax = max(gVals);

        mMin = ceil(gMin/(2*pi));
        mMax = floor(gMax/(2*pi));

        for m = mMin:mMax

            target = 2*pi*m;
            H = gVals - target;

            [gridRes,imin] = min(abs(H));

            if gridRes < resBest
                resBest = gridRes;
                alphaBest = alphaGrid(imin);
                branchBest = branch;
                targetBest = target;
                thetaBest = KuramotoRing.thetaFromAlphaBranch(alphaBest,K,n,S,branch);
            end

            if gridRes < opts.resTol
                found = true;
                return
            end

            for ia = 1:numel(alphaGrid)-1

                h1 = H(ia);
                h2 = H(ia+1);

                if ~isfinite(h1) || ~isfinite(h2)
                    continue
                end

                if h1 == 0 || h1*h2 < 0
                    try
                        fun = @(a) KuramotoRing.closureSum(a,K,n,S,branch) - target;
                        alphaRoot = fzero(fun,[alphaGrid(ia),alphaGrid(ia+1)]);
                        res = abs(fun(alphaRoot));

                        if res < resBest
                            resBest = res;
                            alphaBest = alphaRoot;
                            branchBest = branch;
                            targetBest = target;
                            thetaBest = KuramotoRing.thetaFromAlphaBranch(alphaBest,K,n,S,branch);
                        end

                        if res < opts.resTol
                            found = true;
                            return
                        end
                    catch
                    end
                end
            end
        end
    end
end

function g = closureSum(alpha,K,n,S,branch)

    x = n*(alpha - S(:))/K;

    if any(abs(x) > 1 + 1e-12)
        g = nan;
        return
    end

    x = max(-1,min(1,x));

    delta = asin(x);
    idx = branch(:) == 1;
    delta(idx) = pi - delta(idx);

    g = sum(delta);
end

function theta = thetaFromAlphaBranch(alpha,K,n,S,branch)

    x = n*(alpha - S(:))/K;
    x = max(-1,min(1,x));

    delta = asin(x);

    idx = branch(:) == 1;
    delta(idx) = pi - delta(idx);

    theta = zeros(n,1);

    for i = 1:n-1
        theta(i+1) = theta(i) + delta(i);
    end

    theta = theta - mean(theta);
    theta = KuramotoRing.wrapToPi(theta);
end

function branchData = computeLockedBranchCloud( ...
    n,omega,A,Kgrid,numRandomStarts,rootTol,uniqueTol,observableName)

    rng(1);

    optsFS = optimset('Display','off', ...
                      'TolFun',1e-12, ...
                      'TolX',1e-12, ...
                      'MaxIter',1000, ...
                      'MaxFunEvals',50000);

    branchData = struct('K',{},'theta',{},'Y',{},'isStable',{},'maxEig',{});
    prevRoots = [];

    for ik = 1:numel(Kgrid)

        K = Kgrid(ik);
        guesses = [];

        if ~isempty(prevRoots)
            guesses = [guesses prevRoots]; 
        end

        guesses = [guesses, 2*pi*(rand(n-1,numRandomStarts)-0.5)]; 

        rootsK = [];

        for g = 1:size(guesses,2)

            u0 = guesses(:,g);

            try
                fun = @(u) KuramotoRing.lockedResidualReduced(u,omega,A,K);
                [u,~,exitflag] = fsolve(fun,u0,optsFS);
                res = norm(fun(u),Inf);

                if exitflag > 0 && res < rootTol

                    theta = [u(:); 0];
                    theta = KuramotoRing.normalizeTheta(theta);

                    if KuramotoRing.isNewRoot(theta,rootsK,uniqueTol)

                        rootsK = [rootsK theta]; 

                        [isStable,lam] = KuramotoRing.lockedStateStability(theta,A,K);
                        Y = KuramotoRing.lockedObservable(theta,A,observableName);

                        branchData(end+1).K = K; 
                        branchData(end).theta = theta;
                        branchData(end).Y = Y;
                        branchData(end).isStable = isStable;
                        branchData(end).maxEig = lam(2);
                    end
                end
            catch
            end
        end

        if isempty(rootsK)
            prevRoots = [];
        else
            prevRoots = rootsK(1:n-1,:);
        end

        if mod(ik,20) == 0 || ik == 1
            fprintf('K = %.4f: found %d roots\n',K,size(rootsK,2));
        end
    end
end

function F = lockedResidualReduced(u,omega,A,K)

    n = length(omega);
    theta = [u(:); 0];

    rhs = KuramotoRing.thetaRHS(theta,omega,A,K);

    F = rhs(1:n-1);
end

function rhs = thetaRHS(theta,omega,A,K)

    n = length(theta);
    rhs = omega(:);

    for i = 1:n

        coupling = 0;

        for j = 1:n
            if A(i,j) ~= 0
                coupling = coupling + A(i,j)*sin(theta(j)-theta(i));
            end
        end

        rhs(i) = rhs(i) + (K/n)*coupling;
    end
end

function [isStable,lam] = lockedStateStability(theta,A,K)

    n = length(theta);
    J = zeros(n,n);

    for i = 1:n

        for j = 1:n
            if i ~= j && A(i,j) ~= 0
                J(i,j) = (K/n)*A(i,j)*cos(theta(j)-theta(i));
            end
        end

        J(i,i) = -sum(J(i,[1:i-1,i+1:n]));
    end

    lam = eig(J);
    lam = sort(real(lam),'descend');

    isStable = lam(2) < -1e-7;
end

function val = lockedObservable(theta,A,observableName)

    switch observableName

        case 'edgeCoherence'
            val = KuramotoRing.edgeCoherence(theta,A);

        case 'R2'
            val = KuramotoRing.classicalR2(theta);

        case 'maxEdgeDiff'
            val = KuramotoRing.maxEdgeDiff(theta);

        otherwise
            error('Unknown observableName: %s',observableName);
    end
end

function val = edgeCoherence(theta,A)

    n = length(theta);
    val = 0;
    m = nnz(A);

    for i = 1:n
        for j = 1:n
            if A(i,j) ~= 0
                val = val + (1 + cos(theta(j)-theta(i)))/2;
            end
        end
    end

    val = val/m;
end

function val = classicalR2(theta)

    n = length(theta);
    val = abs(sum(exp(1i*theta))/n)^2;
end

function val = maxEdgeDiff(theta)

    n = length(theta);
    val = 0;

    for i = 1:n
        ip = mod(i,n)+1;
        val = max(val,abs(KuramotoRing.wrapToPi(theta(ip)-theta(i))));
    end
end

function theta = normalizeTheta(theta)

    theta = theta(:);
    theta = theta - theta(end);
    theta = KuramotoRing.wrapToPi(theta);
end

function tf = isNewRoot(theta,roots,uniqueTol)

    if isempty(roots)
        tf = true;
        return
    end

    tf = true;

    for k = 1:size(roots,2)

        dist = KuramotoRing.phaseDistance(theta,roots(:,k));

        if dist < uniqueTol
            tf = false;
            return
        end
    end
end

function d = phaseDistance(theta1,theta2)

    theta1 = KuramotoRing.normalizeTheta(theta1);
    theta2 = KuramotoRing.normalizeTheta(theta2);

    diff = KuramotoRing.wrapToPi(theta1 - theta2);

    d = norm(diff,Inf);
end

function x = wrapToPi(x)

    x = mod(x + pi,2*pi) - pi;
end

function savePDF(fig,fileName)

    set(fig,'PaperUnits','centimeters');
    set(fig,'Units','centimeters');

    pos = get(fig,'Position');

    set(fig,'PaperSize',[pos(3) pos(4)]);
    set(fig,'PaperPositionMode','manual');
    set(fig,'PaperPosition',[0 0 pos(3) pos(4)]);

    print(fig,'-dpdf',fileName);
    fprintf('Saved figure: %s\n',fileName);
end

end
end