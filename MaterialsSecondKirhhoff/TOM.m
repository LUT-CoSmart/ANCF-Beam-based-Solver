function SS = TOM(F,const,a0)
% TOM second Piola-Kirchhoff stress
%
% The local fibre direction a0 is supplied separately by AddTensors.
% The final entry of const is the bulk/compressibility parameter.

    % TOM material parameters
    muphi = const(1);
    Ephi  = const(2);
    a     = const(3);
    b     = const(4);
    c     = const(5);

    % Bulk/compressibility parameter added by
    % AssigningMaterialParameters
    d = const(end);

    % Use the fibre direction supplied by AddTensors
    a0 = a0(:);

    if norm(a0) == 0
        error('TOM: fibre direction cannot be zero.');
    end

    a0 = a0/norm(a0);

    % Kinematics
    II = eye(3);
    J  = det(F);

    if ~isreal(J) || ~isfinite(J) || J <= 0
        error('TOM: invalid deformation. det(F) = %g.',J);
    end

    J_inv23    = J^(-2/3);
    C          = F.'*F;
    C_dash     = J_inv23*C;
    Cinv       = inv(C);
    C_dash_inv = inv(C_dash);

    % Structural tensor and fibre invariant
    A0    = a0*a0.';
    I4bar = trace(C_dash*A0);

    if ~isreal(I4bar) || ~isfinite(I4bar) || I4bar <= 0
        error('TOM: invalid I4bar = %g.',I4bar);
    end

    % --- kinematics & isochoric split (same as GOH)
    II        = eye(3);
    J         = det(F);
    J_inv23   = J^(-2/3);
    C         = F.'*F;
    C_dash    = J_inv23 * C;
    Cinv      = inv(C);
    C_dash_inv= inv(C_dash);

    % --- structural tensor & invariant
    A0    = a0(:) * a0(:).';               % 3x3
    I4bar = trace(C_dash * A0);            % = a0·Cbar·a0

    % --- TOM piecewise coefficients A(I4bar), B(I4bar), C(I4bar), D(I4bar)
    %     (careful: C(...) named Cc to avoid clash with C=F'*F)
    if I4bar < a^2
        Acoef = 0; Bcoef = 0; Cc = 0; Dcoef = 0;
    elseif I4bar <= c^2
        denom = (b - a) * (c - a);
        Acoef = -a^2 / denom;
        Bcoef =  2*a*log(a) / denom;
        Cc    =  1 / denom;
        Dcoef = -2*a / denom;
    elseif I4bar <= b^2
        Acoef =  c^2 / ((c - a)*(b - c)) - a^2 / ((b - a)*(c - a));
        Bcoef =  2*a*log(a)/((b - a)*(c - a)) - 2*c*log(c)/((c - a)*(b - c));
        Cc    = -1 / ((b - a)*(b - c));
        Dcoef =  2*b / ((b - a)*(b - c));
    else
        Acoef = -1;
        Bcoef =  2*a*log(a)/((b - a)*(c - a)) ...
               + 2*b*log(b)/((b - a)*(b - c)) ...
               - 2*c*log(c)/((c - a)*(b - c));
        Cc    = 0;
        Dcoef = 0;
    end

    % --- energy derivatives w.r.t. C_bar
    % Isochoric matrix part: W_iso = (muphi/2)(I1_bar - 3)  ->  dW/dC_bar = (muphi/2) * I
    PSI_1 = muphi/2;

    % Fiber part: W_f = Ephi * G(I4bar)  ->  dW/dC_bar = Ephi * G'(I4bar) * dI4bar/dC_bar
    % For transversely isotropic case, this reduces to: PSI_4 * A0 with
    %   PSI_4 = Ephi * dG/dI4bar
    % Using your provided derivative (simplified to avoid duplicate D terms):
    sqrtI4 = sqrt(I4bar);
    if ~(I4bar > 0)
        error('TOM: I4bar must be positive. Got %g.', I4bar);
    end
    PSI_4 = Ephi * ( Acoef/(2*I4bar) + Cc/2 + (Bcoef/(2*sqrtI4)) + (Dcoef*log(I4bar))/(4*sqrtI4) );

    % --- deviatoric projection exactly like GOH
    W_C_dash = PSI_1*II + PSI_4*A0;
    DEV      = W_C_dash - (1/3)*trace(W_C_dash * C_dash.') * C_dash_inv;

    % --- final second Piola–Kirchhoff stress
    SS = 2/d*(J - 1)*J*Cinv + 2*J_inv23*DEV;
end
