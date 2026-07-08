function [x_out, v_out, dx, dv] = obtainDiffs(x, v, delta_v, size, tol, window)

  % Always column vector
  x = x(:);
  v = v(:);

  % x array monotonical
  vec = [true; diff(x) ~= 0];
  x = x(vec);
  v = v(vec);

  % Define array size
  if ~exist('size','var')
    size = 1e4;
  end
  x_new = linspace(x(1), x(end), size)';
  v = linearinterp(x, v, x_new);
  x = x_new;

  % Initialization.
  if ~exist('delta_v','var')
    delta_v = 1e-3;
  end
  if ~exist('tol','var')
    tol = 5e-5;
  end
  nsamp = length(v);

  % Define voltage points.
  vmin = min(v) - tol;
  vmax = max(v) + tol;
  vout = (vmin:delta_v:vmax)';
  
  % Approximate differential capacity dz/dv.
  Nv = histcounts(v, length(vout));
  dx = (max(x)-min(x))/nsamp*Nv/delta_v;

  % Interpolate output soc.
  xout = linearinterp(v, x, vout);
  
  % Force the boundaries to agree with the original data-set.
  [~,ind0] = min(xout); xout(ind0) = min(x); 
  [~,ind1] = max(xout); xout(ind1) = max(x);
  [~,ind0] = min(vout); vout(ind0) = min(v);
  [~,ind1] = max(vout); vout(ind1) = max(v);
  
  % Collect output.
  xout = xout(:); 
    x_out = x(:);
  vout = vout(:);
    v_out = linearinterp(xout, vout, x_out);
  dx = dx(:);
    dx = linearinterp(xout, dx, x_out);
  dv = 1./dx(:);

  % Somoothing
  if ~exist('window','var')
    window = 200;
  end
  kernel = ones(window, 1) / window;
  dx = conv(dx, kernel, 'same');
  dv = conv(dv, kernel, 'same');


