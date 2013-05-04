## Copyright (C) 2013 Andraz Vrhovec
## 
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or
## (at your option) any later version.
## 
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU General Public License
## along with Octave; see the file COPYING.  If not, see
## <http://www.gnu.org/licenses/>.

## hessenberg

## Author: Andraz Vrhovec <az@ares>
## Created: 2013-05-02

function [ H, Q ] = hessenberg (A)
% pretvori matriko v zgornje hessenbergovo obliko

	[m, n] = size(A);
	if issparse(A)
		I = speye(size(A));
		v = sparse(zeros(m, 1));
	else
		I = eye(size(A));
		v = zeros(m, 1);
	end

	PP = I;

	for k=1:m-2
		alpha = -sign(A(k+1, k)) * sqrt(sum(A(k+1:end, k).^2));
		r = sqrt((alpha^2-A(k+1, k)*alpha)/2);
		v(1:k) = zeros(k, 1);
		v(k+1) = (A(k+1, k)-alpha)/(2*r);
		v(k+2:m) = A(k+2:m, k)./(2*r);
		P = I - 2*v*v';
		PP = PP * P;
		A = P*A*P;
	end

	H = triu(A, -1);
	Q = PP;
endfunction
