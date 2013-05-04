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

## testi

## Author: Andraz Vrhovec <az@ares>
## Created: 2013-05-02

function [ times, steps ] = testi ()

	times = [];
	steps = [];
	err = [];

	REPS = 10;
	num_mat = 25;
	AAA = cell(num_mat, 1);

	for i = 1:num_mat
		AAA(i) = rand(10*i);
	end

	for i=1:num_mat
		A = AAA{i};
		tic();
		ste = 0;
		for j = 1:REPS
			[lv, lambda, st] = inv_lastni(A, i*4);
			ste = st + ste;
		end
		e = norm(A*lv - lv*lambda);
		time = toc()/REPS;
		times = [times time];
		steps = [steps ste/REPS];
		err = [err e];
	end

	x = 10:10:num_mat*10;

	figure(1);
	plot(x, times);
	figure(2);
	plot(x, steps);
	figure(3);
	plot(x, err);



endfunction
