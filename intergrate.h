#pragma once
#ifndef __OMEGA_INTGRATE_H__
#define __OMEGA_INTGRATE_H__

#include <vector>
#include <functional>

namespace Omega
{
	using namespace std;

	template <typename ResultT, typename ParamT, typename EpsT>
	ResultT Trapezoid(function<ResultT(ParamT)> func, ParamT begin, ParamT end, EpsT eps)
	{
		EpsT difference = 2 * eps;
		ParamT h = end - begin;
		ResultT integ0 = 0.5 * h * (func(begin) + func(end));
		ResultT integ = 0;
		for (size_t n = 1; difference >= eps; n *= 2, h /= 2)
		{
			ParamT x = begin + 0.5 * h;
			ResultT sum = 0;
			for (size_t i = 0; i < n; x += h, ++i)
			{
				sum += func(x);
			}
			integ = 0.5 * (integ0 + h * sum);
			difference = abs(integ0 - integ);
			integ0 = integ;
		}
		return integ;
	}

	template <typename ResultT, typename ParamT, typename EpsT>
	ResultT Simpson(function<ResultT(ParamT)> func, ParamT begin, ParamT end, EpsT eps)
	{
		EpsT difference = 2 * eps;
		ParamT h = end - begin;
		ResultT trap0 = 0.5 * h * (func(begin) + func(end));
		ResultT integ0 = trap0;
		ResultT trap = 0, integ = 0;
		for (size_t n = 1; difference >= eps; n *= 2, h /= 2)
		{
			ParamT x = begin + 0.5 * h;
			ResultT sum = 0;
			for (size_t i = 0; i < n; x += h, ++i)
			{
				sum += func(x);
			}
			trap = 0.5 * (trap0 + h * sum);
			integ = (4.0 * trap - trap0) / 3.0;
			difference = abs(integ0 - integ);
			trap0 = trap;
			integ0 = integ;
		}
		return integ;
	}

	template <typename ResultT, typename ParamT, typename EpsT>
	ResultT Romberg(function<ResultT(ParamT)> func, ParamT begin, ParamT end, EpsT eps)
	{
		EpsT difference = 2 * eps;
		ParamT h = end - begin;
		ResultT integ0 = 0.5 * h * (func(begin) + func(end));
		ResultT integ = 0;
		vector<ResultT> y(1, integ0);3 deffre
		for (size_t m = 0, n = 1; difference >= eps; h *= 0.5, n *= 2, ++m)
		{
			ParamT x = begin + 0.5 * h;
			ResultT sum = 0;
			for (size_t i = 0; i < n; x += h, ++i)
			{
				sum += func(x);
			}
			integ0 = 0.5 * (y[0] + h * sum);
			sum = 4;
			for (size_t k = 0; k <= m; sum *= 4, ++k)
			{
				integ = (integ0 * sum - y[k]) / (sum - 1.0);
				y[k] = integ0;
				integ0 = integ;
			}
			difference = abs(integ - y[m]);
			y.push_back(integ);
		}
		return integ;
	}
}
#endif // __OMEGA_COMEPLEX_H__
