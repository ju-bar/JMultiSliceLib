// file: 'string_format.cpp'
// author: J. Barthel, Forschungszentrum Juelich GmbH, Juelich, Germany
//         ju.barthel@fz-juelich.de
//
// This file contains the implementation of a string formatting routine.
//
/* -----------------------------------------------------------------------

	This file is part of JMSBench1.

	JMSBench1 is free software: you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation, either version 3 of the License, or
	(at your option) any later version.

	JMSBench1 is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details.

	You should have received a copy of the GNU General Public License
	along with JMSBench1.  If not, see <http://www.gnu.org/licenses/>.

----------------------------------------------------------------------- */

// original code copied from
// https://codereview.stackexchange.com/questions/187183/create-a-c-string-using-printf-style-formatting
// 2020 - Jan - 08

#include "string_format.h"
#include <cstdio>
#include <cstdarg>

#if __cplusplus < 201703L
#include <memory>
#endif

std::string format(const char *fmt, ...)
{
	char buf[STRING_BUF_LEN];

	va_list args;
	va_start(args, fmt);
	const auto r = std::vsnprintf(buf, sizeof buf, fmt, args);
	va_end(args);

	if (r < 0)
		// conversion failed
		return {};

	const size_t len = r;
	if (len < sizeof buf)
		// we fit in the buffer
		return { buf, len };

#if __cplusplus >= 201703L
	// C++17: Create a string and write to its underlying array
	std::string s(len, '\0');
	va_start(args, fmt);
	std::vsnprintf(s.data(), len + 1, fmt, args);
	va_end(args);

	return s;
#else
	// C++11 or C++14: We need to allocate scratch memory
	auto vbuf = std::unique_ptr<char[]>(new char[len + 1]);
	va_start(args, fmt);
	std::vsnprintf(vbuf.get(), len + 1, fmt, args);
	va_end(args);

	return { vbuf.get(), len };
#endif
}
