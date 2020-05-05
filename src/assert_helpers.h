/*
 * Copyright 2011, Ben Langmead <langmea@cs.jhu.edu>
 *
 * This file is part of Bowtie 2.
 *
 * Bowtie 2 is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Bowtie 2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Bowtie 2.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef ASSERT_HELPERS_H_
#define ASSERT_HELPERS_H_

#include <stdexcept>
#include <string>
#include <cassert>
#include <iostream>

/**
 * Assertion for release-enabled assertions
 */
class ReleaseAssertException : public std::runtime_error {
public:
	ReleaseAssertException(const std::string& msg = "") : std::runtime_error(msg) {}
};

/**
 * Macros for release-enabled assertions, and helper macros to make
 * all assertion error messages more helpful.
 */
#ifndef NDEBUG
#define ASSERT_ONLY(...) __VA_ARGS__
#else
#define ASSERT_ONLY(...)
#endif

#define rt_assert(b)  \
	if(!(b)) { \
		std::cerr << "rt_assert at " << __FILE__ << ":" << __LINE__ << std::endl; \
		throw ReleaseAssertException(); \
	}
#define rt_assert_msg(b,msg)  \
	if(!(b)) { \
		std::cerr << msg <<  " at " << __FILE__ << ":" << __LINE__ << std::endl; \
		throw ReleaseAssertException(msg); \
	}

#define rt_assert_eq(ex,ac)  \
	if(!((ex) == (ac))) { \
		std::cerr << "rt_assert_eq: expected (" << (ex) << ", 0x" << std::hex << (ex) << std::dec << ") got (" << (ac) << ", 0x" << std::hex << (ac) << std::dec << ")" << std::endl; \
		std::cerr << __FILE__ << ":" << __LINE__ << std::endl; \
		throw ReleaseAssertException(); \
	}
#define rt_assert_eq_msg(ex,ac,msg)  \
	if(!((ex) == (ac))) { \
		std::cerr << "rt_assert_eq: " << msg <<  ": (" << (ex) << ", 0x" << std::hex << (ex) << std::dec << ") got (" << (ac) << ", 0x" << std::hex << (ac) << std::dec << ")" << std::endl; \
		std::cerr << __FILE__ << ":" << __LINE__ << std::endl; \
		throw ReleaseAssertException(msg); \
	}

#ifndef NDEBUG
#define assert_eq(ex,ac)  \
	if(!((ex) == (ac))) { \
		std::cerr << "assert_eq: expected (" << (ex) << ", 0x" << std::hex << (ex) << std::dec << ") got (" << (ac) << ", 0x" << std::hex << (ac) << std::dec << ")" << std::endl; \
		std::cerr << __FILE__ << ":" << __LINE__ << std::endl; \
		assert(0); \
	}
#define assert_eq_msg(ex,ac,msg)  \
	if(!((ex) == (ac))) { \
		std::cerr << "assert_eq: " << msg <<  ": (" << (ex) << ", 0x" << std::hex << (ex) << std::dec << ") got (" << (ac) << ", 0x" << std::hex << (ac) << std::dec << ")" << std::endl; \
		std::cerr << __FILE__ << ":" << __LINE__ << std::endl; \
		assert(0); \
	}
#else
#define assert_eq(ex,ac)
#define assert_eq_msg(ex,ac,msg)
#endif

#define rt_assert_neq(ex,ac)  \
	if(!((ex) != (ac))) { \
		std::cerr << "rt_assert_neq: expected not (" << (ex) << ", 0x" << std::hex << (ex) << std::dec << ") got (" << (ac) << ", 0x" << std::hex << (ac) << std::dec << ")" << std::endl; \
		std::cerr << __FILE__ << ":" << __LINE__ << std::endl; \
		throw ReleaseAssertException(); \
	}
#define rt_assert_neq_msg(ex,ac,msg)  \
	if(!((ex) != (ac))) { \
		std::cerr << "rt_assert_neq: " << msg << ": (" << (ex) << ", 0x" << std::hex << (ex) << std::dec << ") got (" << (ac) << ", 0x" << std::hex << (ac) << std::dec << ")" << std::endl; \
		std::cerr << __FILE__ << ":" << __LINE__ << std::endl; \
		throw ReleaseAssertException(msg); \
	}

#ifndef NDEBUG
#define assert_neq(ex,ac)  \
	if(!((ex) != (ac))) { \
		std::cerr << "assert_neq: expected not (" << (ex) << ", 0x" << std::hex << (ex) << std::dec << ") got (" << (ac) << ", 0x" << std::hex << (ac) << std::dec << ")" << std::endl; \
		std::cerr << __FILE__ << ":" << __LINE__ << std::endl; \
		assert(0); \
	}
#define assert_neq_msg(ex,ac,msg)  \
	if(!((ex) != (ac))) { \
		std::cerr << "assert_neq: " << msg << ": (" << (ex) << ", 0x" << std::hex << (ex) << std::dec << ") got (" << (ac) << ", 0x" << std::hex << (ac) << std::dec << ")" << std::endl; \
		std::cerr << __FILE__ << ":" << __LINE__ << std::endl; \
		assert(0); \
	}
#else
#define assert_neq(ex,ac)
#define assert_neq_msg(ex,ac,msg)
#endif

#define rt_assert_gt(a,b) \
	if(!((a) > (b))) { \
		std::cerr << "rt_assert_gt: expected (" << (a) << ") > (" << (b) << ")" << std::endl; \
		std::cerr << __FILE__ << ":" << __LINE__ << std::endl; \
		throw ReleaseAssertException(); \
	}
#define rt_assert_gt_msg(a,b,msg) \
	if(!((a) > (b))) { \
		std::cerr << "rt_assert_gt: " << msg << ": (" << (a) << ") > (" << (b) << ")" << std::endl; \
		std::cerr << __FILE__ << ":" << __LINE__ << std::endl; \
		throw ReleaseAssertException(msg); \
	}

#ifndef NDEBUG
#define assert_gt(a,b) \
	if(!((a) > (b))) { \
		std::cerr << "assert_gt: expected (" << (a) << ") > (" << (b) << ")" << std::endl; \
		std::cerr << __FILE__ << ":" << __LINE__ << std::endl; \
		assert(0); \
	}
#define assert_gt_msg(a,b,msg) \
	if(!((a) > (b))) { \
		std::cerr << "assert_gt: " << msg << ": (" << (a) << ") > (" << (b) << ")" << std::endl; \
		std::cerr << __FILE__ << ":" << __LINE__ << std::endl; \
		assert(0); \
	}
#else
#define assert_gt(a,b)
#define assert_gt_msg(a,b,msg)
#endif

#define rt_assert_geq(a,b) \
	if(!((a) >= (b))) { \
		std::cerr << "rt_assert_geq: expected (" << (a) << ") >= (" << (b) << ")" << std::endl; \
		std::cerr << __FILE__ << ":" << __LINE__ << std::endl; \
		throw ReleaseAssertException(); \
	}
#define rt_assert_geq_msg(a,b,msg) \
	if(!((a) >= (b))) { \
		std::cerr << "rt_assert_geq: " << msg << ": (" << (a) << ") >= (" << (b) << ")" << std::endl; \
		std::cerr << __FILE__ << ":" << __LINE__ << std::endl; \
		throw ReleaseAssertException(msg); \
	}

#ifndef NDEBUG
#define assert_geq(a,b) \
	if(!((a) >= (b))) { \
		std::cerr << "assert_geq: expected (" << (a) << ") >= (" << (b) << ")" << std::endl; \
		std::cerr << __FILE__ << ":" << __LINE__ << std::endl; \
		assert(0); \
	}
#define assert_geq_msg(a,b,msg) \
	if(!((a) >= (b))) { \
		std::cerr << "assert_geq: " << msg << ": (" << (a) << ") >= (" << (b) << ")" << std::endl; \
		std::cerr << __FILE__ << ":" << __LINE__ << std::endl; \
		assert(0); \
	}
#else
#define assert_geq(a,b)
#define assert_geq_msg(a,b,msg)
#endif

#define rt_assert_lt(a,b) \
	if(!(a < b)) { \
		std::cerr << "rt_assert_lt: expected (" << a << ") < (" << b << ")" << std::endl; \
		std::cerr << __FILE__ << ":" << __LINE__ << std::endl; \
		throw ReleaseAssertException(); \
	}
#define rt_assert_lt_msg(a,b,msg) \
	if(!(a < b)) { \
		std::cerr << "rt_assert_lt: " << msg << ": (" << a << ") < (" << b << ")" << std::endl; \
		std::cerr << __FILE__ << ":" << __LINE__ << std::endl; \
		throw ReleaseAssertException(msg); \
	}

#ifndef NDEBUG
#define assert_lt(a,b) \
	if(!(a < b)) { \
		std::cerr << "assert_lt: expected (" << a << ") < (" << b << ")" << std::endl; \
		std::cerr << __FILE__ << ":" << __LINE__ << std::endl; \
		assert(0); \
	}
#define assert_lt_msg(a,b,msg) \
	if(!(a < b)) { \
		std::cerr << "assert_lt: " << msg << ": (" << a << ") < (" << b << ")" << std::endl; \
		std::cerr << __FILE__ << ":" << __LINE__ << std::endl; \
		assert(0); \
	}
#else
#define assert_lt(a,b)
#define assert_lt_msg(a,b,msg)
#endif

#define rt_assert_leq(a,b) \
	if(!((a) <= (b))) { \
		std::cerr << "rt_assert_leq: expected (" << (a) << ") <= (" << (b) << ")" << std::endl; \
		std::cerr << __FILE__ << ":" << __LINE__ << std::endl; \
		throw ReleaseAssertException(); \
	}
#define rt_assert_leq_msg(a,b,msg) \
	if(!((a) <= (b))) { \
		std::cerr << "rt_assert_leq: " << msg << ": (" << (a) << ") <= (" << (b) << ")" << std::endl; \
		std::cerr << __FILE__ << ":" << __LINE__ << std::endl; \
		throw ReleaseAssertException(msg); \
	}

#ifndef NDEBUG
#define assert_leq(a,b) \
	if(!((a) <= (b))) { \
		std::cerr << "assert_leq: expected (" << (a) << ") <= (" << (b) << ")" << std::endl; \
		std::cerr << __FILE__ << ":" << __LINE__ << std::endl; \
		assert(0); \
	}
#define assert_leq_msg(a,b,msg) \
	if(!((a) <= (b))) { \
		std::cerr << "assert_leq: " << msg << ": (" << (a) << ") <= (" << (b) << ")" << std::endl; \
		std::cerr << __FILE__ << ":" << __LINE__ << std::endl; \
		assert(0); \
	}
#else
#define assert_leq(a,b)
#define assert_leq_msg(a,b,msg)
#endif

#ifndef NDEBUG
#define assert_in(c, s) assert_in2(c, s, __FILE__, __LINE__)
static inline void assert_in2(char c, const char *str, const char *file, int line) {
	const char *s = str;
	while(*s != '\0') {
		if(c == *s) return;
		s++;
	}
	std::cerr << "assert_in: (" << c << ") not in  (" << str << ")" << std::endl;
	std::cerr << file << ":" << line << std::endl;
	assert(0);
}
#else
#define assert_in(c, s)
#endif

#ifndef NDEBUG
#define assert_range(b, e, v) assert_range_helper(b, e, v, __FILE__, __LINE__)
template<typename T>
inline static void assert_range_helper(const T& begin,
                                       const T& end,
                                       const T& val,
                                       const char *file,
                                       int line)
{
	if(val < begin || val > end) {
		std::cerr << "assert_range: (" << val << ") not in  ["
		          << begin << ", " << end << "]" << std::endl;
		std::cerr << file << ":" << line << std::endl;
		assert(0);
	}
}
#else
#define assert_range(b, e, v)
#endif

// define a macro to indicate variables that are only required for asserts
// used to make production build happy, i.e. disable "warning: variable ‘x’ set but not used [-Wunused-but-set-variable]"
#define _unused(x) ((void)x)

#endif /*ASSERT_HELPERS_H_*/

