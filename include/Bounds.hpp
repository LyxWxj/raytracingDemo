

#ifndef __BOUNDS__
#define __BOUNDS__

#include <glm/glm.hpp>
#include <utility>
#include <glm/detail/type_vec.hpp>
#include <float.h>

namespace rt {
	using namespace glm;
	using namespace std;

	struct Bounds {
		glm::vec3 min, max;
		Bounds(glm::vec3 const& min, glm::vec3 const& max) : min(min), max(max) {};
		Bounds(glm::vec3&& min, glm::vec3&& max)noexcept: min(min), max(max) {};
		Bounds() :min(-FLT_MAX), max(FLT_MAX) {};
		glm::vec3 Center()const { return min + (max - min) / 2.f; }
	};

	bool operator&(const Bounds& a, const Bounds& b) {
		return a.min.x <= b.max.x && a.max.x >= b.min.x &&
			a.min.y <= b.max.y && a.max.y >= b.min.y &&
			a.min.z <= b.max.z && a.max.z >= b.min.z;
	}

	vec3 vec3_min(const vec3& a, const vec3& b) {
		return vec3(std::min(a.x, b.x), std::min(a.y, b.y), std::min(a.z, b.z));
	}

	vec3 vec3_max(const vec3& a, const vec3& b) {
		return vec3(std::max(a.x, b.x), std::max(a.y, b.y), std::max(a.z, b.z));
	}

	void mergeBounds(Bounds const& srcBoundA, Bounds const& srcBoundB, Bounds& dstBound) {
		dstBound.min = vec3_min(srcBoundA.min, srcBoundB.min);
		dstBound.max = vec3_max(srcBoundA.max, srcBoundB.max);
	}

	float surfaceArea(Bounds const& bound) {
		vec3 d = bound.max - bound.min;
		return 2.0f * (d.x * d.y + d.x * d.z + d.y * d.z);
	}
}

#endif // !__BOUNDS__