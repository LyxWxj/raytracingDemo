#pragma once

#ifndef __TRIANGLE__
#define __TRIANGLE__

#include <glm/glm.hpp>
#include <vector>

namespace rt {
	enum Axis {X = 0, Y = 1, Z = 2};
	using namespace std;
	using namespace glm;

	struct Material {
		vec3 emissive = vec3(0, 0, 0);  // 作为光源时的发光颜色
		vec3 baseColor = vec3(1, 1, 1);
		float subsurface = 0.0;
		float metallic = 0.0;
		float specular = 0.5;
		float specularTint = 0.0;
		float roughness = 0.5;
		float anisotropic = 0.0;
		float sheen = 0.0;
		float sheenTint = 0.5;
		float clearcoat = 0.0;
		float clearcoatGloss = 1.0;
		float IOR = 1.0;
		float transmission = 0.0;
		float tex_id;
	};

	struct Triangle {
	public:
		vec3 p1, p2, p3;
		vec3 n1, n2, n3;
		Material material;
		Triangle() = default;
		Triangle(vec3 const& a, vec3 const& b, vec3 const& c, vec3 const& d, vec3 const& e, vec3 const& f, Material const& m):
			p1(a), p2(b), p3(c), n1(d), n2(e), n3(f), material(m) {};
		Triangle(vec3&& a, vec3&& b, vec3&& c, vec3&& d, vec3&& e, vec3&& f, Material&& m) noexcept:
			p1(a), p2(b), p3(c) ,n1(d), n2(e), n3(f), material(m){};
		vec3 center() const {
			return (p1 + p2 + p3) / 3.0f;
		}
	};


	template<Axis axis>
	class cmp {
	public:
		bool operator()(const Triangle& t1, const Triangle& t2) {
			return t1.center()[axis] < t2.center()[axis];
		}	
	};

	struct Triangle_encoded {
		vec3 p1, p2, p3;
		vec3 n1, n2, n3;
		vec3 emissive;
		vec3 baseColor;
		vec3 subsurface_metallic_specular; // (subsurface, metallic, specular)
		vec3 specularTint_roughness_anisotropic; // (specularTint, roughness, anisotropic)
		vec3 sheen_sheenTint_clearcoat; // (sheen, sheenTint, clearcoat)
		vec3 clearcoatGloss_IOR_transmission; // (clearcoatGloss, IOR, transmission)
	};

	void encodeTriangle(vector<Triangle> const& triangles, vector<Triangle_encoded>&);
}

void rt::encodeTriangle(std::vector<Triangle> const& triangles, std::vector <Triangle_encoded>& triangles_encoded) {
	size_t nTriangles = triangles.size();
	triangles_encoded.resize(nTriangles);
	for (int i = 0; i < nTriangles; ++i) {
		auto& t = triangles[i];
		auto& m = t.material;

		triangles_encoded[i].p1 = t.p1;
		triangles_encoded[i].p2 = t.p2;
		triangles_encoded[i].p3 = t.p3;
		
		triangles_encoded[i].n1 = t.n1;
		triangles_encoded[i].n2 = t.n2;
		triangles_encoded[i].n3 = t.n3;

		triangles_encoded[i].emissive = m.emissive;
		triangles_encoded[i].baseColor = m.baseColor;	
		triangles_encoded[i].subsurface_metallic_specular = vec3(m.subsurface, m.metallic, m.specular);

		triangles_encoded[i].specularTint_roughness_anisotropic = vec3(m.specularTint, m.roughness, m.anisotropic);
		triangles_encoded[i].sheen_sheenTint_clearcoat = vec3(m.sheen, m.sheenTint, m.clearcoat);
		triangles_encoded[i].clearcoatGloss_IOR_transmission = vec3(m.clearcoatGloss, m.IOR, m.transmission);
	}
}

#endif // !TRIANGLE