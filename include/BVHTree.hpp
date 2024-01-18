#pragma once

#ifndef __BVHTREE__
#define __BVHTREE__

#include "Triangle.hpp"
#include <glm/glm.hpp>
#include <vector>
#include <list>
#include <algorithm>
#include "Bounds.hpp"
#include <string>
#include <glad/glad.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include "tiny_obj_loader.h"
#include "model.hpp"

namespace rt {
	namespace BVH {
		using namespace std;
		using namespace glm;

		struct BVH_Node {
			vec3 AA, BB;
			int left = -1, right = -1;
			int n, index;
		};

		struct BVH_Node_encoded {
			vec3 childs; // (left, right , -1)
			vec3 leafInfo; // (n, index,  -1)
			vec3 AA, BB; // AAΪmin, BBΪmax
		};

		class BVHTree {
		private:
			vector<BVH_Node> nodes;
			vector<Triangle> triangles;
			float T_tri = 1.f;
			float T_aabb = 1.f;
			float INF = 1e9;
		public:
			BVHTree() {
				BVH_Node testNode;
				testNode.left = 255;
				testNode.right = 128;
				testNode.n = 30;
				testNode.AA = vec3(1, 1, 0);
				testNode.BB = vec3(0, 1, 0);
				nodes.push_back(testNode);
			}
			void loadModel(string const& filepath, Material const& material, mat4 const& trans, bool smoothNormal = false);
			int construct(int l, int r, int n);
			vector<BVH_Node>const& getNodes() const { return nodes; }
			vector<Triangle>const& getTriangles() const { return triangles; }
		};
	}
	void encodeBVH(vector<BVH::BVH_Node> const& nodes, vector<BVH::BVH_Node_encoded>& nodes_encoded);
	void sendTriangles2GPU(vector<Triangle_encoded> const&, GLuint&, GLuint&);
	void sendBVHNodes2GPU(vector<BVH::BVH_Node_encoded> const&, GLuint&, GLuint&);
}

void rt::BVH::BVHTree::loadModel(std::string const& filepath, rt::Material const& material, glm::mat4 const& trans, bool smoothNormal) {
	model::loadModel(filepath, triangles, material, trans, smoothNormal);
}


int rt::BVH::BVHTree::construct(int l, int r, int n) {
    if (l > r) return 0;
    nodes.push_back(BVH_Node());
    int id = nodes.size() - 1;
    nodes[id].left = nodes[id].right = nodes[id].n = nodes[id].index = 0;
    nodes[id].AA = vec3(1145141919, 1145141919, 1145141919);
    nodes[id].BB = vec3(-1145141919, -1145141919, -1145141919);
    // 计算 AABB
    for (int i = l; i <= r; i++) {
        // 最小点 AA
        float minx = std::min(triangles[i].p1.x, std::min(triangles[i].p2.x,
            triangles[i].p3.x));
        float miny = std::min(triangles[i].p1.y, std::min(triangles[i].p2.y,
            triangles[i].p3.y));
        float minz = std::min(triangles[i].p1.z, std::min(triangles[i].p2.z,
            triangles[i].p3.z));
        nodes[id].AA.x = std::min(nodes[id].AA.x, minx);
        nodes[id].AA.y = std::min(nodes[id].AA.y, miny);
        nodes[id].AA.z = std::min(nodes[id].AA.z, minz);
        // 最大点 BB
        float maxx = std::max(triangles[i].p1.x, std::max(triangles[i].p2.x,
            triangles[i].p3.x));
        float maxy = std::max(triangles[i].p1.y, std::max(triangles[i].p2.y,
            triangles[i].p3.y));
        float maxz = std::max(triangles[i].p1.z, std::max(triangles[i].p2.z,
            triangles[i].p3.z));
        nodes[id].BB.x = std::max(nodes[id].BB.x, maxx);
        nodes[id].BB.y = std::max(nodes[id].BB.y, maxy);
        nodes[id].BB.z = std::max(nodes[id].BB.z, maxz);
    }
    // 不多于 n 个三角形 返回叶子节点
    if ((r - l + 1) <= n) {
        nodes[id].n = r - l + 1;
        nodes[id].index = l;
        return id;
    }
    // 否则递归建树
    float Cost = INF;
    int Axis = 0;
    int Split = (l + r) / 2;
    for (int axis = 0; axis < 3; axis++) {
        // 分别按 x，y，z 轴排序
        if (axis == 0) std::sort(&triangles[0] + l, &triangles[0] + r + 1,
            cmp<X>());
        if (axis == 1) std::sort(&triangles[0] + l, &triangles[0] + r + 1,
            cmp<Y>());
        if (axis == 2) std::sort(&triangles[0] + l, &triangles[0] + r + 1,
            cmp<Z>());
        // leftMax[i]: [l, i] 中最大的 xyz 值
        // leftMin[i]: [l, i] 中最小的 xyz 值
        std::vector<vec3> leftMax(r - l + 1, vec3(-INF, -INF, -INF));
        std::vector<vec3> leftMin(r - l + 1, vec3(INF, INF, INF));
        // 计算前缀 注意 i-l 以对齐到下标 0
        for (int i = l; i <= r; i++) {
            Triangle& t = triangles[i];
            int bias = (i == l) ? 0 : 1;  // 第一个元素特殊处理
            leftMax[i - l].x = std::max(leftMax[i - l - bias].x, std::max(t.p1.x,
                std::max(t.p2.x, t.p3.x)));
            leftMax[i - l].y = std::max(leftMax[i - l - bias].y, std::max(t.p1.y,
                std::max(t.p2.y, t.p3.y)));
            leftMax[i - l].z = std::max(leftMax[i - l - bias].z, std::max(t.p1.z,
                std::max(t.p2.z, t.p3.z)));
            leftMin[i - l].x = std::min(leftMin[i - l - bias].x, std::min(t.p1.x,
                std::min(t.p2.x, t.p3.x)));
            leftMin[i - l].y = std::min(leftMin[i - l - bias].y, std::min(t.p1.y,
                std::min(t.p2.y, t.p3.y)));
            leftMin[i - l].z = std::min(leftMin[i - l - bias].z, std::min(t.p1.z,
                std::min(t.p2.z, t.p3.z)));
        }
        // rightMax[i]: [i, r] 中最大的 xyz 值
        // rightMin[i]: [i, r] 中最小的 xyz 值
        std::vector<vec3> rightMax(r - l + 1, vec3(-INF, -INF, -INF));
        std::vector<vec3> rightMin(r - l + 1, vec3(INF, INF, INF));
        // 计算后缀 注意 i-l 以对齐到下标 0
        for (int i = r; i >= l; i--) {
            Triangle& t = triangles[i];
            int bias = (i == r) ? 0 : 1;  // 第一个元素特殊处理
            rightMax[i - l].x = std::max(rightMax[i - l + bias].x, std::max(t.p1.x,
                std::max(t.p2.x, t.p3.x)));
            rightMax[i - l].y = std::max(rightMax[i - l + bias].y, std::max(t.p1.y,
                std::max(t.p2.y, t.p3.y)));
            rightMax[i - l].z = std::max(rightMax[i - l + bias].z, std::max(t.p1.z,
                std::max(t.p2.z, t.p3.z)));
            rightMin[i - l].x = std::min(rightMin[i - l + bias].x, std::min(t.p1.x,
                std::min(t.p2.x, t.p3.x)));
            rightMin[i - l].y = std::min(rightMin[i - l + bias].y, std::min(t.p1.y,
                std::min(t.p2.y, t.p3.y)));
            rightMin[i - l].z = std::min(rightMin[i - l + bias].z, std::min(t.p1.z,
                std::min(t.p2.z, t.p3.z)));
        }
        // 遍历寻找分割
        float cost = INF;
        int split = l;
        for (int i = l; i <= r - 1; i++) {
            float lenx, leny, lenz;
            // 左侧 [l, i]
            vec3 leftAA = leftMin[i - l];
            vec3 leftBB = leftMax[i - l];
            lenx = leftBB.x - leftAA.x;
            leny = leftBB.y - leftAA.y;
            lenz = leftBB.z - leftAA.z;
            float leftS = 2.0 * ((lenx * leny) + (lenx * lenz) + (leny * lenz));
            float leftCost = leftS * (i - l + 1);
            // 右侧 [i+1, r]
            vec3 rightAA = rightMin[i + 1 - l];
            vec3 rightBB = rightMax[i + 1 - l];
            lenx = rightBB.x - rightAA.x;
            leny = rightBB.y - rightAA.y;
            lenz = rightBB.z - rightAA.z;
            float rightS = 2.0 * ((lenx * leny) + (lenx * lenz) + (leny *
                lenz));
            float rightCost = rightS * (r - i);
            // 记录每个分割的最小答案
            float totalCost = leftCost + rightCost;
            if (totalCost < cost) {
                cost = totalCost;
                split = i;
            }
        }
        // 记录每个轴的最佳答案
        if (cost < Cost) {
            Cost = cost;
            Axis = axis;
            Split = split;
        }
    }
    // 按最佳轴分割
    if (Axis == 0) std::sort(&triangles[0] + l, &triangles[0] + r + 1, cmp<X>());
    if (Axis == 1) std::sort(&triangles[0] + l, &triangles[0] + r + 1, cmp<Y>());
    if (Axis == 2) std::sort(&triangles[0] + l, &triangles[0] + r + 1, cmp<Z>());
    // 递归
    int left = construct(l, Split, n);
    int right = construct(Split + 1, r, n);
    nodes[id].left = left;
    nodes[id].right = right;
    return id;
}

void rt::encodeBVH(std::vector<rt::BVH::BVH_Node> const& nodes, std::vector<rt::BVH::BVH_Node_encoded>& nodes_encoded) {
	nodes_encoded.resize(nodes.size());
	for (int i = 0; i < nodes.size(); i++) {
		nodes_encoded[i].childs = glm::vec3(nodes[i].left, nodes[i].right, 0);
		nodes_encoded[i].leafInfo = glm::vec3(nodes[i].n, nodes[i].index, 0);
		nodes_encoded[i].AA = nodes[i].AA;
		nodes_encoded[i].BB = nodes[i].BB;
	}
}

void rt::sendTriangles2GPU(std::vector<rt::Triangle_encoded> const& triangles_encoded, GLuint& tbo0, GLuint& trianglesTextureBuffer) {
	glGenBuffers(1, &tbo0);
	glBindBuffer(GL_TEXTURE_BUFFER, tbo0);
	glBufferData(GL_TEXTURE_BUFFER, triangles_encoded.size() *
		sizeof(Triangle_encoded), &triangles_encoded[0], GL_STATIC_DRAW);
	glGenTextures(1, &trianglesTextureBuffer);
	glBindTexture(GL_TEXTURE_BUFFER, trianglesTextureBuffer);
	glTexBuffer(GL_TEXTURE_BUFFER, GL_RGB32F, tbo0);
}

void rt::sendBVHNodes2GPU(std::vector<rt::BVH::BVH_Node_encoded> const& nodes_encoded, GLuint& tbo1, GLuint& nodesTextureBuffer) {
	glGenBuffers(1, &tbo1);
	glBindBuffer(GL_TEXTURE_BUFFER, tbo1);
	glBufferData(GL_TEXTURE_BUFFER, nodes_encoded.size() *
		sizeof(rt::BVH::BVH_Node_encoded), &nodes_encoded[0], GL_STATIC_DRAW);
	glGenTextures(1, &nodesTextureBuffer);
	glBindTexture(GL_TEXTURE_BUFFER, nodesTextureBuffer);
	glTexBuffer(GL_TEXTURE_BUFFER, GL_RGB32F, tbo1);
}

#endif // !__BVHTREE__