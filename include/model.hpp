#pragma once
#ifndef __MODELREADER__
#define __MODELREADER__

#include "tiny_obj_loader.h"
#include "tiny_gltf.h"
#include <string>
#include <vector>
#include <glm/glm.hpp>
#include "Triangle.hpp"
#include "Bounds.hpp"

namespace model{
    using namespace std;
    using namespace rt;

	typedef enum modelType { ERR = -1, OBJ = 0, GLTF = 1 } Type;

	void make_TrianglesFromOBJ(vector<Triangle>& triangles, vector<tinyobj::material_t> const& materials, vector<tinyobj::shape_t> const& shapes, tinyobj::attrib_t const& attrib, glm::mat4 const& trans, Material m, bool smooth) {
        vec3 maxvec, minvec;
        for (auto const& shape : shapes) {
			for (size_t f = 0; f < shape.mesh.num_face_vertices.size(); ++f) {
				// 三角形
				if (shape.mesh.num_face_vertices[f] == 3) {
					// 顶点索引
					tinyobj::index_t idx0 = shape.mesh.indices[3 * f + 0];
					tinyobj::index_t idx1 = shape.mesh.indices[3 * f + 1];
					tinyobj::index_t idx2 = shape.mesh.indices[3 * f + 2];
					// 顶点坐标
					glm::vec3 v0 = glm::vec3(attrib.vertices[3 * idx0.vertex_index + 0], attrib.vertices[3 * idx0.vertex_index + 1], attrib.vertices[3 * idx0.vertex_index + 2]);
					glm::vec3 v1 = glm::vec3(attrib.vertices[3 * idx1.vertex_index + 0], attrib.vertices[3 * idx1.vertex_index + 1], attrib.vertices[3 * idx1.vertex_index + 2]);
					glm::vec3 v2 = glm::vec3(attrib.vertices[3 * idx2.vertex_index + 0], attrib.vertices[3 * idx2.vertex_index + 1], attrib.vertices[3 * idx2.vertex_index + 2]);
                    v0 = vec3(trans * vec4(v0, 1.f));
                    v1 = vec3(trans * vec4(v1, 1.f));
                    v2 = vec3(trans * vec4(v2, 1.f));
                    // 法线坐标
					glm::vec3 n0, n1, n2;
					if (attrib.normals.size() || (!smooth)) {
						n0 = glm::vec3(attrib.normals[3 * idx0.normal_index + 0], attrib.normals[3 * idx0.normal_index + 1], attrib.normals[3 * idx0.normal_index + 2]);
						n1 = glm::vec3(attrib.normals[3 * idx1.normal_index + 0], attrib.normals[3 * idx1.normal_index + 1], attrib.normals[3 * idx1.normal_index + 2]);
						n2 = glm::vec3(attrib.normals[3 * idx2.normal_index + 0], attrib.normals[3 * idx2.normal_index + 1], attrib.normals[3 * idx2.normal_index + 2]);
					}
					else {
						// n1 = n2 = n0 = glm::cross(v1 - v0, v2 - v0);
					}
					// 材质
					if (materials.size()) {
						tinyobj::material_t material = materials[shape.mesh.material_ids[f]];
						m.baseColor = glm::vec3(material.diffuse[0], material.diffuse[1], material.diffuse[2]);
						m.specular = (material.specular[0] + material.specular[1] + material.specular[2]) / 3.f;
					}
                    maxvec = vec3_max(v0, vec3_max(v1, v2));
                    minvec = vec3_min(v0, vec3_min(v1, v2));

					triangles.push_back(Triangle(v0, v1, v2, n0, n1, n2, m));
				}
			}
		}
        vec3 len = maxvec - minvec;
        float maxaxis = std::max(len.x, std::max(len.y, len.z));
        for (auto& triangle : triangles) {
            auto& v0 = triangle.p1, v1 = triangle.p2, v2 = triangle.p3;
            v0.x /= maxaxis;
            v1 /= maxaxis;
            v2 /= maxaxis;
            /*auto& n0 = triangle.n1, n1 = triangle.n2, n2 = triangle.n3;
            if (smooth) {
                vec3 n = normalize(glm::cross(v1 - v0, v2 - v1));
                n0 = n1 = n2 = n;
            } else {
                n0 = normalize(n0);
                n1 = normalize(n1);
                n2 = normalize(n2);
            }*/
        }
	}

	void myloadOBJ(string const& filepath, vector<Triangle>& triangles, rt::Material const& material, glm::mat4 const& trans, bool smoothNormal = false) {
        // 顶点位置，索引
        std::vector<vec3> vertices;
        std::vector<GLuint> indices;
        // 打开文件流
        std::ifstream fin(filepath);
        std::string line;
        if (!fin.is_open()) {
            std::cout << "文件 " << filepath << " 打开失败" << std::endl;
            exit(-1);
        }
        // 计算 AABB 盒，归一化模型大小
        float maxx = -11451419.19;
        float maxy = -11451419.19;
        float maxz = -11451419.19;
        float minx = 11451419.19;
        float miny = 11451419.19;
        float minz = 11451419.19;
        // 按行读取
        while (std::getline(fin, line)) {
            std::istringstream sin(line);   // 以一行的数据作为 string stream 解析并且读取
            std::string type;
            GLfloat x, y, z;
            int v0, v1, v2;
            int vn0, vn1, vn2;
            int vt0, vt1, vt2;
            char slash;
            // 统计斜杆数目，用不同格式读取
            int slashCnt = 0;
            for (int i = 0; i < line.length(); i++) {
                if (line[i] == '/') slashCnt++;
            }
            // 读取obj文件
            sin >> type;
            if (type == "v") {
                sin >> x >> y >> z;
                vertices.push_back(vec3(x, y, z));
                maxx = std::max(maxx, x); maxy = std::max(maxx, y); maxz = std::max(maxx, z);
                minx = std::min(minx, x); miny = std::min(minx, y); minz = std::min(minx, z);
            }
            if (type == "f") {
                if (slashCnt == 6) {
                    sin >> v0 >> slash >> vt0 >> slash >> vn0;
                    sin >> v1 >> slash >> vt1 >> slash >> vn1;
                    sin >> v2 >> slash >> vt2 >> slash >> vn2;
                }
                else if (slashCnt == 3) {
                    sin >> v0 >> slash >> vt0;
                    sin >> v1 >> slash >> vt1;
                    sin >> v2 >> slash >> vt2;
                }
                else {
                    sin >> v0 >> v1 >> v2;
                }
                indices.push_back(v0 - 1);
                indices.push_back(v1 - 1);
                indices.push_back(v2 - 1);
            }
        }
        // 模型大小归一化
        float lenx = maxx - minx;
        float leny = maxy - miny;
        float lenz = maxz - minz;
        float maxaxis = std::max(lenx, std::max(leny, lenz));
        for (auto& v : vertices) {
            v.x /= maxaxis;
            v.y /= maxaxis;
            v.z /= maxaxis;
        }
        // 通过矩阵进行坐标变换
        for (auto& v : vertices) {
            vec4 vv = vec4(v.x, v.y, v.z, 1);
            vv = trans * vv;
            v = vec3(vv.x, vv.y, vv.z);
        }
        // 生成法线
        std::vector<vec3> normals(vertices.size(), vec3(0, 0, 0));
        for (int i = 0; i < indices.size(); i += 3) {
            vec3 p1 = vertices[indices[i]];
            vec3 p2 = vertices[indices[i + 1]];
            vec3 p3 = vertices[indices[i + 2]];
            vec3 n = normalize(cross(p2 - p1, p3 - p1));
            normals[indices[i]] += n;
            normals[indices[i + 1]] += n;
            normals[indices[i + 2]] += n;
        }
        // 构建 Triangle 对象数组
        int offset = triangles.size();  // 增量更新
        triangles.resize(offset + indices.size() / 3);
        for (int i = 0; i < indices.size(); i += 3) {
            Triangle& t = triangles[offset + i / 3];
            // 传顶点属性
            t.p1 = vertices[indices[i]];
            t.p2 = vertices[indices[i + 1]];
            t.p3 = vertices[indices[i + 2]];
            if (!smoothNormal) {
                vec3 n = normalize(cross(t.p2 - t.p1, t.p3 - t.p1));
                t.n1 = n; t.n2 = n; t.n3 = n;
            }
            else {
                t.n1 = normalize(normals[indices[i]]);
                t.n2 = normalize(normals[indices[i + 1]]);
                t.n3 = normalize(normals[indices[i + 2]]);
            }
            // 传材质
            t.material = material;
        }
	}

	void loadOBJ(string const& filepath, vector<Triangle>& triangles, rt::Material const& material, glm::mat4 const& trans, bool smoothNormal) {
		vector<tinyobj::material_t> materials;
		vector<tinyobj::shape_t> shapes;
		// 解析obj文件和mtl文件
		tinyobj::attrib_t attrib;
		std::string err;
		std::string warn;
		std::string dir = filepath.substr(0, filepath.find_last_of('/') + 1);
		bool ret = tinyobj::LoadObj(&attrib, &shapes, &materials, &warn, &err, filepath.c_str(), dir.c_str());
		if (!ret) {
			cerr << "Failed to load obj file: " << err << endl;
			exit(1);
		}
		make_TrianglesFromOBJ(triangles, materials, shapes, attrib, trans, material, smoothNormal);
	}

    void loadGLTF(string const& filepath, vector<Triangle>& triangles, rt::Material const& material, glm::mat4 const& trans) {
        //// LoadModel From GLTF into trangles with material by tiny_gltf.cc
        //tinygltf::Model model;
        //tinygltf::TinyGLTF loader;
        //std::string err;
        //std::string warn;
        //bool ret = loader.LoadASCIIFromFile(&model, &err, &warn, filepath);
        //if (!warn.empty()) {
        //    std::cout << "Warn: " << warn << std::endl;
        //}
        //if (!err.empty()) {
        //    std::cout << "Err: " << err << std::endl;
        //}
        //if (!ret) {
        //    std::cout << "Failed to load glTF: " << filepath << std::endl;
        //    exit(-1);
        //}
        //// 读取模型
        //for (auto const& mesh : model.meshes) {
        //    // 将mesh中的每个primitive转换为三角形
        //    for (auto const& primitive : mesh.primitives) {
        //        // 顶点坐标
        //        auto const& posAccessor = model.accessors[primitive.attributes.find("POSITION")->second];
        //        auto const& posBufferView = model.bufferViews[posAccessor.bufferView];
        //        auto const& posBuffer = model.buffers[posBufferView.buffer];
        //        auto const& posData = posBuffer.data;
        //        // 法线坐标
        //        auto const& norAccessor = model.accessors[primitive.attributes.find("NORMAL")->second];
        //        auto const& norBufferView = model.bufferViews[norAccessor.bufferView];
        //        auto const& norBuffer = model.buffers[norBufferView.buffer];
        //        auto const& norData = norBuffer.data;
        //        // 纹理坐标
        //        auto const& texAccessor = model.accessors[primitive.attributes.find("TEXCOORD_0")->second];
        //        auto const& texBufferView = model.bufferViews[texAccessor.bufferView];
        //        auto const& texBuffer = model.buffers[texBufferView.buffer];
        //        auto const& texData = texBuffer.data;
        //        // 索引
        //        auto const& idxAccessor = model.accessors[primitive.indices];
        //        auto const& idxBufferView = model.bufferViews[idxAccessor.bufferView];
        //        auto const& idxBuffer = model.buffers[idxBufferView.buffer];
        //        auto const& idxData = idxBuffer.data;
        //        // 材质
        //        auto const& material = model.materials[primitive.material];
        //        // 读取顶点
        //        for (int i = 0; i < posAccessor.count; i++) {
        //            // 顶点坐标
        //            float x = *(float*)&posData[posAccessor.byteOffset + posBufferView.byteOffset + i * posAccessor.ByteStride(posBufferView)];
        //            float y = *(float*)&posData[posAccessor.byteOffset + posBufferView.byteOffset + i * posAccessor.ByteStride(posBufferView) + 4];
        //            float z = *(float*)&posData[posAccessor.byteOffset + posBufferView.byteOffset + i * posAccessor.ByteStride(posBufferView) + 8];
        //            
        //        }

        //    }
        // }
    }
    modelType ModelType(string const& extension) {
        if (extension == ".obj" || extension == ".OBJ") return OBJ;
        else if (extension == ".gltf" || extension == ".GLTF") return GLTF;
        else return ERR;
    }

    void loadModel(string const& filepath, vector<Triangle>& triangles, rt::Material const& material, glm::mat4 const& trans, bool smoothNormal = false) {
        string extension = filepath.substr(filepath.find_last_of('.'));
        Type etype = ModelType(extension);
        assert(etype != ERR);
        switch (etype) {
        case OBJ:
			myloadOBJ(filepath, triangles, material, trans, smoothNormal);
            break;
        /*case GLTF:
			loadGLTF(filepath, triangles, material, trans);
            break;*/
        default:
            break;
        }
    }
}

#endif __MODELREADER__
