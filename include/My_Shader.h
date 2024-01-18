#pragma once

#ifndef MY_SHADER
#define MY_SHADER

#include "glad/glad.h"
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include "glm/glm.hpp"

class My_Shader
{
private:
	void checkCompileErrors(GLuint shader, const std::string& type);
public:
	//程序ID
	unsigned int ID;

	// Constructor:
	My_Shader(const std::string& vertexPath, const std::string& fragmentPath, const std::string& geoPath="");
	
	// 使用/激活程序
	void use();

	// uniform工具函数
	void setBool(const std::string& name, bool value) const;
	void setInt(const std::string& name, int value) const;
    void setuInt(const std::string& name, unsigned int value) const;
	void setFloat(const std::string& name, float value) const;
	void setVec2(const std::string& name, float x, float y) const;
	void setVec2(const std::string& name, const glm::vec2& value) const;
	void setVec3(const std::string& name, const glm::vec3& value) const;
	void setVec3(const std::string& name, float x, float y, float z) const;
	void setVec4(const std::string& name, float x, float y, float z, float w) const;
	void setVec4(const std::string& name, const glm::vec4& value) const;
	void setMat2(const std::string& name, const glm::mat2& mat) const;
	void setMat3(const std::string& name, const glm::mat3& mat) const;
	void setMat4(const std::string& name, const glm::mat4& mat) const;

};

#endif // MY_SHADER

