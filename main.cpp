#include "Bounds.hpp"
#include "BVHTree.hpp"
#include "HDRreader.hpp"
#include "Triangle.hpp"
#include "RenderPass.hpp""
#include <vector>
#include <string>
#include <GLFW/glfw3.h>
#include "My_shader.h"
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include "Camera_Movement.h"
#include <iostream>
#include <iomanip>
#include <time.h>
#include "Init.h"

glm::mat4 getTransformMatrix(glm::vec3 const& rotateCtrl, glm::vec3 const& translateCtrl, glm::vec3 const& scaleCtrl);

GLuint tbo0, trianglesTextureBuffer;
GLuint tbo1, nodesTextureBuffer;
GLuint lastFrame;
GLuint hdrMap;
GLuint hdrCache;
int hdrResolution;

double lastX = 0.0, lastY = 0.0;

clock_t t1, t2;
double dt, fps;
unsigned int frameCounter = 0;

float upAngle = 0.0;
float rotatAngle = 0.0;
float r = 3.0;
int maxBounce = 1;

float lastupAngle = 0.0;
float lastrotatAngle = 0.0;
float lastr = 3.0;

extern Camera camera;

int width = 600, height = 600;

rt::RenderPass pass1(width, height);
rt::RenderPass pass2(width, height);
rt::RenderPass pass3(width, height);


int main() {

    GLFWwindow* window = init(width, height, "Ray Tracing");

    rt::Material m;
    m.roughness = 0.5;
    m.specular = 1.0;
    m.metallic = 1.0;
    m.clearcoat = 1.0;
    m.clearcoatGloss = 0.0;
    m.baseColor = glm::vec3(1, 0.73, 0.25);
    rt::BVH::BVHTree BVHtree;
    BVHtree.loadModel("resources//stanford bunny.obj", m, getTransformMatrix(glm::vec3(0, 0, 0), glm::vec3(0.3, -1.6, 0), glm::vec3(1.5)), true);
    m.roughness = 0.01;
    m.metallic = 0.1;
    m.specular = 1.0;
    m.baseColor = glm::vec3(1, 1, 1);
    BVHtree.loadModel("resources//quad.obj", m,  getTransformMatrix(glm::vec3(0, 0, 0), glm::vec3(0, -1.4, 0), glm::vec3(18.83, 0.01, 18.83)));
    m.baseColor = glm::vec3(1, 1, 1);
    m.emissive = glm::vec3(20, 20, 20);
    BVHtree.loadModel("resources//sphere.obj", m, getTransformMatrix(glm::vec3(0, 0, 0), glm::vec3(0.0, 0.9, -0.0), glm::vec3(1, 1, 1)));
    BVHtree.construct(0, BVHtree.getTriangles().size() - 1, 8);
    std::cout << "模型读取完成: 共 " << BVHtree.getTriangles().size() << " 个三角形" << std::endl;
    std::vector<rt::Triangle_encoded> Triangle_encoded;
    rt::encodeTriangle(BVHtree.getTriangles(), Triangle_encoded);
    std::vector<rt::BVH::BVH_Node_encoded> BVHNode_encoded;
    rt::encodeBVH(BVHtree.getNodes(), BVHNode_encoded);
    rt::sendTriangles2GPU(Triangle_encoded, tbo0, trianglesTextureBuffer);
    rt::sendBVHNodes2GPU(BVHNode_encoded, tbo1, nodesTextureBuffer);

    HDRLoaderResult hdrRes;
    hdrMap = rt::loadHDR("HDR//chinese_garden_2k.hdr", hdrRes);
    float* cache = rt::calculateHdrCache(hdrRes.cols, hdrRes.width, hdrRes.height);
    hdrCache = rt::getTextureRGB32F(hdrRes.width, hdrRes.height);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB32F, hdrRes.width, hdrRes.height, 0, GL_RGB, GL_FLOAT, cache);
    hdrResolution = hdrRes.width;

    My_Shader shaderPass1("shaders//rt.vs", "shaders//pass1new.fs");
    pass1.program = shaderPass1.ID;
    pass1.colorAttachments.push_back(rt::getTextureRGB32F(pass1.width, pass1.height));
    pass1.colorAttachments.push_back(rt::getTextureRGB32F(pass1.width, pass1.height));
    pass1.colorAttachments.push_back(rt::getTextureRGB32F(pass1.width, pass1.height));
    pass1.bindData();

    shaderPass1.use();
    shaderPass1.setInt("nTriangles", (int)BVHtree.getTriangles().size());
    shaderPass1.setInt("nNodes", (int)BVHtree.getNodes().size());
    shaderPass1.setInt("width", height);
    shaderPass1.setInt("height", width);

    My_Shader shaderPass2("shaders//rt.vs", "shaders//pass2.fs");
    pass2.program = shaderPass2.ID;
    lastFrame = rt::getTextureRGB32F(pass2.width, pass2.height);
    pass2.colorAttachments.push_back(lastFrame);
    pass2.bindData();

    My_Shader shaderPass3("shaders//rt.vs", "shaders//pass3.fs");
    pass3.program = shaderPass3.ID;
    pass3.bindData(true);

    glEnable(GL_DEPTH_TEST);
    extern Camera camera;
    clock_t t1 = clock(), t2;
    double dt, fps;
    int frameCounter = 0;
    std::cout<<BVHtree.getTriangles().size()<<std::endl;
    std::cout<<BVHtree.getNodes().size()<<std::endl;       
    ImGui::StyleColorsLight();
    while (!glfwWindowShouldClose(window)) {
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        glClearColor(0.0, 0.0, 0.0, 1.0);

        t2 = clock();
        dt = (double)(t2 - t1) / CLOCKS_PER_SEC;
        fps = 1.0 / dt;
        std::cout << "\r";
        std::cout << std::fixed << std::setprecision(2) << "FPS : " << fps << "    迭代次数: " << frameCounter;
        t1 = t2;
        glm::vec3 eye = glm::vec3(-sin(glm::radians(rotatAngle)) * cos(glm::radians(upAngle)),
            sin(glm::radians(upAngle)), cos(glm::radians(rotatAngle)) * cos(glm::radians(upAngle)));
        eye.x *= r; eye.y *= r; eye.z *= r;
        glm::mat4 cameraRotate = glm::lookAt(eye, glm::vec3(0, 0, 0), glm::vec3(0, 1, 0));  // 相机注视着原点
        cameraRotate = inverse(cameraRotate);   // lookat 的逆矩阵将光线方向进行转换
        shaderPass1.use();
        shaderPass1.setInt("maxBounce", maxBounce);
        shaderPass1.setVec3("eye", eye);
        shaderPass1.setMat4("cameraRotate", cameraRotate);
        shaderPass1.setuInt("frameCounter", frameCounter++);
        shaderPass1.setInt("hdrResolution", hdrResolution);
        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_BUFFER, trianglesTextureBuffer);
        shaderPass1.setInt("triangles", 0);
        glActiveTexture(GL_TEXTURE1);
        glBindTexture(GL_TEXTURE_BUFFER, nodesTextureBuffer);
        shaderPass1.setInt("nodes", 1);
        glActiveTexture(GL_TEXTURE2);
        glBindTexture(GL_TEXTURE_2D, lastFrame);
        shaderPass1.setInt("lastFrame", 2);
        glActiveTexture(GL_TEXTURE3);
        glBindTexture(GL_TEXTURE_2D, hdrMap);
        shaderPass1.setInt("hdrMap", 3);
        glActiveTexture(GL_TEXTURE4);
        glBindTexture(GL_TEXTURE_2D, hdrCache);
        glUniform1i(glGetUniformLocation(pass1.program, "hdrCache"), 4);
        
        pass1.draw();
        pass2.draw(pass1.colorAttachments);
        pass3.draw(pass2.colorAttachments);

        /*-----------------------------I M G U I------------------------------------------------------*/
        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();

        ImGui::Begin("Control Pannel");

        ImGui::DragInt("maxBounce", &maxBounce, 0.05, 1, 10);
        ImGui::DragFloat("UpAngle", &upAngle, 0.5, -90, 89);
        ImGui::DragFloat("RotateAngle", &rotatAngle, 0.5, -180, 179);
        ImGui::DragFloat("r", &r, 0.05, 0.1, 10);

        ImGui::End();
        ImGui::Render();
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
        if (lastupAngle != upAngle || lastrotatAngle != rotatAngle || lastr != r) {
			frameCounter = 0;
		}
        lastupAngle = upAngle;
        lastrotatAngle = rotatAngle;
        lastr = r;
        glfwSwapBuffers(window);
        glfwPollEvents();
    }
    glfwTerminate();
    delete[] cache;
    return 0;
}

glm::mat4 getTransformMatrix(glm::vec3 const& rotateCtrl, glm::vec3 const& translateCtrl, glm::vec3 const& scaleCtrl) {
    glm::mat4 unit(    // 单位矩阵
        glm::vec4(1, 0, 0, 0),
        glm::vec4(0, 1, 0, 0),
        glm::vec4(0, 0, 1, 0),
        glm::vec4(0, 0, 0, 1)
    );
    glm::mat4 scale = glm::scale(unit, scaleCtrl);
    glm::mat4 translate = glm::translate(unit, translateCtrl);
    glm::mat4 rotate = unit;
    rotate = glm::rotate(rotate, glm::radians(rotateCtrl.x), glm::vec3(1, 0,
        0));
    rotate = glm::rotate(rotate, glm::radians(rotateCtrl.y), glm::vec3(0, 1,
        0));
    rotate = glm::rotate(rotate, glm::radians(rotateCtrl.z), glm::vec3(0, 0,
        1));
    glm::mat4 model = translate * rotate * scale;
    return model;
}