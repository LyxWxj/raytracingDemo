#pragma once

#ifndef __INIT__
#define __INIT__

#include "GL/glfw3.h"
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <filesystem>
#include "stb_image.h"

unsigned int loadTexture(char const* path)
{
    unsigned int textureID;
    glGenTextures(1, &textureID);

    int width, height, nrComponents;
    unsigned char* data = stbi_load(path, &width, &height, &nrComponents, 0);
    if (data)
    {
        GLenum format;
        if (nrComponents == 1)
            format = GL_RED;
        else if (nrComponents == 3)
            format = GL_RGB;
        else if (nrComponents == 4)
            format = GL_RGBA;

        glBindTexture(GL_TEXTURE_2D, textureID);
        glTexImage2D(GL_TEXTURE_2D, 0, format, width, height, 0, format, GL_UNSIGNED_BYTE, data);
        glGenerateMipmap(GL_TEXTURE_2D);

        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

        stbi_image_free(data);
    }
    else
    {
        std::cout << "Texture failed to load at path: " << path << std::endl;
        stbi_image_free(data);
    }

    return textureID;
}

unsigned int loadCubemap(std::vector<std::string> faces)
{
    unsigned int textureID;
    glGenTextures(1, &textureID);
    glBindTexture(GL_TEXTURE_2D_MULTISAMPLE, textureID);
    glTexImage2DMultisample(GL_TEXTURE_2D_MULTISAMPLE, 4, GL_RGB, 800, 600, GL_TRUE);
    glBindTexture(GL_TEXTURE_2D_MULTISAMPLE, 0);

    int width, height, nrChannels;
    for (unsigned int i = 0; i < faces.size(); i++)
    {
        unsigned char* data = stbi_load(faces[i].c_str(), &width, &height, &nrChannels, 0);
        if (data)
        {
            glTexImage2D(GL_TEXTURE_CUBE_MAP_POSITIVE_X + i,
                0, GL_RGB, width, height, 0, GL_RGB, GL_UNSIGNED_BYTE, data
            );
            stbi_image_free(data);
        }
        else
        {
            std::cout << "Cubemap texture failed to load at path: " << faces[i] << std::endl;
            stbi_image_free(data);
        }
    }
    glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);

    return textureID;
}


#endif 


#ifndef __IMGUIINIT__
#define __IMGUIINIT__

#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <string>
#include <sstream>
#include <imgui/imgui.h>
#include <imgui/imgui_impl_glfw.h>
#include <imgui/imgui_impl_opengl3.h>
#include <imgui/imstb_rectpack.h>
#include <imgui/imgui_internal.h>
#include <imgui/imconfig.h>
#include <imgui/imgui_impl_opengl3_loader.h>
#include <imgui/imstb_truetype.h>


GLFWwindow* init(int width, int height, const char* title) {
    std::string font = "C:\\Windows\\Fonts\\";
    glfwGetError(0);
    glfwInit();
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    // glfwWindowHint(GLFW_TRANSPARENT_FRAMEBUFFER, GLFW_TRUE);
    // 需要在创建窗口前设置标志
    // Create window with graphics context
    GLFWwindow* Window = glfwCreateWindow(width, height, title, NULL, NULL);

    glfwMakeContextCurrent(Window);
    glfwSwapInterval(0);

    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
    {
        std::cout << "Failed to initialize GLAD" << std::endl;
        return nullptr;
    }

    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGuiIO& io = ImGui::GetIO();
    (void)io;   // 避免 unused 变量的警告
    io.AddMouseViewportEvent(0);
    io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;       // Enable Keyboard Controls
    io.ConfigFlags |= ImGuiConfigFlags_NavEnableGamepad;      // Enable Gamepad Controls
    //io.ConfigFlags |= ImGuiConfigFlags_DockingEnable;           // !!! 启用 docking 功能的支持
    //io.ConfigFlags |= ImGuiConfigFlags_ViewportsEnable;   // !!! 启用 viewport 功能的支持

    io.Fonts->AddFontFromFileTTF((font + "consola.ttf").c_str(), 20, NULL, io.Fonts->GetGlyphRangesChineseFull());
    io.Fonts->AddFontFromFileTTF((font + "Candara.ttf").c_str(), 22, NULL, io.Fonts->GetGlyphRangesChineseFull());
    io.Fonts->AddFontFromFileTTF((font + "Arial.ttf").c_str(), 22, NULL, io.Fonts->GetGlyphRangesChineseFull());
    io.Fonts->AddFontFromFileTTF("D:\\imgui\\imgui-master\\misc\\fonts\\Cousine-Regular.ttf", 22, NULL, io.Fonts->GetGlyphRangesChineseFull());
    io.Fonts->AddFontFromFileTTF("D:\\imgui\\imgui-master\\misc\\fonts\\DroidSans.ttf", 22, NULL, io.Fonts->GetGlyphRangesChineseFull());
    io.Fonts->AddFontFromFileTTF("D:\\imgui\\imgui-master\\misc\\fonts\\Karla-Regular.ttf", 22, NULL, io.Fonts->GetGlyphRangesChineseFull());
    io.Fonts->AddFontFromFileTTF("D:\\imgui\\imgui-master\\misc\\fonts\\ProggyClean.ttf", 22, NULL, io.Fonts->GetGlyphRangesChineseFull());
    io.Fonts->AddFontFromFileTTF("D:\\imgui\\imgui-master\\misc\\fonts\\ProggyTiny.ttf", 22, NULL, io.Fonts->GetGlyphRangesChineseFull());
    io.Fonts->AddFontFromFileTTF("D:\\imgui\\imgui-master\\misc\\fonts\\Roboto-Medium.ttf", 22, NULL, io.Fonts->GetGlyphRangesChineseFull());
    io.Fonts->AddFontFromFileTTF((font + "Calibri.ttf").c_str(), 22, NULL, io.Fonts->GetGlyphRangesChineseFull());

    ImGui_ImplGlfw_InitForOpenGL(Window, true);
    ImGui_ImplOpenGL3_Init("#version 330");
    return Window;
}

std::string formatDoubleValue(int val) {
    std::ostringstream oss;
    oss << val;
    return oss.str();
}

std::string formatDoubleValue(double val) {
    std::ostringstream oss;
    oss << val;
    return oss.str();
}

#endif