cmake_minimum_required(VERSION 3.10)
project(raytracingDemo)


add_executable(${PROJECT_NAME} main.cpp )

set(INCLUDE_DIR ${PROJECT_SOURCE_DIR}/include)
set(SRC_DIR ${PROJECT_SOURCE_DIR}/src)

include_directories(
	${INCLUDE_DIR}
    	${INCLUDE_DIR}/GL
	${INCLUDE_DIR}/imgui
)

target_link_libraries(raytracingDemo ${PROJECT_SOURCE_DIR}/lib/glfw3.lib) 
target_link_libraries(raytracingDemo ${PROJECT_SOURCE_DIR}/lib/tinyobjloader.lib) 

target_sources(raytracingDemo PUBLIC ${SRC_DIR}/glad.c)
target_sources(raytracingDemo PUBLIC ${SRC_DIR}/hdrloader.cpp)
target_sources(raytracingDemo PUBLIC ${SRC_DIR}/My_Shader.cpp)
target_sources(raytracingDemo PUBLIC ${SRC_DIR}/tiny_obj_loader.cpp)
target_sources(raytracingDemo PUBLIC ${SRC_DIR}/stb_image.cpp)
target_sources(raytracingDemo PUBLIC ${SRC_DIR}/tiny_gltf.cc)
target_sources(raytracingDemo PUBLIC ${SRC_DIR}/imgui/imgui.cpp)
target_sources(raytracingDemo PUBLIC ${SRC_DIR}/imgui/imgui_demo.cpp)
target_sources(raytracingDemo PUBLIC ${SRC_DIR}/imgui/imgui_draw.cpp)
target_sources(raytracingDemo PUBLIC ${SRC_DIR}/imgui/imgui_impl_glfw.cpp)
target_sources(raytracingDemo PUBLIC ${SRC_DIR}/imgui/imgui_impl_opengl3.cpp)
target_sources(raytracingDemo PUBLIC ${SRC_DIR}/imgui/imgui_tables.cpp)
target_sources(raytracingDemo PUBLIC ${SRC_DIR}/imgui/imgui_widgets.cpp)




