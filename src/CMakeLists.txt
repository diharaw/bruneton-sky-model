cmake_minimum_required(VERSION 3.8 FATAL_ERROR)

find_program(CLANG_FORMAT_EXE NAMES "clang-format" DOC "Path to clang-format executable")

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)

file(GLOB_RECURSE ATMOSPHERE_HEADERS ${PROJECT_SOURCE_DIR}/src/*.h)
file(GLOB_RECURSE ATMOSPHERE_SOURCES ${PROJECT_SOURCE_DIR}/src/*.cpp)

if (EMSCRIPTEN)
    set(CMAKE_EXECUTABLE_SUFFIX ".html")
endif()

if(APPLE)
    add_executable(AtmosphericScatttering MACOSX_BUNDLE ${ATMOSPHERE_HEADERS} ${ATMOSPHERE_SOURCES})
    set(MACOSX_BUNDLE_BUNDLE_NAME "com.dihara.atmosphere") 
else()
    add_executable(AtmosphericScatttering ${ATMOSPHERE_HEADERS} ${ATMOSPHERE_SOURCES}) 
endif()

if(CLANG_FORMAT_EXE)
    add_custom_target(clang-format-project-files COMMAND ${CLANG_FORMAT_EXE} -i -style=file ${ATMOSPHERE_HEADERS} ${ATMOSPHERE_SOURCES})
endif()

target_link_libraries(AtmosphericScatttering dwSampleFramework)