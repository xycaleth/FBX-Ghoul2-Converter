cmake_minimum_required(VERSION 2.8)
project("FBX-G2 Converter")

set(SOURCE_FILES
	"src/EntryPoint.cpp"
	"src/GLA.cpp"
	"src/IO.cpp"
	)

set(HEADER_FILES
	"src/GLA.h"
	"src/GLM.h"
	"src/IO.h"
	)

if(NOT FBX_ROOT)
	if(APPLE)
		set(FBX_ROOT "/Applications/Autodesk/FBX\ SDK/2014.1")
	elseif(WIN32)
		set(FBX_ROOT "$ENV{ProgramFiles}\\Autodesk\\FBX\\FBX SDK\\2014.1")
	endif()
endif()

set(INCLUDE_PATH
	${INCUDE_PATH}
	${FBX_ROOT}/include
)

if(APPLE)
	set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} "-stdlib=libstdc++ -std=c++11 -Wno-reserved-user-defined-literal")
	set(LINK_LIBRARIES
		${LINK_LIBRARIES}
		/System/Library/Frameworks/Carbon.framework
		debug ${FBX_ROOT}/lib/gcc4/ub/debug/libfbxsdk.a
		optimized ${FBX_ROOT}/lib/gcc4/ub/release/libfbxsdk.a
	)
elseif(WIN32)
	if(MSVC)
		add_definitions(-D_CRT_SECURE_NO_WARNINGS)
		add_definitions(-D_SCL_SECURE_NO_WARNINGS)
	endif()

	set(LINK_LIBRARIES
		${LINK_LIBRARIES}
		debug ${FBX_ROOT}/lib/vs2012/x86/debug/libfbxsdk-md.lib
		optimized ${FBX_ROOT}/lib/vs2012/x86/release/libfbxsdk-md.lib
	)
endif()

include_directories(${INCLUDE_PATH})

add_executable(fbx2ghoul2 ${SOURCE_FILES} ${HEADER_FILES})
target_link_libraries(fbx2ghoul2 ${LINK_LIBRARIES})
