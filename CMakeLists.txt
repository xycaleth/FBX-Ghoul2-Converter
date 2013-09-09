cmake_minimum_required(VERSION 2.8)
project("FBX-G2 Converter")

set(SOURCE_FILES
	"src/EntryPoint.cpp"
)

if(APPLE)
	set(FBX_ROOT "/Applications/Autodesk/FBX\ SDK/2014.1")
endif()

set(INCLUDE_PATH
	${INCUDE_DIR}
	${FBX_ROOT}/include
)

set(LINK_LIBRARIES
	${LINK_LIBRARIES}
	${FBX_ROOT}/lib/gcc4/ub/debug/libfbxsdk.a
	/System/Library/Frameworks/Carbon.framework
)

include_directories(${INCLUDE_PATH})

add_executable(fbx2ghoul2 ${SOURCE_FILES})
target_link_libraries(fbx2ghoul2 ${LINK_LIBRARIES})
