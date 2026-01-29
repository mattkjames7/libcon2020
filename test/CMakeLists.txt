include(FetchContent)

FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/releases/download/v1.17.0/googletest-1.17.0.tar.gz
  DOWNLOAD_EXTRACT_TIMESTAMP TRUE
)

FetchContent_Declare(
  nlohmann_json
  URL https://github.com/nlohmann/json/archive/refs/tags/v3.12.0.tar.gz
  DOWNLOAD_EXTRACT_TIMESTAMP TRUE
)

# Prevent GoogleTest from being installed
set(INSTALL_GTEST OFF CACHE BOOL "Disable installation of GoogleTest" FORCE)
set(BUILD_GMOCK OFF CACHE BOOL "Disable GoogleMock" FORCE)


set(_OLD_BSL "${BUILD_SHARED_LIBS}")
set(BUILD_SHARED_LIBS OFF CACHE BOOL "" FORCE)

FetchContent_MakeAvailable(googletest)
FetchContent_MakeAvailable(nlohmann_json)

# Restore parent setting
set(BUILD_SHARED_LIBS "${_OLD_BSL}" CACHE BOOL "" FORCE)

find_package(ZLIB REQUIRED)

set(TEST_INCLUDE_DIR "${CMAKE_SOURCE_DIR}/include")

set(TEST_DATA_FILES
  "${CMAKE_CURRENT_SOURCE_DIR}/test_data.json.gz"
  "${CMAKE_CURRENT_SOURCE_DIR}/bessel_data.json.gz"
)

file(COPY ${TEST_DATA_FILES} DESTINATION "${CMAKE_CURRENT_BINARY_DIR}")

add_executable(save_test_data
  save_test_data.cc
)
target_include_directories(save_test_data PRIVATE "${TEST_INCLUDE_DIR}")
target_link_libraries(save_test_data PRIVATE con2020::con2020 nlohmann_json::nlohmann_json ZLIB::ZLIB)
if(UNIX AND NOT APPLE)
  target_link_libraries(save_test_data PRIVATE m)
endif()

add_custom_target(save
  COMMAND save_test_data
  DEPENDS save_test_data
  WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
)

add_executable(model_output_tests
  unit/model_output.cc
)
target_include_directories(model_output_tests PRIVATE "${TEST_INCLUDE_DIR}")
target_link_libraries(model_output_tests PRIVATE con2020::con2020 GTest::gtest GTest::gtest_main nlohmann_json::nlohmann_json ZLIB::ZLIB)
if(UNIX AND NOT APPLE)
  target_link_libraries(model_output_tests PRIVATE m)
endif()

add_executable(bessel_tests
  unit/bessel.cc
)
target_include_directories(bessel_tests PRIVATE "${TEST_INCLUDE_DIR}")
target_link_libraries(bessel_tests PRIVATE con2020::con2020 GTest::gtest GTest::gtest_main nlohmann_json::nlohmann_json ZLIB::ZLIB)
if(UNIX AND NOT APPLE)
  target_link_libraries(bessel_tests PRIVATE m)
endif()

set(_test_rpath "")
if(APPLE)
  set(_test_rpath "@loader_path/../lib")
elseif(UNIX)
  set(_test_rpath "\$ORIGIN/../lib")
endif()

if(_test_rpath)
  set_target_properties(save_test_data model_output_tests bessel_tests PROPERTIES
    BUILD_RPATH "${_test_rpath}"
  )
endif()

# On Windows the test executable won't find DLLs located in other directories
# (RPATH doesn't apply). If we built shared libs, copy the spline DLL
# next to the test executable so it can be loaded at runtime in CI.
if(BUILD_SHARED_LIBS AND WIN32)
    # Copy the spline DLL to the test exe dir
    add_custom_command(TARGET model_output_tests POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
            $<TARGET_FILE:libcon2020>
            $<TARGET_FILE_DIR:model_output_tests>
    )
    add_custom_command(TARGET bessel_tests POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
            $<TARGET_FILE:libcon2020>
            $<TARGET_FILE_DIR:bessel_tests>
    )
    # Copy GoogleTest DLLs to the test exe dir (they're in bin/)
    add_custom_command(TARGET model_output_tests POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
            $<TARGET_FILE:gtest>
            $<TARGET_FILE_DIR:model_output_tests>
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
            $<TARGET_FILE:gtest_main>
            $<TARGET_FILE_DIR:model_output_tests>
    )
    add_custom_command(TARGET bessel_tests POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
            $<TARGET_FILE:gtest>
            $<TARGET_FILE_DIR:bessel_tests>
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
            $<TARGET_FILE:gtest_main>
            $<TARGET_FILE_DIR:bessel_tests>
    )
    
    # On Windows with shared libs, skip automatic test discovery since DLLs
    # aren't available during CMake configuration. Register the test manually instead.
    add_test(NAME model_output_tests COMMAND model_output_tests WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")
    add_test(NAME bessel_tests COMMAND bessel_tests WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")
else()
    include(GoogleTest)
    gtest_discover_tests(model_output_tests WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")
    gtest_discover_tests(bessel_tests WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")
endif()