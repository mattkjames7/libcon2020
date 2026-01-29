include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/releases/download/v1.17.0/googletest-1.17.0.tar.gz
  DOWNLOAD_EXTRACT_TIMESTAMP TRUE
)

# Prevent GoogleTest from being installed
set(INSTALL_GTEST OFF CACHE BOOL "Disable installation of GoogleTest" FORCE)
set(BUILD_GMOCK OFF CACHE BOOL "Disable GoogleMock" FORCE)


set(_OLD_BSL "${BUILD_SHARED_LIBS}")
set(BUILD_SHARED_LIBS OFF CACHE BOOL "" FORCE)

FetchContent_MakeAvailable(googletest)

# Restore parent setting
set(BUILD_SHARED_LIBS "${_OLD_BSL}" CACHE BOOL "" FORCE)

