if(BUILD_SHARED_LIBS)
    add_library(libcon2020 SHARED)
else()
    add_library(libcon2020 STATIC)
endif()

add_library(con2020::con2020 ALIAS libcon2020)

set_target_properties(libcon2020 PROPERTIES OUTPUT_NAME con2020)
set_target_properties(libcon2020 PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

target_sources(libcon2020 PRIVATE
    bessel.cc
    con2020.cc
    flux.cc
    libcon2020.cc
    lmic.cc
    polyeval.cc
    scalarpot.cc
    smoothd.cc
    trap.cc
    trimstring.cc
)

install(TARGETS libcon2020
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    RUNTIME DESTINATION bin
)
