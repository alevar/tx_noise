cmake_minimum_required(VERSION 3.10)
project(get_results)

set (CMAKE_STATIC_LINKER_FLAGS "-Wl,--as-needed.-lcurl")

include(ExternalProject)

# thanks for the macro: David Grayson
macro(use_cxx11)
    if (CMAKE_VERSION VERSION_LESS "3.1")
        if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
            set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++11")
        endif ()
    else ()
        set (CMAKE_CXX_STANDARD 11)
    endif ()
endmacro(use_cxx11)

#use_cxx11()

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++11")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")

include_directories(src)

add_executable(get_results
        src/arg_parse.cpp
        src/arg_parse.h
        src/codons.cpp
        src/codons.h
        src/GBase.cpp
        src/GBase.h
        src/gdna.cpp
        src/gdna.h
        src/GFaSeqGet.cpp
        src/GFaSeqGet.h
        src/GHash.hh
        src/GList.hh
        src/GVec.hh
        src/gff.h
        src/gff.cpp
        src/tokenize.cpp
        src/FastaTools.cpp
        src/FastaTools.h
        src/QuantRes.cpp
        src/QuantRes.h
        main.cpp)

target_link_libraries(get_results -lbz2 -llzma -pthread -lcrypto -lcurl -lz)