cmake_minimum_required(VERSION 3.16)

project(Ditto LANGUAGES CXX)

if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX "${CMAKE_SOURCE_DIR}" CACHE PATH
      "Installation prefix" FORCE)
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

find_package(ROOT CONFIG REQUIRED Core RIO Tree Hist Physics EG)
find_package(CURL REQUIRED)

set(PYTHIA_ROOT "$ENV{PYTHIA_ROOT}" CACHE PATH "Pythia8 installation prefix")
if(NOT PYTHIA_ROOT)
  message(FATAL_ERROR "PYTHIA_ROOT must point to the Pythia8 installation prefix")
endif()

find_path(PYTHIA8_INCLUDE_DIR
  NAMES Pythia8/Pythia.h
  PATHS "${PYTHIA_ROOT}/include"
  NO_DEFAULT_PATH
)
find_library(PYTHIA8_LIBRARY
  NAMES pythia8
  PATHS "${PYTHIA_ROOT}/lib" "${PYTHIA_ROOT}/lib64"
  NO_DEFAULT_PATH
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Pythia8
  REQUIRED_VARS PYTHIA8_INCLUDE_DIR PYTHIA8_LIBRARY
)

add_library(Ditto::Pythia8 UNKNOWN IMPORTED)
set_target_properties(Ditto::Pythia8 PROPERTIES
  IMPORTED_LOCATION "${PYTHIA8_LIBRARY}"
  INTERFACE_INCLUDE_DIRECTORIES "${PYTHIA8_INCLUDE_DIR}"
)
set(Ditto_PYTHIA_TARGET Ditto::Pythia8)

add_library(Ditto SHARED
  Generation/Ditto.cxx
  Tuning/DittoTune.cxx
  Tuning/DittoTuneDownload.cxx
)

add_library(Ditto::Ditto ALIAS Ditto)

ROOT_GENERATE_DICTIONARY(DittoDict
  ${CMAKE_CURRENT_SOURCE_DIR}/Tuning/DittoTune.h
  LINKDEF ${CMAKE_CURRENT_SOURCE_DIR}/Tuning/LinkDef.h
  MODULE Ditto
)

target_include_directories(Ditto
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Generation>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Tuning>
    $<INSTALL_INTERFACE:include>
)

target_link_libraries(Ditto
  PUBLIC
    ROOT::Core
    ROOT::RIO
    ROOT::Tree
    ROOT::Hist
    ROOT::Physics
    ROOT::EG
    ${Ditto_PYTHIA_TARGET}
  PRIVATE
    CURL::libcurl
)

add_library(DittoTuner SHARED Tuning/DittoPythiaTuner.cxx
                              Tuning/DittoTuneAccumulator.cxx)
add_library(Ditto::DittoTuner ALIAS DittoTuner)

target_include_directories(DittoTuner
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Generation>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Tuning>
    $<INSTALL_INTERFACE:include>
)

target_link_libraries(DittoTuner
  PUBLIC
    Ditto::Ditto
    ${Ditto_PYTHIA_TARGET}
)

install(TARGETS Ditto DittoTuner
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  RUNTIME DESTINATION bin
)

install(DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/"
  DESTINATION lib
  FILES_MATCHING
    PATTERN "*.pcm"
    PATTERN "*.rootmap"
)

install(FILES
  Generation/Ditto.h
  Tuning/DittoTune.h
  Tuning/DittoPythiaTuner.h
  Tuning/DittoTuneAccumulator.h
  Tuning/DittoTuneDownload.h
  DESTINATION include/Ditto
)
