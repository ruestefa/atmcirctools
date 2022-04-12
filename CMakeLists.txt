cmake_minimum_required(VERSION 3.22.0)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if (NOT DEFINED($ENV{CMAKE_PYTHON_VERSION}))
    set(CMAKE_PYTHON_VERSION 3.9)
endif()
if (DEFINED($ENV{CMAKE_PREFIX_PATH}))
    set(CMAKE_PREFIX_PATH $ENV{CMAKE_PREFIX_PATH})
elseif (NOT DEFINED CMAKE_PREFIX_PATH)
    message(
        FATAL_ERROR
        "Must set CMAKE_PREFIX_PATH to the root of your Python environment.\n"
    )
endif()

set(Python_FIND_STRATEGY LOCATION)
set(Python_ROOT ${CMAKE_PREFIX_PATH})

# Manually set numpy include dir to ensure that `find_package(NumPy, ...)` below finds the correct
# installation (may otherwise, e.g., chose system installation over active conda environment)
set(NumPy_INCLUDE_DIR ${CMAKE_PREFIX_PATH}/lib/python${CMAKE_PYTHON_VERSION}/site-packages/numpy/core/include)

message("-- Setup:")
message("   CMAKE_PYTHON_VERSION=${CMAKE_PYTHON_VERSION}")
message("   CMAKE_PREFIX_PATH=${CMAKE_PREFIX_PATH}")
message("   Python_FIND_STRATEGY=${Python_FIND_STRATEGY}")
message("   Python_ROOT=${Python_ROOT}")
message("   NumPy_INCLUDE_DIR=${NumPy_INCLUDE_DIR}")
message("--")

project(
    dpv_th
    VERSION 0.0.0
    DESCRIPTION "Compute gradient of PV on TH surfaces."
    LANGUAGES C Fortran
)

# Safety net
if(PROJECT_SOURCE_DIR STREQUAL PROJECT_BINARY_DIR)
    message(
        FATAL_ERROR
        "In-source builds not allowed. Please make a new directory (called a build directory) and run CMake from there.\n"
    )
endif()

# Grab Python
find_package(Python3 ${CMAKE_PYTHON_VERSION} REQUIRED COMPONENTS Interpreter Development NumPy)

# Ensure scikit-build modules
if (NOT SKBUILD)
    # Kanged --> https://github.com/Kitware/torch_liberator/blob/master/CMakeLists.txt
    # If skbuild is not the driver; include its utilities in CMAKE_MODULE_PATH
    execute_process(
        COMMAND "${Python3_EXECUTABLE}"
        -c "import os, skbuild; print(os.path.dirname(skbuild.__file__))"
        OUTPUT_VARIABLE SKBLD_DIR
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    set(SKBLD_CMAKE_DIR "${SKBLD_DIR}/resources/cmake")
    list(APPEND CMAKE_MODULE_PATH ${SKBLD_CMAKE_DIR})
endif()

# `scikit-build`-style includes
find_package(PythonExtensions REQUIRED) # for ${PYTHON_EXTENSION_MODULE_SUFFIX}
find_package(PythonExtensions REQUIRED) # for ${PYTHON_EXTENSION_MODULE_SUFFIX}
message(Python_INCLUDE_DIRS=${Python_INCLUDE_DIRS})
message(Python_SITELIB=${Python_SITELIB})
find_package(NumPy REQUIRED MODULE) # for ${NumPy_INCLUDE_DIRS}
find_package(F2PY REQUIRED) # for ${F2PY_INCLUDE_DIR}

#-------------------------------------------------------------------------------
# Extension module form _deriv_f77.f
#-------------------------------------------------------------------------------

# Prepping the module
set(f2py_mod_name_f77 "_deriv_f77")
set(f2py_mod_dir_f77 "src/dpv_th/deriv")
set(f_src_file_f77 "${CMAKE_SOURCE_DIR}/${f2py_mod_dir_f77}/${f2py_mod_name_f77}.f")
set(f2py_cmod_file_f77 "${f2py_mod_dir_f77}/${f2py_mod_name_f77}module.c")
set(ext_mod_file_f77 "${f2py_mod_name_f77}${PYTHON_EXTENSION_MODULE_SUFFIX}")

# Target for enforcing dependencies
add_custom_target(${f2py_mod_name_f77} ALL DEPENDS "${f_src_file_f77}")

# Custom command for generating .c
add_custom_command(
    OUTPUT "${f2py_cmod_file_f77}"
    COMMAND ${F2PY_EXECUTABLE}
        -m ${f2py_mod_name_f77}
        ${f_src_file_f77}
        --build-dir ${f2py_mod_dir_f77}
        --lower
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    DEPENDS ${f_src_file_f77}
)

add_library(
    ${ext_mod_file_f77} MODULE
    "${f2py_cmod_file_f77}"
    "${F2PY_INCLUDE_DIR}/fortranobject.c"
    "${f_src_file_f77}"
)

target_include_directories(
    ${ext_mod_file_f77} PUBLIC
    ${F2PY_INCLUDE_DIRS}
    ${PYTHON_INCLUDE_DIRS}
)
set_target_properties(${ext_mod_file_f77} PROPERTIES SUFFIX "")
set_target_properties(${ext_mod_file_f77} PROPERTIES PREFIX "")

# Linker fixes
if (UNIX)
    if (APPLE)
        set_target_properties(
            ${ext_mod_file_f77} PROPERTIES
            LINK_FLAGS "-Wl,-dylib,-undefined,dynamic_lookup"
        )
    else()
        set_target_properties(
            ${ext_mod_file_f77} PROPERTIES
            LINK_FLAGS "-Wl,--allow-shlib-undefined"
        )
    endif()
endif()

if (SKBUILD)
    install(TARGETS ${ext_mod_file_f77} DESTINATION ${f2py_mod_dir_f77})
else()
    install(TARGETS ${ext_mod_file_f77} DESTINATION ${CMAKE_SOURCE_DIR}/${f2py_mod_dir_f77})
endif()

#-------------------------------------------------------------------------------
# Extension module form _deriv_f90.f90
#-------------------------------------------------------------------------------

# Prepping the module
set(f2py_mod_name_f90 "_deriv_f90")
set(f2py_mod_dir_f90 "src/dpv_th/deriv")
set(f_src_file_f90 "${CMAKE_SOURCE_DIR}/${f2py_mod_dir_f90}/${f2py_mod_name_f90}.f90")
set(f2py_cmod_file_f90 "${f2py_mod_dir_f90}/${f2py_mod_name_f90}module.c")
set(ext_mod_file_f90 "${f2py_mod_name_f90}${PYTHON_EXTENSION_MODULE_SUFFIX}")

# Target for enforcing dependencies
add_custom_target(${f2py_mod_name_f90} ALL DEPENDS "${f_src_file_f90}")

# Custom command for generating .c
add_custom_command(
    OUTPUT "${f2py_cmod_file_f90}"
    COMMAND ${F2PY_EXECUTABLE}
        -m ${f2py_mod_name_f90}
        ${f_src_file_f90}
        --build-dir ${f2py_mod_dir_f90}
        --lower
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    DEPENDS ${f_src_file_f90}
)

add_library(
    ${ext_mod_file_f90} MODULE
    "${f2py_cmod_file_f90}"
    "${F2PY_INCLUDE_DIR}/fortranobject.c"
    "${f_src_file_f90}"
)

target_include_directories(
    ${ext_mod_file_f90} PUBLIC
    ${F2PY_INCLUDE_DIRS}
    ${PYTHON_INCLUDE_DIRS}
)
set_target_properties(${ext_mod_file_f90} PROPERTIES SUFFIX "")
set_target_properties(${ext_mod_file_f90} PROPERTIES PREFIX "")

# Linker fixes
if (UNIX)
    if (APPLE)
        set_target_properties(
            ${ext_mod_file_f90} PROPERTIES
            LINK_FLAGS "-Wl,-dylib,-undefined,dynamic_lookup"
        )
    else()
        set_target_properties(
            ${ext_mod_file_f90} PROPERTIES
            LINK_FLAGS "-Wl,--allow-shlib-undefined"
        )
    endif()
endif()

if (SKBUILD)
    install(TARGETS ${ext_mod_file_f90} DESTINATION ${f2py_mod_dir_f90})
else()
    install(TARGETS ${ext_mod_file_f90} DESTINATION ${CMAKE_SOURCE_DIR}/${f2py_mod_dir_f90})
endif()
