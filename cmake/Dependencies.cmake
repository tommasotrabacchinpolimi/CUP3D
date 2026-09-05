# Look for HDF5 and GSL, all required dependencies.
# If not found, but if the CUP_AUTO_INSTALL_DEPENDENCIES is set, the
# ./install_dependencies.sh script will be run automatically.

set(_CUP_HDF5_ROOT "${DEP_BUILD_DIR}/hdf5-1.10.1-parallel")
set(_CUP_GSL_ROOT "${DEP_BUILD_DIR}/gsl-2.1")

set(_missing_dep)
set(_dep_install_flags)

# Update search dirs (immediately and after running ./install_dependencies.sh).
macro(_update_dirs)
    if (EXISTS "${_CUP_HDF5_ROOT}")
        set(HDF5_ROOT "${_CUP_HDF5_ROOT}" CACHE PATH "Prefix for HDF5 installation" FORCE)
    endif()
    if (EXISTS "${_CUP_GSL_ROOT}")
        set(GSL_ROOT_DIR "${_CUP_GSL_ROOT}" CACHE PATH "Prefix for GSL installation" FORCE)
        # Clear any stale cached results from previous CMake runs.
        unset(GSL_INCLUDE_DIR CACHE)
        unset(GSL_LIBRARY CACHE)
        unset(GSL_CBLAS_LIBRARY CACHE)
        unset(GSL_VERSION CACHE)
    endif()
endmacro()

# Helper for adding a package to the list of missing dependency if it is not found.
macro(_find_package package flag)
    find_package(${package})
    if (NOT ${package}_FOUND)
        list(APPEND _missing_dep ${package})
        list(APPEND _dep_install_flags "${flag}")
    endif()
endmacro()

_update_dirs()

set(HDF5_PREFER_PARALLEL ON)
set(HDF5_NO_FIND_PACKAGE_CONFIG_FILE ON)
_find_package(HDF5 "--hdf5")

_find_package(GSL "--gsl")

if (_missing_dep)
    string(JOIN " " _dep_install_flags_str ${_dep_install_flags})
    if (NOT CUP_AUTO_INSTALL_DEPENDENCIES)
        string(JOIN ", " _missing_dep_str ${_missing_dep})
        message(FATAL_ERROR
                "One or more dependencies not found: ${_missing_dep_str}. "
                "Dependencies can be installed by running\n"
                "    ./install_dependencies.sh ${_dep_install_flags_str}\n"
                "from the repository root, or by rerunning cmake with "
                "-DCUP_AUTO_INSTALL_DEPENDENCIES=ON")
    endif()
    message(STATUS "Installing dependencies: ${_dep_install_flags_str}")
    execute_process(
        COMMAND "./install_dependencies.sh" ${_dep_install_flags}
        WORKING_DIRECTORY "${ROOT_DIR}"
        ERROR_VARIABLE error
        RESULT_VARIABLE error_code)
    if (error_code)
        message(FATAL_ERROR "Installing dependencies ${_dep_install_flags_str} failed:\n${error}")
    endif()
    message(STATUS "./install_dependencies.sh ${_dep_install_flags_str} completed successfully.")
    _update_dirs()
    foreach (_dep IN LISTS _missing_dep)
        if (_dep STREQUAL "HDF5")
            find_package(HDF5 REQUIRED)
        elseif (_dep STREQUAL "GSL")
            find_package(GSL REQUIRED)
        endif()
    endforeach()
endif()

if (NOT HDF5_FOUND)
    find_package(HDF5 REQUIRED)
endif()
if (NOT GSL_FOUND)
    find_package(GSL REQUIRED)
endif()
